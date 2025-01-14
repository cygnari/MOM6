module MOM_conv_self_attr_load

use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_grid,            only : ocean_grid_type, get_global_grid_size
use MOM_file_parser,     only : get_param
use MOM_coms_infra,      only : sum_across_PEs, max_across_PEs, num_PEs, PE_here, broadcast, sync_PEs
use MOM_coms_infra,      only : send_to_PE, recv_from_PE
use MOM_coms,            only : reproducing_sum
use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_domains,         only : pass_var

implicit none ; private

public sal_conv_init, sal_conv_eval

#include <MOM_memory.h>

type, private :: cube_panel
    ! data type to represent a panel of the cubed sphere in the tree structure
    integer :: level = 0 ! level 0 is the base level
    logical :: is_leaf = .true. ! all leaves start out as leaves, become not leaf when refined
    integer :: id
    integer :: parent_panel = -1
    integer :: child_panel_1 = -1
    integer :: child_panel_2 = -1
    integer :: child_panel_3 = -1
    integer :: child_panel_4 = -1
    integer :: face 
    real :: min_xi ! xi and eta are the angle coordinates on the face of the cube
    real :: mid_xi ! based on the equiangular gnomonic cubed sphere
    real :: max_xi
    real :: min_eta
    real :: mid_eta
    real :: max_eta
    real :: radius ! distance from center of panel to corner
    integer, allocatable :: points_inside(:) ! indices of contained points
    integer, allocatable :: relabeled_points_inside(:)
    integer :: panel_point_count = 0

    contains
        procedure :: contains_point
end type cube_panel

type, private :: interaction_pair 
    ! data type to represent a pp/pc interaction between a point and a cube panel
    integer :: index_target
    integer :: index_source 
    integer :: interact_type ! 0 for PP, 1 for PC
end type interaction_pair

type, public :: SAL_conv_type ; private
    ! type to contain all the needed data structures
    ! including the cubed sphere tree
    ! and communication things
    logical :: reprod_sum !< True if use reproducible global sums
    real, allocatable :: e_xs(:), e_ys(:), e_zs(:) ! x/y/z coordinates of unowned points needed for particle-particle interactions
    type(interaction_pair), allocatable :: pp_interactions(:), pc_interactions(:) ! interaction lists
    type(cube_panel), allocatable :: tree_struct(:)
    integer, allocatable :: points_panels(:,:) ! points_panels(lev, i)=k means that point i is contained in panel k at level lev
    integer :: p ! total ranks
    integer :: id ! rank
    integer, allocatable :: points_to_give_proc(:) ! communication patterns
    integer, allocatable :: points_to_get_proc(:) ! second index is rank
    integer, allocatable :: points_to_give(:,:), points_to_get(:,:), unowned_source_points(:)
    integer :: unowned_sources
    integer :: interp_degree 
    integer, allocatable :: indexsg(:), indexeg(:), pcg(:)
    integer, allocatable :: two_d_to_1d(:,:) ! maps between 1d and 2d
    integer, allocatable :: one_d_to_2d_i(:), one_d_to_2d_j(:)
    integer :: own_ocean_points
    integer :: total_ocean_points
end type SAL_conv_type

integer :: id_clock_SAL   !< CPU clock for self-attraction and loading
integer :: id_clock_SAL_pc
integer :: id_clock_SAL_pp

contains

integer function face_from_xyz(x, y, z) 
    ! takes (x, y, z) coordinates of a point and returns which face of the cubed sphere it corresponds to
    real, intent(in) :: x, y, z
    real :: ax, ay, az
    integer :: face
    ax = abs(x)
    ay = abs(y)
    az = abs(z)
    face = 0
    if ((ax >= ay) .and. (ax >= az)) then
        if (x >= 0) then
        face = 1
        else 
        face = 3
        end if
    else if ((ay >= ax) .and. (ay >= az)) then
        if (y >= 0) then
        face = 2
        else
        face = 4
        end if
    else
        if (z >= 0) then
        face = 5
        else
        face = 6
        end if
    end if
    face_from_xyz = face
    return
end function face_from_xyz

subroutine xieta_from_xyz(x, y, z, xi, eta, in_face)
    ! computes xi eta angle coordinates on the cubed sphere from x, y, z points and optional face
    ! based on Ronchi et al, 1996 cubed sphere
    real, intent(in) :: x, y, z
    integer, intent(in), optional :: in_face
    real, intent(out) :: xi, eta
    integer :: face
    IF (.not. present(in_face)) THEN
        face = face_from_xyz(x, y, z)
    ELSE 
        face = in_face
    END IF
    IF (face == 1) THEN
        xi = ATAN(y/x); eta = ATAN(z/x)
    ELSE IF (face == 2) THEN
        xi = ATAN(-x/y); eta = ATAN(z/y)
    ELSE IF (face == 3) THEN
        xi = ATAN(y/x); eta = ATAN(-z/x)
    ELSE IF (face == 4) THEN
        xi = ATAN(-x/y); eta = ATAN(-z/y)
    ELSE IF (face == 5) THEN
        xi = ATAN(y/z); eta = ATAN(-x/z)
    ELSE IF (face == 6) THEN
        xi = ATAN(-y/z); eta = ATAN(-x/z)
    END IF
end subroutine xieta_from_xyz

SUBROUTINE xyz_from_xieta(x, y, z, xi, eta, face) 
    ! converts xi eta on the cubed sphere to x y z
    ! based on Ronchi et al, 1996 cubed sphere
    real, intent(in) :: xi, eta
    integer, intent(in) :: face
    real, intent(out) :: x, y, z
    real :: ax, ay
    ax = TAN(xi); ay = TAN(eta)
    IF (face == 1) THEN
        x = 1.0 / SQRT(1+ax*ax+ay*ay)
        y = ax*x
        z = ay*x
    ELSE IF (face == 2) THEN
        y = 1.0 / SQRT(1+ax*ax+ay*ay)
        x = -ax*y
        z = ay*y
    ELSE IF (face == 3) THEN
        x = -1.0 / SQRT(1+ax*ax+ay*ay)
        y = ax*x
        z = -ay*x
    ELSE IF (face == 4) THEN
        y = -1.0 / SQRT(1+ax*ax+ay*ay)
        x = -ax*y
        z = -ay*y
    ELSE IF (face == 5) THEN
        z = 1.0 / SQRT(1+ax*ax+ay*ay)
        x = -ay*z
        y = ax*z
    ELSE 
        z = -1.0 / SQRT(1+ax*ax+ay*ay)
        x = -ay*z
        y = -ax*z
    END IF
end subroutine xyz_from_xieta

logical function contains_point(self, x, y, z) result(contains)
    class(cube_panel), intent(in) :: self
    real, intent(in) :: x, y, z
    integer :: face
    real :: xi, eta
    face = face_from_xyz(x, y, z)
    if (face == self%face) then
        call xieta_from_xyz(x, y, z, xi, eta, face)
        if ((xi >= self%min_xi) .and. (xi < self%max_xi) .and. (eta >= self%min_eta) &
                        .and. (eta < self%max_eta)) then
            contains = .true.
        else
            contains = .false.
        end if
    else 
        contains = .false.
    end if
end function contains_point

subroutine tree_traversal(G, tree_panels, xg, yg, zg, cluster_thresh, point_count)
    ! constructs cubed sphere tree of points
    type(ocean_grid_type), intent(inout) :: G ! ocean grid
    type(cube_panel), allocatable, intent(out) :: tree_panels(:)
    type(cube_panel), allocatable :: tree_panels_temp(:)
    integer, intent(in) :: cluster_thresh, point_count
    real, intent(in) :: xg(:), yg(:), zg(:)
    real :: pi, xval, yval, zval, min_xi, mid_xi, max_xi, min_eta, mid_eta, max_eta, xi, eta
    integer, allocatable :: curr_loc(:), temp(:)
    integer :: face, i, panel_count, j, count, index, k
    integer :: which_panel, total
    real :: x1, x2, x3, y1, y2, y3, d1, d2, d3, d4

    pi = 4.D0*DATAN(1.D0)
    allocate(tree_panels_temp(max(6, point_count)))
    allocate(curr_loc(6))

    ! initialize the six top level cube panels
    do i = 1, 6
        tree_panels_temp(i)%id = i
        tree_panels_temp(i)%face = i
        tree_panels_temp(i)%min_xi = -pi/4.D0
        tree_panels_temp(i)%mid_xi = 0.D0
        tree_panels_temp(i)%max_xi = pi/4.D0
        tree_panels_temp(i)%min_eta = -pi/4.D0
        tree_panels_temp(i)%mid_eta = 0.D0
        tree_panels_temp(i)%max_eta = pi/4.D0
        allocate(tree_panels_temp(i)%points_inside(point_count))
        tree_panels_temp(i)%panel_point_count = 0
        curr_loc(i) = 0
    enddo

    do i = 1, point_count
        xval = xg(i)
        yval = yg(i)
        zval = zg(i)
        face = face_from_xyz(xval, yval, zval)
        curr_loc(face) = curr_loc(face) + 1
        tree_panels_temp(face)%panel_point_count = tree_panels_temp(face)%panel_point_count + 1
        tree_panels_temp(face)%points_inside(curr_loc(face)) = i
    enddo

    do i = 1, 6
        allocate(temp(curr_loc(i)))
        do j = 1, tree_panels_temp(i)%panel_point_count
            temp(j) = tree_panels_temp(i)%points_inside(j)
        enddo
        call move_alloc(from=temp, to=tree_panels_temp(i)%points_inside)
        curr_loc(i) = 0
    enddo

    i = 1
    panel_count = 6
    do while (i <= panel_count)
        ! first check if panel needs to be divided
        count = tree_panels_temp(i)%panel_point_count
        IF ((count >= cluster_thresh) .and. (tree_panels_temp(i)%is_leaf)) THEN
            ! panel is a leaf and has many points => refine
            tree_panels_temp(i)%is_leaf = .false.
            min_xi = tree_panels_temp(i)%min_xi
            mid_xi = tree_panels_temp(i)%mid_xi
            max_xi = tree_panels_temp(i)%max_xi
            min_eta = tree_panels_temp(i)%min_eta
            mid_eta = tree_panels_temp(i)%mid_eta
            max_eta = tree_panels_temp(i)%max_eta
            ! four new panels are index panel_count+1 through panel_count+4
            DO j = 1, 4
                index = panel_count+j
                tree_panels_temp(panel_count+j)%parent_panel = i
                tree_panels_temp(panel_count+j)%level = tree_panels_temp(i)%level + 1
                tree_panels_temp(panel_count+j)%face = tree_panels_temp(i)%face
                tree_panels_temp(panel_count+j)%id = panel_count+j
                allocate(tree_panels_temp(index)%points_inside(count))
            END DO
            ! coordinates of four new panels in angle space
            tree_panels_temp(panel_count+1)%min_xi = min_xi; tree_panels_temp(panel_count+1)%min_eta = min_eta
            tree_panels_temp(panel_count+1)%max_xi = mid_xi; tree_panels_temp(panel_count+1)%max_eta = mid_eta
            tree_panels_temp(panel_count+1)%mid_xi = 0.5*(min_xi+mid_xi)
            tree_panels_temp(panel_count+1)%mid_eta = 0.5*(min_eta+mid_eta)

            tree_panels_temp(panel_count+2)%min_xi = mid_xi; tree_panels_temp(panel_count+2)%min_eta = min_eta
            tree_panels_temp(panel_count+2)%max_xi = max_xi; tree_panels_temp(panel_count+2)%max_eta = mid_eta
            tree_panels_temp(panel_count+2)%mid_xi = 0.5*(mid_xi+max_xi)
            tree_panels_temp(panel_count+2)%mid_eta = 0.5*(min_eta+mid_eta)

            tree_panels_temp(panel_count+3)%min_xi = mid_xi; tree_panels_temp(panel_count+3)%min_eta = mid_eta
            tree_panels_temp(panel_count+3)%max_xi = max_xi; tree_panels_temp(panel_count+3)%max_eta = max_eta
            tree_panels_temp(panel_count+3)%mid_xi = 0.5*(mid_xi+max_xi)
            tree_panels_temp(panel_count+3)%mid_eta = 0.5*(mid_eta+max_eta)

            tree_panels_temp(panel_count+4)%min_xi = min_xi; tree_panels_temp(panel_count+4)%min_eta = mid_eta
            tree_panels_temp(panel_count+4)%max_xi = mid_xi; tree_panels_temp(panel_count+4)%max_eta = max_eta
            tree_panels_temp(panel_count+4)%mid_xi = 0.5*(min_xi+mid_xi)
            tree_panels_temp(panel_count+4)%mid_eta = 0.5*(mid_eta+max_eta)

            tree_panels_temp(i)%child_panel_1 = panel_count+1
            tree_panels_temp(i)%child_panel_2 = panel_count+2
            tree_panels_temp(i)%child_panel_3 = panel_count+3
            tree_panels_temp(i)%child_panel_4 = panel_count+4

            curr_loc(:) = 0
            DO j = 1, tree_panels_temp(i)%panel_point_count
                ! loop through points contained in parent panel, assign to sub panels
                index = tree_panels_temp(i)%points_inside(j)
                xval = xg(index); yval = yg(index); zval = zg(index)
                call xieta_from_xyz(xval, yval, zval, xi, eta, tree_panels_temp(i)%face)
                IF (xi < mid_xi) THEN
                    IF (eta < mid_eta) THEN
                        which_panel = 1
                    ELSE 
                        which_panel = 4
                    END IF
                ELSE 
                    IF (eta < mid_eta) THEN
                        which_panel = 2
                    ELSE 
                        which_panel = 3
                    END IF
                END IF
                curr_loc(which_panel) = curr_loc(which_panel) + 1
                tree_panels_temp(panel_count+which_panel)%panel_point_count = & 
                                            tree_panels_temp(panel_count+which_panel)%panel_point_count + 1
                tree_panels_temp(panel_count+which_panel)%points_inside(curr_loc(which_panel)) = index
            END DO

            ! resize points_inside array
            DO j = 1, 4
                allocate(temp(curr_loc(j)))
                DO k = 1, tree_panels_temp(panel_count+j)%panel_point_count
                    temp(k) = tree_panels_temp(panel_count+j)%points_inside(k)
                END DO
                call move_alloc(from=temp, to=tree_panels_temp(panel_count+j)%points_inside)
                curr_loc(j) = 0
            END DO

            ! sanity check to make sure all points are assigned
            total = tree_panels_temp(panel_count+1)%panel_point_count
            total = total + tree_panels_temp(panel_count+2)%panel_point_count
            total = total + tree_panels_temp(panel_count+3)%panel_point_count
            total = total + tree_panels_temp(panel_count+4)%panel_point_count
            IF (total /= tree_panels_temp(i)%panel_point_count) THEN
                print *, 'Error in refining tree at parent panel ', i
            END IF
            panel_count = panel_count + 4
        END IF
        i = i + 1
    enddo

    ! tree_panels_temp is the wrong size so move everything over to the correctly sized tree_panels
    allocate(tree_panels(panel_count))
    DO i = 1, panel_count
        tree_panels(i) = tree_panels_temp(i)
        tree_panels(i)%relabeled_points_inside = tree_panels(i)%points_inside
        min_xi = tree_panels(i)%min_xi
        mid_xi = tree_panels(i)%mid_xi
        max_xi = tree_panels(i)%max_xi
        min_eta = tree_panels(i)%min_eta
        mid_eta = tree_panels(i)%mid_eta
        max_eta = tree_panels(i)%max_eta
        ! compute the furthest distance from panel center to vertex for each panel
        call xyz_from_xieta(x1, x2, x3, mid_xi, mid_eta, tree_panels(i)%face)
        call xyz_from_xieta(y1, y2, y3, min_xi, min_eta, tree_panels(i)%face)
        d1 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
        call xyz_from_xieta(y1, y2, y3, min_xi, max_eta, tree_panels(i)%face)
        d2 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
        call xyz_from_xieta(y1, y2, y3, max_xi, max_eta, tree_panels(i)%face)
        d3 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
        call xyz_from_xieta(y1, y2, y3, max_xi, min_eta, tree_panels(i)%face)
        d4 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
        tree_panels(i)%radius = MAX(d1, d2, d3, d4)
    END DO
end subroutine tree_traversal

subroutine assign_points_to_panels(G, tree_panels, x, y, z, points_panels, levs)
    ! finds the panels containing the points in the computational domain, for the purposes of later computing proxy source potentials
    type(ocean_grid_type), intent(in) :: G
    type(cube_panel), intent(in) :: tree_panels(:)
    real, intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: levs
    integer, intent(inout) :: points_panels(:,:)
    integer :: level, i, j, k, isc, iec, jsc, jec, ic, jc
    real :: xco, yco, zco

    isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec
    ic = iec-isc+1; jc = jec-jsc+1

    DO i = 1, size(x)
        level = 1
        j = 1
        jloop: DO ! do loop over tree panels
            IF (j == -1) THEN
                exit jloop
            ELSE IF (tree_panels(j)%contains_point(x(i), y(i), z(i))) THEN
                ! point i is contained in panel j
                points_panels(level, i) = j
                level = level + 1
                j = tree_panels(j)%child_panel_1
            ELSE
                j = j + 1
            END IF
        END DO jloop
    END DO
end subroutine assign_points_to_panels

subroutine interaction_list_compute(G, pp_interactions, pc_interactions, tree_panels, x, y, z, theta, point_count)
    type(ocean_grid_type), intent(in) :: G
    type(interaction_pair), intent(out), allocatable :: pp_interactions(:), pc_interactions(:)
    type(cube_panel), intent(in) :: tree_panels(:)
    real, intent(in) :: x(:), y(:), z(:), theta
    integer, intent(in) :: point_count
    integer, allocatable :: source_index(:)
    type(interaction_pair), allocatable :: interaction_lists_temp(:)
    integer :: interaction_count, pp_count, pc_count, i, j, tt_count, k, curr_loc, i_s, c_s, l
    integer :: isc, iec, jsc, jec, ic, jc, index, tree_traverse_count, id
    real :: xco, yco, zco, xs, ys, zs, dist, separation, x1s, x2s, x3s, x1t, x2t, x3t
    logical :: well_separated

    isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec
    allocate(source_index(size(tree_panels)))
    allocate(interaction_lists_temp(100*point_count))

    interaction_count = 0
    pp_count = 0
    pc_count = 0

    do i = 1, point_count
        tt_count = 6
        do k = 1, 6
            source_index(k) = k
        enddo
        curr_loc = 1
        xco = x(i)
        yco = y(i)
        zco = z(i)
        panelloop: do while (curr_loc <= tt_count) ! go through source panels
            i_s = source_index(curr_loc)
            c_s = tree_panels(i_s)%panel_point_count
            if (c_s > 0) then ! cluster has points
                call xyz_from_xieta(xs, ys, zs, tree_panels(i_s)%mid_xi, tree_panels(i_s)%mid_eta, tree_panels(i_s)%face)
                dist = ACOS(MIN(MAX(xco*xs+yco*ys+zco*zs, -1.0), 1.0))
                separation = tree_panels(i_s)%radius/dist
                if ((dist > 0) .and. (separation < theta)) then
                    interaction_count = interaction_count + 1
                    pc_count = pc_count + 1
                    interaction_lists_temp(interaction_count)%index_target = i
                    interaction_lists_temp(interaction_count)%index_source = i_s
                    interaction_lists_temp(interaction_count)%interact_type = 1
                else  ! not well separated
                    if (tree_panels(i_s)%is_leaf) then ! leaf, do pp interaction
                        interaction_count = interaction_count + 1
                        pp_count = pp_count + 1
                        interaction_lists_temp(interaction_count)%index_target = i
                        interaction_lists_temp(interaction_count)%index_source = i_s
                        interaction_lists_temp(interaction_count)%interact_type = 0
                    else ! refine source panel
                        source_index(tt_count+1) = tree_panels(i_s)%child_panel_1
                        source_index(tt_count+2) = tree_panels(i_s)%child_panel_2
                        source_index(tt_count+3) = tree_panels(i_s)%child_panel_3
                        source_index(tt_count+4) = tree_panels(i_s)%child_panel_4
                        tt_count = tt_count + 4
                    end if
                end if
            end if
            curr_loc = curr_loc + 1
        enddo panelloop
    enddo

    allocate(pc_interactions(pc_count))
    allocate(pp_interactions(pp_count))

    k = 0
    l = 0
    do i = 1, interaction_count
        if (interaction_lists_temp(i)%interact_type == 1) then
            k = k + 1
            pc_interactions(k) = interaction_lists_temp(i)
        else
            l = l + 1
            pp_interactions(l) = interaction_lists_temp(i)
        end if
    enddo
end subroutine interaction_list_compute

subroutine calculate_communications(sal_ct, xg, yg, zg, G)
    type(sal_conv_type), intent(inout) :: sal_ct
    real, intent(in) :: xg(:), yg(:), zg(:)
    type(ocean_grid_type), intent(in) :: G
    integer :: p, id, unowned_source_count, pp_count, i, i_s, j, i_sp, j_sp, k, max_p, pl, count
    integer :: isg, ieg, jsg, jeg, isc, iec, jsc, jec, i_off, j_off, i_sp2, j_sp2, isdg, iedg, jsdg, jedg
    integer :: idgo, jdgo, isd, jsd, own_points
    integer, allocatable :: proc_start_i(:), proc_start_j(:), proc_end_i(:), proc_end_j(:)
    integer, allocatable :: unowned_temp_i(:), unowned_temp_j(:), unowned_sources_i(:), unowned_sources_j(:)
    integer, allocatable :: points_needed_from_proc(:), points_from_proc_i(:,:), proc_loc(:), temp_locs(:)
    integer, allocatable :: points_from_proc_j(:,:), points_needed_from_procs(:,:), points_to_give_proc(:)
    integer, allocatable :: points_to_give_proc_i(:,:), points_to_give_proc_j(:,:), unowned_temp(:), unowned_sources(:)
    integer, allocatable :: points_from_proc(:,:), points_to_give_proc_index(:,:)
    integer, allocatable :: point_counts_to_communicate(:), points_from_proc_2d(:,:,:), points_to_give_proc_2d(:,:,:) 
    real, allocatable :: e_xs(:), e_ys(:), e_zs(:)
    logical :: found

    p = sal_ct%p
    id = sal_ct%id

    isc = sal_ct%indexsg(id+1); iec = sal_ct%indexeg(id+1)
    own_points = sal_ct%pcg(id+1)

    ! find unowned sources
    allocate(unowned_temp(size(xg)), source=-1)
    allocate(points_needed_from_proc(p), source=0)
    allocate(proc_loc(size(xg)), source=-1)
    pp_count = size(sal_ct%pp_interactions)
    unowned_source_count = 0
    do i = 1, pp_count
        i_s = sal_ct%pp_interactions(i)%index_source
        do j = 1, sal_ct%tree_struct(i_s)%panel_point_count
            ! loop over points in source panel, check if owned
            i_sp = sal_ct%tree_struct(i_s)%points_inside(j)
            ! check if point is in the data domain
            if ((i_sp < isc) .or. (i_sp > iec)) then ! outside of computational domain
                found = .false.
                kloop: do k = 1, unowned_source_count
                    if (unowned_temp(k) == i_sp) then
                        found = .true.
                    end if
                enddo kloop
                if (.not. found) then
                    ! new unowned point
                    unowned_source_count = unowned_source_count + 1
                    unowned_temp(unowned_source_count) = i_sp
                    ! find which processor owns i_sp
                    kloop2: do k = 1, p
                        if ((i_sp >= sal_ct%indexsg(k)) .and. (i_sp <= sal_ct%indexeg(k))) then
                            points_needed_from_proc(k) = points_needed_from_proc(k) + 1
                            proc_loc(unowned_source_count) = k
                        end if
                    enddo kloop2
                end if
            end if
        enddo
    enddo

    ! print *, id, " has ", unowned_source_count, " unowned sources"

    allocate(unowned_sources(unowned_source_count))
    allocate(e_xs(unowned_source_count))
    allocate(e_ys(unowned_source_count))
    allocate(e_zs(unowned_source_count))

    sal_ct%unowned_sources = unowned_source_count

    max_p = 0
    do i = 1, p
        max_p = max(max_p, points_needed_from_proc(i))
    enddo

    allocate(points_from_proc(max_p, p), source=-1)
    allocate(points_from_proc_2d(1, max_p, p), source=-1)
    allocate(temp_locs(p), source=0)

    ! points needed from each proc in global indices
    do i = 1, unowned_source_count
        pl = proc_loc(i)
        temp_locs(pl) = temp_locs(pl) + 1
        points_from_proc(temp_locs(pl), pl) = unowned_temp(i)
        points_from_proc_2d(1, temp_locs(pl), pl) = unowned_temp(i)
    enddo

    count = 0
    do i = 1, p ! rearrange the unowned source points so continuous by owner
        do j = 1, max_p
            i_sp = points_from_proc(j, i)
            if (i_sp .ne. -1) then
                count = count + 1
                unowned_sources(count) = i_sp
                e_xs(count) = xg(i_sp)
                e_ys(count) = yg(i_sp)
                e_zs(count) = zg(i_sp)
            end if
        enddo
    enddo

    sal_ct%e_xs = e_xs
    sal_ct%e_ys = e_ys
    sal_ct%e_zs = e_zs
    allocate(points_to_give_proc(p), source=0)

    do i = 1,p
        call send_to_PE(points_needed_from_proc(i), i-1)
        ! print *, id, ' needs ', points_needed_from_proc(i), ' from proc ', i-1
    enddo

    do i = 1,p
        call recv_from_PE(points_to_give_proc(i), i-1)
        ! print *, id, ' sends ', points_to_give_proc(i), ' to proc ', i-1
    enddo

    call sync_PEs()

    max_p = 0
    do i = 1, p
        max_p = max(max_p, points_to_give_proc(i))
    enddo

    allocate(points_to_give_proc_index(max_p, p), source=-1)
    allocate(points_to_give_proc_2d(1, max_p, p), source=-1)

    ! do i=0, p-1
    !     call sync_PEs()
    !     if (id == i) then
    !         print *, id, ' needs points: '
    !         do j = 1, points_needed_from_proc(p-id)
    !             print *, points_from_proc_2d(1,j,p-id)
    !         enddo
    !         print *, '-------------------------'
    !     endif
    ! enddo

    do i=0, p-1 ! send point indices
        if (points_needed_from_proc(i+1)>0) then
            call send_to_PE(points_from_proc_2d(:,:,i+1), points_needed_from_proc(i+1), i)
        endif
    enddo

    do i=0, p-1 ! receive point indices
        if (points_to_give_proc(i+1)>0) then
            call recv_from_PE(points_to_give_proc_2d(:,:,i+1), points_to_give_proc(i+1), i)
        endif
    enddo

    ! do i=0, p-1
    !     call sync_PEs()
    !     if (id == i) then
    !         print *, id, ' sends points: '
    !         do j = 1, points_to_give_proc(p-id)
    !             print *, points_to_give_proc_2d(1,j,p-id)
    !         enddo
    !         print *, '-------------------------'
    !     endif
    ! enddo

    call sync_PEs()
    do i = 1,p
        if (points_to_give_proc(i)>0) then
            do j = 1, points_to_give_proc(i)
                if (points_to_give_proc_2d(1,j,i)>=0) then
                    points_to_give_proc_index(j,i) = points_to_give_proc_2d(1,j,i)-isc+1 ! shift to local indices 
                endif
            enddo
        endif
    enddo

    ! do i=0, p-1
    !     call sync_PEs()
    !     if (id == i) then
    !         print *, id, ' sends points: '
    !         do j = 1, points_to_give_proc(p-id)
    !             print *, points_to_give_proc_index(j,p-id)
    !         enddo
    !         print *, '-------------------------'
    !     endif
    ! enddo

    sal_ct%points_to_give = points_to_give_proc_index
    sal_ct%points_to_give_proc = points_to_give_proc
    sal_ct%points_to_get = points_from_proc
    sal_ct%points_to_get_proc = points_needed_from_proc

    ! relabel sources in tree for locally owned points
    do i = 1, size(sal_ct%tree_struct)
        do j = 1, sal_ct%tree_struct(i)%panel_point_count
            i_sp = sal_ct%tree_struct(i)%points_inside(j)
            ! contained in data domain
            if ((i_sp >= isc) .and. (i_sp <= iec)) then
                ! owned point, relabel
                sal_ct%tree_struct(i)%relabeled_points_inside(j) = i_sp - isc + 1
            end if
        enddo
    enddo

    ! relabel sources in tree for unowned points needed for pp interactions
    count = 0
    do i = 1, unowned_source_count
        i_sp = unowned_sources(i)
        j = 1
        ! loop over tree panels, relabel points i_sp, j_sp => i+own_points
        treeloop: do
            if (j == -1) then
                exit treeloop
            else 
                found = .false.
                kloop3: do k = 1, sal_ct%tree_struct(j)%panel_point_count
                    i_sp2 = sal_ct%tree_struct(j)%points_inside(k)
                    if (i_sp == i_sp2) then
                        found = .true.
                        sal_ct%tree_struct(j)%relabeled_points_inside(k) = i+own_points
                        exit kloop3
                    end if
                enddo kloop3
                if (found) then
                    j = sal_ct%tree_struct(j)%child_panel_1
                else
                    j = j + 1
                end if
            end if
        enddo treeloop
    enddo
end subroutine calculate_communications

subroutine sal_conv_init(sal_ct, G)
    ! initialize all the sal convolution things
    ! call only once
    ! does tree traversal, interaction list computation, sets up communication patterns
    type(SAL_conv_type), intent(out) :: sal_ct
    type(ocean_grid_type), intent(inout) :: G ! ocean grid
    character(len=12) :: mdl = "MOM_sal_conv" ! This module's name.
    integer :: proc_count, isc, iec, jsc, jec, isg, ieg, jsg, jeg, imax, jmax, ic, jc, i, j, ig_off, jg_off
    integer :: max_level, proc_rank, i_off, j_off, isd, ied, jsd, jed, rank, itc, jtc, count, pointcount
    integer, allocatable :: points_panels(:,:,:), indexsg(:), indexeg(:), pcg(:)
    real, allocatable :: xg(:,:), yg(:,:), zg(:,:), xc(:,:), yc(:,:), zc(:,:), xt(:,:), yt(:,:), zt(:,:)
    real, allocatable :: xg1d(:), yg1d(:), zg1d(:), xc1d(:), yc1d(:), zc1d(:), xt1d(:), yt1d(:), zt1d(:)
    real :: lat, lon, colat, x, y, z, pi

    ! fix this
    ! call get_param(param_file, mdl, "SHT_REPRODUCING_SUM", sal_ct%reprod_sum, &
    !                "If true, use reproducing sums (invariant to PE layout) in inverse transform "// &
    !                "of proxy source weights. Otherwise use a simple sum of floating point numbers. ", &
    !                default=.False.)

    sal_ct%p = num_PEs() ! number of ranks
    sal_ct%id = PE_here() ! current rank

    pi = 4.D0*DATAN(1.D0)

    isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec
    isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
    ! ig_off = G%idg_offset+isg-isd
    ! jg_off = G%jdg_offset+jsg-jsd
    ! isg = isc+ig_off; ieg = iec+ig_off; jsg = jsc+jg_off; jeg = jec+jg_off
    ! call get_global_grid_size(G, imax, jmax) ! total size in i/k directions

    ic=iec-isc+1; jc=jec-jsc+1

    count = 0
    allocate(sal_ct%two_d_to_1d(ied, jed), source=-1)
    do j = jsc, jec ! count number of ocean points in computational domain
        do i = isc, iec
            if (G%mask2dT(i, j) > 0.1) then
                count = count + 1
                sal_ct%two_d_to_1d(i, j) = count
            endif
        enddo
    enddo
    allocate(sal_ct%one_d_to_2d_i(count))
    allocate(sal_ct%one_d_to_2d_j(count))
    allocate(xc1d(count), source=0.0)
    allocate(yc1d(count), source=0.0)
    allocate(zc1d(count), source=0.0)
    sal_ct%own_ocean_points = count

    count = 0
    do j = jsc, jec
        do i = isc, iec
            if (G%mask2dT(i, j) > 0.1) then
                count = count + 1
                lat = G%geoLatT(i, j) * pi/180.0
                lon = G%geoLonT(i, j) * pi/180.0
                colat = 0.5*pi-lat
                x = sin(colat)*cos(lon)
                y = sin(colat)*sin(lon)
                z = cos(colat)

                xc1d(count) = x
                yc1d(count) = y
                zc1d(count) = z
                sal_ct%one_d_to_2d_i(count) = i
                sal_ct%one_d_to_2d_j(count) = j
            end if
        enddo
    enddo

    allocate(pcg(sal_ct%p), source=0)
    allocate(indexsg(sal_ct%p), source=0)
    allocate(indexeg(sal_ct%p), source=0)
    pcg(sal_ct%id+1) = count
    call sum_across_PEs(pcg, sal_ct%p)

    indexsg(1) = 1
    do i = 1, sal_ct%p-1
        indexeg(i)=indexsg(i)-1+pcg(i)
        indexsg(i+1)=indexeg(i)+1
    enddo
    indexeg(sal_ct%p) = indexsg(sal_ct%p)+pcg(sal_ct%p)-1

    sal_ct%indexsg = indexsg
    sal_ct%indexeg = indexeg
    sal_ct%pcg = pcg

    pointcount = 0
    do i = 1, sal_ct%p
        pointcount = pointcount + pcg(i)
    enddo
    sal_ct%total_ocean_points = pointcount

    allocate(xg1d(pointcount), source=0.0)
    allocate(yg1d(pointcount), source=0.0)
    allocate(zg1d(pointcount), source=0.0)
    xg1d(indexsg(sal_ct%id+1):indexeg(sal_ct%id+1)) = xc1d(:)
    yg1d(indexsg(sal_ct%id+1):indexeg(sal_ct%id+1)) = yc1d(:)
    zg1d(indexsg(sal_ct%id+1):indexeg(sal_ct%id+1)) = zc1d(:)
    call sum_across_PEs(xg1d, pointcount)
    call sum_across_PEs(yg1d, pointcount)
    call sum_across_PEs(zg1d, pointcount)
    ! xg/yg/zg is now a copy of all the points from all the processors
    call tree_traversal(G, sal_ct%tree_struct, xg1d, yg1d, zg1d, 30, pointcount) ! constructs cubed sphere tree
    max_level = sal_ct%tree_struct(size(sal_ct%tree_struct))%level

    allocate(sal_ct%points_panels(max_level+1, ic*jc), source=-1)
    ! finds which panels contain the computational domain points
    call assign_points_to_panels(G, sal_ct%tree_struct, xc1d, yc1d, zc1d, sal_ct%points_panels, max_level) 

    ! compute the interaction lists for the target points in the target domain
    call interaction_list_compute(G, sal_ct%pp_interactions, sal_ct%pc_interactions, sal_ct%tree_struct, xc1d, yc1d, zc1d, &
                                    0.7, sal_ct%own_ocean_points)
    call sync_PEs()

    ! compute communication patterns 
    call calculate_communications(sal_ct, xg1d, yg1d, zg1d, G)

    id_clock_SAL = cpu_clock_id('(Ocean SAL)', grain=CLOCK_MODULE)
    id_clock_SAL_pc = cpu_clock_id('(Ocean SAL PC interactions)', grain=CLOCK_MODULE)
    id_clock_SAL_pp = cpu_clock_id('(Ocean SAL PP interactions)', grain=CLOCK_MODULE)
    sal_ct%interp_degree=3
end subroutine sal_conv_init

subroutine ssh_pp_communications(sal_ct, G, eta, e_ssh)
    ! does the necessary communication of sshs to perform PP interactions
    type(SAL_conv_type), intent(in) :: sal_ct
    type(ocean_grid_type), intent(in) :: G
    real, intent(in) :: eta(:,:)
    real, intent(inout) :: e_ssh(:)
    integer :: max_give, max_get, p, id, i, j, i_s, j_s, i_off, j_off, count, index
    real, allocatable :: points_to_give(:,:,:), points_received(:,:,:)
    real :: area, rad2

    p = sal_ct%p; id = sal_ct%id
    i_off = G%isg - G%isc; j_off = G%jsg - G%jsc
    rad2 = G%Rad_Earth ** 2

    max_give = 0; max_get = 0
    do i = 1, p
        max_give = max(max_give, sal_ct%points_to_give_proc(i))
        max_get = max(max_get, sal_ct%points_to_get_proc(i))
    enddo

    allocate(points_to_give(1, max_give, p), source=0.0)
    allocate(points_received(1, max_get, p), source=0.0)

    do i = 1, p
        do j = 1, max_give
            index = sal_ct%points_to_give(j, i)
            if (index .ne. -1) then
                i_s = sal_ct%one_d_to_2d_i(index)
                j_s = sal_ct%one_d_to_2d_j(index)
                area = G%areaT(i_s, j_s)/rad2
                points_to_give(1, j, i) = eta(i_s, j_s)*area
            end if
        enddo
    enddo

    do i=0, p-1 ! send points
        if (sal_ct%points_to_give_proc(i+1)>0) then
            call send_to_PE(points_to_give(:,:,i+1), sal_ct%points_to_give_proc(i+1), i)
        endif
    enddo

    do i=0, p-1 ! receive points
        if (sal_ct%points_to_get_proc(i+1)>0) then
            call recv_from_PE(points_received(:,:,i+1), sal_ct%points_to_get_proc(i+1), i)
        endif
    enddo

    call sync_PEs()

    count = 0
    do i = 1, p
        do j = 1, sal_ct%points_to_get_proc(i)
            count = count + 1
            e_ssh(count) = points_received(1, j, i)
        enddo
    enddo
end subroutine ssh_pp_communications

subroutine bli_coeffs(coeffs, degree)
    real, allocatable, intent(out) :: coeffs(:)
    integer, intent(in) :: degree
    integer :: i

    allocate(coeffs(degree+1))

    if (degree == 0) then
        coeffs(1) = 1.0
    else
        do i = 1, degree+1
            if (i == 1) then
                coeffs(i) = 0.5
            else if (i == degree + 1) then
                coeffs(i) = 0.5 * ((-1) ** (i+1))
            else
                coeffs(i) = (-1.0) ** (i+1)
            end if
        enddo
    end if
end subroutine bli_coeffs

subroutine bli_interp_points_shift(interp_points, min_x, max_x, degree)
    real, allocatable, intent(out) :: interp_points(:)
    real, intent(in) :: min_x, max_x
    integer, intent(in) :: degree
    real :: x_range, x_shift, pi
    integer :: i

    pi = 4.0*DATAN(1.0)

    allocate(interp_points(degree+1))
    x_range = 0.5*(max_x - min_x)
    x_shift = 0.5*(max_x + min_x)

    if (degree == 0) then
        interp_points(1) = x_shift
    else 
        do i = 1, degree + 1
            interp_points(i) = COS(pi/degree*(i-1))*x_range + x_shift
        enddo
    end if
end subroutine bli_interp_points_shift

subroutine interp_vals_bli(vals, xi, eta, min_xi, max_xi, min_eta, max_eta, degree)
    real, allocatable, intent(out) :: vals(:,:)
    real, intent(in) :: xi, eta, min_xi, max_xi, min_eta, max_eta
    integer, intent(in) :: degree
    real, allocatable :: bli_xi_vals(:), bli_eta_vals(:), bli_coeff_vals(:), xi_func_vals(:), eta_func_vals(:)
    logical :: found_xi, found_eta
    integer :: i, j
    real :: denom_xi, denom_eta, val

    allocate(vals(degree+1, degree+1))
    allocate(xi_func_vals(degree+1), source=0.0_8)
    allocate(eta_func_vals(degree+1), source=0.0_8)

    call bli_interp_points_shift(bli_xi_vals, min_xi, max_xi, degree)
    call bli_interp_points_shift(bli_eta_vals, min_eta, max_eta, degree)
    call bli_coeffs(bli_coeff_vals, degree)

    found_xi = .false.
    found_eta = .false.

    ! first check if xi/eta are an interpolation point
    do i = 1, degree+1
        if (ABS(xi - bli_xi_vals(i)) < 1.0e-16) then
            found_xi = .true.
            xi_func_vals(i) = 1.0
        end if
        if (ABS(eta - bli_eta_vals(i)) < 1.0e-16) then
            found_eta = .true.
            eta_func_vals(i) = 1.0
        end if
    enddo

    ! xi/eta are not interpolation point, compute all the BLI basis values
    if (.not. found_xi) then
        denom_xi = 0.0
        do i = 1, degree+1
            val = bli_coeff_vals(i) / (xi - bli_xi_vals(i))
            xi_func_vals(i) = val
            denom_xi = denom_xi + val
        enddo
        do i = 1, degree + 1
            xi_func_vals(i) = xi_func_vals(i) / denom_xi
        enddo
    end if

    if (.not. found_eta) then
        denom_eta = 0.0
        do i = 1, degree+1
            val = bli_coeff_vals(i) / (eta - bli_eta_vals(i))
            eta_func_vals(i) = val
            denom_eta = denom_eta + val
        enddo
        do i = 1, degree + 1
            eta_func_vals(i) = eta_func_vals(i) / denom_eta
        enddo
    end if

    ! compute the outer product
    do j = 1, degree + 1
        do i = 1, degree + 1
            vals(i, j) = xi_func_vals(i) * eta_func_vals(j)
        enddo
    enddo
end subroutine interp_vals_bli

subroutine proxy_source_compute(sal_ct, G, sshs, proxy_source_weights) ! add reproducing sum
    type(SAL_conv_type), intent(in) :: sal_ct
    real, intent(in) :: sshs(:,:)
    real, intent(inout) :: proxy_source_weights(:)
    type(ocean_grid_type), intent(in) :: G
    integer :: isc, iec, jsc, jec, i, j, max_levels, k, i_t, shift, offset, l, m, ic, jc, ii, ij
    real :: x, y, z, a, ssh, lat, lon, colat, pi, r2, min_xi, max_xi, min_eta, max_eta, xi, eta
    real, allocatable:: basis_vals(:,:), proxy_source_weights_sep(:,:,:)
    integer, allocatable :: points_in_panel(:), pos_in_array(:)
    real :: sum_tot

    pi = 4.0*DATAN(1.0)
    isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec
    ic = iec-isc+1; jc = jec-jsc+1
    r2 = G%Rad_Earth ** 2

    do i = 1, sal_ct%own_ocean_points
        ii = sal_ct%one_d_to_2d_i(i)
        ij = sal_ct%one_d_to_2d_j(i)
        lat = G%geoLatT(ii, ij)*pi/180.0
        lon = G%geoLonT(ii, ij)*pi/180.0
        colat = 0.5*pi-lat
        x = sin(colat)*cos(lon)
        y = sin(colat)*sin(lon)
        z = cos(colat)
        a = G%areaT(ii, ij)/r2
        ssh = sshs(ii, ij)
        panelloop2: do k = 1, size(sal_ct%points_panels(:,i))
            i_t = sal_ct%points_panels(k, i)
            if (i_t == -1) then
                exit panelloop2
            else 
                shift = (i_t-1)*(sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
                min_xi = sal_ct%tree_struct(i_t)%min_xi
                max_xi = sal_ct%tree_struct(i_t)%max_xi
                min_eta = sal_ct%tree_struct(i_t)%min_eta
                max_eta = sal_ct%tree_struct(i_t)%max_eta
                call xieta_from_xyz(x, y, z, xi, eta, sal_ct%tree_struct(i_t)%face)
                call interp_vals_bli(basis_vals, xi, eta, min_xi, max_xi, min_eta, max_eta, sal_ct%interp_degree)
                offset = 0
                do l = 1, sal_ct%interp_degree+1
                    do m = 1, sal_ct%interp_degree+1
                        offset=offset+1
                        proxy_source_weights(shift+offset)=proxy_source_weights(shift+offset)+basis_vals(m, l)*ssh*a
                    enddo
                enddo
            end if
        enddo panelloop2
    enddo

    call sum_across_PEs(proxy_source_weights, size(proxy_source_weights))
end subroutine proxy_source_compute

subroutine sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal, sal_x, sal_y) ! explore impact of eps and cons more
    real, intent(in) :: tx, ty, tz, sx, sy, sz
    real, intent(out) :: sal_x, sal_y, sal
    real :: g, mp, sqrtp, cons, sqp, p1, p2, x32, val, mp2, cons1, eps

    cons = -7.029770573725803e-9/3.0 ! modify this
    ! cons = 0.0
    cons1 = -0.0447867/40.0 ! modify this
    ! cons1 = 0.0
    eps=1e-4

    sal_x = 0.0
    sal_y = 0.0
    sal = 0.0
    IF ((abs(tz - 1.0) > 1e-15) .and. (abs(tz+1.0) > 1e-15)) THEN
        g = max(min(tx*sx+ty*sy+tz*sz, 1.0), -1.0) ! floating point check
        mp = 2.0-2.0*g
        sqp = sqrt(mp)
        p1 = (1.0-6.21196)/(sqp*mp+eps)
        p2 = (2.7+6.0)*(2*g+sqp) / (2.0*(g*g-1.0)+eps)
        val = (p1+p2)*cons ! modify this
        x32 = tz*tz
        mp2 = sqrt(1.0-x32)
        sal_y = (sz*(1.0-x32)-tz*(tx*sx+ty*sy))/mp2*val
        sal_x = (tx*sy-ty*sx)/mp2*val
        sal = cons1*((1.0-6.21196)/(sqp+eps)+(2.7+6.0)*log(2*sqp+mp+eps))
    END IF
end subroutine sal_grad_gfunc

subroutine pc_interaction_compute(sal_ct, G, proxy_source_weights, sal, sal_x, sal_y)
    type(SAL_conv_type), intent(in) :: sal_ct
    type(ocean_grid_type), intent(in) :: G
    real, intent(in) :: proxy_source_weights(:)
    real, intent(inout) :: sal_x(:,:), sal_y(:,:), sal(:,:)
    integer :: proxy_count, i, i_s, i_ti, i_tj, j, k, offset, i_t, isc, jsc, id
    real, allocatable :: source_proxy_weights(:), cheb_xi(:), cheb_eta(:)
    real :: min_xi, max_xi, min_eta, max_eta, x, y, z, lat, lon, xi, eta, colat, pi
    real :: cx, cy, cz, sal_grad_x, sal_grad_y, sal_val

    pi = 4.0*DATAN(1.0)
    isc = G%isc; jsc = G%jsc

    id = PE_here()

    proxy_count = (sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
    allocate(source_proxy_weights(proxy_count), source=0.0)
    do i = 1, size(sal_ct%pc_interactions)
        i_s = sal_ct%pc_interactions(i)%index_source
        i_t = sal_ct%pc_interactions(i)%index_target
        i_ti = sal_ct%one_d_to_2d_i(i_t)
        i_tj = sal_ct%one_d_to_2d_j(i_t)
        source_proxy_weights = proxy_source_weights(proxy_count*(i_s-1)+1:proxy_count*i_s)
        min_xi = sal_ct%tree_struct(i_s)%min_xi
        max_xi = sal_ct%tree_struct(i_s)%max_xi
        min_eta = sal_ct%tree_struct(i_s)%min_eta
        max_eta = sal_ct%tree_struct(i_s)%max_eta
        call bli_interp_points_shift(cheb_xi, min_xi, max_xi, sal_ct%interp_degree)
        call bli_interp_points_shift(cheb_eta, min_eta, max_eta, sal_ct%interp_degree)
        lat = G%geoLatT(i_ti, i_tj)*pi/180.0; lon = G%geoLonT(i_ti, i_tj)*pi/180.0
        colat = 0.5*pi-lat
        x = sin(colat)*cos(lon)
        y = sin(colat)*sin(lon)
        z = cos(colat)
        offset = 0
        do k = 1, sal_ct%interp_degree+1
            eta = cheb_eta(k)
            do j = 1, sal_ct%interp_degree+1
                xi = cheb_xi(j)
                call xyz_from_xieta(cx, cy, cz, xi, eta, sal_ct%tree_struct(i_s)%face)
                call sal_grad_gfunc(x, y, z, cx, cy, cz, sal_val, sal_grad_x, sal_grad_y)
                offset = offset+1
                sal(i_ti, i_tj) = sal(i_ti, i_tj) + sal_val*source_proxy_weights(offset)
                sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + sal_grad_x*source_proxy_weights(offset)
                sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + sal_grad_y*source_proxy_weights(offset)
            enddo
        enddo
    enddo
end subroutine pc_interaction_compute

subroutine pp_interaction_compute(sal_ct, G, eta, sal, sal_x, sal_y, e_ssh)
    type(SAL_conv_type), intent(in) :: sal_ct
    type(ocean_grid_type), intent(in) :: G
    real, intent(inout) :: sal_x(:,:), sal_y(:,:), sal(:,:)
    real, intent(in) :: e_ssh(:), eta(:,:)
    integer :: i, i_s, i_ti, i_tj, j, i_si, i_sj, i_sp, i_t, isc, jsc
    real :: pi, lat, lon, colat, x, y, z, sx, sy, sz, ssh, r2, sal_grad_x, sal_grad_y, sal_val
    integer :: id, count, ii, ij, iec, jec
    real, allocatable :: sshs(:), a_s(:)

    pi = 4.0*DATAN(1.0)
    r2 = G%Rad_Earth ** 2
    isc = G%isc; jsc = G%jsc; iec = G%iec; jec = G%jec

    do i = 1, size(sal_ct%pp_interactions)
        i_s = sal_ct%pp_interactions(i)%index_source
        i_t = sal_ct%pp_interactions(i)%index_target
        i_ti = sal_ct%one_d_to_2d_i(i_t)
        i_tj = sal_ct%one_d_to_2d_j(i_t)
        lat = G%geoLatT(i_ti, i_tj)*pi/180.0
        lon = G%geoLonT(i_ti, i_tj)*pi/180.0
        colat = 0.5*pi-lat
        x = sin(colat)*cos(lon)
        y = sin(colat)*sin(lon)
        z = cos(colat)
        do j = 1, sal_ct%tree_struct(i_s)%panel_point_count
            i_sp = sal_ct%tree_struct(i_s)%relabeled_points_inside(j)
            if (i_sp > sal_ct%own_ocean_points) then ! unowned source point
                i_si = i_sp - sal_ct%own_ocean_points
                ssh = e_ssh(i_si)
                sx = sal_ct%e_xs(i_si)
                sy = sal_ct%e_ys(i_si)
                sz = sal_ct%e_zs(i_si)
            else
                i_si = sal_ct%one_d_to_2d_i(i_sp)
                i_sj = sal_ct%one_d_to_2d_j(i_sp)
                lat = G%geoLatT(i_si, i_sj)
                lon = G%geoLonT(i_si, i_sj)
                colat = 0.5*pi-lat
                sx = sin(colat)*cos(lon)
                sy = sin(colat)*sin(lon)
                sz = cos(colat)
                ssh = eta(i_si, i_sj)*G%areaT(i_si, i_sj)/r2
            end if
            call sal_grad_gfunc(x, y, z, sx, sy, sz, sal_val, sal_grad_x, sal_grad_y)
            sal(i_ti, i_tj) = sal(i_ti, i_tj)+sal_val*ssh
            sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + sal_grad_x*ssh
            sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + sal_grad_y*ssh
        enddo
    enddo
end subroutine pp_interaction_compute

subroutine sal_conv_eval(sal_ct, G, eta, e_sal, sal_x, sal_y)
    type(SAL_conv_type), intent(in) :: sal_ct ! conv SAL data struct
    type(ocean_grid_type), intent(in) :: G ! ocean grid
    real, intent(in) :: eta(:,:) ! ssh
    real, intent(inout) :: e_sal(:,:), sal_x(:,:), sal_y(:,:) ! x,y components of SAL potential gradient
    real, allocatable :: e_ssh(:), proxy_source_weights(:)
    integer :: source_size

    call cpu_clock_begin(id_clock_SAL)

    sal_x(:,:) = 0.0
    sal_y(:,:) = 0.0
    e_sal(:,:) = 0.0

    call cpu_clock_begin(id_clock_SAL_pc)

    ! compute proxy source weights for computational domain
    source_size = (sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)*size(sal_ct%tree_struct)
    allocate(proxy_source_weights(source_size), source=0.0)
    call proxy_source_compute(sal_ct, G, eta, proxy_source_weights)

    ! compute PC interactions for target domain
    call pc_interaction_compute(sal_ct, G, proxy_source_weights, e_sal, sal_x, sal_y)

    call cpu_clock_end(id_clock_SAL_pc)
    call cpu_clock_begin(id_clock_SAL_pp)

    ! do SSH communication needed for PP interactions
    allocate(e_ssh(sal_ct%unowned_sources), source=0.0)
    call ssh_pp_communications(sal_ct, G, eta, e_ssh)

    ! compute PP interactions for target domain
    call pp_interaction_compute(sal_ct, G, eta, e_sal, sal_x, sal_y, e_ssh)

    call cpu_clock_end(id_clock_SAL_pp)

    call sync_PEs()

    call pass_var(sal_x, G%domain) ! halo update 
    call pass_var(sal_y, G%domain) ! halo update 

    call cpu_clock_end(id_clock_SAL)
end subroutine sal_conv_eval

! subroutine sal_conv_end(sal_ct)
!     type(SAL_conv_type), intent(inout) :: sal_ct
!     deallocate(sal_ct)
! end subroutine sal_conv_end

end module MOM_conv_self_attr_load
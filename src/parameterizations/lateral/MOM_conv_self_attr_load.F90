module MOM_conv_self_attr_load

use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_grid,            only : ocean_grid_type, get_global_grid_size
use MOM_file_parser,     only : read_param, get_param, log_version, param_file_type
use MOM_coms_infra,      only : sum_across_PEs, max_across_PEs, num_PEs, PE_here
use MOM_coms_infra,			 only : broadcast, sync_PEs
use MOM_coms_infra,      only : send_to_PE, recv_from_PE
use MOM_coms,            only : reproducing_sum
use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_domains,         only : pass_var

implicit none ; private

public sal_conv_init, sal_conv_eval, sal_conv_end

#include <MOM_memory.h>

!> data type to represent a panel of the cubed sphere in the tree structure
type, private :: cube_panel
    integer :: level = 0 
      !< level 0 is the base level
    logical :: is_leaf = .true. 
      !< all leaves start out as leaves, set to false if refined
    integer :: id 
      !< panel id
    integer :: parent_panel = -1 
      !< id of parent panel
    integer, allocatable :: child_panels(:) 
      !< child panel ids if present
    integer :: child_panel_count = 0 
      !< number of child panels
    integer :: face 
      !< which face of the cubed sphere the panel belongs to
    real :: min_xi 
      !< xi and eta are the angle coordinates on the face of the cube
    real :: mid_xi 
      !< based on the equiangular gnomonic cubed sphere
    real :: max_xi 
    real :: min_eta
    real :: mid_eta
    real :: max_eta
    real :: radius 
      !< distance from center of panel to corner in cubed sphere angle coordinates
    integer, allocatable :: points_inside(:) 
      !< indices of contained points
    integer, allocatable :: relabeled_points_inside(:) 
      !< indices after relabeling 
    integer :: panel_point_count = 0 
      !< number of source/target points inside this panel
    integer :: own_point_count = 0 
      !< number of owned source/target points inside this panel
end type cube_panel

!> data type to represent a pp/pc/cp/cc interaction between a point/panel and a cube panel
type, private :: interaction_pair 
    integer :: index_target
      !< id of the target panel
    integer :: index_source 
      !< id of the source panel
    integer :: interact_type 
      !< 0 for PP, 1 for PC, 2 for CP, 3 for CC
end type interaction_pair

!> The control structure for the MOM_conv_self_attr_load module
type, public :: SAL_Conv_CS ; 
! public
    logical :: reprod_sum 
      !< If true, use reproducing sums
    logical :: use_fmm 
      !< If true, run in FMM mode, incompatible with reprod_sum
    logical :: use_bpa
      !< If true, use bottom pressure anomaly instead of sea surface height to calculate SAL
    real, allocatable :: e_xs(:), e_ys(:), e_zs(:) 
      !< x/y/z coordinates of unowned points needed for particle-particle interactions
    type(interaction_pair), allocatable :: pp_interactions(:), pc_interactions(:) 
      !< particle-particle, particle-cluster, cluster-particle, cluster-cluster interaction lists
    type(interaction_pair), allocatable :: cp_interactions(:), cc_interactions(:)
    type(cube_panel), allocatable :: tree_struct(:)
      !< tree structure of source points/source panels
    type(cube_panel), allocatable :: tree_struct_targets(:) 
      !< tree structure for own target points/target panels, used for FMM
    integer, allocatable :: points_panels(:,:) 
      !< points_panels(lev, i)=k means that point i is contained in panel k at level lev
    integer :: p 
      !< total MPI ranks
    integer :: id 
      !< current rank
    integer, allocatable :: no_points_to_give_proc(:) 
      !< number of points the current MPI rank sends to every other MPI rank to evaluate 
      !! particle-particle and cluster-particle interactions
    integer, allocatable :: no_points_to_get_proc(:) 
      !< number of points the current MPI rank receives from every other MPI rank to evaluate
      !! particle-particle and cluster-particle interactions
    integer, allocatable :: points_to_give(:,:)
      !< point indices for SSHs that need to be sent to other MPI ranks
    integer, allocatable :: points_to_get(:,:) 
      !< point indices for SSHs that need to be received from other MPI ranks
    integer, allocatable :: unowned_source_points(:)
      !< point indices of source points that are not in the MPI rank's computational domain
    integer :: unowned_sources
      !< number of source points not in the MPI rank's computational domain
    integer :: interp_degree 
      !< interpolation degree used for the barycentric Lagrange interpolation in the CSFMM
    integer, allocatable :: indexsg(:)
      !< end index of each MPI rank's points
    integer, allocatable :: indexeg(:) 
      !< start index of each MPI rank's points
    integer, allocatable :: pcg(:)
      !< each MPI rank's computational domain point count
    integer, allocatable :: two_d_to_1d(:,:) 
      !< index (i,j) returns the corresponding 1d index of a point
    integer, allocatable :: one_d_to_2d_i(:)
      !< index returns the corresponding horizontal i index value of a point
    integer, allocatable :: one_d_to_2d_j(:)
      !< index returns the corresponding horizontal j index value of a point
    integer :: own_ocean_points 
      !< number of ocean points in a MPI rank's computational domain
    integer :: total_ocean_points
      !< total number of ocean points
    integer, allocatable :: point_leaf_panel(:) 
      !< leaf panel index that contains a specific point
    integer :: cluster_thresh
      !< cluster threshold for the CSFMM to use cluster approximations
    logical :: use_sal_conv = .false.
      !< If true, use the convolution/CSFMM SAL
    logical :: use_farfield = .true. 
      !< If false, excluse the particle-cluster, cluster-particle, and cluster-cluster interactions
    real :: scaling_constant 
      !< If using SSH, (3 rho_w)/(4 pi rho_e)
      !< If using BPA, 3 / (4 pi rho_e grav)
end type SAL_Conv_CS

integer :: id_clock_SAL   !< CPU clock for self-attraction and loading
integer :: id_clock_SAL_tc !< CPU clock when used in tree code mode
integer :: id_clock_SAL_fmm !< CPU clock when used in FMM mode
integer :: id_clock_SAL_upward_pass !< CPU clock for the upward pass computation
integer :: id_clock_SAL_downward_pass !< CPU clock for the downward pass computation
integer :: id_clock_SAL_ssh_comm !< CPU clock for MPI communication
integer :: id_clock_SAL_cc_comp !< CPU clock for cluster/cluster computation
integer :: id_clock_SAL_pc_comp !< CPU clock for particle/cluster computation
integer :: id_clock_SAL_cp_comp !< CPU clock for cluster/particle computation
integer :: id_clock_SAL_pp_comp !< CPU clock for particle/particle computation

contains

!> takes (x, y, z) coordinates of a point and returns which face of the cubed sphere it corresponds to
integer function face_from_xyz(x, y, z) 
  real, intent(in) :: x, y, z
  real :: ax, ay, az !< absolute values of x, y, and z
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
    endif
  else if ((ay >= ax) .and. (ay >= az)) then
    if (y >= 0) then
      face = 2
    else
      face = 4
    endif
  else
    if (z >= 0) then
      face = 5
    else
      face = 6
    endif
  endif
  face_from_xyz = face
  return
end function face_from_xyz

!> computes (xi, eta) angle coordinates on the cubed sphere from (x, y, z) points and optional face
subroutine xieta_from_xyz(x, y, z, xi, eta, in_face) !> based on Ronchi et al, 1996 cubed sphere
  real, intent(in) :: x, y, z
  integer, intent(in), optional :: in_face
  real, intent(out) :: xi, eta
  integer :: face
  if (.not. present(in_face)) then
    face = face_from_xyz(x, y, z)
  else 
    face = in_face
  endif
  if (face == 1) then
    xi = atan(y/x); eta = atan(z/x)
  else if (face == 2) then
    xi = atan(-x/y); eta = atan(z/y)
  else if (face == 3) then
    xi = atan(y/x); eta = atan(-z/x)
  else if (face == 4) then
    xi = atan(-x/y); eta = atan(-z/y)
  else if (face == 5) then
    xi = atan(y/z); eta = atan(-x/z)
  else if (face == 6) then
    xi = atan(-y/z); eta = atan(-x/z)
  endif
end subroutine xieta_from_xyz

!> converts (xi, eta) on the cubed sphere to (x, y, z) coordinates
subroutine xyz_from_xieta(x, y, z, xi, eta, face) !> based on Ronchi et al, 1996 cubed sphere
  real, intent(in) :: xi, eta
  integer, intent(in) :: face
  real, intent(out) :: x, y, z
  real :: ax, ay, isqrt
  ax = tan(xi); ay = tan(eta)
  isqrt = 1.0 / sqrt(1.0 + ax*ax + ay*ay)
  if (face == 1) then
    x = isqrt
    y = ax*x
    z = ay*x
  else if (face == 2) then
    y = isqrt
    x = -ax*y
    z = ay*y
  else if (face == 3) then
    x = -isqrt
    y = ax*x
    z = -ay*x
  else if (face == 4) then
    y = -isqrt
    x = -ax*y
    z = -ay*y
  else if (face == 5) then
    z = isqrt
    x = -ay*z
    y = ax*z
  else 
    z = -isqrt
    x = -ay*z
    y = -ax*z
  endif
end subroutine xyz_from_xieta

!> constructs the cubed sphere tree of points from point locations (xg, yg, zg)
subroutine tree_traversal(G, tree_panels, xg, yg, zg, cluster_thresh, point_count, base_panels)
  type(ocean_grid_type), intent(inout) :: G !< ocean grid
  type(cube_panel), allocatable, intent(out) :: tree_panels(:) !< array of tree panels 
      !! representing the tree
  type(cube_panel), allocatable :: tree_panels_temp(:) !< temp array of tree panels
      !! as the size of the tree is not initially known, allocated to be larger than needed
  integer, intent(in) :: cluster_thresh !< point count threshold for refining a panel
  integer, intent(in) :: point_count !< number of points in the global grid
  real, intent(in) :: xg(:), yg(:), zg(:) !< global (x,y,z) coordinate arrays
  integer, intent(out) :: base_panels !< number of top level panels
  real :: pi, xval, yval, zval, min_xi, mid_xi, max_xi, min_eta, mid_eta, max_eta, xi, eta
  integer, allocatable :: curr_point_assigned_count(:) 
      !< number of points currently assigned to a panel
  ! integer, allocatable :: temp(:)
  integer, allocatable :: point_panel(:) !< which panel a point belongs to
  integer, allocatable :: panel_points(:) !< number of points in a panel
  integer, allocatable :: face_id(:) !< id of each face, 1 through 6
  integer, allocatable :: panel_id(:) !< id of each panel
  integer :: face, i, panel_count, j, count, index, k
  integer :: which_panel, total, loc, kids
  real :: x1, x2, x3, y1, y2, y3, d1, d2, d3, d4
  real, allocatable :: point_xi(:), point_eta(:) !< (xi, eta) coordinates of each point

  pi = 4.D0*datan(1.D0)
  allocate(tree_panels_temp(max(6, point_count)))
  allocate(point_panel(point_count), source=0)
  allocate(panel_points(6), source=0)
  allocate(face_id(6), source=-1)
  allocate(point_xi(point_count))
  allocate(point_eta(point_count))

  ! find the top level panels of all the points
  do i = 1, point_count
    xval = xg(i)
    yval = yg(i)
    zval = zg(i)
    face = face_from_xyz(xval, yval, zval)
    point_panel(i) = face
    panel_points(face) = panel_points(face) + 1
    call xieta_from_xyz(xval, yval, zval, xi, eta, face)
    point_xi(i) = xi
    point_eta(i) = eta
  enddo

  ! initialize the six top level cube panels
  loc = 0
  do i = 1, 6
    if (panel_points(i) > 0) then
      loc = loc + 1
      tree_panels_temp(loc)%id = loc
      tree_panels_temp(loc)%face = i
      tree_panels_temp(loc)%min_xi = -pi/4.D0
      tree_panels_temp(loc)%mid_xi = 0.D0
      tree_panels_temp(loc)%max_xi = pi/4.D0
      tree_panels_temp(loc)%min_eta = -pi/4.D0
      tree_panels_temp(loc)%mid_eta = 0.D0
      tree_panels_temp(loc)%max_eta = pi/4.D0
      allocate(tree_panels_temp(loc)%points_inside(panel_points(i)))
      tree_panels_temp(loc)%panel_point_count = panel_points(i)
      face_id(i) = loc
    endif
  enddo

  panel_count = loc

  allocate(curr_point_assigned_count(6), source=0)

  do i = 1, point_count
    face = point_panel(i)
    loc = face_id(face)
    curr_point_assigned_count(loc) = curr_point_assigned_count(loc) + 1
    tree_panels_temp(loc)%points_inside(curr_point_assigned_count(loc)) = i
  enddo

  ! shrink panels
  do i = 1, panel_count
    min_xi = 3.0
    max_xi = -3.0
    min_eta = 3.0
    max_eta = -3.0
    do j = 1, tree_panels_temp(i)%panel_point_count
      index = tree_panels_temp(i)%points_inside(j)
      min_xi = min(min_xi, point_xi(index))
      max_xi = max(max_xi, point_xi(index))
      min_eta = min(min_eta, point_eta(index))
      max_eta = max(max_eta, point_eta(index))
    enddo
    tree_panels_temp(i)%min_xi = min_xi-1e-15
    tree_panels_temp(i)%max_xi = max_xi+1e-15
    tree_panels_temp(i)%mid_xi = 0.5*(min_xi + max_xi)
    tree_panels_temp(i)%min_eta = min_eta-1e-15
    tree_panels_temp(i)%max_eta = max_eta+1e-15
    tree_panels_temp(i)%mid_eta = 0.5*(min_eta + max_eta)
  enddo

  base_panels = panel_count

  i = 1
  deallocate(panel_points)
  allocate(panel_points(4), source=0)
  allocate(panel_id(4), source=-1)
  do while (i <= panel_count)
    ! first check if panel needs to be divided
    count = tree_panels_temp(i)%panel_point_count
    if ((count >= cluster_thresh) .and. (tree_panels_temp(i)%is_leaf)) then
      ! panel is a leaf and has many points => refine
      tree_panels_temp(i)%is_leaf = .false.
      min_xi = tree_panels_temp(i)%min_xi
      mid_xi = tree_panels_temp(i)%mid_xi
      max_xi = tree_panels_temp(i)%max_xi
      min_eta = tree_panels_temp(i)%min_eta
      mid_eta = tree_panels_temp(i)%mid_eta
      max_eta = tree_panels_temp(i)%max_eta
      panel_points(:) = 0
      kids = 0
      panel_id(:) = -1
      loc = 0

      do j = 1, tree_panels_temp(i)%panel_point_count
        index = tree_panels_temp(i)%points_inside(j)
        xval = xg(index); yval = yg(index); zval = zg(index)
        call xieta_from_xyz(xval, yval, zval, xi, eta, tree_panels_temp(i)%face)
        if (xi < mid_xi) then
          if (eta < mid_eta) then
            which_panel = 1
          else 
            which_panel = 4
          endif
        else 
          if (eta < mid_eta) then
            which_panel = 2
          else 
            which_panel = 3
          endif
        endif
        panel_points(which_panel) = panel_points(which_panel) + 1
        point_panel(j) = which_panel
      enddo

      do j = 1, 4
        if (panel_points(j) > 0) then
          kids = kids + 1
        endif
      enddo

      ! up to four new panels are index panel_count+1 through panel_count+4
      do j = 1, 4
        if (panel_points(j) > 0) then
          loc = loc + 1
          index = panel_count + loc
          panel_id(j) = loc
          count = panel_points(j)
          tree_panels_temp(index)%parent_panel = i
          tree_panels_temp(index)%level = tree_panels_temp(i)%level + 1
          tree_panels_temp(index)%face = tree_panels_temp(i)%face
          tree_panels_temp(index)%id = index
          allocate(tree_panels_temp(index)%points_inside(count))
          tree_panels_temp(index)%panel_point_count = count
          ! coordinates of four new panels in angle space
          if (j == 1) then
            tree_panels_temp(index)%min_xi = min_xi; tree_panels_temp(index)%min_eta = min_eta
            tree_panels_temp(index)%max_xi = mid_xi; tree_panels_temp(index)%max_eta = mid_eta
            tree_panels_temp(index)%mid_xi = 0.5*(min_xi+mid_xi)
            tree_panels_temp(index)%mid_eta = 0.5*(min_eta+mid_eta)
          endif
          if (j == 2) then
            tree_panels_temp(index)%min_xi = mid_xi; tree_panels_temp(index)%min_eta = min_eta
            tree_panels_temp(index)%max_xi = max_xi; tree_panels_temp(index)%max_eta = mid_eta
            tree_panels_temp(index)%mid_xi = 0.5*(mid_xi+max_xi)
            tree_panels_temp(index)%mid_eta = 0.5*(min_eta+mid_eta)
          endif
          if (j == 3) then
            tree_panels_temp(index)%min_xi = mid_xi; tree_panels_temp(index)%min_eta = mid_eta
            tree_panels_temp(index)%max_xi = max_xi; tree_panels_temp(index)%max_eta = max_eta
            tree_panels_temp(index)%mid_xi = 0.5*(mid_xi+max_xi)
            tree_panels_temp(index)%mid_eta = 0.5*(mid_eta+max_eta)
          endif
          if (j == 4) then
            tree_panels_temp(index)%min_xi = min_xi; tree_panels_temp(index)%min_eta = mid_eta
            tree_panels_temp(index)%max_xi = mid_xi; tree_panels_temp(index)%max_eta = max_eta
            tree_panels_temp(index)%mid_xi = 0.5*(min_xi+mid_xi)
            tree_panels_temp(index)%mid_eta = 0.5*(mid_eta+max_eta)
          endif
        endif
      enddo

      tree_panels_temp(i)%child_panel_count = kids
      allocate(tree_panels_temp(i)%child_panels(kids), source=-1)
      do j = 1, kids
        tree_panels_temp(i)%child_panels(j) = panel_count+j
      enddo

      curr_point_assigned_count(:) = 0
      do j = 1, tree_panels_temp(i)%panel_point_count      
        which_panel = point_panel(j)
        loc = panel_id(which_panel)         
        curr_point_assigned_count(loc) = curr_point_assigned_count(loc) + 1
        tree_panels_temp(panel_count+loc)%points_inside(curr_point_assigned_count(loc)) = tree_panels_temp(i)%points_inside(j)
      enddo

      do j = 1, kids ! shrink panels
        min_xi = 3.0
        max_xi = -3.0
        min_eta = 3.0
        max_eta = -3.0
        do k = 1, tree_panels_temp(panel_count+loc)%panel_point_count
          index = tree_panels_temp(panel_count+loc)%points_inside(k)
          min_xi = min(min_xi, point_xi(index))
          max_xi = max(max_xi, point_xi(index))
          min_eta = min(min_eta, point_eta(index))
          max_eta = max(max_eta, point_eta(index))
        enddo
        tree_panels_temp(panel_count+loc)%min_xi = min_xi-1e-15
        tree_panels_temp(panel_count+loc)%max_xi = max_xi+1e-15
        tree_panels_temp(panel_count+loc)%mid_xi = 0.5*(min_xi + max_xi)
        tree_panels_temp(panel_count+loc)%min_eta = min_eta-1e-15
        tree_panels_temp(panel_count+loc)%max_eta = max_eta+1e-15
        tree_panels_temp(panel_count+loc)%mid_eta = 0.5*(min_eta + max_eta)
      enddo

      ! sanity check to make sure all points are assigned
      total = 0
      do j = 1, kids
        total = total + tree_panels_temp(panel_count+j)%panel_point_count
      enddo
      if (total /= tree_panels_temp(i)%panel_point_count) THEN
        print *, 'Error in refining tree at parent panel ', i
      endif
      panel_count = panel_count + tree_panels_temp(i)%child_panel_count
    endif
    i = i + 1
  enddo

  ! tree_panels_temp is the wrong size so move everything over to the correctly sized tree_panels
  allocate(tree_panels(panel_count))
  do i = 1, panel_count
      tree_panels(i) = tree_panels_temp(i)
      tree_panels(i)%relabeled_points_inside = tree_panels(i)%points_inside
  enddo

  do i = 1, panel_count
    ! compute the furthest distance from panel center to vertex for each panel
    min_xi = tree_panels(i)%min_xi
    mid_xi = tree_panels(i)%mid_xi
    max_xi = tree_panels(i)%max_xi
    min_eta = tree_panels(i)%min_eta
    mid_eta = tree_panels(i)%mid_eta
    max_eta = tree_panels(i)%max_eta
    call xyz_from_xieta(x1, x2, x3, mid_xi, mid_eta, tree_panels(i)%face)
    call xyz_from_xieta(y1, y2, y3, min_xi, min_eta, tree_panels(i)%face)
    d1 = acos(max(min(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
    call xyz_from_xieta(y1, y2, y3, min_xi, max_eta, tree_panels(i)%face)
    d2 = acos(max(min(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
    call xyz_from_xieta(y1, y2, y3, max_xi, max_eta, tree_panels(i)%face)
    d3 = acos(max(min(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
    call xyz_from_xieta(y1, y2, y3, max_xi, min_eta, tree_panels(i)%face)
    d4 = acos(max(min(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
    tree_panels(i)%radius = max(d1, d2, d3, d4)
  enddo
end subroutine tree_traversal

!> assigns points to tree panels, making use of the tree structure
subroutine assign_points_to_panels(tree_panels, points_panels, point_leaf_panel, start_index, count)
  type(cube_panel), intent(inout) :: tree_panels(:) !< array of tree panels
  integer, intent(inout) :: points_panels(:,:) !< panel containing a point at each level
  integer, intent(out), allocatable :: point_leaf_panel(:) !< leaf panel containing each point
  integer, intent(in) :: start_index, count
  integer :: i, index, j, level

  allocate(point_leaf_panel(count))
  do i = 1, count
    index = i + start_index - 1
    level = 1
    j = 1
    jloop: do
      if (any(tree_panels(j)%points_inside == index)) then
        points_panels(level, i) = j
        point_leaf_panel(i) = j
        level = level + 1
        tree_panels(j)%own_point_count = tree_panels(j)%own_point_count + 1
        if (tree_panels(j)%is_leaf) then 
          exit jloop
        else
          j = tree_panels(j)%child_panels(1)
        endif
      else
        j = j + 1
      endif
    enddo jloop
  enddo
end subroutine assign_points_to_panels

!> computes the interaction list when running in tree code mode
subroutine interaction_list_compute_tc(pp_interactions, pc_interactions, tree_panels, x, y, z, &
      theta, cluster_thresh, point_count, base_panels_source) 
  type(interaction_pair), intent(out), allocatable :: pp_interactions(:), pc_interactions(:)
    !< list of particle/particle and particle/cluster interactions
  type(cube_panel), intent(in) :: tree_panels(:)
  real, intent(in) :: x(:), y(:), z(:), theta
  integer, intent(in) :: point_count, cluster_thresh, base_panels_source
  integer, allocatable :: source_index(:)
  type(interaction_pair), allocatable :: interaction_lists_temp(:)
  integer :: interaction_count, pp_count, pc_count, i, j, tt_count, k, curr_loc, i_s, c_s, l
  integer :: isc, iec, jsc, jec, ic, jc, index, tree_traverse_count, id
  real :: xco, yco, zco, xs, ys, zs, dist, separation, x1s, x2s, x3s, x1t, x2t, x3t
  logical :: well_separated

  allocate(source_index(size(tree_panels)))
  allocate(interaction_lists_temp(100*point_count))

  interaction_count = 0
  pp_count = 0
  pc_count = 0

  do i = 1, point_count
    tt_count = 0
    do k = 1, base_panels_source
      source_index(k) = k
      tt_count = tt_count + 1
    enddo
    curr_loc = 1
    xco = x(i)
    yco = y(i)
    zco = z(i)
    panelloop: do while (curr_loc <= tt_count) ! go through source panels
      i_s = source_index(curr_loc)
      c_s = tree_panels(i_s)%panel_point_count
      if (c_s > 0) then ! cluster has points
        call xyz_from_xieta(xs, ys, zs, tree_panels(i_s)%mid_xi, tree_panels(i_s)%mid_eta, &
            tree_panels(i_s)%face)
        dist = acos(min(max(xco*xs+yco*ys+zco*zs, -1.0), 1.0))
        separation = tree_panels(i_s)%radius/dist
        if ((dist > 0) .and. (separation < theta)) then
          interaction_count = interaction_count + 1
          if (c_s > cluster_thresh) then
            pc_count = pc_count + 1
            interaction_lists_temp(interaction_count)%index_target = i
            interaction_lists_temp(interaction_count)%index_source = i_s
            interaction_lists_temp(interaction_count)%interact_type = 1
          else 
            pp_count = pp_count + 1
            interaction_lists_temp(interaction_count)%index_target = i
            interaction_lists_temp(interaction_count)%index_source = i_s
            interaction_lists_temp(interaction_count)%interact_type = 0
          endif
        else  ! not well separated
          if (tree_panels(i_s)%is_leaf) then ! leaf, do pp interaction
            interaction_count = interaction_count + 1
            pp_count = pp_count + 1
            interaction_lists_temp(interaction_count)%index_target = i
            interaction_lists_temp(interaction_count)%index_source = i_s
            interaction_lists_temp(interaction_count)%interact_type = 0
          else ! refine source panel
            do j = 1, tree_panels(i_s)%child_panel_count
              source_index(tt_count+j) = tree_panels(i_s)%child_panels(j)
            enddo
            tt_count = tt_count + tree_panels(i_s)%child_panel_count
          endif
        endif
      endif
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
    endif
  enddo
end subroutine interaction_list_compute_tc

!> computes the interaction list in fmm mode
subroutine interaction_list_compute_fmm(pp_ints, pc_ints, cp_ints, cc_ints, source_tree, &
      target_tree, theta, cluster_thresh, ppc, base_panels_source, base_panels_target) 
  type(interaction_pair), allocatable, intent(out) :: pp_ints(:), pc_ints(:)
  type(interaction_pair), allocatable, intent(out) :: cp_ints(:), cc_ints(:)
    !< arrays for the four types of interactions
  type(cube_panel), intent(in) :: source_tree(:), target_tree(:)
  real, intent(in) :: theta
  integer, intent(in) :: cluster_thresh, ppc, base_panels_source, base_panels_target
  integer :: i, j, tree_traverse_count, curr_loc, pp_count, pc_count, cp_count, cc_count
  integer :: i_t, i_s, c_t, c_s, int_count, id
  integer, allocatable :: source_index(:), target_index(:), loc(:)
  type(interaction_pair), allocatable :: interaction_lists_temp(:)
  real :: x1t, x2t, x3t, x1s, x2s, x3s, dist, separation

  allocate(target_index(size(source_tree)*32))
  allocate(source_index(size(source_tree)*32))
  allocate(interaction_lists_temp(size(source_tree)*32))
  tree_traverse_count = 0
  curr_loc = 1
  pp_count = 0
  pc_count = 0
  cp_count = 0
  cc_count = 0
  int_count = 0

  id = PE_here()

  if (ppc > 0) then
    do i = 1, base_panels_source
      do j = 1, base_panels_target
        tree_traverse_count = tree_traverse_count + 1
        source_index(tree_traverse_count) = i
        target_index(tree_traverse_count) = j
      enddo
    enddo

    do while (curr_loc <= tree_traverse_count)
      i_t = target_index(curr_loc)
      i_s = source_index(curr_loc)
      c_t = target_tree(i_t)%panel_point_count
      c_s = source_tree(i_s)%panel_point_count
      if ((c_t > 0) .and. (c_s > 0)) then
        call xyz_from_xieta(x1t, x2t, x3t, target_tree(i_t)%mid_xi, target_tree(i_t)%mid_eta, &
            target_tree(i_t)%face)
        call xyz_from_xieta(x1s, x2s, x3s, source_tree(i_s)%mid_xi, source_tree(i_s)%mid_eta, &
            source_tree(i_s)%face)
        dist = acos(min(max(x1t*x1s+x2t*x2s+x3t*x3s, -1.0_8), 1.0_8))
        separation = (target_tree(i_t)%radius+source_tree(i_s)%radius)/dist
        if ((dist > 0) .and. (separation < theta)) then
          ! two panels are well separated
          int_count = int_count + 1
          interaction_lists_temp(int_count)%index_target = i_t
          interaction_lists_temp(int_count)%index_source = i_s
          if (c_s > cluster_thresh) then
            if (c_t > cluster_thresh) then ! CC 
              cc_count = cc_count + 1
              interaction_lists_temp(int_count)%interact_type = 3
            else ! PC
              pc_count = pc_count + 1
              interaction_lists_temp(int_count)%interact_type = 1
            endif
          else
            if (c_t > cluster_thresh) then ! CP
              cp_count = cp_count + 1
              interaction_lists_temp(int_count)%interact_type = 2
            else ! PP
              pp_count = pp_count + 1
              interaction_lists_temp(int_count)%interact_type = 0
            endif
          endif
        else
          ! two points are not well separated
          if ((c_t < cluster_thresh) .and. (c_s < cluster_thresh)) then ! both have few points, pp interaction
            int_count = int_count + 1
            pp_count = pp_count + 1
            interaction_lists_temp(int_count)%index_target = i_t
            interaction_lists_temp(int_count)%index_source = i_s
            interaction_lists_temp(int_count)%interact_type = 0
          else if (source_tree(i_s)%is_leaf .and. target_tree(i_t)%is_leaf) then
            ! both panels are leaves, pp interaction
            int_count = int_count + 1
            pp_count = pp_count + 1
            interaction_lists_temp(int_count)%index_target = i_t
            interaction_lists_temp(int_count)%index_source = i_s
            interaction_lists_temp(int_count)%interact_type = 0
          else if (target_tree(i_t)%is_leaf) then
            ! target panel is leaf, refine source
            do i = 1, source_tree(i_s)%child_panel_count
              target_index(tree_traverse_count+i) = i_t
              source_index(tree_traverse_count+i) = source_tree(i_s)%child_panels(i)
            enddo
            tree_traverse_count = tree_traverse_count + source_tree(i_s)%child_panel_count
          else if (source_tree(i_s)%is_leaf) then
            ! source panel is leaf, refine target
            do i = 1, target_tree(i_t)%child_panel_count
              target_index(tree_traverse_count+i) = target_tree(i_t)%child_panels(i)
              source_index(tree_traverse_count+i) = i_s
            enddo
            tree_traverse_count = tree_traverse_count + target_tree(i_t)%child_panel_count
          else
            ! neither is a leaf, refine the panel with more points
            if (c_s > c_t) then
              ! source has more points, refine source
              do i = 1, source_tree(i_s)%child_panel_count
                target_index(tree_traverse_count+i) = i_t
                source_index(tree_traverse_count+i) = source_tree(i_s)%child_panels(i)
              enddo
              tree_traverse_count = tree_traverse_count + source_tree(i_s)%child_panel_count
            else
              ! target has more points, refine target
              do i = 1, target_tree(i_t)%child_panel_count
                target_index(tree_traverse_count+i) = target_tree(i_t)%child_panels(i)
                source_index(tree_traverse_count+i) = i_s
              enddo
              tree_traverse_count = tree_traverse_count + target_tree(i_t)%child_panel_count
            endif
          endif
        endif
      endif
      curr_loc = curr_loc + 1
    enddo
  endif

  allocate(pp_ints(pp_count))
  allocate(pc_ints(pc_count))
  allocate(cp_ints(cp_count))
  allocate(cc_ints(cc_count))
  allocate(loc(4), source=0)
  
  do i = 1, int_count
    if (interaction_lists_temp(i)%interact_type == 0) then
      loc(1) = loc(1) + 1
      pp_ints(loc(1)) = interaction_lists_temp(i)
    else if (interaction_lists_temp(i)%interact_type == 1) then
      loc(2) = loc(2) + 1
      pc_ints(loc(2)) = interaction_lists_temp(i)
    else if (interaction_lists_temp(i)%interact_type == 2) then
      loc(3) = loc(3) + 1
      cp_ints(loc(3)) = interaction_lists_temp(i)
    else
      loc(4) = loc(4) + 1
      cc_ints(loc(4)) = interaction_lists_temp(i)
    endif
  enddo
end subroutine interaction_list_compute_fmm

!> computes the communication patterns for particle/particle and cluster/particle interactions
subroutine calculate_communications(sal_ct, xg, yg, zg, G)
  type(SAL_Conv_CS), intent(inout) :: sal_ct
  real, intent(in) :: xg(:), yg(:), zg(:)
  type(ocean_grid_type), intent(in) :: G
  integer :: p, id, unowned_source_count, pp_count, i, i_s, j, i_sp, j_sp, k, max_p, pl, count
  integer :: isg, ieg, jsg, jeg, isc, iec, jsc, jec
  integer :: i_off, j_off, i_sp2, j_sp2, isdg, iedg, jsdg, jedg
  integer :: idgo, jdgo, isd, jsd, own_points, cp_count
  ! integer, allocatable :: proc_start_i(:)
  ! integer, allocatable :: proc_start_j(:)
  ! integer, allocatable :: proc_end_i(:)
  ! integer, allocatable :: proc_end_j(:)
  ! integer, allocatable :: unowned_temp_i(:)
  ! integer, allocatable :: unowned_temp_j(:)
  ! integer, allocatable :: unowned_sources_i(:)
  ! integer, allocatable :: unowned_sources_j(:)
  integer, allocatable :: points_needed_from_proc(:) 
    !< the number of SSHs needed from every other MPI rank
  ! integer, allocatable :: points_from_proc_i(:,:)
  integer, allocatable :: proc_loc(:)
    !< which MPI rank owns a given point index for unowned source points
  integer, allocatable :: temp_locs(:)
    !< running total of points needed from each MPI rank
  ! integer, allocatable :: points_from_proc_j(:,:)
  ! integer, allocatable :: points_needed_from_procs(:,:)
  integer, allocatable :: points_to_give_proc(:)
    !< the number of SSHs to send to every other MPI rank
  integer, allocatable :: points_to_give_proc_i(:,:)
    !< indices of points/SSHs to send to other MPI ranks
  ! integer, allocatable :: points_to_give_proc_j(:,:)
  integer, allocatable :: unowned_temp(:)
    !< temp array for unowned source points
  integer, allocatable :: unowned_sources(:)
    !< indices of unowned source points
  integer, allocatable :: points_from_proc(:,:)
    !< point indices of SSHs to send to other MPI ranks
  integer, allocatable :: points_to_give_proc_index(:,:)
    !< local indices of points/SSHs to send to other MPI ranks
  ! integer, allocatable :: point_counts_to_communicate(:)
  integer, allocatable :: points_from_proc_2d(:,:,:)
    !< indices of points/SSHs being received from other MPI ranks
  integer, allocatable :: points_to_give_proc_2d(:,:,:) 
    !< indices of points/SSHs being sent to other MPI ranks
  real, allocatable :: e_xs(:), e_ys(:), e_zs(:)
    !< (x,y,z) coordinates of unowned source points
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
          endif
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
            endif
          enddo kloop2
        endif
      endif
    enddo
  enddo

  if (sal_ct%use_fmm) then
    cp_count = size(sal_ct%cp_interactions)
    do i = 1, cp_count
      i_s = sal_ct%cp_interactions(i)%index_source
      do j = 1, sal_ct%tree_struct(i_s)%panel_point_count
        i_sp = sal_ct%tree_struct(i_s)%points_inside(j)
        if ((i_sp < isc) .or. (i_sp > iec)) then 
          found = .false.
          kloop4: do k = 1, unowned_source_count
            if (unowned_temp(k) == i_sp) then
              found = .true.
            endif
          enddo kloop4
          if (.not. found) then
            unowned_source_count = unowned_source_count + 1
            unowned_temp(unowned_source_count) = i_sp
            kloop5: do k = 1, p
              if ((i_sp >= sal_ct%indexsg(k)) .and. (i_sp <= sal_ct%indexeg(k))) then
                points_needed_from_proc(k) = points_needed_from_proc(k) + 1
                proc_loc(unowned_source_count) = k
              end if
            enddo kloop5
          endif
        endif
      enddo
    enddo
  endif

  allocate(unowned_sources(unowned_source_count))
  allocate(e_xs(unowned_source_count+sal_ct%own_ocean_points))
  allocate(e_ys(unowned_source_count+sal_ct%own_ocean_points))
  allocate(e_zs(unowned_source_count+sal_ct%own_ocean_points))

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
        e_xs(count+sal_ct%own_ocean_points) = xg(i_sp)
        e_ys(count+sal_ct%own_ocean_points) = yg(i_sp)
        e_zs(count+sal_ct%own_ocean_points) = zg(i_sp)
      end if
    enddo
  enddo

  sal_ct%e_xs = e_xs
  sal_ct%e_ys = e_ys
  sal_ct%e_zs = e_zs
  allocate(points_to_give_proc(p), source=0)

  do i = 1,p
    call send_to_PE(points_needed_from_proc(i), i-1)
  enddo

  do i = 1,p
    call recv_from_PE(points_to_give_proc(i), i-1)
  enddo

  call sync_PEs()

  max_p = 0
  do i = 1, p
    max_p = max(max_p, points_to_give_proc(i))
  enddo

  allocate(points_to_give_proc_index(max_p, p), source=-1)
  allocate(points_to_give_proc_2d(1, max_p, p), source=-1)

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

  sal_ct%points_to_give = points_to_give_proc_index
  sal_ct%no_points_to_give_proc = points_to_give_proc
  sal_ct%points_to_get = points_from_proc
  sal_ct%no_points_to_get_proc = points_needed_from_proc

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
          if (sal_ct%tree_struct(j)%is_leaf) then
            exit treeloop
          else
            j = sal_ct%tree_struct(j)%child_panels(1)
          endif
        else
          j = j + 1
        endif
      endif
    enddo treeloop
  enddo
end subroutine calculate_communications

!> This subroutine initializes the convolution SAL control structure
subroutine sal_conv_init(sal_ct, G, param_file)
  type(SAL_Conv_CS), intent(out) :: sal_ct
  type(ocean_grid_type), intent(inout) :: G ! ocean grid
  type(param_file_type), intent(in) :: param_file
  integer :: proc_count, isc, iec, jsc, jec, isg, ieg, jsg, jeg, imax, jmax
  integer :: ic, jc, i, j, ig_off, jg_off
  integer :: max_level, proc_rank, i_off, j_off, isd, ied, jsd, jed, rank, itc, jtc
  integer :: count, pointcount, cluster_thresh
  integer, allocatable :: points_panels(:,:,:), indexsg(:), indexeg(:), pcg(:)
  real, allocatable :: xg(:,:), yg(:,:), zg(:,:)
    !< (x,y,z) coordinates of the global grid
  real, allocatable :: xc(:,:), yc(:,:), zc(:,:)
    !< (x,y,z) coordinates of the computational domain for each MPI rank
  ! real, allocatable :: xt(:,:), yt(:,:), zt(:,:)
  real, allocatable :: xg1d(:), yg1d(:), zg1d(:)
    !< (x,y,z) coordinates of the global grid in a 1d array
  real, allocatable :: xc1d(:), yc1d(:), zc1d(:)
    !< (x,y,z) coordinates of the computational domain in a 1d array
  ! real, allocatable :: xt1d(:), yt1d(:), zt1d(:)
  real :: lat, lon, colat, x, y, z, pi, theta
  integer :: base_panels_source, base_panels_target
  logical :: sht

# include "version_variable.h"
  character(len=40) :: mdl = "MOM_conv_self_attr_load" ! This module's name.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CONV_SAL_THETA", theta, &
                  "Theta parameter for Convolution SAL Tree Code", units="m m-1", default=0.9)
  call get_param(param_file, mdl, "CONV_SAL_CLUSTER_SIZE", sal_ct%cluster_thresh, &
                  "Cluster threshold size for Convolution SAL Tree Code", default=50)
  call get_param(param_file, mdl, "CONV_SAL_INTERP_DEGREE", sal_ct%interp_degree, &
                  "Interpolation degree for Convolution SAL Tree Code", default=2)
  call get_param(param_file, mdl, "CONV_SAL_REPRODUCING", sal_ct%reprod_sum, &
                  "Convolution SAL Tree Code Reproducing sum mode", default=.false.)
  call get_param(param_file, mdl, "SAL_CONVOLUTION", sal_ct%use_sal_conv, &
                  "Whether or not to use SAL convolution", default=.false.)
  call get_param(param_file, mdl, "SAL_HARMONICS", sht, "SHT", default=.false., do_not_log=.True.)
  call get_param(param_file, mdl, "SAL_USE_BPA", sal_ct%use_bpa, &
                 "If true, use bottom pressure anomaly to calculate self-attraction and "// &
                 "loading (SAL). Otherwise sea surface height anomaly is used, which is "// &
                 "only accurate for uniform density fluid.", default=.false.)
  ! if (sht) then
  !   call get_param(param_file, mdl, "CONV_SAL_TRUNC_THRESH", sal_ct%trunc_thresh, &
  !                   "Truncation threshold for Conv SAL Greens Function", units="m m-1", &
  !                   default=0.9, do_not_log=.True.)
  ! else
  !   sal_ct%trunc_thresh = -2.0
  ! endif
  sal_ct%use_farfield = .not. sht
  if (.not. sal_ct%reprod_sum) then
    call get_param(param_file, mdl, "CONV_SAL_FMM", sal_ct%use_fmm, &
                "Convolution SAL FMM mode or not, FMM is incompatible with reproducing sum", &
                default=.true.)
  else 
    sal_ct%use_fmm = .false.
  endif

  if (sal_ct%use_bpa) then 
    sal_ct%scaling_constant = 3.148e-13 ! [M-1 L T2 ~> m3 J-1]
  else
    sal_ct%scaling_constant = 3.19535026078445591e-9 ! [L-1]
  endif

  if (sal_ct%use_sal_conv) then
    sal_ct%p = num_PEs() ! number of ranks
    sal_ct%id = PE_here() ! current rank

    pi = 4.D0*DATAN(1.D0)

    isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec
    isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed

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
        endif
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
    call sum_across_PEs(xg1d, pointcount) ! potentially a problem at high resolution, memory issues
    call sum_across_PEs(yg1d, pointcount)
    call sum_across_PEs(zg1d, pointcount)

    ! xg/yg/zg is now a copy of all the points from all the processors
    call tree_traversal(G, sal_ct%tree_struct, xg1d, yg1d, zg1d, sal_ct%cluster_thresh, &
                        pointcount, base_panels_source) ! constructs cubed sphere tree
    max_level = sal_ct%tree_struct(size(sal_ct%tree_struct))%level

    if (sal_ct%use_fmm) then
      call tree_traversal(G, sal_ct%tree_struct_targets, xc1d, yc1d, zc1d, sal_ct%cluster_thresh, &
                          count, base_panels_target) ! constructs tree of targets
    endif

    allocate(sal_ct%points_panels(max_level+1, ic*jc), source=-1)
    ! finds which panels contain the computational domain points
    call assign_points_to_panels(sal_ct%tree_struct, sal_ct%points_panels, sal_ct%point_leaf_panel, &
                                  sal_ct%indexsg(sal_ct%id+1), sal_ct%own_ocean_points)
    ! compute the interaction lists for the target points in the target domain
    if (sal_ct%use_fmm) then
      call interaction_list_compute_fmm(sal_ct%pp_interactions, sal_ct%pc_interactions, &
                                        sal_ct%cp_interactions, sal_ct%cc_interactions, &
                                        sal_ct%tree_struct, sal_ct%tree_struct_targets, &
                                        theta, sal_ct%cluster_thresh, sal_ct%own_ocean_points, &
                                        base_panels_source, base_panels_target)
    else
      call interaction_list_compute_tc(sal_ct%pp_interactions, sal_ct%pc_interactions, &
                                        sal_ct%tree_struct, xc1d, yc1d, zc1d, theta, &
                                        sal_ct%cluster_thresh, sal_ct%own_ocean_points, &
                                        base_panels_source)
    endif
    call sync_PEs()
    ! compute communication patterns 
    call calculate_communications(sal_ct, xg1d, yg1d, zg1d, G)

    do i = 1, sal_ct%own_ocean_points
      sal_ct%e_xs(i) = xc1d(i)
      sal_ct%e_ys(i) = yc1d(i)
      sal_ct%e_zs(i) = zc1d(i)
    enddo

    id_clock_SAL = cpu_clock_id('(Ocean SAL)', grain=CLOCK_MODULE)
    id_clock_SAL_ssh_comm = cpu_clock_id('(Ocean SAL SSH communication)', grain=CLOCK_MODULE)
    
    if (sal_ct%use_fmm) then
      id_clock_SAL_fmm = cpu_clock_id('(Ocean SAL FMM)', grain=CLOCK_MODULE)
      id_clock_SAL_upward_pass = cpu_clock_id('(Ocean SAL Upward Pass)', grain=CLOCK_MODULE)
      id_clock_SAL_downward_pass = cpu_clock_id('(Ocean SAL downward pass)', grain=CLOCK_MODULE)
      id_clock_SAL_cc_comp = cpu_clock_id('(Ocean SAL CC interaction comp)', grain=CLOCK_MODULE)
      id_clock_SAL_cp_comp = cpu_clock_id('(Ocean SAL CP interaction comp)', grain=CLOCK_MODULE)
    else
      id_clock_SAL_tc = cpu_clock_id('(Ocean SAL Tree Code)', grain=CLOCK_MODULE)
      id_clock_SAL_upward_pass = cpu_clock_id('(Ocean SAL Upward Pass)', grain=CLOCK_MODULE)
    endif
    id_clock_SAL_pc_comp = cpu_clock_id('(Ocean SAL PC interaction comp)', grain=CLOCK_MODULE)
    id_clock_SAL_pp_comp = cpu_clock_id('(Ocean SAL PP interaction comp)', grain=CLOCK_MODULE)
  endif
end subroutine sal_conv_init

!> does the necessary MPI communications for particle/particle and particle/cluster interactions
subroutine ssh_communications(sal_ct, G, eta, e_ssh)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(in) :: eta(:,:) !< sea surface heights or bottom pressure anomalies
  real, intent(inout) :: e_ssh(:) !< data points received from other MPI ranks
  integer :: max_give, max_get, p, id, i, j, i_s, j_s, i_off, j_off, count, index
  real, allocatable :: points_to_give(:,:,:)
    !< SSH/BPAs to send to other MPI ranks
  real, allocatable :: points_received(:,:,:)
    !< SSH/BPAs received from other MPI ranks
  real :: area, r2

  p = sal_ct%p; id = sal_ct%id
  i_off = G%isg - G%isc; j_off = G%jsg - G%jsc
  r2 = G%Rad_Earth_L * G%Rad_Earth_L 

  call cpu_clock_begin(id_clock_SAL_ssh_comm)

  max_give = 0; max_get = 0
  do i = 1, p
    max_give = max(max_give, sal_ct%no_points_to_give_proc(i))
    max_get = max(max_get, sal_ct%no_points_to_get_proc(i))
  enddo

  allocate(points_to_give(1, max_give, p), source=0.0)
  allocate(points_received(1, max_get, p), source=0.0)

  do i = 1, p
    do j = 1, max_give
      index = sal_ct%points_to_give(j, i)
      if (index .ne. -1) then
        i_s = sal_ct%one_d_to_2d_i(index)
        j_s = sal_ct%one_d_to_2d_j(index)
        area = G%areaT(i_s, j_s)/r2
        points_to_give(1, j, i) = eta(i_s, j_s)*area
      endif
    enddo
  enddo

  do i=0, p-1 ! send points
    if (sal_ct%no_points_to_give_proc(i+1)>0) then
      call send_to_PE(points_to_give(:,:,i+1), sal_ct%no_points_to_give_proc(i+1), i)
    endif
  enddo

  do i=0, p-1 ! receive points
    if (sal_ct%no_points_to_get_proc(i+1)>0) then
      call recv_from_PE(points_received(:,:,i+1), sal_ct%no_points_to_get_proc(i+1), i)
    endif
  enddo

  call sync_PEs()

  count = 0
  do i = 1, p
    do j = 1, sal_ct%no_points_to_get_proc(i)
      count = count + 1
      e_ssh(count+sal_ct%own_ocean_points) = points_received(1, j, i)
    enddo
  enddo

  do i = 1, sal_ct%own_ocean_points
    i_s = sal_ct%one_d_to_2d_i(i)
    j_s = sal_ct%one_d_to_2d_j(i)
    area = G%areaT(i_s, j_s)/r2
    e_ssh(i) = eta(i_s,j_s)*area
  enddo

  call cpu_clock_end(id_clock_SAL_ssh_comm)
end subroutine ssh_communications

!> computes the interpolation weights for barycentric Lagrange interpolation 
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
      endif
    enddo
  endif
end subroutine bli_coeffs

!> computes the shifted Chebyshev nodes in the interval [min_x, max_x]
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
      interp_points(i) = cos(pi/degree*(i-1))*x_range + x_shift
    enddo
  endif
end subroutine bli_interp_points_shift

!> returns the barycentric Lagrange basis values for a point (xi, eta)
subroutine interp_vals_bli(vals, xi, eta, min_xi, max_xi, min_eta, max_eta, degree)
  real, allocatable, intent(out) :: vals(:,:)
  real, intent(in) :: xi, eta, min_xi, max_xi, min_eta, max_eta
  integer, intent(in) :: degree
  real, allocatable :: bli_xi_vals(:), bli_eta_vals(:), bli_coeff_vals(:)
  real, allocatable :: xi_func_vals(:), eta_func_vals(:)
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
    if (abs(xi - bli_xi_vals(i)) < 1.0e-15) then
      found_xi = .true.
      xi_func_vals(i) = 1.0
    endif
    if (abs(eta - bli_eta_vals(i)) < 1.0e-15) then
      found_eta = .true.
      eta_func_vals(i) = 1.0
    endif
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
  endif

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
  endif

  ! compute the outer product
  do j = 1, degree + 1
    do i = 1, degree + 1
      vals(i, j) = xi_func_vals(i) * eta_func_vals(j)
    enddo
  enddo
end subroutine interp_vals_bli

!> computes the upward pass without using reproducing sums
subroutine proxy_source_compute_nonrep(sal_ct, G, e_ssh, proxy_source_weights) 
  type(SAL_Conv_CS), intent(in) :: sal_ct
  real, intent(in) :: e_ssh(:)
  real, intent(inout) :: proxy_source_weights(:)
  type(ocean_grid_type), intent(in) :: G
  real :: pi, r2, lat, lon, colat, x, y, z, a, ssh, min_xi, max_xi, min_eta, max_eta, xi, eta
  real :: min_xi_p, max_xi_p, min_eta_p, max_eta_p, val
  integer :: isc, iec, jsc, jec, ic, jc, i, ii, ij, leaf_i, shift, j, k
  integer :: parent_i, shift1, offset, l, m, offset1
  real, allocatable :: basis_vals(:,:), xi_vals(:), eta_vals(:)

  pi = 4.0*DATAN(1.0)
  isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec
  ic = iec-isc+1; jc = jec-jsc+1
  r2 = G%Rad_Earth_L * G%Rad_Earth_L 

  call cpu_clock_begin(id_clock_SAL_upward_pass)

  ! calculate proxy source potentials for leaf panels
  do i = 1, sal_ct%own_ocean_points
    x = sal_ct%e_xs(i)
    y = sal_ct%e_ys(i)
    z = sal_ct%e_zs(i)
    ssh = e_ssh(i)
    leaf_i = sal_ct%point_leaf_panel(i)
    shift = (leaf_i-1)*(sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
    min_xi = sal_ct%tree_struct(leaf_i)%min_xi
    max_xi = sal_ct%tree_struct(leaf_i)%max_xi
    min_eta = sal_ct%tree_struct(leaf_i)%min_eta
    max_eta = sal_ct%tree_struct(leaf_i)%max_eta
    call xieta_from_xyz(x, y, z, xi, eta, sal_ct%tree_struct(leaf_i)%face)
    call interp_vals_bli(basis_vals, xi, eta, min_xi, max_xi, min_eta, max_eta, &
                          sal_ct%interp_degree)
    offset = 0
    do j = 1, sal_ct%interp_degree+1
      do k= 1, sal_ct%interp_degree+1
        offset=offset+1
        proxy_source_weights(shift+offset)=proxy_source_weights(shift+offset)+basis_vals(k, j)*ssh
      enddo
    enddo
  enddo

  ! upward pass to interpolate from child to parent panel proxy source potentials
  do i = size(sal_ct%tree_struct), 1, -1
    parent_i = sal_ct%tree_struct(i)%parent_panel
    if ((parent_i > 0) .and. (sal_ct%tree_struct(i)%own_point_count > 0)) then
        min_xi = sal_ct%tree_struct(i)%min_xi
      max_xi = sal_ct%tree_struct(i)%max_xi
      min_eta = sal_ct%tree_struct(i)%min_eta
      max_eta = sal_ct%tree_struct(i)%max_eta
      call bli_interp_points_shift(xi_vals, min_xi, max_xi, sal_ct%interp_degree)
      call bli_interp_points_shift(eta_vals, min_eta, max_eta, sal_ct%interp_degree)
      ! loop over proxy source points in child panel
      min_xi_p = sal_ct%tree_struct(parent_i)%min_xi
      max_xi_p = sal_ct%tree_struct(parent_i)%max_xi
      min_eta_p = sal_ct%tree_struct(parent_i)%min_eta
      max_eta_p = sal_ct%tree_struct(parent_i)%max_eta
      shift = (parent_i-1)*(sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
      shift1 = (i-1)*(sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
      offset1 = 0
      do j = 1, sal_ct%interp_degree+1
        do k = 1, sal_ct%interp_degree+1
          call interp_vals_bli(basis_vals, xi_vals(k), eta_vals(j), min_xi_p, max_xi_p, min_eta_p, &
                                  max_eta_p, sal_ct%interp_degree)
          offset = 0
          offset1 = offset1+1
          val = proxy_source_weights(shift1+offset1)
          do l = 1, sal_ct%interp_degree+1
            do m = 1, sal_ct%interp_degree+1
              offset = offset+1
              proxy_source_weights(shift+offset) = proxy_source_weights(shift+offset)+ &
                                                    basis_vals(m,l)*val
            enddo
          enddo
        enddo
      enddo
    endif
  enddo

  call sum_across_PEs(proxy_source_weights, size(proxy_source_weights))
  call cpu_clock_end(id_clock_SAL_upward_pass)
end subroutine proxy_source_compute_nonrep

! perform the upward pass using reproducing sums
subroutine proxy_source_compute_reprod(sal_ct, G, e_ssh, proxy_source_weights) 
  type(SAL_Conv_CS), intent(in) :: sal_ct
  real, intent(in) :: e_ssh(:)
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
  r2 = G%Rad_Earth_L * G%Rad_Earth_L 

  call cpu_clock_begin(id_clock_SAL_upward_pass)
  allocate(proxy_source_weights_sep(1, sal_ct%own_ocean_points, size(proxy_source_weights)), &
            source=0.0)

  do i = 1, sal_ct%own_ocean_points
    ii = sal_ct%one_d_to_2d_i(i)
    ij = sal_ct%one_d_to_2d_j(i)
    x = sal_ct%e_xs(i)
    y = sal_ct%e_ys(i)
    z = sal_ct%e_zs(i)
    ssh = e_ssh(i)
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
        call interp_vals_bli(basis_vals, xi, eta, min_xi, max_xi, min_eta, max_eta, &
                              sal_ct%interp_degree)
        offset = 0
        do l = 1, sal_ct%interp_degree+1
          do m = 1, sal_ct%interp_degree+1
            offset=offset+1
            proxy_source_weights_sep(1, i, shift+offset) = basis_vals(m, l)*ssh
          enddo
        enddo
      endif
    enddo panelloop2
  enddo
  sum_tot = reproducing_sum(proxy_source_weights_sep(:,:,1:size(proxy_source_weights)), &
                            sums=proxy_source_weights(:))
  call cpu_clock_end(id_clock_SAL_upward_pass)
end subroutine proxy_source_compute_reprod

!> Computes the SAL gradient Green's function for two points (tx, ty, tz) and (sx, sy, sz)
subroutine sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal_x, sal_y) 
  real, intent(in) :: tx, ty, tz, sx, sy, sz
  real, intent(out) :: sal_x, sal_y
  real :: g, mp, sqrtp, sqp, p1, p2, x32m, mp2iv, eps, coeff
  real :: sum, pn1, pn2, pn
  integer :: i

  eps=1e-6

  sal_x = 0.0
  sal_y = 0.0
  if ((abs(tz - 1.0) > 1e-15) .and. (abs(tz+1.0) > 1e-15)) then
    g = max(min(tx*sx+ty*sy+tz*sz, 1.0), -1.0) 
      !< prevent floating point errors, dot product of x and y
    mp = 2.0-2.0*g ! 2-2 x.y
    sqp = sqrt(mp) ! sqrt(2-2x.y)
    p1 = (1.0+6.21196)/(sqp*mp+eps) ! (1-b0)/((2-2x.y)^(3/2))
    p2 = (2.7+6.12)*(2.0*g+sqp) / (2.0*(g*g-1.0)+eps) ! (a1-b1)(2x.y+sqrt(2-2x.y))/(2((x.y)^2-1))
    x32m = 1.0-tz*tz ! 1-(x3)^2
    mp2iv = (p1+p2)/sqrt(x32m) ! 3 rho_w/(4pi rho_e R_E sqrt(1-(x3)^2))
    sal_y = (sz*x32m-tz*(tx*sx+ty*sy))*mp2iv
    sal_x = (tx*sy-ty*sx)*mp2iv
  endif
end subroutine sal_grad_gfunc

!> computes a particle-cluster interaction, not reproducing
subroutine pc_interaction_compute_nonreprod(sal_ct, G, proxy_source_weights, sal_x, sal_y)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(in) :: proxy_source_weights(:)
  real, intent(inout) :: sal_x(:,:), sal_y(:,:)
  integer :: proxy_count, i, i_s, i_ti, i_tj, j, k, offset, i_t, isc, jsc
  real, allocatable :: source_proxy_weights(:), cheb_xi(:), cheb_eta(:)
  real :: min_xi, max_xi, min_eta, max_eta, x, y, z, lat, lon, xi, eta, colat, pi
  real :: cx, cy, cz, sal_grad_x, sal_grad_y, sal_val

  pi = 4.0*DATAN(1.0)
  isc = G%isc; jsc = G%jsc

  call cpu_clock_begin(id_clock_SAL_pc_comp)
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
    x = sal_ct%e_xs(i_t)
    y = sal_ct%e_ys(i_t)
    z = sal_ct%e_zs(i_t)
    offset = 0
    do k = 1, sal_ct%interp_degree+1
      eta = cheb_eta(k)
      do j = 1, sal_ct%interp_degree+1
        xi = cheb_xi(j)
        sal_grad_x = 0.0
        sal_grad_y = 0.0
        call xyz_from_xieta(cx, cy, cz, xi, eta, sal_ct%tree_struct(i_s)%face)
        call sal_grad_gfunc(x, y, z, cx, cy, cz, sal_grad_x, sal_grad_y)
        offset = offset+1
        sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + sal_grad_x*source_proxy_weights(offset)
        sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + sal_grad_y*source_proxy_weights(offset)
      enddo
    enddo
  enddo
  call cpu_clock_end(id_clock_SAL_pc_comp)
end subroutine pc_interaction_compute_nonreprod

!> computes a particle-cluster interaction, reproducing
subroutine pc_interaction_compute_reprod(sal_ct, G, proxy_source_weights, sal_x, sal_y)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(in) :: proxy_source_weights(:)
  real, intent(inout) :: sal_x(:,:), sal_y(:,:)
  integer :: proxy_count, i, i_s, i_ti, i_tj, j, k, offset, i_t, isc, jsc
  real, allocatable :: source_proxy_weights(:), cheb_xi(:), cheb_eta(:)
  real :: min_xi, max_xi, min_eta, max_eta, x, y, z, lat, lon, xi, eta, colat, pi
  real :: cx, cy, cz, sal_grad_x, sal_grad_y, sal_val
  real, allocatable :: sal_x_vals(:,:)
    !< values of sal_grad_x * BLI proxy weights
  real, allocatable :: sal_y_vals(:,:)
    !< values of sal_grad_y * BLI proxy weights

  pi = 4.0*DATAN(1.0)
  isc = G%isc; jsc = G%jsc

  call cpu_clock_begin(id_clock_SAL_pc_comp)
  proxy_count = (sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
  allocate(source_proxy_weights(proxy_count), source=0.0)
  allocate(sal_x_vals(sal_ct%interp_degree+1, sal_ct%interp_degree+1))
  allocate(sal_y_vals(sal_ct%interp_degree+1, sal_ct%interp_degree+1))
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
    x = sal_ct%e_xs(i_t)
    y = sal_ct%e_ys(i_t)
    z = sal_ct%e_zs(i_t)
    offset = 0
    do k = 1, sal_ct%interp_degree+1
      eta = cheb_eta(k)
      do j = 1, sal_ct%interp_degree+1
        xi = cheb_xi(j)
        sal_grad_x = 0.0
        sal_grad_y = 0.0
        call xyz_from_xieta(cx, cy, cz, xi, eta, sal_ct%tree_struct(i_s)%face)
        call sal_grad_gfunc(x, y, z, cx, cy, cz, sal_grad_x, sal_grad_y)
        offset = offset+1
        sal_x_vals(j,k) = sal_grad_x * source_proxy_weights(offset)
        sal_y_vals(j,k) = sal_grad_y * source_proxy_weights(offset)
      enddo
    enddo
    sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + reproducing_sum(sal_x_vals, only_on_PE=.true.)
    sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + reproducing_sum(sal_y_vals, only_on_PE=.true.)
  enddo
  call cpu_clock_end(id_clock_SAL_pc_comp)
end subroutine pc_interaction_compute_reprod

!> particle-cluster interaction computation in FMM mode
subroutine pc_interaction_compute_fmm(sal_ct, G, proxy_source_weights, sal_x, sal_y)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(in) :: proxy_source_weights(:)
  real, intent(inout) :: sal_x(:,:), sal_y(:,:)
  integer :: proxy_count, i, i_s, i_ti, i_tj, j, k, offset, i_t, isc, jsc, l, i_tk
  real, allocatable :: source_proxy_weights(:), cheb_xi(:), cheb_eta(:)
  real :: min_xi, max_xi, min_eta, max_eta, x, y, z, lat, lon, xi, eta, colat, pi
  real :: cx, cy, cz, sal_grad_x, sal_grad_y, sal_val

  pi = 4.0*DATAN(1.0)
  isc = G%isc; jsc = G%jsc

  call cpu_clock_begin(id_clock_SAL_pc_comp)
  proxy_count = (sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
  allocate(source_proxy_weights(proxy_count), source=0.0)
  do i = 1, size(sal_ct%pc_interactions)
    i_s = sal_ct%pc_interactions(i)%index_source
    i_t = sal_ct%pc_interactions(i)%index_target
    source_proxy_weights = proxy_source_weights(proxy_count*(i_s-1)+1:proxy_count*i_s)
    min_xi = sal_ct%tree_struct(i_s)%min_xi
    max_xi = sal_ct%tree_struct(i_s)%max_xi
    min_eta = sal_ct%tree_struct(i_s)%min_eta
    max_eta = sal_ct%tree_struct(i_s)%max_eta
    call bli_interp_points_shift(cheb_xi, min_xi, max_xi, sal_ct%interp_degree)
    call bli_interp_points_shift(cheb_eta, min_eta, max_eta, sal_ct%interp_degree)
    do l = 1, sal_ct%tree_struct_targets(i_t)%panel_point_count
      i_tk = sal_ct%tree_struct_targets(i_t)%points_inside(l)
      i_ti = sal_ct%one_d_to_2d_i(i_tk)
      i_tj = sal_ct%one_d_to_2d_j(i_tk)
      x = sal_ct%e_xs(i_tk)
      y = sal_ct%e_ys(i_tk)
      z = sal_ct%e_zs(i_tk)
      offset = 0
      do k = 1, sal_ct%interp_degree+1
        eta = cheb_eta(k)
        do j = 1, sal_ct%interp_degree+1
          xi = cheb_xi(j)
          call xyz_from_xieta(cx, cy, cz, xi, eta, sal_ct%tree_struct(i_s)%face)
          call sal_grad_gfunc(x, y, z, cx, cy, cz, sal_grad_x, sal_grad_y)
          offset = offset+1
          sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + sal_grad_x*source_proxy_weights(offset) 
          sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + sal_grad_y*source_proxy_weights(offset)
        enddo
      enddo
    enddo 
  enddo
  call cpu_clock_end(id_clock_SAL_pc_comp)
end subroutine pc_interaction_compute_fmm

!> cluster cluster interaction for FMM mode
subroutine cc_interaction_compute(sal_ct, proxy_source_weights, proxy_target_weights_x, &
                                  proxy_target_weights_y)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  real, intent(in) :: proxy_source_weights(:)
  real, intent(inout) :: proxy_target_weights_x(:,:,:), proxy_target_weights_y(:,:,:)
  real :: min_x_s, max_x_s, min_e_s, max_e_s, min_x_t, max_x_t, min_e_t, max_e_t
  real :: e_t, x_t, e_s, x_s
  real :: t_x, t_y, t_z, s_x, s_y, s_z, sal_grad_x, sal_grad_y
  integer :: i, i_s, i_t, j, shift_s, offset_s, k, l, m
  real, allocatable :: cheb_x_s(:), cheb_e_s(:), cheb_x_t(:), cheb_e_t(:)

  call cpu_clock_begin(id_clock_SAL_cc_comp)

  do i = 1, size(sal_ct%cc_interactions)
    i_s = sal_ct%cc_interactions(i)%index_source
    i_t = sal_ct%cc_interactions(i)%index_target    
    min_x_s = sal_ct%tree_struct(i_s)%min_xi
    max_x_s = sal_ct%tree_struct(i_s)%max_xi
    min_e_s = sal_ct%tree_struct(i_s)%min_eta
    max_e_s = sal_ct%tree_struct(i_s)%max_eta
    min_x_t = sal_ct%tree_struct_targets(i_t)%min_xi
    max_x_t = sal_ct%tree_struct_targets(i_t)%max_xi
    min_e_t = sal_ct%tree_struct_targets(i_t)%min_eta
    max_e_t = sal_ct%tree_struct_targets(i_t)%max_eta
    call bli_interp_points_shift(cheb_x_s, min_x_s, max_x_s, sal_ct%interp_degree)
    call bli_interp_points_shift(cheb_e_s, min_e_s, max_e_s, sal_ct%interp_degree)
    call bli_interp_points_shift(cheb_x_t, min_x_t, max_x_t, sal_ct%interp_degree)
    call bli_interp_points_shift(cheb_e_t, min_e_t, max_e_t, sal_ct%interp_degree)
    shift_s = (i_s-1)*(sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)
    do j = 1, sal_ct%interp_degree+1 ! target eta loop
      e_t = cheb_e_t(j)
      do m = 1, sal_ct%interp_degree+1 ! target xi loop
        x_t = cheb_x_t(m)
        call xyz_from_xieta(t_x, t_y, t_z, x_t, e_t, sal_ct%tree_struct_targets(i_t)%face)
        offset_s = 0
        do l = 1, sal_ct%interp_degree+1 ! source eta loop
          e_s = cheb_e_s(l)
          do k = 1, sal_ct%interp_degree+1 ! source xi loop
            x_s = cheb_x_s(k)
            call xyz_from_xieta(s_x, s_y, s_z, x_s, e_s, sal_ct%tree_struct(i_s)%face)
            offset_s = offset_s + 1
            call sal_grad_gfunc(t_x, t_y, t_z, s_x, s_y, s_z, sal_grad_x, sal_grad_y)
            proxy_target_weights_x(m,j,i_t) = proxy_target_weights_x(m,j,i_t) + &
                                    proxy_source_weights(shift_s+offset_s)*sal_grad_x
            proxy_target_weights_y(m,j,i_t) = proxy_target_weights_y(m,j,i_t) + &
                                    proxy_source_weights(shift_s+offset_s)*sal_grad_y
          enddo
        enddo
      enddo
    enddo
  enddo
  call cpu_clock_end(id_clock_SAL_cc_comp)
end subroutine cc_interaction_compute

!> cluster particle interaction for FMM mode
subroutine cp_interaction_compute(sal_ct, G, eta, e_ssh, proxy_target_weights_x, proxy_target_weights_y)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(in) :: e_ssh(:), eta(:,:)
  real, intent(inout) :: proxy_target_weights_x(:,:,:), proxy_target_weights_y(:,:,:)
  real :: pi, r2, lat, lon, colat, sx, sy, sz, min_eta, max_eta, min_xi, max_xi, ssh, eta_t, xi_t
  real :: sal_grad_x, sal_grad_y, tx, ty, tz
  integer :: i, i_s, i_t, j, i_sp, k, l, i_si, i_sj
  real, allocatable :: cheb_xi(:), cheb_eta(:)

  pi = 4.0*DATAN(1.0)
  r2 = G%Rad_Earth_L * G%Rad_Earth_L 

  call cpu_clock_begin(id_clock_SAL_cp_comp)

  do i = 1, size(sal_ct%cp_interactions)
    i_s = sal_ct%cp_interactions(i)%index_source
    i_t = sal_ct%cp_interactions(i)%index_target
    min_xi = sal_ct%tree_struct_targets(i_t)%min_xi
    max_xi = sal_ct%tree_struct_targets(i_t)%max_xi
    min_eta = sal_ct%tree_struct_targets(i_t)%min_eta
    max_eta = sal_ct%tree_struct_targets(i_t)%max_eta
    call bli_interp_points_shift(cheb_xi, min_xi, max_xi, sal_ct%interp_degree)
    call bli_interp_points_shift(cheb_eta, min_eta, max_eta, sal_ct%interp_degree)
    do j = 1, sal_ct%tree_struct(i_s)%panel_point_count
      i_sp = sal_ct%tree_struct(i_s)%relabeled_points_inside(j)
      sx = sal_ct%e_xs(i_sp)
      sy = sal_ct%e_ys(i_sp)
      sz = sal_ct%e_zs(i_sp)
      ssh = e_ssh(i_sp)
      do k = 1, sal_ct%interp_degree+1 ! target eta loop
        eta_t = cheb_eta(k)
        do l = 1, sal_ct%interp_degree+1 ! target xi loop
          xi_t = cheb_xi(l)
          call xyz_from_xieta(tx, ty, tz, xi_t, eta_t, sal_ct%tree_struct_targets(i_t)%face)
          call sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal_grad_x, sal_grad_y)
          proxy_target_weights_x(l,k,i_t) = proxy_target_weights_x(l,k,i_t)+sal_grad_x*ssh
          proxy_target_weights_y(l,k,i_t) = proxy_target_weights_y(l,k,i_t)+sal_grad_y*ssh
        enddo
      enddo
    enddo
  enddo
  call cpu_clock_end(id_clock_SAL_cp_comp)
end subroutine cp_interaction_compute

!> particle particle interaction for tree code mode, without reproducing sums
subroutine pp_interaction_compute_nonreprod(sal_ct, G, sal_x, sal_y, e_ssh)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(inout) :: sal_x(:,:), sal_y(:,:)
  real, intent(in) :: e_ssh(:)
  integer :: i, i_s, i_ti, i_tj, j, i_si, i_sj, i_sp, i_t, id, count, ii, ij
  real :: pi, lat, lon, colat, x, y, z, sx, sy, sz, ssh, r2, sal_grad_x, sal_grad_y, sal_val
  real, allocatable :: sshs(:), a_s(:)

  pi = 4.0*DATAN(1.0)
  r2 = G%Rad_Earth_L * G%Rad_Earth_L 

  call cpu_clock_begin(id_clock_SAL_pp_comp)

  do i = 1, size(sal_ct%pp_interactions)
    i_s = sal_ct%pp_interactions(i)%index_source
    i_t = sal_ct%pp_interactions(i)%index_target
    i_ti = sal_ct%one_d_to_2d_i(i_t)
    i_tj = sal_ct%one_d_to_2d_j(i_t)
    x = sal_ct%e_xs(i_t)
    y = sal_ct%e_ys(i_t)
    z = sal_ct%e_zs(i_t)
    do j = 1, sal_ct%tree_struct(i_s)%panel_point_count
      i_sp = sal_ct%tree_struct(i_s)%relabeled_points_inside(j)
      sx = sal_ct%e_xs(i_sp)
      sy = sal_ct%e_ys(i_sp)
      sz = sal_ct%e_zs(i_sp)
      ssh = e_ssh(i_sp)
      call sal_grad_gfunc(x, y, z, sx, sy, sz, sal_grad_x, sal_grad_y)
      sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + sal_grad_x*ssh
      sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + sal_grad_y*ssh
    enddo
  enddo
  call cpu_clock_end(id_clock_SAL_pp_comp)
end subroutine pp_interaction_compute_nonreprod

!> particle particle interaction for tree code mode, with reproducing sums
subroutine pp_interaction_compute_reprod(sal_ct, G, sal_x, sal_y, e_ssh)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(inout) :: sal_x(:,:), sal_y(:,:)
  real, intent(in) :: e_ssh(:)
  integer :: i, i_s, i_ti, i_tj, j, i_si, i_sj, i_sp, i_t, id, count, ii, ij
  real :: pi, lat, lon, colat, x, y, z, sx, sy, sz, ssh, r2, sal_grad_x, sal_grad_y, sal_val
  real, allocatable :: sshs(:), a_s(:), sal_x_vals(:,:), sal_y_vals(:,:)

  pi = 4.0*DATAN(1.0)
  r2 = G%Rad_Earth_L * G%Rad_Earth_L 

  call cpu_clock_begin(id_clock_SAL_pp_comp)

  do i = 1, size(sal_ct%pp_interactions)
    i_s = sal_ct%pp_interactions(i)%index_source
    i_t = sal_ct%pp_interactions(i)%index_target
    i_ti = sal_ct%one_d_to_2d_i(i_t)
    i_tj = sal_ct%one_d_to_2d_j(i_t)
    x = sal_ct%e_xs(i_t)
    y = sal_ct%e_ys(i_t)
    z = sal_ct%e_zs(i_t)
    allocate(sal_x_vals(1, sal_ct%tree_struct(i_s)%panel_point_count))
    allocate(sal_y_vals(1, sal_ct%tree_struct(i_s)%panel_point_count))
    do j = 1, sal_ct%tree_struct(i_s)%panel_point_count
      i_sp = sal_ct%tree_struct(i_s)%relabeled_points_inside(j)
      sx = sal_ct%e_xs(i_sp)
      sy = sal_ct%e_ys(i_sp)
      sz = sal_ct%e_zs(i_sp)
      ssh = e_ssh(i_sp)
      call sal_grad_gfunc(x, y, z, sx, sy, sz, sal_grad_x, sal_grad_y)
      sal_x_vals(1, j) = sal_grad_x * ssh
      sal_y_vals(1, j) = sal_grad_y * ssh
    enddo
    sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + reproducing_sum(sal_x_vals, only_on_PE=.true.)
    sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + reproducing_sum(sal_y_vals, only_on_PE=.true.)
    deallocate(sal_x_vals)
    deallocate(sal_y_vals)
  enddo
  call cpu_clock_end(id_clock_SAL_pp_comp)
end subroutine pp_interaction_compute_reprod

!> particle particle interaction in FMM mode
subroutine pp_interaction_compute_fmm(sal_ct, G, sal_x, sal_y, e_ssh)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(inout) :: sal_x(:,:), sal_y(:,:)
  real, intent(in) :: e_ssh(:)
  integer :: i, i_s, i_ti, i_tj, j, i_si, i_sj, i_sp, i_t, id, count, ii, ij, i_tp, k, sc
  real :: pi, lat, lon, colat, x, y, z, sx, sy, sz, ssh, r2, sal_grad_x, sal_grad_y, sal_val
  real, allocatable :: sxs(:), sys(:), szs(:), ssshs(:)

  pi = 4.0*DATAN(1.0)
  r2 = G%Rad_Earth_L * G%Rad_Earth_L 

  call cpu_clock_begin(id_clock_SAL_pp_comp)

  do i = 1, size(sal_ct%pp_interactions)
    i_s = sal_ct%pp_interactions(i)%index_source
    i_t = sal_ct%pp_interactions(i)%index_target
    sc = sal_ct%tree_struct(i_s)%panel_point_count
    allocate(sxs(sc))
    allocate(sys(sc))
    allocate(szs(sc))
    allocate(ssshs(sc))
    do j = 1, sc
      i_sp = sal_ct%tree_struct(i_s)%relabeled_points_inside(j)
      sxs(j) = sal_ct%e_xs(i_sp)
      sys(j) = sal_ct%e_ys(i_sp)
      szs(j) = sal_ct%e_zs(i_sp)
      ssshs(j) = e_ssh(i_sp)
    enddo
    do k = 1, sal_ct%tree_struct_targets(i_t)%panel_point_count
      i_tp = sal_ct%tree_struct_targets(i_t)%points_inside(k)
      i_ti = sal_ct%one_d_to_2d_i(i_tp)
      i_tj = sal_ct%one_d_to_2d_j(i_tp)
      x = sal_ct%e_xs(i_tp)
      y = sal_ct%e_ys(i_tp)
      z = sal_ct%e_zs(i_tp)
      do j = 1, sc             
        call sal_grad_gfunc(x, y, z, sxs(j), sys(j), szs(j), sal_grad_x, sal_grad_y)
        sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + sal_grad_x*ssshs(j)
        sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + sal_grad_y*ssshs(j)
      enddo
    enddo
    deallocate(sxs)
    deallocate(sys)
    deallocate(szs)
    deallocate(ssshs)
  enddo
  call cpu_clock_end(id_clock_SAL_pp_comp)
end subroutine pp_interaction_compute_fmm

!> performs the downward pass when run in FMM mode
subroutine fmm_downward_pass(sal_ct, G, sal_x, sal_y, proxy_target_weights_x, &
                              proxy_target_weights_y)
  type(SAL_Conv_CS), intent(in) :: sal_ct
  type(ocean_grid_type), intent(in) :: G
  real, intent(inout) :: sal_x(:,:), sal_y(:,:)
  real, intent(inout) :: proxy_target_weights_x(:,:,:), proxy_target_weights_y(:,:,:)
  real :: pi, min_xi_p, max_xi_p, min_eta_p, max_eta_p, min_xi_c, max_xi_c, min_eta_c, max_eta_c
  real :: xi, eta, tx, ty, tz
  real :: lat, lon, colat, min_xi, max_xi, min_eta, max_eta
  integer :: i, i_c, k, l, m, n, i_t, i_ti, i_tj, j
  real, allocatable :: cheb_xi_c(:), cheb_eta_c(:), basis_vals(:,:) 

  pi = 4.0*DATAN(1.0)
  call cpu_clock_begin(id_clock_SAL_downward_pass)
  ! interpolate from parent panels to leaf panels
  iloop: do i = 1, size(sal_ct%tree_struct_targets)
    min_xi_p = sal_ct%tree_struct_targets(i)%min_xi
    max_xi_p = sal_ct%tree_struct_targets(i)%max_xi
    min_eta_p = sal_ct%tree_struct_targets(i)%min_eta
    max_eta_p = sal_ct%tree_struct_targets(i)%max_eta
    if (.not. sal_ct%tree_struct_targets(i)%is_leaf) then 
      ! interpolate from parent panel to child panel
      do j = 1, sal_ct%tree_struct_targets(i)%child_panel_count
        i_c = sal_ct%tree_struct_targets(i)%child_panels(j)
        min_xi_c = sal_ct%tree_struct_targets(i_c)%min_xi
        max_xi_c = sal_ct%tree_struct_targets(i_c)%max_xi
        min_eta_c = sal_ct%tree_struct_targets(i_c)%min_eta
        max_eta_c = sal_ct%tree_struct_targets(i_c)%max_eta
        call bli_interp_points_shift(cheb_xi_c, min_xi_c, max_xi_c, sal_ct%interp_degree)
        call bli_interp_points_shift(cheb_eta_c, min_eta_c, max_eta_c, sal_ct%interp_degree)
        do k = 1, sal_ct%interp_degree+1 ! child panel eta loop
          eta = cheb_eta_c(k)
          do l = 1, sal_ct%interp_degree+1 ! child panel xi loop
            xi = cheb_xi_c(l)
            call interp_vals_bli(basis_vals, xi, eta, min_xi_p, max_xi_p, min_eta_p, max_eta_p, &
                                  sal_ct%interp_degree)
            do m = 1, sal_ct%interp_degree+1
              do n = 1, sal_ct%interp_degree+1
                proxy_target_weights_x(l,k,i_c) = proxy_target_weights_x(l,k,i_c) + &
                                                  basis_vals(n,m)*proxy_target_weights_x(n,m,i)
                proxy_target_weights_y(l,k,i_c) = proxy_target_weights_y(l,k,i_c) + &
                                                  basis_vals(n,m)*proxy_target_weights_y(n,m,i)
              enddo
            enddo
          enddo
        enddo
      enddo
    else 
    ! interpolate from leaf panel to target points inside
      do j = 1, sal_ct%tree_struct_targets(i)%panel_point_count
        i_t = sal_ct%tree_struct_targets(i)%points_inside(j)
        i_ti = sal_ct%one_d_to_2d_i(i_t)
        i_tj = sal_ct%one_d_to_2d_j(i_t)
        tx = sal_ct%e_xs(i_t)
        ty = sal_ct%e_ys(i_t)
        tz = sal_ct%e_zs(i_t)
        call xieta_from_xyz(tx, ty, tz, xi, eta)
        call interp_vals_bli(basis_vals, xi, eta, min_xi_p, max_xi_p, min_eta_p, max_eta_p, &
                              sal_ct%interp_degree)
        do k = 1, sal_ct%interp_degree+1 ! eta loop
          do l = 1, sal_ct%interp_degree+1 ! xi loop
            sal_x(i_ti, i_tj) = sal_x(i_ti, i_tj) + proxy_target_weights_x(l,k,i)*basis_vals(l,k)
            sal_y(i_ti, i_tj) = sal_y(i_ti, i_tj) + proxy_target_weights_y(l,k,i)*basis_vals(l,k)
          enddo
        enddo
      enddo
    endif
  enddo iloop
  call cpu_clock_end(id_clock_SAL_downward_pass)
end subroutine fmm_downward_pass

!> evaluates the convolution/CSFMM SAL during runtime
subroutine sal_conv_eval(sal_ct, G, eta, sal_x, sal_y)
  type(SAL_Conv_CS), intent(in) :: sal_ct ! conv SAL data struct
  type(ocean_grid_type), intent(in) :: G ! ocean grid
  real, intent(in) :: eta(:,:) ! ssh or bpa
  real, intent(inout) :: sal_x(:,:), sal_y(:,:) ! x,y components of SAL potential gradient
  real, allocatable :: e_ssh(:), proxy_source_weights(:), proxy_target_weights_x(:,:,:)
  real, allocatable :: proxy_target_weights_y(:,:,:)
  integer :: source_size, target_weights, i, id, is, ie, js, je, j

  sal_x(:,:) = 0.0
  sal_y(:,:) = 0.0
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  if (sal_ct%use_sal_conv) then
    call cpu_clock_begin(id_clock_SAL)
    
    id = PE_here()

    if (sal_ct%use_fmm) then ! fmm with upward and downward pass
      call cpu_clock_begin(id_clock_SAL_fmm)
      ! do SSH communication needed for interactions
      allocate(e_ssh(sal_ct%unowned_sources+sal_ct%own_ocean_points), source=0.0)
      call ssh_communications(sal_ct, G, eta, e_ssh)
      if (sal_ct%use_farfield) then
        ! compute proxy source weights
        source_size = (sal_ct%interp_degree+1) * (sal_ct%interp_degree+1) * &
                          size(sal_ct%tree_struct)
        target_weights = (sal_ct%interp_degree+1) * (sal_ct%interp_degree+1)* &
                          size(sal_ct%tree_struct_targets)
        allocate(proxy_source_weights(source_size), source=0.0)
        allocate(proxy_target_weights_x(sal_ct%interp_degree+1, sal_ct%interp_degree+1, &
                    size(sal_ct%tree_struct_targets)), source=0.0)
        allocate(proxy_target_weights_y(sal_ct%interp_degree+1, sal_ct%interp_degree+1, &
                    size(sal_ct%tree_struct_targets)), source=0.0)

        ! upward pass
        call proxy_source_compute_nonrep(sal_ct, G, e_ssh, proxy_source_weights)
        ! perform PC and CC interactions
        call cc_interaction_compute(sal_ct, proxy_source_weights, proxy_target_weights_x, &
                                    proxy_target_weights_y)
        call pc_interaction_compute_fmm(sal_ct, G, proxy_source_weights, sal_x, sal_y)
        
        call cp_interaction_compute(sal_ct, G, eta, e_ssh, proxy_target_weights_x, &
                                    proxy_target_weights_y)
        ! downward pass
        call fmm_downward_pass(sal_ct, G, sal_x, sal_y, proxy_target_weights_x, &
                                proxy_target_weights_y)
      endif
      
      ! compute PP interactions for target domain
      call pp_interaction_compute_fmm(sal_ct, G, sal_x, sal_y, e_ssh)
      call cpu_clock_end(id_clock_SAL_fmm)
    else ! standard tree code
      call cpu_clock_begin(id_clock_SAL_tc)
      ! do SSH communication needed for interactions
      allocate(e_ssh(sal_ct%unowned_sources+sal_ct%own_ocean_points), source=0.0)
      call ssh_communications(sal_ct, G, eta, e_ssh)
      if (sal_ct%use_farfield) then
        ! compute proxy source weights for computational domain
        source_size = (sal_ct%interp_degree+1)*(sal_ct%interp_degree+1)*size(sal_ct%tree_struct)
        allocate(proxy_source_weights(source_size), source=0.0)
        if (sal_ct%reprod_sum) then
          call proxy_source_compute_reprod(sal_ct, G, e_ssh, proxy_source_weights)
          call pc_interaction_compute_reprod(sal_ct, G, proxy_source_weights, sal_x, sal_y)
        else
          call proxy_source_compute_nonrep(sal_ct, G, e_ssh, proxy_source_weights)
          call pc_interaction_compute_nonreprod(sal_ct, G, proxy_source_weights, sal_x, sal_y)
        endif          
      endif

      ! compute PP interactions for target domain
      if (sal_ct%reprod_sum) then
        call pp_interaction_compute_reprod(sal_ct, G, sal_x, sal_y, e_ssh)
      else
        call pp_interaction_compute_nonreprod(sal_ct, G, sal_x, sal_y, e_ssh)
      endif
      
      call cpu_clock_end(id_clock_SAL_tc)
    endif

    call sync_PEs()

    call cpu_clock_end(id_clock_SAL)
  endif
  do j = js, je
    do i = is, ie
      sal_x(i,j) = sal_x(i,j) * sal_ct%scaling_constant
      sal_y(i,j) = sal_y(i,j) * sal_ct%scaling_constant
    enddo
  enddo
  call pass_var(sal_x, G%domain) ! halo update 
  call pass_var(sal_y, G%domain) ! halo update 
end subroutine sal_conv_eval

!> deallocates all the memory used in the convolution SAL control structure
subroutine sal_conv_end(sal_ct)
  type(SAL_Conv_CS), intent(inout) :: sal_ct
  integer :: i
  if (sal_ct%use_sal_conv) then
    if (allocated(sal_ct%e_xs)) deallocate(sal_ct%e_xs)
    if (allocated(sal_ct%e_ys)) deallocate(sal_ct%e_ys)
    if (allocated(sal_ct%e_zs)) deallocate(sal_ct%e_zs)
    if (allocated(sal_ct%pp_interactions)) deallocate(sal_ct%pp_interactions)
    if (allocated(sal_ct%pc_interactions)) deallocate(sal_ct%pc_interactions)
    if (allocated(sal_ct%points_panels)) deallocate(sal_ct%points_panels)

    do i = 1, size(sal_ct%tree_struct)
      if (allocated(sal_ct%tree_struct(i)%child_panels)) &
        deallocate(sal_ct%tree_struct(i)%child_panels)
      if (allocated(sal_ct%tree_struct(i)%points_inside)) &
        deallocate(sal_ct%tree_struct(i)%points_inside)
      if (allocated(sal_ct%tree_struct(i)%relabeled_points_inside)) &
        deallocate(sal_ct%tree_struct(i)%relabeled_points_inside)
    enddo

    if (allocated(sal_ct%tree_struct)) deallocate(sal_ct%tree_struct)
    if (allocated(sal_ct%no_points_to_give_proc)) deallocate(sal_ct%no_points_to_give_proc)
    if (allocated(sal_ct%no_points_to_get_proc)) deallocate(sal_ct%no_points_to_get_proc)
    if (allocated(sal_ct%points_to_give)) deallocate(sal_ct%points_to_give)
    if (allocated(sal_ct%points_to_get)) deallocate(sal_ct%points_to_get)
    if (allocated(sal_ct%unowned_source_points)) deallocate(sal_ct%unowned_source_points)
    if (allocated(sal_ct%indexsg)) deallocate(sal_ct%indexsg)
    if (allocated(sal_ct%indexeg)) deallocate(sal_ct%indexeg)
    if (allocated(sal_ct%pcg)) deallocate(sal_ct%pcg)
    if (allocated(sal_ct%two_d_to_1d)) deallocate(sal_ct%two_d_to_1d)
    if (allocated(sal_ct%one_d_to_2d_i)) deallocate(sal_ct%one_d_to_2d_i)
    if (allocated(sal_ct%one_d_to_2d_j)) deallocate(sal_ct%one_d_to_2d_j)
    if (allocated(sal_ct%point_leaf_panel)) deallocate(sal_ct%point_leaf_panel)

    if (sal_ct%use_fmm) then
      if (allocated(sal_ct%cp_interactions)) deallocate(sal_ct%cp_interactions)
      if (allocated(sal_ct%cc_interactions)) deallocate(sal_ct%cc_interactions)

      do i = 1, size(sal_ct%tree_struct_targets)
        if (allocated(sal_ct%tree_struct_targets(i)%child_panels)) &
          deallocate(sal_ct%tree_struct_targets(i)%child_panels)
        if (allocated(sal_ct%tree_struct_targets(i)%points_inside)) &
          deallocate(sal_ct%tree_struct_targets(i)%points_inside)
        if (allocated(sal_ct%tree_struct_targets(i)%relabeled_points_inside)) &
          deallocate(sal_ct%tree_struct_targets(i)%relabeled_points_inside)
      enddo

      if (allocated(sal_ct%tree_struct_targets)) deallocate(sal_ct%tree_struct_targets)
    endif
  endif
end subroutine sal_conv_end

!> \namespace self_attr_load
!!
!! \section section_conv_SAL Convolution based Self Attraction and Loading
!!
!! This module contains methods to calculate self-attraction and loading (SAL) using a convolution 
!! as a function of sea surface height or bottom pressure anomaly. SAL is primarily used for 
!! fast evolving processes like tides or storm surges, but the effect applies to all motions.
!!
!! If <code>SAL_CONVOLUTION</code> is true, a convolution is used to calculate SAL. 
!! Subroutines in module MOM_conv_self_attr_load are called and the convolution is computed either using a 
!! tree code or the CSFMM, based on <code>CONV_SAL_FMM</code>. The algorithm is based on the CSFMM
!! developed by University of Michigan
!! [\cite Chen2025 and \cite Chen2026].
!! 
!! References
!! 
!! Chen, Anthony, and Robert Krasny. "A Cubed Sphere Fast Multipole Method." 
!! arXiv preprint arXiv:2508.13550 (2025).
!!
!! Chen, Anthony, et al. "Convolution Based Self Attraction and Loading." 
!! arXiv preprint arXiv:2602.01416 (2026).

end module MOM_conv_self_attr_load
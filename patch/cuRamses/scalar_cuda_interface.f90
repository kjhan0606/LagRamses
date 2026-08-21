!==========================================================================
! ISO_C_BINDING interface to scalar_cuda_kernels.cu — GPU Newton-GS sweeps
! for the nGR scalar-field solvers (f(R), symmetron, dilaton, nDGP,
! Galileon). Model ids must match the SCAL_MODEL_* defines in the .cu file.
!==========================================================================
module scalar_cuda_interface
  use iso_c_binding
  implicit none

  integer(c_int), parameter :: SCAL_MODEL_FR        = 0
  integer(c_int), parameter :: SCAL_MODEL_SYMMETRON = 1
  integer(c_int), parameter :: SCAL_MODEL_DILATON   = 2
  integer(c_int), parameter :: SCAL_MODEL_NDGP      = 3
  integer(c_int), parameter :: SCAL_MODEL_GALILEON  = 4

  interface

     subroutine cuda_scal_upload_c(field, rho, ncell, igrid, face, edge, &
          & bnd_slot, bnd_live, bnd_val, ngrid, nbnd, noff) &
          & bind(C, name='cuda_scal_upload')
       import :: c_double, c_int, c_long_long
       real(c_double), dimension(*), intent(in) :: field, rho, bnd_val
       integer(c_long_long), value :: ncell
       integer(c_int), dimension(*), intent(in) :: igrid, face, edge
       integer(c_int), dimension(*), intent(in) :: bnd_slot, bnd_live
       integer(c_int), value :: ngrid, nbnd, noff
     end subroutine cuda_scal_upload_c

     function cuda_scal_is_ready_c() result(ready) &
          & bind(C, name='cuda_scal_is_ready')
       import :: c_int
       integer(c_int) :: ready
     end function cuda_scal_is_ready_c

     subroutine cuda_scal_sweep_c(model, params, ngridmax, ncoarse, &
          & block_size, child_count, tracker, res_max, src_max) &
          & bind(C, name='cuda_scal_sweep')
       import :: c_double, c_int
       integer(c_int), value :: model, ngridmax, ncoarse
       integer(c_int), value :: block_size, child_count, tracker
       real(c_double), dimension(*), intent(in) :: params
       real(c_double), intent(out) :: res_max, src_max
     end subroutine cuda_scal_sweep_c

     subroutine cuda_scal_download_c(field, ncell) &
          & bind(C, name='cuda_scal_download')
       import :: c_double, c_long_long
       real(c_double), dimension(*), intent(out) :: field
       integer(c_long_long), value :: ncell
     end subroutine cuda_scal_download_c

     subroutine cuda_scal_halo_setup_c(emit, n_emit, recv, n_recv) &
          & bind(C, name='cuda_scal_halo_setup')
       import :: c_int
       integer(c_int), dimension(*), intent(in) :: emit, recv
       integer(c_int), value :: n_emit, n_recv
     end subroutine cuda_scal_halo_setup_c

     subroutine cuda_scal_halo_gather_c(buf, n) &
          & bind(C, name='cuda_scal_halo_gather')
       import :: c_double, c_int
       real(c_double), dimension(*), intent(out) :: buf
       integer(c_int), value :: n
     end subroutine cuda_scal_halo_gather_c

     subroutine cuda_scal_halo_scatter_c(buf, n) &
          & bind(C, name='cuda_scal_halo_scatter')
       import :: c_double, c_int
       real(c_double), dimension(*), intent(in) :: buf
       integer(c_int), value :: n
     end subroutine cuda_scal_halo_scatter_c

     subroutine cuda_scal_release_c() bind(C, name='cuda_scal_release')
     end subroutine cuda_scal_release_c

     subroutine cuda_scal_finalize_c() bind(C, name='cuda_scal_finalize')
     end subroutine cuda_scal_finalize_c

  end interface

end module scalar_cuda_interface

!==========================================================================
! Host-side working arrays for the GPU scalar solve (packed grid tables,
! coarse-fine boundary blocks, halo cell lists). Filled per level solve
! by scalar_gpu_begin in force_fine.
!==========================================================================
module scalar_gpu_commons
  use amr_parameters, only: dp
  implicit none

  ! Sentinel matching SCAL_SENTINEL in the .cu (values >= this mean
  ! "use the centre cell value", the zero-gradient fallback)
  real(dp), parameter :: scal_gpu_sentinel = 1.0d300

  logical :: scal_gpu_active = .false.   ! GPU path live for current solve
  integer :: scal_gpu_level  = -1

  ! Packed per-grid tables (row-per-grid layout for the GPU)
  integer, allocatable, dimension(:) :: sgpu_face   ! (6*ngrid)
  integer, allocatable, dimension(:) :: sgpu_edge   ! (12*ngrid)
  integer, allocatable, dimension(:) :: sgpu_bnd_slot ! (ngrid)

  ! Coarse-fine boundary blocks: (noff per cell, 8 cells, nbnd grids)
  integer,  allocatable, dimension(:) :: sgpu_bnd_live
  real(dp), allocatable, dimension(:) :: sgpu_bnd_val
  integer :: sgpu_nbnd = 0
  integer :: sgpu_noff = 6

  ! Halo cell lists (same enumeration as make_virtual_fine_dp packing)
  integer,  allocatable, dimension(:) :: sgpu_emit_cells, sgpu_recv_cells
  real(dp), allocatable, dimension(:) :: sgpu_emit_buf,  sgpu_recv_buf
  integer :: sgpu_n_emit = 0
  integer :: sgpu_n_recv = 0

  ! Offset table (must match c_sc_off in scalar_cuda_kernels.cu):
  ! 1..6 faces (-x,+x,-y,+y,-z,+z), 7..10 xy, 11..14 xz, 15..18 yz
  integer, parameter, dimension(3,18) :: sgpu_off = reshape( (/ &
       & -1,0,0,  1,0,0,  0,-1,0,  0,1,0,  0,0,-1,  0,0,1, &
       &  1,1,0,  1,-1,0, -1,1,0, -1,-1,0, &
       &  1,0,1,  1,0,-1, -1,0,1, -1,0,-1, &
       &  0,1,1,  0,1,-1,  0,-1,1,  0,-1,-1 /), (/3,18/) )

end module scalar_gpu_commons

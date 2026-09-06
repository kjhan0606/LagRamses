module snrt_cuda_multigroup_interface
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  implicit none

  interface
     function snrt_cuda_available() bind(C, name='snrt_cuda_available_c') result(count)
       import :: c_int
       integer(c_int) :: count
     end function snrt_cuda_available

     function snrt_cuda_multigroup_rt_step(state, direction, neighbor, &
          optical_depth, neutral_hydrogen, absorbed, absorbed_group, &
          ncell, ndirection, ngroup, cdt_over_dx) &
          bind(C, name='snrt_cuda_multigroup_rt_step_c') &
          result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: direction(*), optical_depth(*), neutral_hydrogen(*)
       integer(c_int), intent(in) :: neighbor(*)
       real(c_float), intent(out) :: absorbed(*), absorbed_group(*)
       integer(c_int), value :: ncell, ndirection, ngroup
       real(c_float), value :: cdt_over_dx
       integer(c_int) :: ierr
     end function snrt_cuda_multigroup_rt_step

     function snrt_cuda_multigroup_rt_step_owned(state, direction, neighbor, &
          optical_depth, neutral_hydrogen, absorbed, absorbed_group, &
          nowned, nwork, ndirection, ngroup, cdt_over_dx) &
          bind(C, name='snrt_cuda_multigroup_rt_step_owned_c') result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: direction(*), optical_depth(*), neutral_hydrogen(*)
       integer(c_int), intent(in) :: neighbor(*)
       real(c_float), intent(out) :: absorbed(*), absorbed_group(*)
       integer(c_int), value :: nowned, nwork, ndirection, ngroup
       real(c_float), value :: cdt_over_dx
       integer(c_int) :: ierr
     end function snrt_cuda_multigroup_rt_step_owned

     function snrt_cuda_multigroup_rt_step_species(state, direction, neighbor, &
          optical_depth, optical_depth_species, available_species, absorbed, &
          absorbed_group, nowned, nwork, ndirection, ngroup, cdt_over_dx) &
          bind(C, name='snrt_cuda_multigroup_rt_step_species_c') result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: direction(*), optical_depth(*)
       real(c_float), intent(in) :: optical_depth_species(*)
       real(c_float), intent(inout) :: available_species(*)
       integer(c_int), intent(in) :: neighbor(*)
       real(c_float), intent(out) :: absorbed(*), absorbed_group(*)
       integer(c_int), value :: nowned, nwork, ndirection, ngroup
       real(c_float), value :: cdt_over_dx
       integer(c_int) :: ierr
     end function snrt_cuda_multigroup_rt_step_species

     ! DUST-7 ABI.  Species arrays use the contiguous Fortran caller layout
     ! (cell,group,species), seen by CUDA as (species,group,cell).  Dust has
     ! no inventory array: it is a non-depleting absorber.  `absorbed_group`
     ! is the assigned total, not raw removal; `raw_group` and `returned_group`
     ! make the photon closure independently testable.
     function snrt_cuda_multigroup_rt_step_species_dust(state, direction, neighbor, &
          optical_depth, optical_depth_species, optical_depth_dust, available_species, &
          absorbed_hhe_species, absorbed_dust_group, returned_group, raw_group, &
          absorbed_group, absorbed, nowned, nwork, ndirection, ngroup, cdt_over_dx) &
          bind(C, name='snrt_cuda_multigroup_rt_step_species_dust_c') result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: direction(*), optical_depth(*)
       real(c_float), intent(in) :: optical_depth_species(*), optical_depth_dust(*)
       real(c_float), intent(inout) :: available_species(*)
       real(c_float), intent(out) :: absorbed_hhe_species(*)
       real(c_float), intent(out) :: absorbed_dust_group(*), returned_group(*), raw_group(*)
       real(c_float), intent(out) :: absorbed_group(*), absorbed(*)
       integer(c_int), intent(in) :: neighbor(*)
       integer(c_int), value :: nowned, nwork, ndirection, ngroup
       real(c_float), value :: cdt_over_dx
       integer(c_int) :: ierr
     end function snrt_cuda_multigroup_rt_step_species_dust
  end interface

end module snrt_cuda_multigroup_interface

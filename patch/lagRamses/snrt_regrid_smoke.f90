! Exercises production state/IR/regrid objects, not a surrogate remap.
program snrt_regrid_smoke
  use amr_parameters, only: dp,ngridmax,amr_block_size,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: ncoarse,myid,ncpu,active,father,headl,next
  use hydro_commons, only: uold
  use snrt_state
  use snrt_regrid
  use snrt_dust_contract
  use snrt_dust_live
  use mpi_mod
  implicit none
  real(dp) :: p(snrt_checkpoint_cell_width),q(snrt_checkpoint_cell_width),saved(snrt_checkpoint_cell_width)
  real(dp), allocatable :: ir(:),out_ir(:),saved_ir(:)
  real(dp) :: expected_ions,observed_ions
  integer :: ierr,info,j,slot,old_slots
  call MPI_INIT(info)
  ncoarse=1; ngridmax=4; amr_block_size=1; myid=1; ncpu=1
  allocate(uold(1+8*ngridmax,nvar)); uold=0; uold(:,1)=1
  call snrt_dust_contract_load_from_environment(ierr)
  if(ierr/=0.or.snrt_dust_contract_version/=3.or..not.snrt_dust_contract_runtime_allowed)stop 1
  allocate(ir(snrt_dust_contract_number_ir*snrt_ndirection))
  allocate(out_ir(size(ir)),saved_ir(size(ir)))
  p=1d-4; p(1:4)=[1d0,.25d0,.2d0,.3d0]; ir=1d-23
  call snrt_state_restore_cell(1,p,ierr)
  if(ierr/=0)stop 2
  call snrt_dust_live_restore(1,ir,ierr)
  if(ierr/=0)stop 3
  call snrt_state_pack_cell(1,p,ierr)
  call snrt_regrid_refine(1,1,ierr)
  if(ierr/=0)stop 4
  do j=1,8
     call snrt_state_pack_cell(1+j,q,ierr)
     if(ierr/=0.or.any(q/=p))stop 5
     call snrt_dust_live_pack(1+j,out_ir,ierr)
     if(ierr/=0.or.any(out_ir/=ir))stop 6
  end do
  expected_ions=0
  do j=1,8
     uold(1+j,1)=j
     q=p; q(2:4)=[real(j,dp)/16,.1d0,.2d0]; q(5:)=p(5:)*j
     call snrt_state_restore_cell(1+j,q,ierr)
     if(ierr/=0)stop 7
     call snrt_dust_live_restore(1+j,ir*j,ierr)
     if(ierr/=0)stop 8
     expected_ions=expected_ions+j*q(2)
  end do
  call snrt_regrid_coarsen(1,1,ierr)
  if(ierr/=0)stop 9
  call snrt_state_pack_cell(1,q,ierr)
  observed_ions=q(2)*36
  if(abs(observed_ions-expected_ions)>1d-13)stop 10
  if(maxval(abs(q(5:)/p(5:)-4.5d0))>1d-6)stop 11
  if(abs(q(3)-.1d0)>1d-15.or.abs(q(4)-.2d0)>1d-15)stop 12
  call snrt_dust_live_pack(1,out_ir,ierr)
  if(ierr/=0.or.maxval(abs(out_ir/ir-4.5d0))>1d-14)stop 13
  saved=q; saved_ir=out_ir
  do j=1,8
     call snrt_state_pack_cell(1+j,q,ierr)
     if(ierr/=0.or.any(q(2:)/=0))stop 14
     call snrt_dust_live_pack(1+j,out_ir,ierr)
     if(ierr/=0.or.any(out_ir/=0))stop 15
  end do
  ! The recycled grid under an unrelated, uninitialized parent stays empty.
  old_slots=snrt_nslot
  call snrt_regrid_refine(10,1,ierr)
  if(ierr/=0.or.snrt_nslot/=old_slots)stop 16
  call snrt_state_pack_cell(2,q,ierr)
  if(any(q(2:)/=0))stop 17
  ! Invalid child state/density must leave the retained parent untouched.
  call snrt_regrid_refine(1,1,ierr)
  slot=snrt_state_get_slot(2); snrt_intensity(1,1,slot)=-1
  call snrt_regrid_coarsen(1,1,ierr)
  if(ierr==0)stop 18
  call snrt_state_pack_cell(1,q,ierr)
  if(any(q/=saved))stop 19
  call snrt_regrid_refine(1,1,ierr)
  uold(2,1)=-1
  call snrt_regrid_coarsen(1,1,ierr)
  if(ierr==0)stop 20
  call snrt_state_pack_cell(1,q,ierr)
  call snrt_dust_live_pack(1,out_ir,ierr)
  if(any(q/=saved).or.any(out_ir/=saved_ir))stop 21
  ! RAMSES capacity growth must extend, not discard, the cell-to-slot map.
  ngridmax=8
  call snrt_state_restore_cell(65,p,ierr)
  if(ierr/=0.or.snrt_state_get_slot(65)==0)stop 22
  call snrt_state_pack_cell(1,q,ierr)
  if(any(q/=saved))stop 23
  ! Completed level-2 oct (grid 3) restricts to cell 10 in level-1 grid 2.
  ! Unlike deletion, this synchronization must retain every child payload.
  allocate(active(2),father(ngridmax),headl(1,2),next(ngridmax))
  father=0; headl=0; next=0
  active(2)%ngrid=1; allocate(active(2)%igrid(1)); active(2)%igrid(1)=3
  father(3)=10; headl(1,1)=2
  call snrt_state_restore_cell(10,p,ierr)
  call snrt_dust_live_restore(10,ir,ierr)
  call snrt_regrid_refine(10,3,ierr)
  if(ierr/=0)stop 24
  do j=1,8
     q=p; q(2)=real(j,dp)/16; q(5:)=p(5:)*j
     call snrt_state_restore_cell(17+j,q,ierr)
     call snrt_dust_live_restore(17+j,ir*j,ierr)
     uold(17+j,1)=j
  end do
  call snrt_regrid_upload(2,ierr)
  if(ierr/=0)stop 25
  call snrt_state_pack_cell(10,q,ierr)
  if(abs(q(2)*36-expected_ions)>1d-13)stop 26
  call snrt_dust_live_pack(10,out_ir,ierr)
  if(maxval(abs(out_ir/ir-4.5d0))>1d-14)stop 27
  do j=1,8
     call snrt_state_pack_cell(17+j,q,ierr)
     if(abs(q(2)-real(j,dp)/16)>1d-15)stop 28
     call snrt_dust_live_pack(17+j,out_ir,ierr)
     if(maxval(abs(out_ir/ir-j))>1d-14)stop 29
  end do
  write(*,'(A)')'SNRT_NATIVE_LEVEL_UPLOAD_PASS parent_updated=1 children_retained=1'
  write(*,'(A)')'SNRT_NATIVE_REGRID_PASS refine=1 restrict=1 rho_weighted_ions=1 retired_clear=1 reject_before_write=1'
  call MPI_FINALIZE(info)
end program

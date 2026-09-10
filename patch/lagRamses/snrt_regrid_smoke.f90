! Exercises production state/IR/regrid objects, not a surrogate remap.
program snrt_regrid_smoke
  use amr_parameters, only: dp,ngridmax,amr_block_size,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: ncoarse,myid,ncpu,active,father,headl,next,son
  use hydro_commons, only: uold
  use snrt_state
  use snrt_regrid
  use snrt_dust_contract
  use snrt_dust_live
  use snrt_dust_ir, only: dust_ir_diagnostics,dust_ok
  use snrt_runtime_backend, only: snrt_backend_initialize
  use snrt_spectral_contract, only: snrt_spectral_contract_load_from_environment,snrt_group_mean_energy_ev
#ifdef HDF5
  use snrt_hdf5, only: snrt_hdf5_write,snrt_hdf5_read
#endif
  use snrt_amr_topology, only: snrt_face_kind,snrt_face_cell, &
       SNRT_FACE_COARSE_TO_FINE,SNRT_FACE_FINE_TO_COARSE
  use mpi_mod
  implicit none
  real(dp) :: p(snrt_checkpoint_cell_width),q(snrt_checkpoint_cell_width),saved(snrt_checkpoint_cell_width)
  real(dp), allocatable :: ir(:),out_ir(:),saved_ir(:)
  real(dp) :: expected_ions,observed_ions
  integer :: ierr,info,j,slot,old_slots,g
  call MPI_INIT(info)
  call snrt_backend_initialize(ierr)
  if(ierr/=0)stop 54
  ncoarse=1; ngridmax=4; amr_block_size=1; myid=1; ncpu=1
  allocate(uold(1+8*ngridmax,nvar)); uold=0; uold(:,1)=1
  call snrt_dust_contract_load_from_environment(ierr)
  if(ierr/=0.or.snrt_dust_contract_version/=3.or..not.snrt_dust_contract_runtime_allowed)stop 1
  call snrt_spectral_contract_load_from_environment(ierr)
  allocate(ir(snrt_dust_contract_number_ir*snrt_ndirection))
  allocate(out_ir(size(ir)),saved_ir(size(ir)))
  p=0.125_dp; p(1:4)=[1d0,.25d0,.2d0,.3d0]; ir=1d-23
  ! The old standalone number-only fixture needs no spectral environment.
  ! When a contract is supplied, exercise both signs of the correction.
  p(5+snrt_checkpoint_number_width:)=0.0_dp
  if(ierr==0)then
  do g=1,snrt_ngroups
     p(5+snrt_checkpoint_number_width+(g-1)*snrt_ndirection: &
          4+snrt_checkpoint_number_width+g*snrt_ndirection)=(-1.0_dp)**g*snrt_group_mean_energy_ev(g)/32.0_dp
  end do
  end if
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
  if(maxval(abs(q(5:4+snrt_checkpoint_number_width)/p(5:4+snrt_checkpoint_number_width)-4.5d0))>1d-6)stop 11
  if(maxval(abs(q(5+snrt_checkpoint_number_width:)-4.5_dp*p(5+snrt_checkpoint_number_width:)))> &
       1d-6*max(maxval(abs(4.5_dp*p(5+snrt_checkpoint_number_width:))),tiny(1.0_dp)))stop 47
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
  call coarse_ir_checks()
  call defrag_checks()
#ifdef HDF5
  call hdf5_checks()
#endif
  call MPI_FINALIZE(info)
contains
#ifdef HDF5
  subroutine hdf5_checks()
    use ramses_hdf5_io
    use amr_commons, only: numbl,varcpu_restart,varcpu_ngrid_file,varcpu_grid_file_idx
    use amr_parameters, only: nlevelmax,i8b
    use snrt_agn_efficiency, only: snrt_agn_rt_requested
    use snrt_thermochemistry, only: snrt_secondary_tables_load_from_environment
    character(len=1024) :: directory,path
    character(len=32) :: mode
    integer(HID_T) :: grp
    integer :: status,width,old_width,case_id,j,k,first_shift,base
    logical :: with_ir
    real(dp), allocatable :: disk(:),legacy(:),expected(:,:),expected_ir(:,:)
    call get_command_argument(1,directory)
    if(len_trim(directory)==0)return
    call snrt_spectral_contract_load_from_environment(status)
    if(status/=0)stop 47
    call snrt_secondary_tables_load_from_environment(status)
    if(status/=0)stop 55
    call get_command_argument(2,mode)
    with_ir=trim(mode)/='primary'
    if(.not.snrt_agn_rt_requested())stop 48
    path=trim(directory)//'/radiation.h5'
    nlevelmax=1
    allocate(numbl(1,1)); numbl=2
    headl(1,1)=1; next=0; next(1)=2
    width=snrt_checkpoint_cell_width
    if(with_ir)then
       width=width+size(ir)
    else
       snrt_dust_contract_version=2
    end if
    old_width=width-snrt_checkpoint_number_width
    first_shift=5+snrt_checkpoint_number_width
    allocate(expected(snrt_checkpoint_cell_width,16),expected_ir(size(ir),16),disk(16*width),legacy(16*old_width))
    do j=1,16
       q=p; q(5:)=p(5:)*j
       call snrt_state_restore_cell(1+j,q,status)
       if(status/=0)stop 49
       if(with_ir)call snrt_dust_live_restore(1+j,ir*j,status)
       call snrt_state_pack_cell(1+j,expected(:,j),status)
       if(with_ir)call snrt_dust_live_pack(1+j,expected_ir(:,j),status)
    end do
    call hdf5_create_parallel(trim(path),MPI_COMM_WORLD)
    call snrt_hdf5_write()
    call hdf5_close_file()
    ! Read the exact writer payload, then exercise current and old layouts
    ! through both same-rank and file-grid-map restart paths.
    call hdf5_open_parallel(trim(path),MPI_COMM_WORLD)
    call hdf5_open_group('/snrt',grp)
    call hdf5_read_dataset_1d_dp_checked(grp,'level_1',disk,size(disk),0_i8b,int(size(disk),i8b),status)
    if(status/=0)stop 50
    call hdf5_close_group(grp)
    call hdf5_close_file()
    allocate(varcpu_ngrid_file(1),varcpu_grid_file_idx(ngridmax))
    varcpu_ngrid_file=2; varcpu_grid_file_idx=0
    ! Destination grids 3 and 4 are different from the original grids 1/2.
    varcpu_grid_file_idx(3)=1; varcpu_grid_file_idx(4)=2
    do case_id=1,4
       varcpu_restart=mod(case_id,2)==0
       base=1
       if(varcpu_restart)base=3
       headl(1,1)=base; next=0; next(base)=base+1
       if(case_id==3)then
          do j=1,16
             legacy((j-1)*old_width+1:(j-1)*old_width+first_shift-1)= &
                  disk((j-1)*width+1:(j-1)*width+first_shift-1)
             legacy((j-1)*old_width+first_shift:j*old_width)= &
                  disk((j-1)*width+snrt_checkpoint_cell_width+1:j*width)
          end do
          call hdf5_open_rw(path)
          call hdf5_open_group('/snrt',grp)
          call h5adelete_f(grp,'format_version',status)
          call h5adelete_f(grp,'cell_width',status)
          call hdf5_write_attr_int(grp,'format_version',merge(2,1,with_ir))
          call hdf5_write_attr_int(grp,'cell_width',old_width)
          call h5ldelete_f(grp,'level_1',status)
          call hdf5_write_dataset_1d_dp(grp,'level_1',legacy,size(legacy),0_i8b,int(size(legacy),i8b))
          call hdf5_close_group(grp)
          call hdf5_close_file()
       end if
       if(case_id==1.and.len_trim(mode)>0.and.trim(mode)/='primary')then
          ! Corruption is in the last cell: an early valid cell must never
          ! be published during the checkpoint's first validation pass.
          call hdf5_open_rw(path)
          call hdf5_open_group('/snrt',grp)
          select case(trim(mode))
          case('bad-shift')
             disk(15*width+first_shift)=-huge(0.0_dp)
          case('bad-ir')
             disk(15*width+snrt_checkpoint_cell_width+1)=-1.0_dp
          case('bad-width')
             call h5adelete_f(grp,'cell_width',status)
             call hdf5_write_attr_int(grp,'cell_width',old_width)
          case('bad-version')
             call h5adelete_f(grp,'format_version',status)
             call hdf5_write_attr_int(grp,'format_version',2)
          case default
             stop 51
          end select
          call h5ldelete_f(grp,'level_1',status)
          call hdf5_write_dataset_1d_dp(grp,'level_1',disk,size(disk),0_i8b,int(size(disk),i8b))
          call hdf5_close_group(grp)
          call hdf5_close_file()
       end if
       do j=1,16
          k=1+(base-1)*8+j
          call snrt_state_restore_cell(k,p,status)
          if(with_ir)call snrt_dust_live_restore(k,0.0_dp*ir,status)
       end do
       call hdf5_open_parallel(trim(path),MPI_COMM_WORLD)
       call snrt_hdf5_read()
       call hdf5_close_file()
       do j=1,16
          k=1+(base-1)*8+j
          call snrt_state_pack_cell(k,q,status)
          if(case_id>=3)expected(first_shift:,j)=0.0_dp
          if(status/=0.or.any(q/=expected(:,j)))stop 52
          if(with_ir)then
             call snrt_dust_live_pack(k,out_ir,status)
             if(status/=0.or.any(out_ir/=expected_ir(:,j)))stop 53
          end if
       end do
    end do
    write(*,'(A)')'SNRT_NATIVE_HDF5_SHIFT_PASS signed=1 legacy_zero=1 IR_preserved=1 file_grid_map=1'
  end subroutine hdf5_checks

  subroutine hdf5_open_rw(path)
    use ramses_hdf5_io
    character(len=*), intent(in) :: path
    integer(HID_T) :: access
    integer :: status
    call h5open_f(status)
    call h5pcreate_f(H5P_FILE_ACCESS_F,access,status)
    call h5pset_fapl_mpio_f(access,MPI_COMM_WORLD,MPI_INFO_NULL,status)
    call h5fopen_f(trim(path),H5F_ACC_RDWR_F,hdf5_file_id,status,access_prp=access)
    if(status/=0)stop 56
    call h5pclose_f(access,status)
  end subroutine hdf5_open_rw
#endif
  subroutine defrag_checks()
    real(dp) :: before_p(snrt_checkpoint_cell_width,2),before_ir(size(ir),2)
    integer :: old_cells(2),new_cells(2),map(8),status,j
    integer :: bad(snrt_nslot)
    old_cells=[2,18]; new_cells=[10,2]
    do j=1,2
       call snrt_state_pack_cell(old_cells(j),before_p(:,j),status)
       call snrt_dust_live_pack(old_cells(j),before_ir(:,j),status)
    end do
    ! Cyclic overlapping renumbering must move identities, not overwrite
    ! payloads in place. Grid 9 (the old capacity-growth slot) is retired.
    map=0; map(1:3)=[2,3,1]
    call snrt_regrid_defrag(map,status)
    if(status/=0)stop 43
    do j=1,2
       call snrt_state_pack_cell(new_cells(j),q,status)
       call snrt_dust_live_pack(new_cells(j),out_ir,status)
       if(any(q/=before_p(:,j)).or.any(out_ir/=before_ir(:,j)))stop 44
    end do
    bad=2
    call snrt_state_rebind_cells(bad,status)
    if(status==0)stop 45
    call snrt_state_pack_cell(10,q,status)
    if(any(q/=before_p(:,1)))stop 46
    write(*,'(A)')'SNRT_NATIVE_DEFRAG_PASS cyclic_rebind=1 primary_IR_exact=1 duplicate_rejected=1'
  end subroutine

  subroutine coarse_ir_checks()
    type(dust_live_coarse_trial) :: coarse
    type(dust_ir_diagnostics) :: diag
    real(dp), allocatable :: trial(:,:,:)
    real(dp) :: rays(3,snrt_ndirection),weights(snrt_ndirection),material(4),temperature(4)
    real(dp) :: expected,coarse_expected,total
    integer :: cells(4),slots(4),links(6,4),d,j,status
    ! Four fine faces (area 4) meet one coarse face (area 4). All other
    ! faces are fine-owned/blocked so the native transport has no escape.
    allocate(son(65),snrt_face_kind(6,4),snrt_face_cell(6,4))
    son=1; son(2)=0; headl(1,1)=1; next=0
    cells=[18,19,20,21]; links=0
    snrt_face_kind=SNRT_FACE_COARSE_TO_FINE; snrt_face_cell=0
    snrt_face_kind(1,:)=SNRT_FACE_FINE_TO_COARSE; snrt_face_cell(1,:)=2
    rays=0; rays(1,1:snrt_ndirection/2)=1; rays(1,snrt_ndirection/2+1:)=-1; weights=1d0/snrt_ndirection
    call snrt_dust_live_restore(2,2*ir,status)
    do j=1,4
       slots(j)=snrt_state_get_slot(cells(j))
       call snrt_dust_live_restore(cells(j),ir,status)
    end do
    ! dt/dx=2.4 forces three native substeps. The coarse donor must be
    ! refreshed between them, not repeatedly sampled from the old field.
    call snrt_dust_live_stage(2,cells,slots,links,rays,weights,1d0,2.4d0,1d0, &
         [0d0,0d0,0d0,0d0],[0d0,0d0,0d0,0d0],[0d0,0d0,0d0,0d0],[1d0,1d0,1d0,1d0], &
         trial,material,temperature,diag,status,coarse)
    if(status/=dust_ok)then
       write(*,*)'COARSE_IR_STAGE_ERROR',status
       stop 30
    end if
    if(diag%escaped_erg/=0.or.size(coarse%slots)/=1)stop 31
    do d=1,snrt_ndirection
       if(d<=snrt_ndirection/2)then
          coarse_expected=2*.6d0**3; expected=1+1.6d0*(1+.6d0+.6d0**2)
       else
          coarse_expected=2+.4d0*(1+.2d0+.2d0**2); expected=.2d0**3
       end if
       if(maxval(abs(trial(:,d,:)/ir(1)-expected))>2d-14)stop 32
       if(maxval(abs(coarse%energy(:,d,1)/ir(1)-coarse_expected))>2d-14)stop 33
       total=sum(trial(1,d,:))+8*coarse%energy(1,d,1)
       if(abs(total/ir(1)-20)>2d-14)stop 34
    end do
    call snrt_dust_live_pack(2,out_ir,status)
    if(any(out_ir/=2*ir))stop 35
    do j=1,4
       call snrt_dust_live_pack(cells(j),out_ir,status)
       if(any(out_ir/=ir))stop 36
    end do
    call snrt_dust_live_commit(slots,trial,coarse)
    call snrt_dust_live_pack(2,saved_ir,status)
    ! Failure after coarse trial construction must not change either side.
    call snrt_dust_live_stage(2,cells,slots,links,rays,weights,1d0,2.4d0,1d0, &
         [0d0,0d0,0d0,0d0],[-1d-23,0d0,0d0,0d0],[0d0,0d0,0d0,0d0],[1d0,1d0,1d0,1d0], &
         trial,material,temperature,diag,status,coarse)
    if(status==dust_ok)stop 37
    call snrt_dust_live_pack(2,out_ir,status)
    if(any(out_ir/=saved_ir))stop 38
    do j=1,4
       call snrt_dust_live_pack(cells(j),out_ir,status)
       if(any(out_ir/=reshape(trial(:,:,j),[size(ir)])))stop 39
    end do
    ! A face pointing at a covered/unmapped donor must fail, not inject zero.
    snrt_face_cell(1,:)=3
    call snrt_dust_live_stage(2,cells,slots,links,rays,weights,1d0,2.4d0,1d0, &
         [0d0,0d0,0d0,0d0],[0d0,0d0,0d0,0d0],[0d0,0d0,0d0,0d0],[1d0,1d0,1d0,1d0], &
         trial,material,temperature,diag,status,coarse)
    if(status==dust_ok)stop 40
    ! The later coarse-level advance must not apply the fine-owned flux again.
    deallocate(snrt_face_kind,snrt_face_cell)
    allocate(snrt_face_kind(6,1),snrt_face_cell(6,1))
    snrt_face_kind=SNRT_FACE_COARSE_TO_FINE; snrt_face_cell=0
    call snrt_dust_live_stage(1,[2],[snrt_state_get_slot(2)],links(:,1:1),rays,weights,2d0,2.4d0,1d0, &
         [0d0],[0d0],[0d0],[1d0],trial,material(1:1),temperature(1:1),diag,status,coarse)
    if(status/=dust_ok.or.diag%escaped_erg/=0)stop 41
    if(any(reshape(trial,[size(ir)])/=saved_ir))stop 42
    write(*,'(A)')'SNRT_NATIVE_IR_COARSE_FINE_PASS bidirectional=1 substeps=3 conserved=1 atomic=1 no_double_flux=1'
  end subroutine
end program

program test_agn_coarse_state
  use amr_commons
  use pm_commons
  implicit none

  character(len=256)::output_path

  call get_command_argument(1,output_path)
  if(len_trim(output_path)==0) error stop 'usage: test_agn_coarse_state OUTPUT.jsonl'

  myid=1
  ncpu=1
  nsinkmax=2
  nsink=2
  cosmo=.false.
  aexp=0.5d0
  t=1.25d0
  units_density=1.0d-24
  units_time=3.15576d13
  units_length=3.085677581d21
  agn_coarse_dump=.true.
  agn_coarse_dump_file=trim(output_path)
  selfgrav=.true.
  mad_jet=.true.
  X_floor=1.0d-2

  allocate(msink(nsinkmax),dMBHoverdt(nsinkmax),dMEdoverdt(nsinkmax))
  allocate(dMBH_coarse(nsinkmax),dMEd_coarse(nsinkmax),dMsmbh(nsinkmax))
  allocate(jsink(nsinkmax,ndim),Esave(nsinkmax))
  allocate(c_avgptr(nsinkmax),v_avgptr(nsinkmax),d_avgptr(nsinkmax))
  allocate(spinmag(nsinkmax),bhspin(nsinkmax,ndim),eps_sink(nsinkmax))
  allocate(idsink(nsinkmax))

  idsink=(/101,202/)
  msink=(/1.0d-4,2.0d-4/)
  dMBHoverdt=(/2.0d-8,1.0d-10/)
  dMEdoverdt=(/1.0d-8,2.0d-8/)
  dMBH_coarse=(/2.0d-7,1.0d-9/)
  dMEd_coarse=(/1.0d-7,2.0d-7/)
  dMsmbh=(/5.0d-8,6.0d-10/)
  jsink=0d0
  jsink(1,:)=(/1d0,0d0,0d0/)
  jsink(2,:)=(/0d0,1d0,0d0/)
  bhspin=0d0
  bhspin(1,:)=(/0d0,1d0,0d0/)
  bhspin(2,:)=(/0d0,1d0,0d0/)
  spinmag=(/0.7d0,0.3d0/)
  eps_sink=(/0.1d0,0.08d0/)
  Esave=(/3.5d-7,4.5d-7/)
  d_avgptr=(/2.0d0,3.0d0/)
  c_avgptr=(/0.2d0,0.3d0/)
  v_avgptr=(/0.05d0,0.1d0/)

  nstep_coarse_old=6
  nstep_coarse=7
  call dump_agn_coarse_state
  call dump_agn_coarse_state

  nstep_coarse_old=7
  nstep_coarse=8
  t=1.5d0
  call dump_agn_coarse_state

  write(*,'(A)') 'AGN_COARSE_STATE_HARNESS_PASS'
end program test_agn_coarse_state

! Native regression for the actual 68.65 Myr failing photo cell. The fixture
! stores run-length encoded IEEE754 words, preserving its input bits without
! an HDF5 snapshot or a platform-endian binary fixture. Existing data-bank
! and FS2010 environment variables initialize the production implementation.
program chimes_photo_long_interval_test
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_chimes
  use snrt_thermochemistry,only:snrt_secondary_tables_load_from_environment
  implicit none
  integer :: status,j,nd,nr,nsh,count,nvalues,pos,last,offset
  integer(c_int64_t) :: bits
  character(500) :: path,groups(9),group_dir
  character(2) :: num
  type(c_ptr) :: bank,mol
  real(c_double) :: nh,dt,chat,pumping,old(157),new(157),alpha(128,9),shield(2)
  real(c_double) :: ledger(10),gn(9),ge(9),identity(320),ident2(320)
  real(c_double),allocatable :: values(:),n(:,:),e(:,:),nn(:,:),ee(:,:)
  real(c_double) :: initial(11),final(11),q0,q1,t0,t1,elem_error,photon_error,energy_error
  real(c_double) :: cmb_elements(11),cmb_old(157),cmb_new(157),cmb_photons(9),cmb_ctl(9),cold_t,hot_t
  call get_command_argument(1,path)
  open(10,file=trim(path),form='formatted',status='old',action='read')
  read(10,*)nd,nvalues
  if(nd/=80.or.nvalues/=2755)stop 1
  allocate(values(nvalues),n(nd,9),e(nd,9),nn(nd,9),ee(nd,9))
  pos=1
  do while(pos<=nvalues)
    read(10,*)count,path
    read(path,'(Z16)')bits
    last=pos+count-1
    if(count<1.or.last>nvalues)stop 2
    values(pos:last)=transfer(bits,0d0)
    pos=last+1
  enddo
  close(10)
  nh=values(1);dt=values(2);chat=values(3);pumping=values(4)
  old=values(5:161);offset=161
  n=reshape(values(offset+1:offset+9*nd),shape(n));offset=offset+9*nd
  e=reshape(values(offset+1:offset+9*nd),shape(e));offset=offset+9*nd
  alpha=reshape(values(offset+1:offset+1152),shape(alpha));offset=offset+1152
  shield=values(offset+1:offset+2)
  call snrt_secondary_tables_load_from_environment(status)
  if(status/=0)stop 3
  call get_environment_variable('SNRT_CHIMES_MAIN_DATA',path)
  call get_environment_variable('SNRT_CHIMES_GROUP_DIR',group_dir)
  do j=1,9
    write(num,'(I2.2)')j
    groups(j)=trim(group_dir)//'/group_'//num//'.hdf5'//c_null_char
  enddo
  status=chimes_initialize(trim(path)//c_null_char,9,groups,500)
  if(status/=0)stop 4
  ! CMB epoch changes are serial; invalid inputs must not mutate configuration.
  if(chimes_cmb_temperature()/=2.727d0)stop 20
  status=chimes_set_expansion(.01d0)
  if(status/=0.or.abs(chimes_cmb_temperature()-272.7d0)>1d-12)stop 21
  status=chimes_set_expansion(0d0)
  if(status==0.or.abs(chimes_cmb_temperature()-272.7d0)>1d-12)stop 22
  status=chimes_set_expansion(1d0)
  if(status/=0.or.chimes_cmb_temperature()/=2.727d0)stop 23
  ! Exercise the real Compton term, not only the configuration accessor.
  ! Dilute ionized H suppresses collisional losses; CMB above/below100K
  ! must reverse the thermal response from the identical initial state.
  cmb_elements=0;cmb_elements(1)=1
  status=chimes_neutral(cmb_elements,cmb_old)
  if(status/=0)stop 24
  cmb_old(1)=1;cmb_old(2)=0;cmb_old(3)=1
  cmb_ctl=[1d-8,100d0,20d0,1d12,1d18,0d0,1d0,0d0,1d0]
  status=chimes_cell(cmb_ctl,cmb_elements,cmb_old,[(0d0,j=1,9)],cold_t,cmb_new,cmb_photons)
  if(status/=0)stop 25
  status=chimes_set_expansion(.01d0)
  if(status/=0)stop 26
  status=chimes_cell(cmb_ctl,cmb_elements,cmb_old,[(0d0,j=1,9)],hot_t,cmb_new,cmb_photons)
  if(status/=0.or.hot_t<=100d0.or.cold_t>=100d0)stop 27
  print *, 'CMB_NATIVE_THERMAL_RESPONSE cold,hot=',cold_t,hot_t
  status=chimes_budget(cmb_new,final,q1)
  if(status/=0.or.maxval(abs(final-cmb_elements))>1d-8.or.abs(q1)>1d-8)stop 28
  status=chimes_set_expansion(1d0)
  if(status/=0)stop 29
  call get_environment_variable('SNRT_CHIMES_BAND_TABLE',path)
  status=chimes_band_load(trim(path)//c_null_char,bank,nr,nsh,identity)
  if(status/=0)stop 5
  call get_environment_variable('SNRT_CHIMES_MOLECULAR_TABLE',path)
  status=chimes_molecular_load(trim(path)//c_null_char,mol,ident2)
  if(status/=0)stop 6
  status=chimes_budget(old,initial,q0)
  if(status/=0)stop 7
  call cpu_time(t0)
  status=chimes_band_photo_molecular_groups(bank,mol,nd,nh,dt,chat,alpha,shield,pumping,old,n,e,new,nn,ee, &
      ledger,gn,ge)
  call cpu_time(t1)
  print *, 'long interval status,dt_seconds,cpu_seconds=',status,dt,t1-t0
  if(status/=0)stop 8
  if(any(.not.ieee_is_finite(new)).or.any(new<0))stop 9
  if(any(.not.ieee_is_finite(nn)).or.any(nn<0))stop 10
  if(any(.not.ieee_is_finite(ee)).or.any(ee<0))stop 11
  if(any(.not.ieee_is_finite(ledger)).or.any(ledger<0))stop 12
  status=chimes_budget(new,final,q1)
  if(status/=0)stop 13
  elem_error=maxval(abs(final-initial)/max(initial,1d-20))
  photon_error=abs(sum(n)-sum(nn)-ledger(8)-ledger(10))/sum(n)
  energy_error=abs(sum(e)-sum(ee)-ledger(7)-ledger(9))/sum(e)
  if(elem_error>1d-8.or.abs(q1-q0)>1d-8.or.photon_error>1d-7.or.energy_error>1d-7)stop 14
  if(abs(sum(ledger(1:6))-ledger(7))>1d-8*sum(e))stop 15
  if(abs(sum(gn)-ledger(10))>1d-12*sum(n).or.abs(sum(ge)-ledger(9))>1d-12*sum(e))stop 16
  print *, 'nuclear,charge,photon,energy_errors=',elem_error,abs(q1-q0),photon_error,energy_error
  ! No-capture and failure must still obey the original transactional API.
  status=chimes_band_photo_molecular_groups(bank,mol,nd,nh,0d0,chat,alpha,shield,pumping,old,n,e,new,nn,ee, &
      ledger,gn,ge)
  if(status/=0.or.any(new/=old).or.any(nn/=n).or.any(ee/=e).or.any(ledger/=0))stop 17
  new=-999;nn=-999;ee=-999;ledger=-999;gn=-999;ge=-999
  status=chimes_band_photo_molecular_groups(bank,mol,nd,nh,-dt,chat,alpha,shield,pumping,old,n,e,new,nn,ee, &
      ledger,gn,ge)
  if(status==0.or.any(new/=-999).or.any(nn/=-999).or.any(ee/=-999).or.any(ledger/=-999))stop 18
  if(any(gn/=-999).or.any(ge/=-999))stop 19
  call chimes_molecular_free(mol)
  call chimes_band_free(bank)
  print *, 'CHIMES_PHOTO_LONG_INTERVAL_OK'
end program

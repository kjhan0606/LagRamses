program parsec_source_smoke
  use stellar_enrichment_config
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_interpolation, only: source_metallicity_bracket
  use snrt_spectral_contract
  use snrt_parsec_source
  use snrt_stellar_source
  implicit none
  type(stellar_yield_table_t)::table,trial
  character(len=1024)::path,history
  character(len=32)::mode
  integer::ierr,i,k,j
  real(stellar_dp)::q(9),e(9),qa(9),ea(9),qb(9),eb(9),ql(9),el(9),qh(9),eh(9)
  real(stellar_dp),allocatable::identity(:),again(:)
  real(stellar_dp),allocatable::zcheck(:)
  real(stellar_dp)::zl,zh,w
  call get_command_argument(1,path);call get_command_argument(2,history)
  call get_command_argument(3,mode)
  call set_enrichment_defaults()
  enable_agb=.false.;enable_pisn=.true.;configured_imf_mass_max=600
  configured_channel_mass_min([1,3,5])=14;configured_channel_mass_max([1,3,5])=600
  call snrt_spectral_contract_load_from_environment(ierr)
  call require(ierr==0,'spectral contract')
  call stellar_sed_load(ierr)
  if(mode=='reject')then
     call require(ierr/=0.and..not.stellar_sed_enabled,'invalid v4 file admission')
     print *,'PARSEC_SOURCE_FILE_REJECT_PASS'
     stop
  endif
  call require(ierr==0.and.stellar_sed_has_energy,'native v4 source load')
  call stellar_photon_interval(0d0,1d0,.011d0,1d0,q,ierr,e)
  call require(ierr/=0.and.all(q==0),'no unbound source publication')
  call load_yield_table(trim(path),table,ierr)
  call require(ierr==0,'actual feedback table')
  call prepare_high_mass_history(table,trim(history),ierr)
  call require(ierr==0,'actual prepared feedback v4')
  call parsec_sed_bind(table,ierr)
  if(mode=='mismatch')then
     call require(ierr/=0,'five-Z radiation cannot bind two-Z feedback')
     call stellar_photon_interval(0d0,1d0,.025d0,1d0,q,ierr,e)
     call require(ierr/=0.and.all(q==0).and.all(e==0),'mismatched grid cannot publish')
     print *,'PARSEC_ACTUAL_GRID_MISMATCH_PASS'
     stop
  endif
  call require(ierr==0,'common population admission')
  zcheck=pack(table%hm_z,[.true.,table%hm_z(2:)/=table%hm_z(:size(table%hm_z)-1)])
  call source_metallicity_bracket(zcheck,.011d0,zl,zh,w,ierr)
  call require(ierr==0,'reference query bracket')
  call stellar_sed_report()
  call stellar_sed_identity(identity)
  do i=1,4
     default_imf_id=i-1
     if(i==4)default_imf_id=4
     call parsec_sed_bind(table,ierr)
     call require(ierr==0,'four native individual-star IMF weights')
     call stellar_photon_interval(0d0,30d0,.011d0,10000d0,q,ierr,e)
     call require(ierr==0.and.all(q>0).and.all(e>0),'lifetime integrated radiation')
     call require(all(e>=q*snrt_group_edges_ev(:9)).and.all(e<=q*snrt_group_edges_ev(2:)), 'mean support')
     call stellar_photon_interval(0d0,3.4d0,.011d0,10000d0,qa,ierr,ea)
     call require(ierr==0,'early interval')
     call stellar_photon_interval(3.4d0,30d0,.011d0,10000d0,qb,ierr,eb)
     call require(ierr==0,'late interval')
     call require(maxval(abs((qa+qb)/q-1))<2d-13,'Q telescoping')
     call require(maxval(abs((ea+eb)/e-1))<2d-13,'E telescoping')
     call stellar_photon_interval(0d0,30d0,zl,10000d0,ql,ierr,el)
     call require(ierr==0,'lower Z')
     call stellar_photon_interval(0d0,30d0,zh,10000d0,qh,ierr,eh)
     call require(ierr==0,'upper Z')
     call require(maxval(abs(((1-w)*ql+w*qh)/q-1))<2d-13,'same-age Z photons')
     call require(maxval(abs(((1-w)*el+w*eh)/e-1))<2d-13,'same-age Z energy')
     print *,'IMF_QE_SUM',default_imf_id,sum(q),sum(e)
     if(maxval(table%hm_z)>.014d0)then
        do j=1,size(zcheck)-1
           call stellar_photon_interval(0d0,30d0,zcheck(j),10000d0,ql,ierr,el)
           call require(ierr==0,'added lower node')
           call stellar_photon_interval(0d0,30d0,zcheck(j+1),10000d0,qh,ierr,eh)
           call require(ierr==0,'added upper node')
           call stellar_photon_interval(0d0,30d0,.5d0*(zcheck(j)+zcheck(j+1)),10000d0,q,ierr,e)
           call require(ierr==0.and.all(q>0),'added interior radiation')
           call require(maxval(abs(.5d0*(ql+qh)/q-1))<2d-13,'five-Z same-age Q mixture')
           call require(maxval(abs(.5d0*(el+eh)/e-1))<2d-13,'five-Z same-age E mixture')
        enddo
        call stellar_photon_interval(0d0,30d0,nearest(.03d0,1d0),10000d0,q,ierr,e)
        call require(ierr==0.and.all(q==qh).and.all(e==eh),'same material/radiation node roundoff rule')
        call stellar_photon_interval(0d0,30d0,.030001d0,10000d0,q,ierr,e)
        call require(ierr/=0.and.all(q==0),'no upper-Z extrapolation')
        print *,'N_Z_NATIVE_MIXTURE_AND_BOUNDARIES_PASS',size(zcheck),default_imf_id
     endif
  enddo
  default_imf_id=2
  call parsec_sed_bind(table,ierr)
  call require(ierr==0,'restore actual configured IMF')
  call stellar_sed_identity(again)
  call require(size(again)==size(identity),'identity extent')
  call require(all(again==identity),'identity deterministic after weight rebuild')
  do k=1,7
     trial=table
     select case(k)
     case(1)
        trial%hm_fate(1)=5
     case(2)
        trial%hm_age(1)=trial%hm_age(1)*1.01d0
     case(3)
        trial%hm_mass(1)=13.9d0
     case(4)
        trial%high_mass_version=3
     case(5)
        configured_channel_mass_max(5)=599
     case(6)
        trial%high_mass_identity(1)='different_wind_model'
     case(7)
        trial%hm_z(1)=trial%hm_z(1)+1d-6
     end select
     call parsec_sed_bind(trial,ierr)
     call require(ierr/=0,'mismatched actual feedback rejects')
     call stellar_photon_interval(0d0,1d0,.011d0,1d0,q,ierr,e)
     call require(ierr/=0.and.all(q==0).and.all(e==0),'failed rebind cannot publish radiation')
     configured_channel_mass_max(5)=600
  enddo
  call parsec_sed_bind(table,ierr)
  call require(ierr==0,'valid binding restored')
  call stellar_photon_interval(16d0,1d6,.011d0,1d0,q,ierr,e)
  call require(ierr==0.and.all(q==0).and.all(e==0),'all sources dead, exact plateau')
  call stellar_photon_interval(-1d0,0d0,.011d0,1d0,q,ierr,e)
  call require(ierr==0.and.all(q==0),'no pre-birth radiation')
  if(minval(zcheck)<.001d0)then
     call source_metallicity_bracket(zcheck,.001d0,zl,zh,w,ierr)
     call require(ierr==0,'omitted author branch bracket')
     call stellar_photon_interval(0d0,30d0,zl,10000d0,ql,ierr,el)
     call require(ierr==0,'low-Z lower Q/E')
     call stellar_photon_interval(0d0,30d0,zh,10000d0,qh,ierr,eh)
     call require(ierr==0,'low-Z upper Q/E')
     call stellar_photon_interval(0d0,30d0,.001d0,10000d0,q,ierr,e)
     call require(ierr==0,'Z=.001 mixture')
     call require(maxval(abs(((1-w)*ql+w*qh)/q-1))<2d-13,'low-Z Q mixture')
     call require(maxval(abs(((1-w)*el+w*eh)/e-1))<2d-13,'low-Z E mixture')
  endif
  call stellar_photon_interval(0d0,1d0,.5d0*minval(zcheck),1d0,q,ierr,e)
  call require(ierr/=0.and.all(q==0),'no Z extrapolation')
  call stellar_photon_interval(0d0,1d0,.011d0,1d0,q,ierr)
  call require(ierr/=0,'cannot silently drop photon energy')
  print *,'PARSEC_COMMON_NATIVE_PASS IDENTITY_VALUES=',size(identity)
contains
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok)then
       print *,'FAIL ',label
       error stop 1
    endif
  end subroutine
end program

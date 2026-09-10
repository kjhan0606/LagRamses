! Group photon production (and optional v3 energy) per INITIAL stellar mass.
! The supplied table owns its IMF and metallicity; no SED is inferred
! from stellar feedback energy. Integration is exact for the declared linear
! age interpolation, including intervals that cross multiple age nodes.
module snrt_stellar_source
  use amr_parameters, only: dp
  use snrt_spectral_contract, only: snrt_ngroups, snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_group_edges_sha256, snrt_spectral_contract_status, &
       snrt_spectral_contract_approval_id, snrt_spectral_contract_fraction_semantics, &
       snrt_spectral_contract_runtime_allowed,snrt_band_enabled,snrt_group_edges_ev
  use stellar_enrichment_config, only: default_imf_id, population_model_id, &
       configured_imf_mass_min, configured_imf_mass_max, configured_binary_fraction
  use snrt_parsec_source, only: parsec_sed_load,parsec_photon_interval,parsec_sed_identity
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter :: max_age=128,max_z=32
  real(dp),parameter :: myr_s=31557600d6
  logical,save,public :: stellar_sed_enabled=.false.
  logical,save,public :: stellar_sed_has_energy=.false.
  character(len=32),save :: energy_semantics=''
  real(dp),save :: energy_rates(snrt_ngroups,max_age,max_z)
  logical,save :: resolved=.false.
  integer,save :: load_status=0,na=0,nz=0,imf_id=0
  integer,save :: population_id=-1
  integer,save :: version=0
  character(len=64),save :: population_binding='match_feedback',source_sha256=''
  character(len=96),save :: radiation_population=''
  character(len=32),save :: young_age_policy='',spectral_tail_policy=''
  character(len=1024),save :: node_history_file=''
  real(dp),save :: imf_slopes(2)=0d0,imf_break=0d0,escape_fraction=-1d0
  real(dp),save :: imf_min=-1d0,imf_max=-1d0,binary_fraction=-1d0
  real(dp),save :: ages(max_age),metals(max_z),rates(snrt_ngroups,max_age,max_z)
  public :: stellar_sed_load, stellar_photon_interval, stellar_sed_identity
  public :: stellar_sed_consensus, stellar_sed_report
contains
  subroutine stellar_sed_consensus(ierr)
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    integer,intent(out)::ierr
    integer::count,root_count,info,global_error
    real(dp),allocatable::values(:),root_values(:)
    call stellar_sed_identity(values)
    count=size(values);root_count=count;ierr=0
#ifndef WITHOUTMPI
    call MPI_BCAST(root_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
    if(info/=0.or.root_count/=count)ierr=1
    call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(global_error/=0.or.info/=0)then
       ierr=1
       return
    endif
    root_values=values
    call MPI_BCAST(root_values,count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
    if(info/=0.or.any(root_values/=values))ierr=1
    call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    ierr=max(global_error,abs(info))
#endif
  end subroutine

  subroutine stellar_sed_load(ierr)
    integer,intent(out)::ierr
    character(len=1024)::filename
    character(len=64)::transport_sha256,edges_sha256
    character(len=32)::status,interpolation
    character(len=128)::approval_id
    character(len=32)::fraction_semantics
    integer::ios,unit,length,i
    namelist/snrt_stellar_sed/version,na,nz,imf_id,population_id,imf_min,imf_max,binary_fraction,ages,metals,rates, &
         transport_sha256,edges_sha256,status,interpolation,approval_id,fraction_semantics, &
         population_binding,source_sha256,radiation_population,young_age_policy,spectral_tail_policy, &
         imf_slopes,imf_break,escape_fraction,energy_rates,energy_semantics,node_history_file
    if(resolved)then
       ierr=load_status
       return
    endif
    resolved=.true.
    call get_environment_variable('SNRT_STELLAR_SED',filename,length=length,status=ios)
    ierr=0
    if(ios==1.or.length==0)return
    load_status=1
    if(ios/=0)goto 900
    ages=0d0;metals=0d0;rates=-1d0
    energy_rates=-1d0;energy_semantics=''
    version=0;status='';interpolation='';transport_sha256='';edges_sha256=''
    approval_id='';fraction_semantics=''
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if(ios/=0)goto 900
    read(unit,nml=snrt_stellar_sed,iostat=ios)
    close(unit)
    if(ios/=0.or.version<1.or.version>4)goto 900
    if(version<4)then
       if(na<2.or.na>max_age.or.nz<2.or.nz>max_z)goto 900
       if(trim(interpolation)/='linear_age_linear_Z'.or.node_history_file/='')goto 900
    endif
    select case(version)
    case(1)
       if(population_binding/='match_feedback')goto 900
       ! v1 cannot hide unbound v2 metadata in a legacy restart identity.
       if(source_sha256/=''.or.radiation_population/=''.or.young_age_policy/=''.or.spectral_tail_policy/='')goto 900
       if(any(imf_slopes/=0d0).or.imf_break/=0d0.or.escape_fraction/=-1d0)goto 900
       if(imf_id/=default_imf_id.or.population_id/=population_model_id.or. &
            imf_min/=configured_imf_mass_min.or.imf_max/=configured_imf_mass_max.or. &
            binary_fraction/=configured_binary_fraction)goto 900
    case(2,3)
       ! An explicit comparison model, NOT a new globally selected feedback IMF.
       if(population_binding/='independent_radiation_reference'.or.status/='reference_control')goto 900
       if(radiation_population/='BPASS_v2.2.1_bin-imf135_300')goto 900
       if(imf_id/=-1.or.population_id/=-1.or.binary_fraction/=-1d0)goto 900
       if(imf_min/=0.1d0.or.imf_max/=300d0.or.imf_break/=0.5d0)goto 900
       if(any(imf_slopes/=[-1.30d0,-2.35d0]))goto 900
       if(.not.ieee_is_finite(escape_fraction))goto 900
       if(escape_fraction<0d0.or.escape_fraction>1d0.or.fraction_semantics/='escaped')goto 900
       if(len_trim(source_sha256)/=64)goto 900
       do i=1,64
          if(index('0123456789abcdef',source_sha256(i:i))==0)goto 900
       enddo
       if(young_age_policy/='hold_first_to_zero'.or.spectral_tail_policy/='zero_outside_source_domain')goto 900
       if(ages(1)/=0d0.or.ages(2)/=1d0.or.any(rates(:,1,1:nz)/=rates(:,2,1:nz)))goto 900
    case(4)
       if(population_binding/='match_feedback_high_mass_only'.or.status/='reference_control')goto 900
       if(radiation_population/='PARSEC_v2_nonrot_high_mass_only')goto 900
       if(imf_id/=default_imf_id.or.population_id/=population_model_id.or. &
            imf_min/=configured_imf_mass_min.or.imf_max/=configured_imf_mass_max.or. &
            binary_fraction/=configured_binary_fraction)goto 900
       if(na/=0.or.nz/=0.or.any(rates/=-1d0).or.any(energy_rates/=-1d0))goto 900
       if(any(imf_slopes/=0d0).or.imf_break/=0d0)goto 900
       if(.not.ieee_is_finite(escape_fraction).or.escape_fraction<0.or.escape_fraction>1)goto 900
       if(fraction_semantics/='escaped'.or..not.snrt_band_enabled())goto 900
       if(energy_semantics/='photon_number_and_energy_v1')goto 900
       if(interpolation/='linear_cumulative_age_linear_Z'.or.node_history_file=='')goto 900
       if(young_age_policy/='hold_first_to_zero'.or.spectral_tail_policy/='Q5_Planck_and_track_tail_v1')goto 900
       if(len_trim(source_sha256)/=64)goto 900
       do i=1,64
          if(index('0123456789abcdef',source_sha256(i:i))==0)goto 900
       enddo
    end select
    if(.not.snrt_spectral_contract_runtime_allowed)goto 900
    if(trim(snrt_spectral_contract_status)/=trim(status))goto 900
    if(fraction_semantics/=snrt_spectral_contract_fraction_semantics)goto 900
    select case(trim(status))
    case('reference_control')
       if(len_trim(approval_id)/=0)goto 900
    case('approved_production')
       ! Approval must cover the common stellar/AGN transport closure and
       ! the selected population table, not merely an upstream SED filename.
       if(len_trim(approval_id)==0.or.approval_id/=snrt_spectral_contract_approval_id)goto 900
    case default
       goto 900
    end select
    if(transport_sha256/=snrt_spectral_contract_source_sha256.or. &
         edges_sha256/=snrt_spectral_contract_group_edges_sha256)goto 900
    if(version==4)then
       call parsec_sed_load(trim(node_history_file),ios)
       if(ios/=0)goto 900
       load_status=0;stellar_sed_enabled=.true.;stellar_sed_has_energy=.true.
       ierr=0;return
    endif
    if(any(.not.ieee_is_finite(ages(1:na))).or.any(ages(1:na)<0d0))goto 900
    if(any(ages(2:na)<=ages(1:na-1)))goto 900
    if(any(.not.ieee_is_finite(metals(1:nz))).or.any(metals(1:nz)<0d0).or.any(metals(1:nz)>1d0))goto 900
    if(any(metals(2:nz)<=metals(1:nz-1)))goto 900
    if(any(.not.ieee_is_finite(rates(:,1:na,1:nz))).or.any(rates(:,1:na,1:nz)<0d0))goto 900
    if(version==3)then
       if(.not.snrt_band_enabled().or.energy_semantics/='photon_number_and_energy_v1')goto 900
       if(any(.not.ieee_is_finite(energy_rates(:,1:na,1:nz))).or.any(energy_rates(:,1:na,1:nz)<0d0))goto 900
       if(any(energy_rates(:,1,1:nz)/=energy_rates(:,2,1:nz)))goto 900
       do i=1,snrt_ngroups
          if(any(energy_rates(i,1:na,1:nz)<rates(i,1:na,1:nz)*snrt_group_edges_ev(i)*(1d0-1d-12)).or. &
               any(energy_rates(i,1:na,1:nz)>rates(i,1:na,1:nz)*snrt_group_edges_ev(i+1)*(1d0+1d-12)))goto 900
       enddo
    else
       ! Legacy identity must not silently discard unbound energy fields.
       if(energy_semantics/=''.or.any(energy_rates/=-1d0))goto 900
    endif
    load_status=0
    stellar_sed_enabled=.true.
    stellar_sed_has_energy=version==3
900 ierr=load_status
  end subroutine

  subroutine stellar_photon_interval(age0,age1,metallicity,mass_msun,photons,ierr,energy_ev)
    real(dp),intent(in)::age0,age1,metallicity,mass_msun ! Myr, mass fraction, initial Msun
    real(dp),intent(out)::photons(snrt_ngroups)
    real(dp),optional,intent(out)::energy_ev(snrt_ngroups)
    integer,intent(out)::ierr
    integer::ia,iz
    real(dp)::lo,hi,wz,wa,wb,q0(snrt_ngroups),q1(snrt_ngroups)
    photons=0d0;ierr=1
    if(present(energy_ev))energy_ev=0d0
    if(.not.stellar_sed_enabled)return
    if(stellar_sed_has_energy.neqv.present(energy_ev))return
    if(version==4)then
       call parsec_photon_interval(age0,age1,metallicity,mass_msun,photons,energy_ev,ierr)
       return
    endif
    if(any(.not.ieee_is_finite([age0,age1,metallicity,mass_msun])))return
    if(age1<age0.or.mass_msun<0d0)return
    ! Before formation there is no star. Never extrapolate a tabulated SED.
    if(age1<=0d0)then
       ierr=0
       return
    endif
    if(max(0d0,age0)<ages(1).or.age1>ages(na))return
    if(metallicity<metals(1).or.metallicity>metals(nz))return
    iz=1
    do while(iz<nz-1)
       if(metallicity<=metals(iz+1))exit
       iz=iz+1
    enddo
    ! Preserve exact knots even under reciprocal-optimized production builds:
    ! a tiny contaminating weight can matter for extremely steep hard-X SEDs.
    if(metallicity==metals(iz))then
       wz=0d0
    else if(metallicity==metals(iz+1))then
       wz=1d0
    else
       wz=(metallicity-metals(iz))/(metals(iz+1)-metals(iz))
    endif
    do ia=1,na-1
       lo=max(0d0,age0,ages(ia));hi=min(age1,ages(ia+1))
       if(hi<=lo)cycle
       if(wz==0d0)then
          q0=rates(:,ia,iz);q1=rates(:,ia+1,iz)
       else if(wz==1d0)then
          q0=rates(:,ia,iz+1);q1=rates(:,ia+1,iz+1)
       else
          q0=(1d0-wz)*rates(:,ia,iz)+wz*rates(:,ia,iz+1)
          q1=(1d0-wz)*rates(:,ia+1,iz)+wz*rates(:,ia+1,iz+1)
       endif
       wa=(lo-ages(ia))/(ages(ia+1)-ages(ia))
       wb=(hi-ages(ia))/(ages(ia+1)-ages(ia))
       photons=photons+(hi-lo)*((1d0-0.5d0*(wa+wb))*q0+0.5d0*(wa+wb)*q1)
       if(stellar_sed_has_energy)then
          if(wz==0d0)then
             q0=energy_rates(:,ia,iz);q1=energy_rates(:,ia+1,iz)
          else if(wz==1d0)then
             q0=energy_rates(:,ia,iz+1);q1=energy_rates(:,ia+1,iz+1)
          else
             q0=(1d0-wz)*energy_rates(:,ia,iz)+wz*energy_rates(:,ia,iz+1)
             q1=(1d0-wz)*energy_rates(:,ia+1,iz)+wz*energy_rates(:,ia+1,iz+1)
          endif
          energy_ev=energy_ev+(hi-lo)*((1d0-0.5d0*(wa+wb))*q0+0.5d0*(wa+wb)*q1)
       endif
    enddo
    photons=photons*myr_s*mass_msun
    if(any(.not.ieee_is_finite(photons)).or.any(photons<0d0))return
    if(stellar_sed_has_energy)then
       energy_ev=energy_ev*myr_s*mass_msun
       if(any(.not.ieee_is_finite(energy_ev)).or.any(energy_ev<0d0))return
    endif
    ierr=0
  end subroutine

  subroutine stellar_sed_identity(values)
    real(dp),allocatable,intent(out)::values(:)
    real(dp),allocatable::node_values(:)
    character(len=288)::metadata
    integer::i
    values=[real(na,dp),real(nz,dp),real(imf_id,dp),real(population_id,dp),imf_min,imf_max,binary_fraction, &
         ages(1:na),metals(1:nz), &
         reshape(rates(:,1:na,1:nz),[snrt_ngroups*na*nz])]
    if(version>=2)then
       ! Bind source, population and approximations to MPI/restart, even if
       ! two tables happen to have identical rates. Leave v1 bytes unchanged.
       metadata=population_binding//source_sha256//radiation_population//young_age_policy//spectral_tail_policy
       values=[values,real(version,dp),imf_slopes,imf_break,escape_fraction, &
            (real(iachar(metadata(i:i)),dp),i=1,len(metadata))]
    endif
    if(version==3)values=[values,reshape(energy_rates(:,1:na,1:nz),[snrt_ngroups*na*nz]), &
         (real(iachar(energy_semantics(i:i)),dp),i=1,len(energy_semantics))]
    if(version==4)then
       call parsec_sed_identity(node_values)
       values=[values,(real(iachar(energy_semantics(i:i)),dp),i=1,len(energy_semantics)),node_values]
    endif
  end subroutine

  subroutine stellar_sed_report()
    if(.not.stellar_sed_enabled)return
    write(*,'(A,A)')'SNRT stellar population binding: ',trim(population_binding)
    if(version<2)return
    if(version==4)then
       write(*,'(A)')'SNRT PARSEC HIGH-MASS ONLY: matched 14--600 Msun; full IMF denominator .08--600; not full SSP'
       write(*,'(A)')'SNRT PARSEC: measured Q5 + within-band Planck prior; missing photon tails use full-track Planck'
       write(*,'(A)')'SNRT PARSEC: native cumulative Q/E; no post-terminal radiation; no below-14 Msun source supplied'
       write(*,'(A,A)')'SNRT PARSEC source SHA256: ',source_sha256
       write(*,'(A,ES16.8)')'SNRT PARSEC escaped fraction already applied: ',escape_fraction
       return
    endif
    write(*,'(A,A)')'SNRT radiation population (feedback unchanged): ',trim(radiation_population)
    write(*,'(A,A)')'SNRT stellar source SHA256: ',source_sha256
    write(*,'(A,A,A,A)')'SNRT stellar approximations: ',trim(young_age_policy),'; ',trim(spectral_tail_policy)
    write(*,'(A,ES16.8)')'SNRT stellar escaped fraction already applied to rates: ',escape_fraction
    if(stellar_sed_has_energy)then
       write(*,'(A)')'SNRT stellar comparison: actual BPASS Q/E injection; unchanged state Eref; selected spectral receiver'
    else if(snrt_band_enabled())then
       write(*,'(A)')'SNRT stellar comparison: common nominal injection energies; selected spectral receiver'
    else
       write(*,'(A)')'SNRT stellar comparison only: common grey AGN/stellar transport closure retained'
    endif
  end subroutine
end module

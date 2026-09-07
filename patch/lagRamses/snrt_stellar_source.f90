! Group photon production per INITIAL stellar mass. The supplied table owns
! its IMF, metallicity and common transport-closure choices; no SED is inferred
! from stellar feedback energy. Integration is exact for the declared linear
! age interpolation, including intervals that cross multiple age nodes.
module snrt_stellar_source
  use amr_parameters, only: dp
  use snrt_spectral_contract, only: snrt_ngroups, snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_group_edges_sha256, snrt_spectral_contract_status, &
       snrt_spectral_contract_approval_id, snrt_spectral_contract_fraction_semantics, &
       snrt_spectral_contract_runtime_allowed
  use stellar_enrichment_config, only: default_imf_id, population_model_id, &
       configured_imf_mass_min, configured_imf_mass_max, configured_binary_fraction
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter :: max_age=128,max_z=32
  real(dp),parameter :: myr_s=31557600d6
  logical,save,public :: stellar_sed_enabled=.false.
  logical,save :: resolved=.false.
  integer,save :: load_status=0,na=0,nz=0,imf_id=0
  integer,save :: population_id=-1
  real(dp),save :: imf_min=-1d0,imf_max=-1d0,binary_fraction=-1d0
  real(dp),save :: ages(max_age),metals(max_z),rates(snrt_ngroups,max_age,max_z)
  public :: stellar_sed_load, stellar_photon_interval, stellar_sed_identity
  public :: stellar_sed_consensus
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
    integer::ios,unit,length,version
    namelist/snrt_stellar_sed/version,na,nz,imf_id,population_id,imf_min,imf_max,binary_fraction,ages,metals,rates, &
         transport_sha256,edges_sha256,status,interpolation,approval_id,fraction_semantics
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
    version=0;status='';interpolation='';transport_sha256='';edges_sha256=''
    approval_id='';fraction_semantics=''
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if(ios/=0)goto 900
    read(unit,nml=snrt_stellar_sed,iostat=ios)
    close(unit)
    if(ios/=0.or.version/=1)goto 900
    if(na<2.or.na>max_age.or.nz<2.or.nz>max_z)goto 900
    if(imf_id/=default_imf_id.or.trim(interpolation)/='linear_age_linear_Z')goto 900
    if(population_id/=population_model_id.or.imf_min/=configured_imf_mass_min.or. &
         imf_max/=configured_imf_mass_max.or.binary_fraction/=configured_binary_fraction)goto 900
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
    if(any(.not.ieee_is_finite(ages(1:na))).or.any(ages(1:na)<0d0))goto 900
    if(any(ages(2:na)<=ages(1:na-1)))goto 900
    if(any(.not.ieee_is_finite(metals(1:nz))).or.any(metals(1:nz)<0d0).or.any(metals(1:nz)>1d0))goto 900
    if(any(metals(2:nz)<=metals(1:nz-1)))goto 900
    if(any(.not.ieee_is_finite(rates(:,1:na,1:nz))).or.any(rates(:,1:na,1:nz)<0d0))goto 900
    load_status=0
    stellar_sed_enabled=.true.
900 ierr=load_status
  end subroutine

  subroutine stellar_photon_interval(age0,age1,metallicity,mass_msun,photons,ierr)
    real(dp),intent(in)::age0,age1,metallicity,mass_msun ! Myr, mass fraction, initial Msun
    real(dp),intent(out)::photons(snrt_ngroups)
    integer,intent(out)::ierr
    integer::ia,iz
    real(dp)::lo,hi,wz,wa,wb,q0(snrt_ngroups),q1(snrt_ngroups)
    photons=0d0;ierr=1
    if(.not.stellar_sed_enabled)return
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
    wz=(metallicity-metals(iz))/(metals(iz+1)-metals(iz))
    do ia=1,na-1
       lo=max(0d0,age0,ages(ia));hi=min(age1,ages(ia+1))
       if(hi<=lo)cycle
       q0=(1d0-wz)*rates(:,ia,iz)+wz*rates(:,ia,iz+1)
       q1=(1d0-wz)*rates(:,ia+1,iz)+wz*rates(:,ia+1,iz+1)
       wa=(lo-ages(ia))/(ages(ia+1)-ages(ia))
       wb=(hi-ages(ia))/(ages(ia+1)-ages(ia))
       photons=photons+(hi-lo)*((1d0-0.5d0*(wa+wb))*q0+0.5d0*(wa+wb)*q1)
    enddo
    photons=photons*myr_s*mass_msun
    if(any(.not.ieee_is_finite(photons)).or.any(photons<0d0))return
    ierr=0
  end subroutine

  subroutine stellar_sed_identity(values)
    real(dp),allocatable,intent(out)::values(:)
    values=[real(na,dp),real(nz,dp),real(imf_id,dp),real(population_id,dp),imf_min,imf_max,binary_fraction, &
         ages(1:na),metals(1:nz), &
         reshape(rates(:,1:na,1:nz),[snrt_ngroups*na*nz])]
  end subroutine
end module

! Immutable per-star radiation histories; no spectral reconstruction in RAMSES.
! HIGH-MASS ONLY. Cumulative Q/E carry the explicit offline atmosphere closure.
module snrt_parsec_source
  use stellar_enrichment_config
  use stellar_enrichment_contract, only: stellar_population_t
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_ssp_sources, only: build_source_mass_edges,calculate_imf_mass_fraction
  use stellar_yield_interpolation, only: source_metallicity_bracket,interpolation_ok
  use snrt_spectral_contract, only: snrt_group_edges_ev
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter :: dp=stellar_dp
  integer,save :: nn=0,nr=0
  logical,save,public :: parsec_sed_enabled=.false.
  logical,save :: bound=.false.
  real(dp),allocatable,save :: mass(:),metal(:),death(:),age(:),release(:,:),weight(:)
  real(dp),allocatable,save :: z_nodes(:)
  integer,allocatable,save :: first(:),count(:),fate(:)
  character(len=64),save :: model=''
  public :: parsec_sed_load,parsec_sed_bind,parsec_photon_interval,parsec_sed_identity
contains
  subroutine parsec_sed_load(path,ierr)
    character(len=*),intent(in)::path
    integer,intent(out)::ierr
    integer::u,ios,i,j,k
    character(len=64)::magic
    real(dp)::delta(21)
    ierr=1
    if(allocated(mass))return ! Load once, not a mutable runtime source.
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(u,'(A)',iostat=ios)magic
    if(ios/=0.or.magic/='SNRT_PARSEC_QE_INTERVAL_V1')goto 900
    read(u,*,iostat=ios)nn,nr
    if(ios/=0.or.nn<2.or.nn>512.or.nr<2*nn.or.nr>500000)goto 900
    read(u,'(A)',iostat=ios)model
    if(ios/=0)goto 900
    select case(model)
    case('parsec2025_w17_hw02_composite_v1','parsec2025_w17_hw02_phase_escape_f22_v1', &
         'parsec2025_w17_hw02_precision_rate_v1')
    case default
       goto 900
    end select
    allocate(mass(nn),metal(nn),death(nn),weight(nn),first(nn),count(nn),fate(nn),age(nr),release(21,nr))
    weight=0;k=0
    do i=1,nn
       read(u,*,iostat=ios)mass(i),metal(i),death(i),fate(i),count(i)
       if(ios/=0)goto 900
       if(any(.not.ieee_is_finite([mass(i),metal(i),death(i)])))goto 900
       if(mass(i)<14.or.mass(i)>600.or.death(i)<=0.or.fate(i)<1.or.fate(i)>5)goto 900
       if(metal(i)<0.or.metal(i)>1)goto 900
       if(count(i)<2.or.count(i)>nr-k)goto 900
       if(i>1)then
          if(metal(i)<metal(i-1))goto 900
          if(metal(i)==metal(i-1).and.mass(i)<=mass(i-1))goto 900
       endif
       first(i)=k+1
       do j=1,count(i)
          k=k+1
          read(u,*,iostat=ios)age(k),release(:,k)
          if(ios/=0)goto 900
          if(.not.ieee_is_finite(age(k)).or.any(.not.ieee_is_finite(release(:,k))))goto 900
          if(age(k)<0.or.any(release(:,k)<0))goto 900
          if(j==1)then
             if(age(k)/=0.or.any(release(:,k)/=0))goto 900
          else
             if(age(k)<=age(k-1))goto 900
             ! On disk each row is the positive integral of the preceding
             ! interval. Never subtract two rounded release endpoints.
             delta=release(:,k)
             if(any(delta<0))goto 900
             if(any(delta(10:18)<delta(1:9)*snrt_group_edges_ev(:9)*(1d0-1d-10)).or. &
                  any(delta(10:18)>delta(1:9)*snrt_group_edges_ev(2:)*(1d0+1d-10)))goto 900
          endif
          if(abs(sum(release(10:18,k))+sum(release(20:21,k))-release(19,k))> &
               1d-10*max(tiny(1d0),release(19,k)))goto 900
       enddo
       if(age(k)/=death(i))goto 900
    enddo
    if(k/=nr)goto 900
    z_nodes=pack(metal,[.true.,metal(2:)/=metal(:nn-1)])
    if(size(z_nodes)<2)goto 900
    ! Only whitespace may follow the final declared source record.
    do
       read(u,'(A)',iostat=ios)magic
       if(ios<0)exit
       if(ios>0.or.len_trim(magic)>0)goto 900
    enddo
    parsec_sed_enabled=.true.;ierr=0
900 close(u)
    if(ierr/=0)then
       if(allocated(mass))deallocate(mass,metal,death,weight,first,count,fate,age,release)
       if(allocated(z_nodes))deallocate(z_nodes)
       nn=0;nr=0;model=''
    endif
  end subroutine

  subroutine parsec_sed_bind(table,ierr)
    type(stellar_yield_table_t),intent(in)::table
    integer,intent(out)::ierr
    type(stellar_population_t)::population
    real(dp),allocatable::edges(:),trial(:)
    real(dp)::l,r,m,z,best,d,fraction,total
    integer::i,b,node,status,iz
    ierr=0
    if(.not.parsec_sed_enabled)return
    ierr=1;bound=.false.
    if(.not.table%high_mass_ready.or.table%high_mass_version/=4)return
    if(.not.allocated(table%hm_fate).or..not.table%high_mass_linear_z)return
    if(size(table%hm_mass)/=nn)return
    if(trim(table%high_mass_identity(1))/=trim(model))return
    if(any(table%hm_mass/=mass).or.any(table%hm_z/=metal).or.any(table%hm_fate/=fate))return
    if(any(abs(table%hm_age*1000-death)>1d-12*death))return
    if(stellar_feedback_mode/='channel_resolved'.or.high_mass_model/='source_consistent')return
    if(population_model_id/=0.or.configured_binary_fraction/=0)return
    if(configured_imf_mass_min/=.08d0.or.configured_imf_mass_max/=600d0)return
    select case(default_imf_id)
    case(0,1,2,4)
    case default
       return
    end select
    if(.not.enable_wind.or..not.enable_snii.or..not.enable_pisn)return
    if(any(configured_channel_mass_min([1,3,5])/=14d0).or. &
         any(configured_channel_mass_max([1,3,5])/=600d0))return
    population%imf_id=default_imf_id
    population%imf_mass_min=configured_imf_mass_min;population%imf_mass_max=configured_imf_mass_max
    call build_source_mass_edges(table,population,64,edges,status)
    if(status/=0)return
    allocate(trial(nn));trial=0
    do iz=1,size(z_nodes)
       z=z_nodes(iz)
       do b=1,size(edges)-1
          l=max(14d0,edges(b));r=min(600d0,edges(b+1))
          if(r<=l)cycle
          m=sqrt(l*r);best=huge(1d0);node=0
          do i=1,nn
             if(metal(i)/=z.or.((m<40d0).neqv.(mass(i)<40d0)))cycle
             d=abs(m-mass(i))
             if(d<best)then
                node=i;best=d
             endif
          enddo
          if(node==0)return
          call calculate_imf_mass_fraction(default_imf_id,.08d0,600d0,l,r,fraction,status)
          if(status/=0)return
          trial(node)=trial(node)+fraction/mass(node)
       enddo
       call calculate_imf_mass_fraction(default_imf_id,.08d0,600d0,14d0,600d0,total,status)
       if(status/=0.or.abs(sum(trial*mass,mask=metal==z)-total)>1d-12)return
    enddo
    weight=trial;bound=.true.;ierr=0
  end subroutine

  subroutine node_interval(node,t0,t1,result)
    integer,intent(in)::node
    real(dp),intent(in)::t0,t1
    real(dp),intent(out)::result(18)
    integer::l,r,m,k,last
    real(dp)::w,lo,hi
    result=0
    if(t1<=0.or.t0>=death(node).or.t1<=t0)return
    l=first(node);r=l+count(node)-1
    last=r
    do while(r-l>1)
       m=(l+r)/2
       if(age(m)<=max(0d0,t0))then
          l=m
       else
          r=m
       endif
    enddo
    do k=l+1,last
       lo=max(0d0,t0,age(k-1));hi=min(t1,age(k))
       if(hi>lo)then
          w=(hi-lo)/(age(k)-age(k-1))
          result=result+w*release(:18,k)
       endif
       if(age(k)>=t1)exit
    enddo
  end subroutine

  subroutine parsec_photon_interval(t0,t1,z,initial_mass,photons,energy,ierr)
    real(dp),intent(in)::t0,t1,z,initial_mass ! Myr, absolute metal fraction, INITIAL Msun.
    real(dp),intent(out)::photons(9),energy(9)
    integer,intent(out)::ierr
    real(dp)::a(18),result(18),wz,w,zl,zh
    integer::i,status
    ierr=1;photons=0;energy=0
    if(.not.parsec_sed_enabled.or..not.bound)return
    if(any(.not.ieee_is_finite([t0,t1,z,initial_mass])))return
    if(t1<t0.or.initial_mass<0)return
    call source_metallicity_bracket(z_nodes,z,zl,zh,wz,status)
    if(status/=interpolation_ok)return
    result=0
    do i=1,nn
       if(metal(i)==zl)then
          w=weight(i)*(1-wz)
       else if(metal(i)==zh)then
          w=weight(i)*wz
       else
          cycle
       endif
       if(w==0)cycle
       call node_interval(i,t0,t1,a)
       result=result+w*a
    enddo
    result=result*initial_mass
    if(any(.not.ieee_is_finite(result)).or.any(result<0))return
    photons=result(:9);energy=result(10:18);ierr=0
  end subroutine

  subroutine parsec_sed_identity(values)
    real(dp),allocatable,intent(out)::values(:)
    integer::i
    values=[4d0,real(nn,dp),real(nr,dp),(real(iachar(model(i:i)),dp),i=1,len(model))]
    if(.not.parsec_sed_enabled)return
    values=[values,mass,metal,death,real(fate,dp),real(count,dp),age,reshape(release,[21*nr]),weight]
  end subroutine
end module

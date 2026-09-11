! Matched LC18 two-parent source convolution. No AMR/deposition globals.
module stellar_radioactive_sources
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_radioactive_decay
  use stellar_enrichment_config, only: stellar_dp,yield_basis_per_star_cumulative
  use stellar_enrichment_contract, only: stellar_population_t
  use stellar_yield_tables, only: stellar_yield_table_t,yield_mass_assignment_piecewise_constant
  use stellar_yield_interpolation, only: source_metallicity_bracket
  use stellar_ssp_sources, only: build_source_mass_edges,calculate_imf_mass_fraction
  implicit none
  private
  character(len=*),parameter :: prompt_model='prompt_t12_le_100yr_baryonic_v1'
  character(len=*),parameter :: nuclear_sha='1585a5eea86c5e17e90307c7e6e786d060049c4039e392a261ff6db977df9859'
  real(stellar_dp),parameter :: gyr_s=1d9*radioactive_year_s
  type,public :: radioactive_companion_t
     logical :: loaded=.false.
     integer :: n_rows=0
     character(len=128) :: source_identity(4)=''
     ! Selected table bytes, selected history bytes, nuclear bytes, inventory.
     ! Bind these AND actual cumulative(:,:) in the caller's MPI/HDF identity.
     character(len=64) :: digest(4)=''
     real(stellar_dp),allocatable :: cumulative(:,:) ! (Al26/Fe60, native row)
     integer,allocatable :: next_wind_row(:) ! Derived from the matched row clocks.
  end type
  type(radioactive_companion_t),save,public :: active_radioactive_companion
  public :: load_radioactive_companion,radioactive_node_interval,integrate_radioactive_channel
  public :: configure_radioactive_source,radioactive_source_identity
contains
  subroutine configure_radioactive_source(table,filename,ierr)
    type(stellar_yield_table_t),intent(in)::table
    character(len=*),intent(in)::filename
    integer,intent(out)::ierr
    type(radioactive_companion_t)::empty
    if(len_trim(filename)==0)then
       active_radioactive_companion=empty;ierr=0;return
    endif
    call load_radioactive_companion(filename,table,active_radioactive_companion,ierr)
  end subroutine

  subroutine radioactive_source_identity(values)
    real(stellar_dp),allocatable,intent(out)::values(:)
    real(stellar_dp)::labels(4*128+4*64+len(radioactive_model))
    integer::i,j,k
    allocate(values(0))
    if(.not.active_radioactive_companion%loaded)return
    k=0
    do i=1,4
       do j=1,128
          k=k+1;labels(k)=iachar(active_radioactive_companion%source_identity(i)(j:j))
       enddo
       do j=1,64
          k=k+1;labels(k)=iachar(active_radioactive_companion%digest(i)(j:j))
       enddo
    enddo
    do j=1,len(radioactive_model)
       k=k+1;labels(k)=iachar(radioactive_model(j:j))
    enddo
    values=[1d0,radioactive_half_life_s,real(active_radioactive_companion%n_rows,stellar_dp),labels, &
         reshape(active_radioactive_companion%cumulative,[2*active_radioactive_companion%n_rows])]
  end subroutine

  subroutine load_radioactive_companion(filename,table,companion,ierr)
    character(len=*),intent(in)::filename
    type(stellar_yield_table_t),intent(in)::table ! after prepare_high_mass_history
    type(radioactive_companion_t),intent(inout)::companion
    integer,intent(out)::ierr
    type(radioactive_companion_t)::trial
    integer::unit,status,r,k,version,row_count,n
    character(len=80)::model_id,source_projection
    character(len=128)::source_identity(4)
    character(len=64)::source_yields_sha256,source_history_sha256,nuclear_sha256,inventory_sha256
    character(len=4096)::line
    real(stellar_dp)::half_life_s(2),row(34),native(32),tol,dt,expected(2),ptol(2)
    namelist /stellar_radioactive_companion/ version,model_id,source_projection,row_count, &
         source_identity,source_yields_sha256,source_history_sha256,nuclear_sha256,inventory_sha256,half_life_s
    ierr=1
    if(.not.table%loaded.or..not.table%high_mass_ready)return
    if(table%high_mass_version/=2)return
    if(table%mass_assignment_mode/=yield_mass_assignment_piecewise_constant)return
    if(index(table%high_mass_identity(2),'LC18:SetR:')/=1)return
    if(table%high_mass_identity(2)/=table%high_mass_identity(3))return
    version=0;row_count=0;model_id='';source_projection='';source_identity='';half_life_s=0
    source_yields_sha256='';source_history_sha256='';nuclear_sha256='';inventory_sha256=''
    open(newunit=unit,file=filename,status='old',action='read',iostat=status)
    if(status/=0)return
    read(unit,nml=stellar_radioactive_companion,iostat=status)
    if(status/=0)then
       close(unit);return
    endif
    if(.not.all(ieee_is_finite(half_life_s)))then
       close(unit);return
    endif
    if(version/=1.or.model_id/=radioactive_model.or.source_projection/=prompt_model.or. &
         row_count/=table%n_rows.or.any(source_identity/=table%high_mass_identity).or. &
         nuclear_sha256/=nuclear_sha.or.any(half_life_s/=radioactive_half_life_s))then
       close(unit);return
    endif
    trial%digest=[source_yields_sha256,source_history_sha256,nuclear_sha256,inventory_sha256]
    do k=1,4
       if(len_trim(trial%digest(k))/=64.or.verify(trim(trial%digest(k)),'0123456789abcdef')/=0)then
          close(unit);return
       endif
    enddo
    allocate(trial%cumulative(2,table%n_rows))
    allocate(trial%next_wind_row(table%n_rows));trial%next_wind_row=0
    do r=1,table%n_rows
       read(unit,'(A)',iostat=status)line
       if(status/=0)exit
       ! A full native-row echo makes coordinate/content mismatch observable
       ! without relying on trusting a declared hash of the selected table.
       read(line,*,iostat=status)row
       if(status/=0)exit
       if(.not.all(ieee_is_finite(row)))exit
       native=[real(table%channel(r),stellar_dp),table%initial_mass(r),table%birth_metallicity(r), &
            table%age_gyr(r),table%returned_mass(r),table%remnant_mass(r),table%energy(r), &
            table%momentum(r,:),table%ejected_mass(r,:),table%net_yield(r,:)]
       row(4)=row(4)*1d-9
       if(any(row(:32)/=native).or.any(row(33:)<0))exit
       tol=64*epsilon(1d0)*max(table%returned_mass(r),tiny(1d0))
       if(row(33)>table%returned_mass(r)-sum(table%ejected_mass(r,:))+tol)exit
       if(row(34)>table%ejected_mass(r,11)+tol)exit
       if(table%channel(r)/=1.and.table%channel(r)/=3)then
          if(any(row(33:)/=0))exit ! Fully decayed AGB/Ia not processed again.
       endif
       trial%cumulative(:,r)=row(33:)
    enddo
    if(r<=table%n_rows)then
       close(unit);return
    endif
    read(unit,'(A)',iostat=status)line
    close(unit)
    if(status>=0)return ! No trailing inventory silently ignored.
    ! Check cumulative clocks and subset increments, not only endpoints.
    do n=1,size(table%hm_mass)
       do r=1,table%n_rows
          if(table%initial_mass(r)/=table%hm_mass(n).or.table%birth_metallicity(r)/=table%hm_z(n))cycle
          if(table%channel(r)==3)then
             k=table%hm_terminal_row(n)
             if(table%age_gyr(r)<table%hm_age(n))then
                if(any(trial%cumulative(:,r)/=0))return
             else
                if(any(trial%cumulative(:,r)/=trial%cumulative(:,k)))return
             endif
          else if(table%channel(r)==1)then
             if(table%age_gyr(r)==0.and.any(trial%cumulative(:,r)/=0))return
             k=table%hm_wind_row(n)
             expected=0
             if(table%returned_mass(k)>0)expected=trial%cumulative(:,k)*(table%returned_mass(r)/table%returned_mass(k))
             ptol=64*epsilon(1d0)*max(expected,trial%cumulative(:,r),tiny(1d0))
             if(any(abs(trial%cumulative(:,r)-expected)>ptol))return
             do k=1,table%n_rows
                if(table%channel(k)/=1.or.table%initial_mass(k)/=table%hm_mass(n).or. &
                     table%birth_metallicity(k)/=table%hm_z(n))cycle
                dt=table%age_gyr(r)-table%age_gyr(k)
                if(dt<0)then
                   if(trial%next_wind_row(r)==0)then
                      trial%next_wind_row(r)=k
                   else if(table%age_gyr(k)<table%age_gyr(trial%next_wind_row(r)))then
                      trial%next_wind_row(r)=k
                   endif
                endif
                if(dt>0.and.any(trial%cumulative(:,r)<trial%cumulative(:,k)))return
                if(dt>0.and.table%age_gyr(k)>=table%hm_age(n))then
                   if(any(trial%cumulative(:,r)/=trial%cumulative(:,k)))return
                endif
             enddo
          endif
       enddo
    enddo
    trial%n_rows=table%n_rows;trial%source_identity=source_identity;trial%loaded=.true.
    companion=trial;ierr=0
  end subroutine

  subroutine radioactive_node_interval(table,companion,node,channel,previous_age_gyr,current_age_gyr, &
       surviving,decayed,ierr)
    type(stellar_yield_table_t),intent(in)::table
    type(radioactive_companion_t),intent(in)::companion
    integer,intent(in)::node,channel
    real(stellar_dp),intent(in)::previous_age_gyr,current_age_gyr
    real(stellar_dp),intent(out)::surviving(2),decayed(2)
    integer,intent(out)::ierr
    real(stellar_dp)::lo,hi,emitted(2),s(2),d(2),total_s(2),total_d(2)
    integer::r,next,status
    ierr=1;surviving=0;decayed=0;total_s=0;total_d=0
    if(.not.companion%loaded.or..not.table%high_mass_ready)return
    if(companion%n_rows/=table%n_rows.or.any(companion%source_identity/=table%high_mass_identity))return
    if(node<1.or.node>size(table%hm_mass).or.(channel/=1.and.channel/=3))return
    if(.not.all(ieee_is_finite([previous_age_gyr,current_age_gyr])))return
    if(previous_age_gyr<0.or.current_age_gyr<previous_age_gyr.or.current_age_gyr>huge(1d0)/gyr_s)return
    if(channel==3)then
       if(previous_age_gyr<table%hm_age(node).and.current_age_gyr>=table%hm_age(node))then
          r=table%hm_terminal_row(node)
          call radioactive_release(companion%cumulative(:,r),0d0, &
               (current_age_gyr-table%hm_age(node))*gyr_s,total_s,total_d,status)
          if(status/=0)return
       endif
    else
       do r=1,table%n_rows
          if(table%channel(r)/=1.or.table%initial_mass(r)/=table%hm_mass(node).or. &
               table%birth_metallicity(r)/=table%hm_z(node))cycle
          if(table%age_gyr(r)>=min(current_age_gyr,table%hm_age(node)))cycle
          next=companion%next_wind_row(r)
          if(next==0)return
          lo=max(previous_age_gyr,table%age_gyr(r));hi=min(current_age_gyr,table%age_gyr(next),table%hm_age(node))
          if(hi<=lo)cycle
          emitted=(companion%cumulative(:,next)-companion%cumulative(:,r))* &
               ((hi-lo)/(table%age_gyr(next)-table%age_gyr(r)))
          call radioactive_release(emitted,(hi-lo)*gyr_s,(current_age_gyr-hi)*gyr_s,s,d,status)
          if(status/=0)return
          total_s=total_s+s;total_d=total_d+d
       enddo
    endif
    if(.not.all(ieee_is_finite([total_s,total_d])))return
    surviving=total_s;decayed=total_d;ierr=0
  end subroutine

  subroutine integrate_radioactive_channel(table,companion,population,channel_id,previous_age_gyr, &
       current_age_gyr,mass_min,mass_max,n_mass_bins,surviving,decayed,ierr)
    ! Same source-cell boundaries, IMF mass fractions, nearest-node selection,
    ! M_query/M_node scaling and same-age linear-Z mixture as native material.
    ! Outputs are Msun per SSP at END of this interval, not cumulative ages.
    type(stellar_yield_table_t),intent(in)::table
    type(radioactive_companion_t),intent(in)::companion
    type(stellar_population_t),intent(in)::population
    integer,intent(in)::channel_id,n_mass_bins
    real(stellar_dp),intent(in)::previous_age_gyr,current_age_gyr,mass_min,mass_max
    real(stellar_dp),intent(out)::surviving(2),decayed(2)
    integer,intent(out)::ierr
    real(stellar_dp),allocatable::edges(:)
    real(stellar_dp)::zl,zh,w,zs(2),zw(2),left,right,mass,fraction,weight,best,distance
    real(stellar_dp)::s(2),d(2),total_s(2),total_d(2)
    integer::status,bin,iz,node,j,nz
    surviving=0;decayed=0;total_s=0;total_d=0;ierr=1
    if(.not.companion%loaded.or..not.table%high_mass_ready)return
    if(channel_id/=1.and.channel_id/=3)return
    if(.not.all(ieee_is_finite([previous_age_gyr,current_age_gyr,mass_min,mass_max, &
         population%initial_mass,population%imf_mass_min,population%imf_mass_max,population%birth_metallicity])))return
    if(previous_age_gyr<0.or.current_age_gyr<previous_age_gyr.or.population%initial_mass<=0)return
    if(population%yield_basis_id/=yield_basis_per_star_cumulative.or.n_mass_bins<1)return
    if(mass_min<minval(table%hm_mass).or.mass_max>maxval(table%hm_mass).or.mass_max<=mass_min)return
    if(mass_min<population%imf_mass_min.or.mass_max>population%imf_mass_max)return
    call build_source_mass_edges(table,population,n_mass_bins,edges,status)
    if(status/=0)return
    call source_metallicity_bracket(table%hm_z,population%birth_metallicity,zl,zh,w,status)
    if(status/=0)return
    nz=1;zs=[zl,zh];zw=[1d0,0d0]
    if(zl/=zh)then
       if(.not.table%high_mass_linear_z)return
       nz=2;zw=[1-w,w]
    endif
    do bin=1,size(edges)-1
       left=max(mass_min,edges(bin));right=min(mass_max,edges(bin+1))
       if(right<=left)cycle
       mass=sqrt(left*right)
       call calculate_imf_mass_fraction(population%imf_id,population%imf_mass_min, &
            population%imf_mass_max,left,right,fraction,status)
       if(status/=0)return
       do iz=1,nz
          best=huge(1d0);node=0
          do j=1,size(table%hm_mass)
             if(table%hm_z(j)/=zs(iz))cycle
             if((mass<40d0).neqv.(table%hm_mass(j)<40d0))cycle
             distance=abs(mass-table%hm_mass(j))
             if(distance<best)then
                best=distance;node=j
             endif
          enddo
          if(node==0)return
          call radioactive_node_interval(table,companion,node,channel_id,previous_age_gyr,current_age_gyr,s,d,status)
          if(status/=0)return
          weight=population%initial_mass*fraction/table%hm_mass(node)*zw(iz)
          total_s=total_s+weight*s;total_d=total_d+weight*d
       enddo
    enddo
    if(.not.all(ieee_is_finite([total_s,total_d])))return
    surviving=total_s;decayed=total_d;ierr=0
  end subroutine
end module

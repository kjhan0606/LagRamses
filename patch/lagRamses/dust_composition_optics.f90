! Draine 2003 sphere optical coefficients for four composition/size masses.
! Explicit d03_transport_v1 receiver. Fixed-mix defaults remain unchanged.
module dust_composition_optics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'dust_d03_optics_data.inc'
  include 'dust_d03_band_data.inc'
  integer,parameter,public :: d03_ng=size(d03_primary_ev),d03_nir=size(d03_ir_ev)
  public :: d03_radius_cm,d03_solid_density,d03_edges,d03_primary_ev,d03_ir_ev
  public :: d03_optics_binding,d03_opacity_basis,d03_mix_opacity,d03_absorption_depth
  integer,parameter,public :: d03_identity_size=5+size(d03_edges)+d03_ng+d03_nir+12*(d03_ng+d03_nir)
  public :: d03_identity,d03_cell_weights
  public :: d03_band_nodes,d03_band_ev,d03_band_abs,d03_band_transport,d03_band_sha256
contains
  function d03_identity() result(v)
    real(real64)::v(d03_identity_size)
    ! Version 2 binds delta-isotropic primary AND secondary IR transport.
    v=[2d0,d03_radius_cm,d03_solid_density,d03_edges,d03_primary_ev,d03_ir_ev, &
         reshape(d03_primary_qabs,[4*d03_ng]),reshape(d03_primary_qsca,[4*d03_ng]), &
         reshape(d03_primary_g,[4*d03_ng]),reshape(d03_ir_qabs,[4*d03_nir]), &
         reshape(d03_ir_qsca,[4*d03_nir]),reshape(d03_ir_g,[4*d03_nir])]
  end function

  subroutine d03_cell_weights(bins,weights,ierr)
    real(real64),intent(in)::bins(4)
    real(real64),intent(out)::weights(4)
    integer,intent(out)::ierr
    real(real64)::largest
    weights=.25d0;ierr=1
    if(any(.not.ieee_is_finite(bins)).or.any(bins<0))return
    largest=maxval(bins)
    if(largest>0)then
       weights=bins/largest;weights=weights/sum(weights)
    endif
    ! Zero-density cells need a positive dummy U/P table, never dust mass.
    ierr=0
  end subroutine
  logical function d03_optics_binding(edges,primary,ir,radius,solid_density) result(ok)
    real(real64),intent(in)::edges(:),primary(:),ir(:),radius(2),solid_density(2)
    ! Exact values, not just matching group counts. A different group contract,
    ! radius, or density requires a newly constructed physical table.
    ok=.false.
    if(size(edges)/=size(d03_edges).or.size(primary)/=d03_ng.or.size(ir)/=d03_nir)return
    if(any(.not.ieee_is_finite(edges)).or.any(.not.ieee_is_finite(primary)))return
    if(any(.not.ieee_is_finite(ir)).or.any(.not.ieee_is_finite(radius)))return
    if(any(.not.ieee_is_finite(solid_density)))return
    ok=all(edges==d03_edges).and.all(primary==d03_primary_ev).and.all(ir==d03_ir_ev).and. &
         all(radius==d03_radius_cm).and.all(solid_density==d03_solid_density)
  end function

  subroutine d03_opacity_basis(normalization,primary_abs,primary_sca,primary_sca_g, &
       ir_abs,ir_sca,ir_sca_g,ierr)
    ! normalization=1 -> cm2/g of dust; reference mass/H -> cm2/reference H.
    ! Q*pi*a2 / (4*pi*a3*rho/3), with the D03 material densities (2.2,3.8).
    ! Keep scattering*g separately, permitting an opacity-weighted g on mix.
    real(real64),intent(in)::normalization
    real(real64),intent(out)::primary_abs(d03_ng,4),primary_sca(d03_ng,4),primary_sca_g(d03_ng,4)
    real(real64),intent(out)::ir_abs(d03_nir,4),ir_sca(d03_nir,4),ir_sca_g(d03_nir,4)
    integer,intent(out)::ierr
    real(real64)::area
    integer::s,j,k
    primary_abs=0;primary_sca=0;primary_sca_g=0;ir_abs=0;ir_sca=0;ir_sca_g=0;ierr=1
    if(.not.ieee_is_finite(normalization).or.normalization<=0)return
    do s=1,2
       do j=1,2
          k=2*(s-1)+j
          area=normalization*.75d0/(d03_solid_density(s)*d03_radius_cm(j))
          primary_abs(:,k)=area*d03_primary_qabs(:,k)
          primary_sca(:,k)=area*d03_primary_qsca(:,k)
          primary_sca_g(:,k)=primary_sca(:,k)*d03_primary_g(:,k)
          ir_abs(:,k)=area*d03_ir_qabs(:,k)
          ir_sca(:,k)=area*d03_ir_qsca(:,k)
          ir_sca_g(:,k)=ir_sca(:,k)*d03_ir_g(:,k)
       enddo
    enddo
    if(any(.not.ieee_is_finite(primary_abs)).or.any(.not.ieee_is_finite(primary_sca)))return
    if(any(.not.ieee_is_finite(ir_abs)).or.any(.not.ieee_is_finite(ir_sca)))return
    ierr=0
  end subroutine

  subroutine d03_mix_opacity(bins,radius,solid_density,normalization,primary_abs,primary_sca,primary_g, &
       ir_abs,ir_sca,ir_g,ierr)
    ! Local four-bin masses/densities/fractions. Only their ratio is used.
    ! At zero dust the physical coefficients are zero, not a fabricated mix.
    real(real64),intent(in)::bins(4),radius(2),solid_density(2),normalization
    real(real64),intent(out)::primary_abs(d03_ng),primary_sca(d03_ng),primary_g(d03_ng)
    real(real64),intent(out)::ir_abs(d03_nir),ir_sca(d03_nir),ir_g(d03_nir)
    integer,intent(out)::ierr
    real(real64)::pa(d03_ng,4),ps(d03_ng,4),pg(d03_ng,4)
    real(real64)::ia(d03_nir,4),isc(d03_nir,4),ig(d03_nir,4),f(4),largest
    primary_abs=0;primary_sca=0;primary_g=0;ir_abs=0;ir_sca=0;ir_g=0;ierr=1
    if(.not.d03_optics_binding(d03_edges,d03_primary_ev,d03_ir_ev,radius,solid_density))return
    if(any(.not.ieee_is_finite(bins)).or.any(bins<0))return
    call d03_opacity_basis(normalization,pa,ps,pg,ia,isc,ig,ierr)
    if(ierr/=0)return
    largest=maxval(bins)
    if(largest==0)return
    ! Normalize after scaling to avoid overflow/underflow in sum(bins).
    f=bins/largest;f=f/sum(f)
    primary_abs=matmul(pa,f);primary_sca=matmul(ps,f);primary_g=matmul(pg,f)
    ir_abs=matmul(ia,f);ir_sca=matmul(isc,f);ir_g=matmul(ig,f)
    where(primary_sca>0)primary_g=primary_g/primary_sca
    where(ir_sca>0)ir_g=ir_g/ir_sca
    ! The caller decides angular transport. Do not turn Qsca into (1-g)*Qsca
    ! here, and do not use either scattering coefficient as absorption/heating.
  end subroutine

  subroutine d03_absorption_depth(bins,radius,solid_density,path,tau,ierr)
    ! bins in physical g/cm3, path in cm. No nH or reference mixture enters.
    real(real64),intent(in)::bins(4),radius(2),solid_density(2),path
    real(real64),intent(out)::tau(d03_ng)
    integer,intent(out)::ierr
    real(real64)::pa(d03_ng,4),ps(d03_ng,4),pg(d03_ng,4)
    real(real64)::ia(d03_nir,4),isc(d03_nir,4),ig(d03_nir,4),trial(d03_ng)
    tau=0;ierr=1
    if(.not.d03_optics_binding(d03_edges,d03_primary_ev,d03_ir_ev,radius,solid_density))return
    if(any(.not.ieee_is_finite(bins)).or.any(bins<0))return
    if(.not.ieee_is_finite(path).or.path<0)return
    call d03_opacity_basis(1d0,pa,ps,pg,ia,isc,ig,ierr)
    if(ierr/=0)return
    ierr=1
    trial=matmul(pa,bins)*path
    if(any(.not.ieee_is_finite(trial)).or.any(trial<0))return
    tau=trial;ierr=0
  end subroutine
end module

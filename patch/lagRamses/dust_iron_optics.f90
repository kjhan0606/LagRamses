! C/silicate/Fe six-component optical bank. Fe is the electric+eddy base
! only: spin-magnetic absorption is not included. Full admission stays false;
! the explicitly named live comparison does not confer full-band approval.
module dust_iron_optics
  use dust_composition_optics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'dust_fe_electric_base_data.inc'
  include 'dust_fe_band_data.inc'
  public :: fe_band_ev,fe_band_abs,fe_band_transport,fe_band_sha256,fe_band_nodes
  real(real64),parameter,public :: fe_six_band_abs(d03_band_nodes,d03_ng,6)=reshape( &
       [reshape(d03_band_abs,[d03_band_nodes*d03_ng*4]),reshape(fe_band_abs,[fe_band_nodes*d03_ng*2])], &
       [d03_band_nodes,d03_ng,6])
  real(real64),parameter,public :: fe_six_band_transport(d03_band_nodes,d03_ng,6)=reshape( &
       [reshape(d03_band_transport,[d03_band_nodes*d03_ng*4]),reshape(fe_band_transport,[fe_band_nodes*d03_ng*2])], &
       [d03_band_nodes,d03_ng,6])
  integer,parameter,public :: fe_ng=size(fe_primary_ev),fe_nir=size(fe_ir_ev),fe_nt=size(fe_temperature)
  logical,parameter,public :: fe_full_optics_admitted=.false.
  public :: fe_electric_basis,fe_six_opacity_basis,fe_base_binding
  integer,parameter,public :: fe_optics_identity_n=5+size(fe_edges)+fe_ng+2*fe_nir+fe_nt+6*fe_ng+6*fe_nir
  public :: fe_optics_identity
  public :: fe_radius_cm,fe_density,fe_edges,fe_primary_ev,fe_ir_ev,fe_ir_weight_ev,fe_temperature
contains
  function fe_optics_identity() result(v)
    real(real64)::v(fe_optics_identity_n)
    ! Exact arrays and conventions, not just a user label or mutable path.
    v=[1d0,0d0,fe_radius_cm,fe_density,fe_edges,fe_primary_ev,fe_ir_ev,fe_ir_weight_ev,fe_temperature, &
         reshape(fe_primary_qabs,[2*fe_ng]),reshape(fe_primary_qsca,[2*fe_ng]),reshape(fe_primary_g,[2*fe_ng]), &
         reshape(fe_ir_qabs,[2*fe_nir]),reshape(fe_ir_qsca,[2*fe_nir]),reshape(fe_ir_g,[2*fe_nir])]
  end function
  logical function fe_base_binding(edges,primary,ir,radius,density) result(ok)
    real(real64),intent(in)::edges(:),primary(:),ir(:),radius(2),density
    ok=.false.
    if(size(edges)/=size(fe_edges).or.size(primary)/=fe_ng.or.size(ir)/=fe_nir)return
    if(any(.not.ieee_is_finite(edges)).or.any(.not.ieee_is_finite(primary)))return
    if(any(.not.ieee_is_finite(ir)).or.any(.not.ieee_is_finite(radius)).or..not.ieee_is_finite(density))return
    ok=all(edges==fe_edges).and.all(primary==fe_primary_ev).and.all(ir==fe_ir_ev).and. &
         all(radius==fe_radius_cm).and.density==fe_density(1)
  end function

  subroutine fe_electric_basis(normalization,pa,ps,psg,ia,isc,isg,ierr)
    ! normalization 1 gives cm2/g. Qsca and Qsca*g are distinct from Qabs;
    ! neither is gas/dust heating. Retain the forward hard-X-ray scattering.
    real(real64),intent(in)::normalization
    real(real64),intent(out)::pa(fe_ng,2),ps(fe_ng,2),psg(fe_ng,2)
    real(real64),intent(out)::ia(fe_nir,2),isc(fe_nir,2),isg(fe_nir,2)
    integer,intent(out)::ierr
    real(real64)::area
    integer::j
    pa=0;ps=0;psg=0;ia=0;isc=0;isg=0;ierr=1
    if(.not.ieee_is_finite(normalization).or.normalization<=0)return
    do j=1,2
       area=normalization*.75d0/(fe_density(1)*fe_radius_cm(j))
       pa(:,j)=area*fe_primary_qabs(:,j);ps(:,j)=area*fe_primary_qsca(:,j)
       psg(:,j)=ps(:,j)*fe_primary_g(:,j)
       ia(:,j)=area*fe_ir_qabs(:,j);isc(:,j)=area*fe_ir_qsca(:,j)
       isg(:,j)=isc(:,j)*fe_ir_g(:,j)
    enddo
    if(any(.not.ieee_is_finite(pa)).or.any(.not.ieee_is_finite(ps)))return
    if(any(.not.ieee_is_finite(ia)).or.any(.not.ieee_is_finite(isc)))return
    ierr=0
  end subroutine

  subroutine fe_six_opacity_basis(normalization,pa,ps,psg,ia,isc,isg,ierr)
    real(real64),intent(in)::normalization
    real(real64),intent(out)::pa(fe_ng,6),ps(fe_ng,6),psg(fe_ng,6)
    real(real64),intent(out)::ia(fe_nir,6),isc(fe_nir,6),isg(fe_nir,6)
    integer,intent(out)::ierr
    ! One actual group contract: do not concatenate equal-length tables with
    ! different edge/mean/IR energies, including the 2--10keV source group.
    pa=0;ps=0;psg=0;ia=0;isc=0;isg=0;ierr=1
    if(.not.fe_base_binding(d03_edges,d03_primary_ev,d03_ir_ev,fe_radius_cm,fe_density(1)))return
    call d03_opacity_basis(normalization,pa(:,1:4),ps(:,1:4),psg(:,1:4),ia(:,1:4),isc(:,1:4),isg(:,1:4),ierr)
    if(ierr/=0)return
    call fe_electric_basis(normalization,pa(:,5:6),ps(:,5:6),psg(:,5:6),ia(:,5:6),isc(:,5:6),isg(:,5:6),ierr)
  end subroutine
end module

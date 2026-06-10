program test_darkcool
  use amr_parameters
  use dark_cooling_mod, only: dark_net_cooling
  implicit none
  real(dp)::T, aex, Latom_lo, Lmol_lo, Lmol_hi, Ltmp
  real(dp)::n_lo, n_hi, x, logT0, logT1
  integer::i, iu
  integer,parameter::nT=200

  aex = 1.0d0
  n_lo = 1.0d0      ! low-density probe [cm^-3]
  n_hi = 1.0d8      ! high-density (LTE) probe [cm^-3]
  logT0 = log10(1.0d0)
  logT1 = log10(1.0d6)

  ! ---------- SM reference (r=1): must reproduce Glover-Abel SM H2 ----------
  adm_alpha    = 7.2973525693d-3
  adm_mp       = 0.9382720813d0
  adm_me_ratio = 0.5109989461d-3 / 0.9382720813d0
  adm_xi       = 0.5d0
  adm_fH2      = 1.0d-3
  call sweep('darkcool_sm.dat')

  ! ---------- Dark example (M,m,alpha)=(40 GeV,40 keV,0.01), paper Fig 1 ----------
  adm_alpha    = 0.01d0
  adm_mp       = 40.0d0
  adm_me_ratio = 40.0d-6 / 40.0d0
  call sweep('darkcool_dark.dat')

  write(*,*) 'done: darkcool_sm.dat, darkcool_dark.dat'

contains

  subroutine sweep(fname)
    character(len=*),intent(in)::fname
    open(newunit=iu,file=fname,status='replace')
    write(iu,'(A)') '# T[K]  Latom_lo  Lmol_lo  Lmol_hi  [erg/s/cm^3]'
    do i=0,nT
       x=logT0+(logT1-logT0)*dble(i)/dble(nT)
       T=10.0d0**x
       adm_mol=.false.; Latom_lo = dark_net_cooling(T,n_lo,aex)
       adm_mol=.true.;  Lmol_lo  = dark_net_cooling(T,n_lo,aex) - Latom_lo
       adm_mol=.false.; Ltmp     = dark_net_cooling(T,n_hi,aex)
       adm_mol=.true.;  Lmol_hi  = dark_net_cooling(T,n_hi,aex) - Ltmp
       write(iu,'(4ES14.5)') T, Latom_lo, Lmol_lo, Lmol_hi
    end do
    close(iu)
  end subroutine sweep

end program test_darkcool

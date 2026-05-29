subroutine diag_check_eint(label, ilev_check)
  !------------------------------------------------------------------------
  ! Scan leaf cells and report min(eint) for debugging eint<0 crashes.
  !
  ! ilev_check = 0  : scan ALL levels (for coarse-level operations)
  ! ilev_check > 0  : scan only that level (for per-level operations)
  !
  ! For coarse-level (ilev_check=0): always print min eint (with MPI)
  ! For per-level (ilev_check>0): only print if eint < 0 found locally
  !------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=*), intent(in) :: label
  integer, intent(in) :: ilev_check

  integer :: ilevel, igrid, ind, iskip, icell, ncache, info
  integer :: lmin, lmax
  real(dp) :: d, u, v, w, etot, eint
  real(dp) :: eint_min_loc, eint_min_glob
  integer :: ilev_worst_loc, ilev_worst_glob
  ! For negative eint reporting
  real(dp) :: d_worst, u_worst, v_worst, w_worst, etot_worst
  integer :: neg_count_loc, neg_count_glob

  ! Diagnostic reads uold; only allocated when hydro=.true.
  ! Several callers in amr_step.jaehyun.f90 (e.g. the unconditional
  ! 'cooling' call at line 833) fire even on poisson-only / DM-only runs.
  if(.not. hydro) return

  if(ilev_check == 0) then
     lmin = levelmin
     lmax = nlevelmax
  else
     lmin = ilev_check
     lmax = ilev_check
  endif

  eint_min_loc = 1d30
  ilev_worst_loc = 0
  neg_count_loc = 0
  d_worst = 0d0; u_worst = 0d0; v_worst = 0d0; w_worst = 0d0; etot_worst = 0d0

  do ilevel = lmin, lmax
     ncache = active(ilevel)%ngrid
     do igrid = 1, ncache
        do ind = 1, twotondim
           iskip = ncoarse + (ind-1)*ngridmax
           icell = iskip + active(ilevel)%igrid(igrid)
           if(son(icell) /= 0) cycle  ! skip non-leaf
           d = uold(icell,1)
           if(d <= 0d0) cycle
           u = uold(icell,2)/d
           v = uold(icell,3)/d
           w = uold(icell,4)/d
           etot = uold(icell,5)
           eint = etot - 0.5d0*d*(u*u + v*v + w*w)
           if(eint < 0d0) neg_count_loc = neg_count_loc + 1
           if(eint < eint_min_loc) then
              eint_min_loc = eint
              ilev_worst_loc = ilevel
              d_worst = d
              u_worst = u
              v_worst = v
              w_worst = w
              etot_worst = etot
           endif
        enddo
     enddo
  enddo

  if(ilev_check == 0) then
     ! Coarse-level: full MPI reduction, always print
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(eint_min_loc, eint_min_glob, 1, &
          MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(neg_count_loc, neg_count_glob, 1, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#else
     eint_min_glob = eint_min_loc
     neg_count_glob = neg_count_loc
#endif
     if(myid == 1) then
        if(eint_min_glob < 0d0) then
           write(*,'(A,A,A,ES11.3,A,I6)') &
                ' *** DIAG ', trim(label), &
                ': eint_min=', eint_min_glob, ' neg_cells=', neg_count_glob
        else
           write(*,'(A,A,A,ES11.3)') &
                ' DIAG ', trim(label), ': eint_min=', eint_min_glob
        endif
     endif
     ! Also print local details on the rank that has the worst cell
     if(eint_min_loc == eint_min_glob .and. eint_min_glob < 0d0) then
        write(*,'(A,I4,A,I2,A,ES11.3,A,ES11.3,A,3(ES10.2))') &
             ' *** DIAG rank=', myid, ' lev=', ilev_worst_loc, &
             ' eint=', eint_min_loc, ' d=', d_worst, &
             ' v=', u_worst, v_worst, w_worst
     endif
  else
     ! Per-level: only print if negative eint found (no MPI to save cost)
     if(neg_count_loc > 0) then
        write(*,'(A,I4,A,A,A,I2,A,ES11.3,A,I6,A,ES11.3,A,3(ES10.2))') &
             ' *** DIAG rank=', myid, ' ', trim(label), &
             ' lev=', ilev_check, &
             ' eint_min=', eint_min_loc, ' neg=', neg_count_loc, &
             ' d=', d_worst, ' v=', u_worst, v_worst, w_worst
     endif
  endif

end subroutine diag_check_eint

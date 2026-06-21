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

subroutine diag_check_nan(label)
  use amr_commons
  use hydro_commons
  use poisson_commons, only: f
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=*), intent(in) :: label

  integer :: ilevel, igrid, ind, iskip, icell, ncache, info, ivar
  integer :: nan_uold_loc, nan_f_loc
  integer :: nan_uold_glob, nan_f_glob
  integer :: dzero_loc, dzero_glob
  integer :: first_lev, first_ivar
  real(dp) :: first_d

  if(.not. hydro) return

  nan_uold_loc = 0
  nan_f_loc = 0
  dzero_loc = 0
  first_lev = 0
  first_ivar = 0
  first_d = 0d0

  do ilevel = levelmin, nlevelmax
     ncache = active(ilevel)%ngrid
     do igrid = 1, ncache
        do ind = 1, twotondim
           iskip = ncoarse + (ind-1)*ngridmax
           icell = iskip + active(ilevel)%igrid(igrid)
           if(uold(icell,1) <= 0d0) dzero_loc = dzero_loc + 1
           do ivar = 1, nvar
              if(uold(icell,ivar) /= uold(icell,ivar)) then
                 nan_uold_loc = nan_uold_loc + 1
                 if(first_lev == 0) then
                    first_lev = ilevel
                    first_ivar = ivar
                    first_d = uold(icell,1)
                 end if
              end if
           end do
           if(poisson .and. allocated(f)) then
              do ivar = 1, ndim
                 if(f(icell,ivar) /= f(icell,ivar)) then
                    nan_f_loc = nan_f_loc + 1
                 end if
              end do
           end if
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nan_uold_loc, nan_uold_glob, 1, &
       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(nan_f_loc, nan_f_glob, 1, &
       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(dzero_loc, dzero_glob, 1, &
       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#else
  nan_uold_glob = nan_uold_loc
  nan_f_glob = nan_f_loc
  dzero_glob = dzero_loc
#endif
  if(myid == 1) then
     write(*,'(A,A,A,I10,A,I10,A,I10)') &
          ' NaN_CHK ', trim(label), &
          ': uold=', nan_uold_glob, ' f=', nan_f_glob, &
          ' d0=', dzero_glob
  end if
  if(nan_uold_loc > 0) then
     write(*,'(A,I4,A,A,A,I10,A,I2,A,I2,A,ES11.3)') &
          ' NaN_CHK rank=', myid, ' ', trim(label), &
          ' uold_nan=', nan_uold_loc, &
          ' first_lev=', first_lev, ' first_ivar=', first_ivar, &
          ' d=', first_d
  end if
  if(dzero_loc > 0 .and. nan_uold_loc == 0) then
     write(*,'(A,I4,A,A,A,I10)') &
          ' NaN_CHK rank=', myid, ' ', trim(label), &
          ' d<=0_cells=', dzero_loc
  end if
end subroutine diag_check_nan

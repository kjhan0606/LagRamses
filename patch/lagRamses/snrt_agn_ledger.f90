module snrt_agn_ledger
  ! Optional runtime evidence for the sink-to-RT source contract.  This is
  ! diagnostic-only: it neither changes sink variables nor emits photons.
  use amr_parameters
  use amr_commons
  use pm_commons
  implicit none

  logical, save :: wrote_header = .false.

contains

  subroutine snrt_agn_ledger_diagnose(ilevel)
    integer, intent(in) :: ilevel
    integer :: unit_id

    if (ilevel /= levelmin .or. myid /= 1 .or. nsink <= 0) return
    if (.not. allocated(dMsmbh) .or. .not. allocated(msink) .or. &
         .not. allocated(tsink)) return

    open(newunit=unit_id, file='snrt_agn_ledger.csv', status='unknown', &
         position='append', action='write')
    if (.not. wrote_header) then
       write(unit_id,'(a)') 'sink_id,tsink,dMsmbh,msink'
       wrote_header = .true.
    end if
    write(unit_id,'(i0,3(a,es24.16))') 1, ',', tsink(1), ',', dMsmbh(1), ',', msink(1)
    close(unit_id)
  end subroutine snrt_agn_ledger_diagnose

end module snrt_agn_ledger

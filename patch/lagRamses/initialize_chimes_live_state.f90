! External entry keeps the shared adaptive loop independent of optional modules.
subroutine initialize_chimes_live_state
  use snrt_chimes_runtime, only: chimes_initialize_hierarchy
  implicit none
  call chimes_initialize_hierarchy()
end subroutine

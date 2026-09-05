! Minimal GNU smoke-only stand-in for the RAMSES amr_parameters module.
! The transaction core imports only the dp kind; the production module graph
! is compiled separately with the real RAMSES amr_parameters module.
module amr_parameters
  integer, parameter :: dp = kind(1.0d0)
end module amr_parameters

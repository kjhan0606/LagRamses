! Build-selected group identity shared by SNRT, CHIMES and Fe UV interfaces.
module snrt_group_profile
  use amr_parameters, only: dp
  implicit none
  private
#include "snrt_group_profile_f90.inc"
end module snrt_group_profile

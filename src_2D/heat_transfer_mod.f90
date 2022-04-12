module heat_transfer_module
  
  ! -----------------------------------------------------------------
  ! This module is used as a wrapper for the modules that contain
  ! the solvers for the heat transfer equations
  ! -----------------------------------------------------------------
  
  use heat_transfer_explicit_module
  use heat_transfer_implicit_module
  use heat_transfer_explicit_no_subcycling
  
end module heat_transfer_module

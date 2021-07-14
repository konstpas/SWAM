module my_amr_fortran_interface_module

  use amrex_amr_module

  implicit none

  private

  public :: amrex_imultifab_swap

contains

  subroutine amrex_imultifab_move (dst, src)
    class(amrex_imultifab), intent(inout) :: dst
    type (amrex_imultifab), intent(inout) :: src
    call amrex_imultifab_destroy(dst)
    dst%owner = src%owner
    dst%p     = src%p
    dst%nc    = src%nc
    dst%ng    = src%ng
    call dst%ba%move(src%ba)
    call dst%dm%move(src%dm)
    src%owner = .false.
    src%p     = c_null_ptr
  end subroutine amrex_imultifab_move

  
  subroutine amrex_imultifab_swap(mf1, mf2)
    type(amrex_imultifab), intent(inout) :: mf1, mf2
    type(amrex_imultifab) :: mftmp
    call amrex_imultifab_move(mftmp, mf1)
    call amrex_imultifab_move(mf1, mf2)
    call amrex_imultifab_move(mf2, mftmp)
    call amrex_imultifab_destroy(mftmp)
    
  end subroutine amrex_imultifab_swap


end module my_amr_fortran_interface_module

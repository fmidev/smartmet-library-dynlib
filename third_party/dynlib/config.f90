! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- configuration variables for all modules
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module config
  use kind
  !
  logical :: grid_cyclic_ew = .false.
  !
  ! Front / Line / Axis detection settings
  real(kind=nr) :: frint_thres =-2.0e-11_nr, &
    &              frspd_thres = 1.5_nr, &
    &                vor_thres = 0.7e-4_nr, &
    &                div_thres =-0.25e-4_nr, &
    &                def_thres = 0.7e-4_nr, &
    &             jetint_thres = 2.5e-4_nr, & 
    &                searchrad = 1.5_nr, & 
    &              smooth_coef = 1.0_nr/6.0_nr, &
    &                   minlen = 1.0e6, &
    &   tfp_mingrad_largescale = 4.5e-5, &
    &   tfp_mingrad_smallscale = 7.5e-5, &
    &                  minsize = 75.0e9, &
    &          thres_min_slope = 1.0e-5, &
    &          thres_max_slope = 1.0e3, &
    &          thres_min_dzdth = 1.0e-4, &
    &          thres_max_dzdth = 1.0e4, &
    &              cyc_minsize = 800.0, &
    &              cyc_maxsize = 4.5e6, &
    &               cyc_maxoro = 1500, &
    &              cyc_mindist = 750, &
    &        cyc_minprominence = 200
  integer(kind=ni) ::  nsmooth = 2_ni, &
    &                 block_dj = 30, &  ! 15 deg latitude at ERAI resolution
    &                 block_di = 15     ! 7.5 deg longitude at ERAI resolution
  !
contains
  !
  !@ Dummy subroutine for f2py to keep this module in the dynfor python library, parameter/variable-only modules
  !@ are skipped for numpy>=1.26.2. The routine does nothing but return a zero.
  subroutine keep_me(res)
    real(kind=nr), intent(out) :: res
    ! -----------------------------------------------------------------
    !
    res = 0.0_nr
  end subroutine
  !
end module

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- some constants
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module consts
  use kind
  !
  character(len=48), parameter :: version  = 'dynlib ####VERREVI####'
  character(len=48), parameter :: verdate  = '####VERDATE####'
  character(len=48), parameter :: compiler = '####FCOMPIL####'
  character(len=48), parameter :: comptime = '####FCOMTIM####'
  character(len=48), parameter :: comphost = '####FCOMHOS####'
  !
  ! Some thermodynamic constants
  real(kind=nr), parameter :: cp = 1004.0_nr     ! specific heat constant of dry air for constant pressure
  real(kind=nr), parameter :: Rl = 287.04_nr     ! Gas constant for dry air
  !real(kind=nr), parameter :: kappa = Rl/cp      ! Exponent of the Exner function
  real(kind=nr), parameter :: Rv = 461.5_nr      ! Gas constant for water vapour
  real(kind=nr), parameter :: Lv = 2.501e6_nr    ! Latent heat of vaporisation
  !
  ! Other constants
  real(kind=nr), parameter :: p0 = 100000.0_nr   ! Reference pressure for T <-> Theta conversion
  real(kind=nr), parameter :: omega_rot = 7.292123516990375e-05 ! angular velocity of the earth's rotation
  real(kind=nr), parameter :: g  = 9.81_nr       ! Gravitiy constant close to Earth surface
  !
  ! Mathematical constants
  real(kind=nr), parameter :: pi = 3.14159265359_nr ! should be clear
  !
  ! Technical stuff
  character, parameter :: cr = char(13_1)        ! Carriage Retirm is ASCI code 13
  real(kind=nr), parameter :: nan = 0.0_nr/0.0_nr ! Not-a-Number, used for missing values. Requires the gfortran compiler option -fno-range-check
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

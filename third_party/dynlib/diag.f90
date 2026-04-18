! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- diagnostical functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Diagnostic functions for dynamic variables. 
!@ 
!@ This module contains diagnostic functions useful for analysing the wind
!@ field. For diagnostics relating to thermodynamic variables, refer to the
!@ :mod:`dynlib.thermodyn` module.
module diag
  use kind
  use config
  use derivatives
  !
  implicit none
contains
  !
  !@ Calculate vorticity::
  !@
  !@     vor = dv/dx - du/dy 
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated vorticity.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`vor_curv`, :meth:`div`, :meth:`def_shear`, :meth:`def_stretch`
  subroutine vor(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dyu(nz,ny,nx),dxv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxv,nx,ny,nz,v,dx,dy)
    call ddy(dyu,nx,ny,nz,u,dx,dy)
    res = dxv - dyu
  end subroutine
  !
  !@ Calculate curvature vorticity::
  !@
  !@     vor_c = dv/ds
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated vorticity.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_stretch_nat`, :meth:`vor`, :meth:`div`, :meth:`def_shear`, :meth:`def_stretch`
  subroutine vor_curv(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: ux(nz,ny,nx), uy(nz,ny,nx), vx(nz,ny,nx), vy(nz,ny,nx), &
       &             ff2(nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    call grad(ux,uy, nx,ny,nz, u, dx,dy)
    call grad(vx,vy, nx,ny,nz, v, dx,dy)
    !
    ff2 = u(:,:,:)**2_ni + v(:,:,:)**2_ni
    res = (u*u*vx + u*v*(vy - ux) - v*v*uy) / ff2
  end subroutine
  !
  !@ Calculate divergence::
  !@
  !@     div = du/dx + dv/dy 
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated divergence.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`vor`, :meth:`def_shear`, :meth:`def_stretch`
  subroutine div(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dxu(nz,ny,nx),dyv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxu,nx,ny,nz,u,dx,dy)
    call ddy(dyv,nx,ny,nz,v,dx,dy)
    res = dxu + dyv
  end subroutine
  !
  !@ Calculate shear deformation::
  !@
  !@     def_shear = du/dy + dv/dx
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated shear deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_stretch`, :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_shear(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dyu(nz,ny,nx),dxv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    call ddx(dxv,nx,ny,nz,v,dx,dy)
    call ddy(dyu,nx,ny,nz,u,dx,dy)
    res = dyu + dxv
  end subroutine
  !
  !@ Calculate stretch deformation::
  !@
  !@    def_stretch =  du/dx - dv/dy
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated stretching deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_shear`, :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_stretch(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dxu(nz,ny,nx),dyv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxu,nx,ny,nz,u,dx,dy)
    call ddy(dyv,nx,ny,nz,v,dx,dy)
    res = dxu - dyv
  end subroutine
  !
  !@ Calculate stretch deformation in natural coordinates::
  !@
  !@    def_stretch =  du/ds - dv/dn
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated stretching deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_shear_nat`, :meth:`def_stretch`,  :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_stretch_nat(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: ux(nz,ny,nx), uy(nz,ny,nx), vx(nz,ny,nx), vy(nz,ny,nx), &
       &             ff2(nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    call grad(ux,uy, nx,ny,nz, u, dx,dy)
    call grad(vx,vy, nx,ny,nz, v, dx,dy)
    !
    ff2 = u(:,:,:)**2_ni + v(:,:,:)**2_ni
    res = ( (u*u - v*v)*(ux - vy) + 2.0_nr*u*v*(vx + uy) ) / ff2
  end subroutine
  !
  !@ Calculate shearing deformation in natural coordinates::
  !@
  !@    def_shear =  du/dn + dv/ds
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated shearing deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`vor_curv`, :meth:`def_stretch_nat`, :meth:`def_shear`, :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_shear_nat(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: ux(nz,ny,nx), uy(nz,ny,nx), vx(nz,ny,nx), vy(nz,ny,nx), &
       &             ff2(nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    call grad(ux,uy, nx,ny,nz, u, dx,dy)
    call grad(vx,vy, nx,ny,nz, v, dx,dy)
    !
    ff2 = u(:,:,:)**2_ni + v(:,:,:)**2_ni
    res = ( (u*u - v*v)*(uy + vy) + 2.0_nr*u*v*(vy - ux) ) / ff2
  end subroutine
  !
  !@ Calculate total deformation
  !@
  !@ Total deformation is a coordinate-system independent measure for the
  !@ strength of deformation::
  !@
  !@    def_total = sqrt(def_stretch^2 + def_shear^2)
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated total deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_angle`, :meth:`def_angle_nat`
  subroutine def_total(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    res = sqrt(sig_sh**2_ni + sig_st**2_ni)
  end subroutine
  !
  !@ Calculate the angle between the x-axis and the axis of dilatation
  !@ 
  !@ This angle goes by many different symbols:
  !@
  !@  * Spensberger and Spengler (2014): ``gamma``
  !@  * Markowski and Richardson (2011): ``alpha``
  !@  * Keyser, Reeder and Reed (1988), Lapeyre, Klein and Hua (1999): ``gamma``
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated deformation angle.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_angle_nat`, :meth:`def_total`
  subroutine def_angle(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)    
    res = 0.5_nr * atan2(sig_sh,sig_st)
  end subroutine
  !
  !@ Calculate the angle between the wind direction and the axis of dilatation
  !@ 
  !@ Yields the equivalent of def_angle in natural coordinates.
  !@
  !@ Parameters
  !@ ----------
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated deformation angle.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_angle`, :meth:`def_total`
  subroutine def_angle_nat(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)    
    res = 0.5_nr * atan2(sig_sh,sig_st)
    !
    ! rotate along wind direction instead of x-axis
    res = res - atan2(v,u)
    where ( res >= pi/2.0_nr )
       res = res - pi
    end where
    where ( res < -pi/2.0_nr )
       res = res + pi
    end where
  end subroutine
  !
  !@ Calculate the local wind shear in natural coordinates
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated wind shear.
  subroutine shear_nat(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: ff(nz,ny,nx), ffx(nz,ny,nx), ffy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ff = sqrt(u(:,:,:)**2_ni + v(:,:,:)**2_ni)
    call grad(ffx, ffy, nx,ny,nz, ff, dx,dy)
    !
    res(:,:,:) = (u(:,:,:)*ffy(:,:,:) - v(:,:,:)*ffx(:,:,:))/ff(:,:,:)
  end subroutine
  !
  !@ Calculate the local wind shear gradient in natural coordinates
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated wind shear gradient.
  subroutine grad_shear_nat(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: shear(nz,ny,nx), shearx(nz,ny,nx), sheary(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call shear_nat(shear, nx,ny,nz, u,v, dx,dy)
    call ddx(shearx, nx,ny,nz, shear, dx,dy)
    call ddy(sheary, nx,ny,nz, shear, dx,dy)
    !
    res(:,:,:) = (u(:,:,:)*sheary(:,:,:) - v(:,:,:)*shearx(:,:,:)) &
            &   /sqrt(u(:,:,:)**2_ni + v(:,:,:)**2_ni)
  end subroutine
  !
  !@ Calculate total deformation and the deformation angle with all shear removed
  !@
  !@ Wind shear is both reflected in vorticity and deformation. This function
  !@ calculates the part of deformation that cannot be interpreted as vorticity.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated total deformation and deformation angle.
  subroutine def_nat_shearless(resabs,resang,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resabs(nz,ny,nx), resang(nz,ny,nx)
    real(kind=nr) :: ff2(nz,ny,nx), shear(nz,ny,nx), ux(nz,ny,nx), uy(nz,ny,nx), &
            &        vx(nz,ny,nx), vy(nz,ny,nx), def_sh(nz,ny,nx), def_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resabs,resang, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call shear_nat(shear, nx,ny,nz, u,v, dx,dy)
    !
    ff2(:,:,:) = u(:,:,:)**2_ni + v(:,:,:)**2_ni
    !
    call grad(ux,uy, nx,ny,nz, u, dx, dy)
    call grad(vx,vy, nx,ny,nz, v, dx, dy)
    !
    ux(:,:,:) = ux(:,:,:) - shear(:,:,:)*(u(:,:,:)*v(:,:,:))/ff2(:,:,:)
    uy(:,:,:) = uy(:,:,:) - shear(:,:,:)*(u(:,:,:)**2_ni)/ff2(:,:,:)
    vx(:,:,:) = vx(:,:,:) - shear(:,:,:)*(v(:,:,:)**2_ni)/ff2(:,:,:)
    vy(:,:,:) = vy(:,:,:) - shear(:,:,:)*(u(:,:,:)*v(:,:,:))/ff2(:,:,:)
    !
    def_sh(:,:,:) = uy(:,:,:) + vx(:,:,:)
    def_st(:,:,:) = ux(:,:,:) - vy(:,:,:)
    !
    resabs = sqrt(def_sh**2_ni + def_st**2_ni)
    resang = 0.5_nr * atan2(def_sh, def_st)
    !
    ! rotate along wind direction instead of x-axis
    resang = resang - atan2(v,u)
    where ( resang >= pi/2.0_nr )
       resang = resang - pi
    end where
    where ( resang < -pi/2.0_nr )
       resang = resang + pi
    end where
  end subroutine
  !
  !@ Calculate 3d deformation
  !@ 
  !@ The generalisation of deformation to 3 dimensions is done using the
  !@ analogy to Lyapunov exponents.
  !@ 
  !@ To be applied on pressure levels.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ w : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Pressure vertical wind velocity.
  !@ rho : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Density of air, used to convert the pressure vertical velocity to a cartensian 
  !@     vertical velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dz : np.ndarray with shape (nz-2,ny,nx) and dtype float64
  !@     The double grid spacing in z-direction to be directly for centered differences.
  !@     ``dz(k,j,i)`` is expected to contain the z-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
  !@ lstretch_only : boolean
  !@     Calculate stretching part of 3D deformation only
  !@ lshear_only : boolean
  !@     Calculate shearing part of 3D deformation only
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z-direction.
  !@ nt : int
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shapes (nt,nz-2,ny,nx), (nt,nz-2,ny,nx,3) and dtype float64
  !@     Calculated 3D eigenvalues and eigenvectors showing the strength and
  !@     orientation of deformation in 3D.
  subroutine def_3d(res_eval,res_evec,nx,ny,nz,nt,u,v,w,rho,dx,dy,dz,lstretch_only,lshear_only)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), w(nt,nz,ny,nx), & 
            &                     rho(nt,nz,ny,nx), dx(ny,nx), dy(ny,nx), dz(nt,nz,ny,nx)
    real(kind=nr), intent(out) :: res_eval(nt,nz,ny,nx,3_ni), &
                                  res_evec(nt,nz,ny,nx,3_ni,3_ni)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    logical, intent(in) :: lstretch_only, lshear_only
    !f2py depend(nx,ny,nz,nt) res, v, w
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: ux(nt,nz,ny,nx), uy(nt,nz,ny,nx), uz(nt,nz,ny,nx), &
            &        vx(nt,nz,ny,nx), vy(nt,nz,ny,nx), vz(nt,nz,ny,nx), &
            &        wx(nt,nz,ny,nx), wy(nt,nz,ny,nx), wz(nt,nz,ny,nx), &
            &        vor_x(nt,nz,ny,nx), vor_y(nt,nz,ny,nx), vor_z(nt,nz,ny,nx), &
            &        div(nt,nz,ny,nx), wm(nt,nz,ny,nx), dzm(nt,nz,ny,nx)
    real(kind=nr) :: jac(3_ni,3_ni), evalr(3_ni), tau(2_ni), &
            &        evec(3_ni,3_ni), work(96_ni), diag(3_ni), subdiag(2_ni)
    integer(kind=ni) :: i,j,k,n,ti, info, ax,zx, minidx,maxidx
    ! -----------------------------------------------------------------
    !
    if ( lshear_only .and. lstretch_only ) then
       write(*,*) 'Either lshear_only or lstretch_only or both must be False'
       stop 1
    end if
    !
    ! Crude transformation from Pa/s to m/s
    wm = -w(:,:,:,:) * rho(:,:,:,:)*g
    !
    ! Transforming vertical grid distances into m
    dzm = -dz(:,:,:,:) * rho(:,:,:,:)*g
    !
    call grad_3d(ux,uy,uz, nx,ny,nz,nt, u, dx,dy,dzm)
    call grad_3d(vx,vy,vz, nx,ny,nz,nt, v, dx,dy,dzm)
    call grad_3d(wx,wy,wz, nx,ny,nz,nt, wm, dx,dy,dzm)
    !
    ! Subtracting divergence and vorticity
    if ( .not. lshear_only ) then
       div = ux + vy + wz
       ux = ux - 1.0_nr/3.0_nr * div
       vy = vy - 1.0_nr/3.0_nr * div
       wz = wz - 1.0_nr/3.0_nr * div
    else
       ux = 0.0_nr
       vy = 0.0_nr
       wz = 0.0_nr
    end if
    !
    if ( .not. lstretch_only ) then
       vor_x = wy - vz
       vor_y = uz - wx
       vor_z = vx - uy
       !
       wy = wy - 1.0_nr/2.0_nr * vor_x
       vz = vz + 1.0_nr/2.0_nr * vor_x
       uz = uz - 1.0_nr/2.0_nr * vor_y
       wx = wx + 1.0_nr/2.0_nr * vor_y
       vx = vx - 1.0_nr/2.0_nr * vor_z
       uy = uy + 1.0_nr/2.0_nr * vor_z
    else
       wy = 0.0_nr
       vz = 0.0_nr
       uz = 0.0_nr
       wx = 0.0_nr
       vx = 0.0_nr
       uy = 0.0_nr
    end if
    !
    if ( grid_cyclic_ew ) then
       ax = 1_ni
       zx = nx
    else 
       ax = 2_ni
       zx = nx-1_ni
    end if
    !
    ! Find the eigenvectors of the remaining symmetric zero-trace matrixes
    do i = ax,zx
       write(*,'(I5,A4,I5,A)', advance='no') i, 'of', nx, cr
       !
       do j = 2_ni,ny-1_ni
          do k = 2_ni,nz-1_ni
             do ti = 1_ni,nt
                ! Eigenvalues of symmetric matrixes
                jac = reshape( (/ ux(ti,k,j,i), vx(ti,k,j,i), wx(ti,k,j,i), &
                    &             uy(ti,k,j,i), vy(ti,k,j,i), wy(ti,k,j,i), &
                    &             uz(ti,k,j,i), vz(ti,k,j,i), wz(ti,k,j,i) /), (/ 3_ni, 3_ni /) )
                call dsytrd('U', 3_ni, jac, 3_ni, diag, subdiag, tau, work, 96_ni, info)
                call dorgtr('U', 3_ni, jac, 3_ni, tau, work, 96_ni, info)
                call dsteqr('V', 3_ni, diag, subdiag, jac, 3_ni, work, info)
                if (info /= 0_ni ) then
                   diag(:) = nan
                   jac(:,:) = nan
                end if
                ! The diagonal elements are overwritten by the eigenvalues,
                ! and the rotation matrix Q by the eigenvectors
                evalr(:) = diag(:)
                evec(:,:) = jac(:,:)
                !
                ! sort eigenvalues and eigenvectors by eigenvalue
                minidx = min(minloc(evalr, dim=1_ni),1_ni)
                maxidx = min(maxloc(evalr, dim=1_ni),1_ni)
                ! If all values of evalr are equal: Use index order
                if (minidx == 1_ni .and. maxidx == 1_ni) then
                   maxidx = 3_ni
                end if
                ! first: max
                n  = maxidx
                res_eval(ti,k,j,i,1_ni)   = evalr(n)
                res_evec(ti,k,j,i,1_ni,:) = evec(:,n)
                ! second: middle
                n = 6_ni - maxidx - minidx
                res_eval(ti,k,j,i,2_ni)   = evalr(n)
                res_evec(ti,k,j,i,2_ni,:) = evec(:,n)
                ! last and least
                n = minidx
                res_eval(ti,k,j,i,3_ni)   = evalr(n)
                res_evec(ti,k,j,i,3_ni,:) = evec(:,n)
             end do
          end do
       end do
    end do
    !
    res_eval(:,1_ni,:,:,:) = nan
    res_eval(:,nz,:,:,:) = nan
    res_evec(:,1_ni,:,:,:,:) = nan
    res_evec(:,nz,:,:,:,:) = nan
    !
    return
  end subroutine
  !
  !@ Calculate 3d deformation in 3d natural coordinates
  !@ 
  !@ The generalisation of deformation to 3 dimensions is done using the
  !@ analogy to Lyapunov exponents. Before calculating the eigenvectors of 
  !@ the velocity gradient tensor, this tensor is transferred to 3D
  !@ natural coordinates. The new x-direction follows the 3D wind vector, the
  !@ new y-direction perpendicular to x but in the horizontal plane. Finally,
  !@ the new z-direction follows from orthogonality.
  !@ 
  !@ To be applied on pressure levels.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ w : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Pressure vertical wind velocity.
  !@ rho : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Density of air, used to convert the pressure vertical velocity to a cartensian 
  !@     vertical velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dz : np.ndarray with shape (nz-2,ny,nx) and dtype float64
  !@     The double grid spacing in z-direction to be directly for centered differences.
  !@     ``dz(k,j,i)`` is expected to contain the z-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
  !@ lstretch_only : boolean
  !@     Calculate stretching part of 3D deformation only
  !@ lshear_only : boolean
  !@     Calculate shearing part of 3D deformation only
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z-direction.
  !@ nt : int
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shapes (nt,nz-2,ny,nx), (nt,nz-2,ny,nx,3) and dtype float64
  !@     Calculated 3D eigenvalues and eigenvectors showing the strength and
  !@     orientation of deformation in 3D.
  subroutine def_3d_nat(res_eval,res_evec,nx,ny,nz,nt,u,v,w,rho,dx,dy,dz,lstretch_only,lshear_only)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), w(nt,nz,ny,nx), & 
            &                     rho(nt,nz,ny,nx), dx(ny,nx), dy(ny,nx), dz(nt,nz,ny,nx)
    real(kind=nr), intent(out) :: res_eval(nt,nz,ny,nx,3_ni), &
                                  res_evec(nt,nz,ny,nx,3_ni,3_ni)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    logical, intent(in) :: lstretch_only, lshear_only
    !f2py depend(nx,ny,nz,nt) res, v, w
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: ux(nt,nz,ny,nx), uy(nt,nz,ny,nx), uz(nt,nz,ny,nx), &
            &        vx(nt,nz,ny,nx), vy(nt,nz,ny,nx), vz(nt,nz,ny,nx), &
            &        wx(nt,nz,ny,nx), wy(nt,nz,ny,nx), wz(nt,nz,ny,nx), &
            &        vor_x(nt,nz,ny,nx), vor_y(nt,nz,ny,nx), vor_z(nt,nz,ny,nx), &
            &        div(nt,nz,ny,nx), wm(nt,nz,ny,nx), dzm(nt,nz,ny,nx)
    real(kind=nr) :: jac(3_ni,3_ni), rot(3_ni,3_ni), evalr(3_ni), tau(2_ni), &
            &        evec(3_ni,3_ni), work(96_ni), diag(3_ni), subdiag(2_ni), ff2, ff3
    integer(kind=ni) :: i,j,k,n,ti, info, ax,zx, minidx,maxidx
    ! -----------------------------------------------------------------
    !
    if ( lshear_only .and. lstretch_only ) then
       write(*,*) 'Either lshear_only or lstretch_only or both must be False'
       stop 1
    end if
    !
    ! Crude transformation from Pa/s to m/s
    wm = -w(:,:,:,:) * rho(:,:,:,:)*g
    !
    ! Transforming vertical grid distances into m
    dzm = -dz(:,:,:,:) * rho(:,:,:,:)*g
    !
    call grad_3d(ux,uy,uz, nx,ny,nz,nt, u, dx,dy,dzm)
    call grad_3d(vx,vy,vz, nx,ny,nz,nt, v, dx,dy,dzm)
    call grad_3d(wx,wy,wz, nx,ny,nz,nt, wm, dx,dy,dzm)
    !
    if ( grid_cyclic_ew ) then
       ax = 1_ni
       zx = nx
    else 
       ax = 2_ni
       zx = nx-1_ni
    end if
    !
    ! Rotate the Jacobian into local natural coordinates
    do i = ax,zx
       do j = 2_ni,ny-1_ni
          do k = 2_ni,nz-1_ni
             do ti = 1_ni,nt
                ff2 = sqrt(u(ti,k,j,i)**2_ni + v(ti,k,j,i)**2_ni)
                ff3 = sqrt(u(ti,k,j,i)**2_ni + v(ti,k,j,i)**2_ni + w(ti,k,j,i)**2_ni)
                !
                ! Rotation matrix into local natural coordinates
                rot = reshape( (/ u(ti,k,j,i)/ff3, -v(ti,k,j,i)/ff2, -u(ti,k,j,i)/ff2*w(ti,k,j,i)/ff3, &
                    &             v(ti,k,j,i)/ff3,  u(ti,k,j,i)/ff2, -v(ti,k,j,i)/ff2*w(ti,k,j,i)/ff3, &
                    &             w(ti,k,j,i)/ff3,  0.0_nr, ff2/ff3 /), (/ 3_ni, 3_ni /) )
                !
                ! Construct the velocity gradient tensor
                jac = reshape( (/ ux(ti,k,j,i), vx(ti,k,j,i), wx(ti,k,j,i), &
                    &             uy(ti,k,j,i), vy(ti,k,j,i), wy(ti,k,j,i), &
                    &             uz(ti,k,j,i), vz(ti,k,j,i), wz(ti,k,j,i) /), (/ 3_ni, 3_ni /) )
                jac = matmul(rot, jac)
                !
                ! Saving back into the gradient arrays
                ux(ti,k,j,i) = jac(1_ni,1_ni)
                uy(ti,k,j,i) = jac(1_ni,2_ni)
                uz(ti,k,j,i) = jac(1_ni,3_ni)
                vx(ti,k,j,i) = jac(2_ni,1_ni)
                vy(ti,k,j,i) = jac(2_ni,2_ni)
                vz(ti,k,j,i) = jac(2_ni,3_ni)
                wx(ti,k,j,i) = jac(3_ni,1_ni)
                wy(ti,k,j,i) = jac(3_ni,2_ni)
                wz(ti,k,j,i) = jac(3_ni,3_ni)
             end do
          end do
       end do
    end do
    !
    ! Subtracting divergence and vorticity
    if ( .not. lshear_only ) then
       div = ux + vy + wz
       ux = ux - 1.0_nr/3.0_nr * div
       vy = vy - 1.0_nr/3.0_nr * div
       wz = wz - 1.0_nr/3.0_nr * div
    else
       ux = 0.0_nr
       vy = 0.0_nr
       wz = 0.0_nr
    end if
    !
    if ( .not. lstretch_only ) then
       vor_x = wy - vz
       vor_y = uz - wx
       vor_z = vx - uy
       !
       wy = wy - 1.0_nr/2.0_nr * vor_x
       vz = vz + 1.0_nr/2.0_nr * vor_x
       uz = uz - 1.0_nr/2.0_nr * vor_y
       wx = wx + 1.0_nr/2.0_nr * vor_y
       vx = vx - 1.0_nr/2.0_nr * vor_z
       uy = uy + 1.0_nr/2.0_nr * vor_z
    else
       wy = 0.0_nr
       vz = 0.0_nr
       uz = 0.0_nr
       wx = 0.0_nr
       vx = 0.0_nr
       uy = 0.0_nr
    end if
    !
    ! Find the eigenvectors of the remaining symmetric zero-trace matrixes
    do i = ax,zx
       write(*,'(I5,A4,I5,A)', advance='no') i, 'of', nx, cr
       !
       do j = 2_ni,ny-1_ni
          do k = 2_ni,nz-1_ni
             do ti = 1_ni,nt
                ! Construct the velocity gradient tensor
                jac = reshape( (/ ux(ti,k,j,i), vx(ti,k,j,i), wx(ti,k,j,i), &
                    &             uy(ti,k,j,i), vy(ti,k,j,i), wy(ti,k,j,i), &
                    &             uz(ti,k,j,i), vz(ti,k,j,i), wz(ti,k,j,i) /), (/ 3_ni, 3_ni /) )
                !
                ! Eigenvalues of symmetric matrixes
                call dsytrd('U', 3_ni, jac, 3_ni, diag, subdiag, tau, work, 96_ni, info)
                call dorgtr('U', 3_ni, jac, 3_ni, tau, work, 96_ni, info)
                call dsteqr('V', 3_ni, diag, subdiag, jac, 3_ni, work, info)
                if (info /= 0_ni ) then
                   diag(:) = nan
                   jac(:,:) = nan
                end if
                ! The diagonal elements are overwritten by the eigenvalues,
                ! and the rotation matrix Q by the eigenvectors
                evalr(:) = diag(:)
                evec(:,:) = jac(:,:)
                !
                ! Transforming the eigenvectors back to the original coordinates
                ff2 = sqrt(u(ti,k,j,i)**2_ni + v(ti,k,j,i)**2_ni)
                ff3 = sqrt(u(ti,k,j,i)**2_ni + v(ti,k,j,i)**2_ni + w(ti,k,j,i)**2_ni)
                !
                ! Rotation matrix back from local natural coordinates (inverse=transpose of the above!)
                rot = reshape( (/ u(ti,k,j,i)/ff3, v(ti,k,j,i)/ff3, w(ti,k,j,i)/ff3, &
                             &   -v(ti,k,j,i)/ff2, u(ti,k,j,i)/ff2, 0.0_nr, &
                             &   -u(ti,k,j,i)/ff2*w(ti,k,j,i)/ff3, -v(ti,k,j,i)/ff2*w(ti,k,j,i)/ff3, ff2/ff3 /), &
                             & (/ 3_ni, 3_ni /) )
                !
                evec(:,1_ni) = matmul(rot, evec(:,1_ni))
                evec(:,2_ni) = matmul(rot, evec(:,2_ni))
                evec(:,3_ni) = matmul(rot, evec(:,3_ni))
                !
                ! sort eigenvalues and eigenvectors by eigenvalue
                minidx = min(minloc(evalr, dim=1_ni),1_ni)
                maxidx = min(maxloc(evalr, dim=1_ni),1_ni)
                ! If all values of evalr are equal: Use index order
                if (minidx == 1_ni .and. maxidx == 1_ni) then
                   maxidx = 3_ni
                end if
                ! first: max
                n  = maxidx
                res_eval(ti,k,j,i,1_ni)   = evalr(n)
                res_evec(ti,k,j,i,1_ni,:) = evec(:,n)
                ! second: middle
                n = 6_ni - maxidx - minidx
                res_eval(ti,k,j,i,2_ni)   = evalr(n)
                res_evec(ti,k,j,i,2_ni,:) = evec(:,n)
                ! last and least
                n = minidx
                res_eval(ti,k,j,i,3_ni)   = evalr(n)
                res_evec(ti,k,j,i,3_ni,:) = evec(:,n)
             end do
          end do
       end do
    end do
    !
    res_eval(:,1_ni,:,:,:) = nan
    res_eval(:,nz,:,:,:) = nan
    res_evec(:,1_ni,:,:,:,:) = nan
    res_evec(:,nz,:,:,:,:) = nan
    !
    return
  end subroutine
  !
  !@ Calculate angle between the x-axis and isolines of a given field
  !@ 
  !@ The direction of isolines is calculated by::
  !@
  !@     k times grad(dat)
  !@
  !@ Keyser, Reeder and Reed (1988) call this angle "alpha".
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input field for the isolines.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated isoline angle.
  subroutine isoline_angle(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: gradx(nz,ny,nx), grady(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(gradx,nx,ny,nz,dat,dx,dy)
    call ddy(grady,nx,ny,nz,dat,dx,dy)  
    res = atan2(-gradx,grady)
  end subroutine
  !
  !@ Calculate angle between the axis of dilatation and isolines of a given field
  !@ 
  !@ This angle goes by different names in the literature:
  !@
  !@  * Keyser, Reeder and Reed (1988), Markowski and Richardson (2011): ``beta``
  !@  * Lapeyre, Klein and Hua (1999): ``theta + phi + pi/4``
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input field for the isolines.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated isoline angle.
  subroutine isoline_to_deformation_angle(res,nx,ny,nz,u,v,dat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: alpha(nz,ny,nx), delta(nz,ny,nx) 
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v, dat
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! beta = delta - alpha (Keyser Reeder Reed) 
    !  alpha = angle between x axis and iso-dat line (direction k X grad dat)
    !  delta = angle between x axis and dilatation axis (dynlib: def_angle)
    call isoline_angle(alpha,nx,ny,nz,dat,dx,dy)
    call def_angle(delta,nx,ny,nz,u,v,dx,dy)
    res = delta - alpha
  end subroutine
  !
  !@ Calculate the frontogenesis function
  !@
  !@ Calculates both the streching and stirring rates of horizontal gradients in a 
  !@ given field. The stretching rate is defined as::
  !@
  !@      1/|grad(dat)| * d/dt(|grad(dat)|)
  !@
  !@ and called gamma in Lapeyre, Klein and Hua (1999). This quantity alone is often
  !@ referred to as the [Petterssen] frontogenesis function 
  !@ (e.g. in Markowski and Richardson 2011). 
  !@ 
  !@ Keyser, Reeder and Reed (1988) generalise the Petterssen frontogenesis function to
  !@ its vector form, introducing the stirring rate of the gradient as its second component.
  !@ The stirring rate is defined as::
  !@
  !@      grad(dat)/|grad(dat)|^2 * (k \times d/dt(grad(dat)))
  !@ 
  !@ In the nomenclature of Lapeyre, Klein and Hua, this quantity is simply::
  !@
  !@      d(theta)/dt
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input field for the isolines.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The streching and stirring component of the generalised frontogenesis function.
  subroutine frontogenesis(resstretch,resstir,nx,ny,nz,u,v,dat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resstretch(nz,ny,nx), resstir(nz,ny,nx)
    real(kind=nr) :: bet(nz,ny,nx),totdef(nz,ny,nx),divergence(nz,ny,nx),vorticity(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resstretch,resstir,v,dat
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! (Lapeyre Klein Hua):
    !  stretch = gamma = 1/rho * d(rho)/dt (rho=|gradPV|) = -sigma/2 * sin(2(theta+phi))
    ! (Keyser Reeder Reed):
    !  d(gradPV)/dt = (Fn,Fs), coordinates rotated to dat-contours (Keyser Reeder Reed) 
    !  stretch = -1/|gradPV| * Fn 
    !       Fn = 0.5*|gradPV|(D-E*cos(2*beta))
    !  stir = Fs/|gradPV|
    !    Fs = 0.5*|gradPV|(vort+E*sin(2*beta))
    call isoline_to_deformation_angle(bet,nx,ny,nz,u,v,dat,dx,dy)
    call def_total(totdef,nx,ny,nz,u,v,dx,dy)
    call div(divergence,nx,ny,nz,u,v,dx,dy)
    call vor(vorticity,nx,ny,nz,u,v,dx,dy)
    !
    resstretch = -0.5_nr * (divergence - totdef*cos(2_nr*bet))
    resstir = 0.5_nr * (vorticity + totdef*sin(2_nr*bet))
  end subroutine
  !
  !@ Calculate terms in the frontogenesis function
  !@
  !@ Calculates the effect of different processes affecting the horizontal gradient of a 
  !@ given field. 
  !@  
  !@ 1. Divergence 
  !@ 2. Horizontal deformation
  !@ 3. Tilting/ vertical deformation
  !@ 4. Differential heating
  !@
  !@ The terms (1) and (2) relate to the divergence and deformation terms calculated in the 
  !@ subroutine :meth:`frontogenesis`.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ w : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dat : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Input field for the isolines.
  !@ heat : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Lagrangian tendency of the input field, for example the diabatic heating rate 
  !@     if dat is potential temperature.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dz : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     The double grid spacing in z-direction to be directly for centered differences.
  !@     ``dz(t,k,j,i)`` is expected to contain the z-distance between ``(t,k+1,j,i)`` and ``(t,k-1,j,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 4-tuple of np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Contributions from the four terms in the order given above.
  subroutine frontogenesis_contributors(tdiv,tdef,ttilt,theat, nx,ny,nz,nt, u,v,w,dat,heat, dx,dy,dz)
    real(kind=nr), intent(in)  :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), w(nt,nz,ny,nx), &
                     &            dat(nt,nz,ny,nx), heat(nt,nz,ny,nx),  &
                     &            dx(ny,nx), dy(ny,nx), dz(nt,nz,ny,nx)
    real(kind=nr), intent(out) :: tdiv(nt,nz,ny,nx), tdef(nt,nz,ny,nx), ttilt(nt,nz,ny,nx), theat(nt,nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !f2py depend(nx,ny,nz,nt) v,w,dat,heat, dz
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: datx(nt,nz,ny,nx), daty(nt,nz,ny,nx), datz(nt,nz,ny,nx), &
         &           absgrad(nt,nz,ny,nx), div(nt,nz,ny,nx), &
         &           def_str(nt,nz,ny,nx), def_she(nt,nz,ny,nx), &
         &           ux(nt,nz,ny,nx), uy(nt,nz,ny,nx), &
         &           vx(nt,nz,ny,nx), vy(nt,nz,ny,nx), &
         &           wx(nt,nz,ny,nx), wy(nt,nz,ny,nx), &
         &           heatx(nt,nz,ny,nx), heaty(nt,nz,ny,nx)
    integer(kind=ni) :: n
    ! -----------------------------------------------------------------
    !
    ! Used in all terms
    call grad_3d(datx,daty,datz, nx,ny,nz,nt, dat, dx,dy,dz)
    absgrad = sqrt(datx**2_ni + daty**2_ni)
    !
    do n = 1_ni,nt
       ! For horizontal kinematic terms (1) and (2)
       call grad(ux(n,:,:,:),uy(n,:,:,:), nx,ny,nz, u(n,:,:,:), dx,dy)
       call grad(vx(n,:,:,:),vy(n,:,:,:), nx,ny,nz, v(n,:,:,:), dx,dy)
       !
       ! Vertical kinematic term (3) 
       call grad(wx(n,:,:,:),wy(n,:,:,:), nx,ny,nz, w(n,:,:,:), dx,dy)
       !
       ! Diabatic effects (4)
       call grad(heatx(n,:,:,:),heaty(n,:,:,:), nx,ny,nz, heat(n,:,:,:), dx,dy)
    end do
    !
    div = ux + vy
    def_str = ux - vy
    def_she = vx + uy
    !
    tdiv  = -0.5_nr*absgrad * div
    tdef  = -1.0_nr/absgrad * (0.5_nr*def_str*(datx**2_ni - daty**2_ni) + def_she*datx*daty)
    ttilt = -1.0_nr/absgrad * (datx*datz*wx + daty*datz*wy)
    theat =  1.0_nr/absgrad * (datx*heatx + daty*heaty)
    !
  end subroutine
  !
  !@ Calculate the Okubo-Weiss parameter
  !@ 
  !@ The parameter is defined as::
  !@
  !@     1/4 (total deformation^2 - vorticity^2)
  !@
  !@ the square of the eigenvalues in Okubo's paper (assuming small divergence).
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated Okubo-Weiss parameter.
  subroutine okuboweiss(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig(nz,ny,nx), omega(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_total(sig,nx,ny,nz,u,v,dx,dy)
    call vor(omega,nx,ny,nz,u,v,dx,dy)
    res = 0.25 * (sig**2 - omega**2)
  end subroutine
  !
  !@ Calculate geostrophic velocity
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ z : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Geopotential on a pressure surface or equivalent potentials on other
  !@     surfaces.
  !@ lat : np.ndarray with shape (ny) and dtype float64
  !@     Latitude of the grid point.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated geostrophic U- and V-velocity components.
  subroutine uv_geo_from_pot(resx,resy,nx,ny,nz,z,lat,dx,dy)
    use consts!, only: pi, omega_rot
    !
    real(kind=nr), intent(in)  :: z(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    real(kind=nr) :: zx(nz,ny,nx), zy(nz,ny,nx),f(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call grad(zx,zy,nx,ny,nz,z,dx,dy)
    !
    forall(k = 1_ni:nz, i = 1_ni:nx)
       f(k,:,i) = 2.0_nr * omega_rot * sin(lat*pi/180._nr)
    end forall
    where (f == 0._nr) f = 9.E99_nr !avoid singularity calculating v_g at equator
    !
    resx = -zy/f
    resy =  zx/f
  end subroutine
  !
  !@ Calculate the geopotential from montgomery potential, potential temperature and pressure
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ mont : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Montgomery potential on an isentropic surface.
  !@ theta : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Potential temperature.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure on the isentropic surface.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated geopotential.
  subroutine geop_from_montgp(res,nx,ny,nz,mont,theta,p,dx,dy)
    use consts!, only: cp, Rl, p0
    !
    real(kind=nr), intent(in)  :: mont(nz,ny,nx), theta(nz,ny,nx), p(nz,ny,nx), &
                                & dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, theta, p
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          do k=1_ni,nz
             res(k,j,i) = mont(k,j,i) - cp*theta(k,j,i)*(p(k,j,i)/p0)**(Rl/cp)
          end do
       end do
    end do
  end subroutine
  !
  !@ Horizontal slope and 3D-gradient of geopotential on isentropic levels, calculated on pressure levels
  !@
  !@ The slope and the vertical gradient are set to NaN where they exceed 
  !@ a maximum threshold, defined in the config module.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ z : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Geopotential at pressure levels.
  !@ th: np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Potential temperature at pressure levels
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dp : np.ndarray with shape (nz-2,ny,nx) and dtype float64
  !@     The double grid spacing in p-direction to be directly for centered differences.
  !@     ``dth(k,j,i)`` is expected to contain the th-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
  !@
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Number of pressure levels.
  !@ nt : int 
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Horizontal slope of the geopotential.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     x-derivative of the geopotential.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     y-derivative of the geopotential.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     th-derivative of the geopotential.
  subroutine pres_geop_slope(slope, dzdx,dzdy,dzdp, dthdx_p,dthdy_p, dthdp, dzdth, nx,ny,nz,nt, z, th, dx,dy,dp)
    use consts !, only: nan
    use config !, only: thres_max_slope, thres_max_dzdth
    !
    real(kind=nr), intent(in) :: z(nt,nz,ny,nx), th(nt,nz,ny,nx), dx(ny,nx), dy(ny,nx), dp(nt,nz,ny,nx) !dp(nz-2_ni,ny,nx)
    real(kind=nr), intent(out) :: slope(nt,nz,ny,nx), dzdx(nt,nz,ny,nx), dzdy(nt,nz,ny,nx), dzdp(nt,nz,ny,nx)
    real(kind=nr), intent(out) :: dthdx_p(nt,nz,ny,nx), dthdy_p(nt,nz,ny,nx), dthdp(nt,nz,ny,nx), dzdth(nt,nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !local variables 
    real(kind=nr)    :: dth(nt,nz,ny,nx), dzdx_p(nt,nz,ny,nx), dzdy_p(nt,nz,ny,nx)
    integer(kind=ni) :: i
    !f2py depend(nt,nz,ny,nx) :: slope,dzdx,dzdy,dzdth, dthdx,dthdy,dthdp,dzdth, N2
    ! -----------------------------------------------------------------
    !
    ! Gradient of z on pressure levels
    call grad_3d(dzdx_p,dzdy_p,dzdp, nx,ny,nz,nt, z, dx,dy,dp)

    ! Gradient of theta on pressure levels
    call grad_3d(dthdx_p,dthdy_p,dthdp, nx,ny,nz,nt, th, dx,dy,dp)
    ! Dzdth:
    do i=1_ni,nt
       call ddz_on_q(dzdth(i,:,:,:),nx,ny,nz,z(i,:,:,:),th(i,:,:,:),dx,dy)
    end do
    !
    ! Mask unreasonable results for dzdth    
    !where (dzdth > thres_max_dzdth .or. dzdth < thres_min_dzdth)
    !    dzdth = nan
    !end where
    where (1.0_nr/dzdth < thres_min_dzdth)
        dzdth = nan
    end where  
    !
    ! Transform the derivatives to be calculated on isentropic surfaces
    dzdx  = (dzdx_p-dthdx_p*dzdth)
    dzdy  = (dzdy_p-dthdy_p*dzdth)
    !
    ! Calculate the slope 
    slope = sqrt(dzdx**2_ni + dzdy**2_ni)
    !
    ! Mask unreasonable results
    where (slope > thres_max_slope .or. slope < thres_min_slope)
        slope = nan
    end where
    !
  end subroutine
  !
  !@ Horizontal slope and 3D-gradient of geopotential on isentropic levels
  !@
  !@ The slope and the vertical gradient are set to NaN where they exceed 
  !@ a maximum threshold, defined in the config module.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ z : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Geopotential at isentropic levels.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dth : np.ndarray with shape (nz-2,ny,nx) and dtype float64
  !@     The double grid spacing in th-direction to be directly for centered differences.
  !@     ``dth(k,j,i)`` is expected to contain the th-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
  !@
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Number of isentropic levels.
  !@ nt : int 
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Horizontal slope of the geopotential.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     x-derivative of the geopotential.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     y-derivative of the geopotential.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     th-derivative of the geopotential.
  subroutine isen_geop_slope(slope, dzdx,dzdy,dzdth, nx,ny,nz,nt, z, dx,dy,dth)
    use consts !, only: nan
    use config !, only: thres_max_slope, thres_max_dzdth
    !
    real(kind=nr), intent(in) :: z(nt,nz,ny,nx), dx(ny,nx), dy(ny,nx), dth(nz-2_ni,ny,nx)
    real(kind=nr), intent(out) :: slope(nz,ny,nx), dzdx(nz,ny,nx), dzdy(nz,ny,nx), dzdth(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !f2py depend(nz,ny,nx) :: slope,dzdx,dzdy,dzdth
    ! -----------------------------------------------------------------
    !
    ! Gradient and slope
    call grad_3d(dzdx,dzdy,dzdth, nx,ny,nz,nt, z, dx,dy,dth)
    slope = sqrt(dzdx**2_ni + dzdy**2_ni)
    !
    ! Alternative calculation: dzdth = -1/(g * rho) * 1 / dthdp
    ! Problem: We can currently not evaluate dthdp on th-levels
    !call deriv(dzdth, th,'z',p,xmin,ymin,dlon,dlat,nx,ny,nz,mdv) 
    !dzdth = -1/(g * rho * dzdth)
    ! TODO: Create a vertical derivative routine simultaneously doing the vertical interpolation.
    !
    ! Mask unreasonable results
    where (slope > thres_max_slope .or. slope < thres_min_slope)
        slope = nan
    end where
    !
    where (dzdth > thres_max_dzdth .or. dzdth > thres_min_dzdth)
        dzdth = nan
    end where
    !
  end subroutine
  !
  !@ Calculate density from temperature and pressure
  !@
  !@ Uses the ideal gas law for dry air.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ T : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Density.
  subroutine rho_from_T_p(rho, nx,ny,nz, T,p)
    use consts !, only: nan
    !
    real(kind=nr), intent(in) :: T(nz,ny,nx), p(nz,ny,nx)
    real(kind=nr), intent(out) :: rho(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    !f2py depend(nz,ny,nx) :: slope,dzdx,dzdy,dzdth
    ! -----------------------------------------------------------------
    !
    rho = p / (Rl * T)
    !
  end subroutine
  !
  !@ Calculate density from temperature and pressure
  !@
  !@ Uses the ideal gas law for dry air.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ T : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Density.
  subroutine w_from_omega(w, nx,ny,nz, omega,rho)
    use consts !, only: nan
    !
    real(kind=nr), intent(in) :: omega(nz,ny,nx), rho(nz,ny,nx)
    real(kind=nr), intent(out) :: w(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    !f2py depend(nz,ny,nx) :: slope,dzdx,dzdy,dzdth
    ! -----------------------------------------------------------------
    !
    w = -omega / (g * rho)
    !
  end subroutine
  !
  !@ Estimate the radius of the local streamline through each grid point
  !@
  !@ R = ds / tan(dalpha), where ds is a distance increment in flow direction 
  !@ and dalpha is the change in wind direction along ds.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     x-wind velocity component.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     y-wind velocity component.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Local stream line radius.
  subroutine streamline_radius(r, nx,ny,nz, u, v, dx, dy)
    use consts
    !
    real(kind=nr), intent(in) :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: r(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    !f2py depend(nz,ny,nx) :: v, r
    !
    real(kind=nr) :: dd(nz,ny,nx), ddx(nz,ny,nx), ddy(nz,ny,nx), dds, ds, ff, &
       &             ones(ny,nx)
    integer(kind=ni) :: i,j,k
    ! -----------------------------------------------------------------
    !
    ! 1. Calculate wind direction
    do i = 1_ni,nx
       do j = 1_ni,ny
          ones(j,i) = 1.0_nr
          do k = 1_ni,nz
             dd(k,j,i) = atan2(-u(k,j,i), v(k,j,i))
          end do
       end do
    end do
    !
    ! 2a. Calculate gradient of wind direction in x and y directions
    call grad(ddx,ddy, nx,ny,nz, dd, ones,ones)
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k = 1_ni,nz
             ! 2b. Take periodicity into account
             if ( ddx(k,j,i) >= pi ) then
                ddx(k,j,i) = ddx(k,j,i) - 2_ni*pi
             end if
             if ( ddx(k,j,i) < -pi ) then
                ddx(k,j,i) = ddx(k,j,i) + 2_ni*pi
             end if
             if ( ddy(k,j,i) >= pi ) then
                ddy(k,j,i) = ddy(k,j,i) - 2_ni*pi
             end if
             if ( ddy(k,j,i) < -pi ) then
                ddy(k,j,i) = ddy(k,j,i) + 2_ni*pi
             end if
             !
             ! 3. Project gradient of wind direction along wind direction
             ff = sqrt(u(k,j,i)**2_ni + v(k,j,i)**2_ni)
             dds = (u(k,j,i) * ddx(k,j,i) + v(k,j,i) * ddy(k,j,i))/ff
             !
             ! 4. Calculate stream line radius
             ds = (abs(u(k,j,i)) * dx(j,i) + abs(v(k,j,i)) * dy(j,i))/ff
             r(k,j,i) = ds / tan(dds)
          end do
       end do
    end do
    !
  end subroutine
  !
  !@ Integrate hydrostasy to calculate from given surface value
  !@
  !@ Uses ideal gas law and assumes hydrostasy, such that d(phi)/d(p) = - RT/p.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ t : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     3D temperature field.
  !@ p : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     3D pressure field.
  !@ ps : np.ndarray with shape (nt,ny,nx) and dtype float64
  !@     Surface pressure.
  !@ phis : np.ndarray with shape (nt,ny,nx) and dtype float64
  !@     Surface geopotential.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z-direction.
  !@ nt : int
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Geopotential.
  subroutine z_from_hydrostasy(phi, nx,ny,nz,nt, t, p, ps, phis)
    use consts
    !
    real(kind=nr), intent(in) :: t(nt,nz,ny,nx), p(nt,nz,ny,nx), ps(nt,ny,nx), phis(ny,nx)
    real(kind=nr), intent(out) :: phi(nt,nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !f2py depend(nt,nz,ny,nx) :: p, phi
    !f2py depend(nt,ny,nx) :: ps
    !f2py depend(ny,nx) :: phis
    !
    real(kind=nr) :: pstag(ny,nx)
    integer(kind=ni) :: n,k
    ! -----------------------------------------------------------------
    !
    do n = 1_ni,nt
       ! Integrate from surface to first model level
       phi(n,nz,:,:) = phis(:,:) + Rl*t(n,nz,:,:)*(log(ps(n,:,:)) - log(p(n,nz,:,:)))
       ! Integrate upwards
       do k = nz-1_ni,1_ni,-1_ni
          pstag(:,:) = (p(n,k+1_ni,:,:) + p(n,k,:,:))/2.0_nr
          ! Assuming constant temperatures in the respective lower and upper halfs of the grid cells
          phi(n,k,:,:) = phi(n,k+1_ni,:,:) + Rl*t(n,k+1_ni,:,:)*(log(p(n,k+1_ni,:,:)) - log(pstag(:,:))) & 
                              &            + Rl*t(n,k,:,:)*(log(pstag(:,:)) - log(p(n,k,:,:)))
       end do
    end do
    !
  end subroutine
  !
  !@ Calculate the F-diagnostic of Parfitt, Czaja and Seao (2017)
  !@
  !@ F := vorticity * abs(temperature gradient) / (fcor * 0.45 K/100km)
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ t : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     3D temperature field.
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
  !@ lat : np.ndarray with shape (ny) and dtype float64
  !@     Latitude of the grid point.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The Parfitt F diagnostic
  subroutine parfitt_f(res, nx,ny,nz, t, u, v, lat, dx, dy)
    use consts
    !
    real(kind=nr), intent(in) :: t(nz,ny,nx), u(nz,ny,nx), v(nz,ny,nx), lat(ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    !f2py depend(nz,ny,nx) :: u,v,res
    !f2py depend(ny,nx) :: lat,dx,dy
    !
    real(kind=nr) :: zeta(nz,ny,nx), tx(nz,ny,nx), ty(nz,ny,nx), fcor
    real(kind=nr), parameter :: t_grad_scale = 0.45e-5
    integer(kind=ni) :: i,j,k
    ! -----------------------------------------------------------------
    !
    call vor(zeta, nx,ny,nz, u, v, dx,dy)
    call grad(tx, ty, nx,ny,nz, t, dx,dy)
    !
    do i = 1_ni,nx
       do j = 1_ni,ny
          fcor = 2.0_nr * omega_rot * sin(pi/180.0_nr * lat(j,i))
          do k = 1_ni,nz
             res(k,j,i) = zeta(k,j,i) * sqrt(tx(k,j,i)**2 + ty(k,j,i)**2) / (t_grad_scale * fcor)
          end do
       end do
    end do
    !
  end subroutine
  !
  !*****************************************************************************************
  !> author: Jacob Williams
  !  date: 7/13/2014
  !
  !  Great circle distance on a spherical body, using the Vincenty algorithm.
  !
  !# References
  !  * T. Vincenty, "[Direct and Inverse Solutions of Geodesics on the Ellipsoid
  !    with Application of Nested Equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)",
  !    Survey Review XXII. 176, April 1975.

  function great_circle_distance(r,long1,lat1,long2,lat2) result(d)

    implicit none

    real(kind=nr)            :: d        !! great circle distance from 1 to 2 [km]
    real(kind=nr),intent(in) :: r        !! radius of the body [km]
    real(kind=nr),intent(in) :: long1    !! longitude of first site [rad]
    real(kind=nr),intent(in) :: lat1     !! latitude of the first site [rad]
    real(kind=nr),intent(in) :: long2    !! longitude of the second site [rad]
    real(kind=nr),intent(in) :: lat2     !! latitude of the second site [rad]

    real(kind=nr) :: c1,s1,c2,s2,dlon,clon,slon

    !Compute aux variables:
    c1    = cos(lat1)
    s1    = sin(lat1)
    c2    = cos(lat2)
    s2    = sin(lat2)
    dlon  = long1-long2
    clon  = cos(dlon)
    slon  = sin(dlon)

    d = r*atan2( sqrt((c2*slon)**2+(c1*s2-s1*c2*clon)**2), (s1*s2+c1*c2*clon) )

  end function great_circle_distance

  function approx_distance(r,long1,lat1,long2,lat2) result(d)
    implicit none
    real(kind=nr)            :: d        !! great circle distance from 1 to 2 [km]
    real(kind=nr),intent(in) :: r        !! radius of the body [km]
    real(kind=nr),intent(in) :: long1    !! longitude of first site [rad]
    real(kind=nr),intent(in) :: lat1     !! latitude of the first site [rad]
    real(kind=nr),intent(in) :: long2    !! longitude of the second site [rad]
    real(kind=nr),intent(in) :: lat2     !! latitude of the second site [rad]

    real(kind=nr) :: x,y
    x = (long2 - long1) * cos( 0.5*(lat2+lat1) )
    y = lat2 - lat1
    d = r * sqrt( x*x + y*y )
  end function approx_distance
  !**************************************
  subroutine calc_weights(weights, nx,ny,np, x, y)
    use consts
    !
    real(kind=nr), intent(in) :: x(ny,nx), y(ny,nx)
    real(kind=nr), intent(out) :: weights(ny,nx) !,np,np
    integer(kind=ni), intent(in) :: np
    integer(kind=ni)           :: nx,ny
    !f2py depend(ny,nx) :: x,y
    !f2py depend(ny,nx) :: weights
    !
    integer(kind=ni) :: i,j,k,l
    real(kind=nr) :: dists(17,17), sigma_dist=50.0
    ! -----------------------------------------------------------------
    !
    do j = 1_ni,ny
       print *, j
       do i = 1_ni,nx
          do l=1_ni,17
             do k=1_ni,17
                if(i + k -8 >= 1 .AND. i + k - 8 <= 720 .AND. j + l - 8  >= 1 .AND. j + l - 8  <= 360) then
                   dists(k,l) = approx_distance(6371.0_nr,x(j,i)*pi/180.0_nr,y(j,i)*pi/180.0_nr, &
                                    x(j + l -8,i + k - 8)*pi/180.0_nr, y(j + l -8,i + k - 8)*pi/180.0_nr)
                end if
             end do
          end do
          weights(:,:) = exp(-dists**2.0/(sigma_dist)**2.0*1.0/2.0)     
	  weights(:,:) = weights(:,:)/sum(weights)

       end do
    end do

    !	  if((lonidx + lonidx2 -8 >= 0) & (lonidx + lonidx2 - 8 < 720) & (latidx + latidx2 - 8  >= 0) & (latidx + latidx2 - 8  < 360)):
    !	      dists[latidx2,lonidx2] = great_circle((ugrid.y[latidx,lonidx],ugrid.x[latidx,lonidx]), 
    !	      (ugrid.y[latidx + latidx2 - 8,lonidx + lonidx2 -8],
    !      ugrid.x[latidx + latidx2 - 8,lonidx + lonidx2 -8])).kilometers 
    !
  end subroutine
  !
end module

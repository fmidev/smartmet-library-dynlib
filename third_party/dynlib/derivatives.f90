! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- partial derivatives
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Discretised partial derivatives 
!@ 
!@ To be used in all dynlib routines that require the calculation of
!@ derivatives.
module derivatives
  use kind
  use config
  use consts !, only: nan
  !
  implicit none
contains
  !
  !@ Calculates partial x derivative: ddatdx = partial(dat)/partial(x)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     x-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddx_o4`, :meth:`ddx_on_q`
  subroutine ddx(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 1_ni:ny, i = 2_ni:nx-1_ni)
       res(k,j,i) = (dat(k,j,i+1_ni)-dat(k,j,i-1_ni))/dx(j,i)
    end forall
    if (grid_cyclic_ew) then
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,1_ni) = (dat(k,j,  2_ni)-dat(k,j     ,nx))/dx(j,1_ni)
          res(k,j,nx  ) = (dat(k,j  ,1_ni)-dat(k,j,nx-1_ni))/dx(j,nx)
       end forall
    else 
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,1_ni) = nan
          res(k,j,nx  ) = nan
       end forall
    end if
  end subroutine
  !
  !@ Calculates partial x derivative, keeping datq constant:
  !@
  !@    ddatdx = partial(dat)/partial(x) - partial(dat)/partial(z)*partial(z)/partial(datq)*partial(datq)/dat(x)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ datq : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array for the vertical coordinate to be kept constant
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
  !@     x-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddx`, :meth:`ddx_o4`
  subroutine ddx_on_q(res,nx,ny,nz,dat,datq,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), datq(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) datq, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! Main part of the field
    do k = 2_ni,nz-1_ni
       do j = 1_ni,ny
          do i = 2_ni,nx-1_ni
             res(k,j,i) = (dat(k,j,i+1_ni) - dat(k,j,i-1_ni))/dx(j,i) &
                  &     - (dat(k+1_ni,j,i) - dat(k-1_ni,j,i))/(datq(k+1_ni,j,i)-datq(k-1_ni,j,i)) &
                  &      *(datq(k,j,i+1_ni) - datq(k,j,i-1_ni))
          enddo
       enddo
    enddo
    ! East-west boundaries (might be periodic)
    if (grid_cyclic_ew) then
       do k = 2_ni,nz-1_ni
          do j = 1_ni,ny
             res(k,j,1_ni) = (dat(k,j,2_ni) - dat(k,j,nx))/dx(j,1_ni) &
                  &        - (dat(k+1_ni,j,i) - dat(k-1_ni,j,i))/(datq(k+1_ni,j,i)-datq(k-1_ni,j,i)) &
                  &         *(datq(k,j,2_ni) - datq(k,j,nx))
             res(k,j,nx  ) = (dat(k,j,1_ni) - dat(k,j,nx-1_ni))/dx(j,nx) &
                  &        - (dat(k+1_ni,j,i) - dat(k-1_ni,j,i))/(datq(k+1_ni,j,i)-datq(k-1_ni,j,i)) &
                  &         *(datq(k,j,1_ni) - datq(k,j,nx-1_ni))
          end do
       end do
    else 
       do k = 1_ni,nz
          do j = 1_ni,ny
             res(k,j,1_ni) = nan
             res(k,j,nx  ) = nan
          end do
       end do
    end if
    ! top/bottom boundaries
    do j = 1_ni,ny
       do i = 1_ni,nx
          res(1_ni,j,i) = nan
          res(nz  ,j,i) = nan
       end do
    end do
    !
  end subroutine
  !
  !@ Calculates partial x derivative: ddatdx = partial(dat)/partial(x)
  !@
  !@ The routine uses 4th-order centered differences. Returns NaN on first and last 
  !@ two longitudes for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     x-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddx`, :meth:`ddx_on_q`
  subroutine ddx_o4(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 1_ni:ny, i = 3_ni:nx-2_ni)
       res(k,j,i) = (-dat(k,j,i+2_ni)+8_ni*dat(k,j,i+1_ni)-8_ni*dat(k,j,i-1_ni)+dat(k,j,i-2_ni))/(6_ni*dx(j,i))
    end forall
    if (grid_cyclic_ew) then
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,nx-1_ni) = (-dat(k,j,1_ni)+8_ni*dat(k,j,nx)-8_ni*dat(k,j,nx-2_ni)+dat(k,j,nx-3_ni))/(6_ni*dx(j,1_ni))
          res(k,j,nx) = (-dat(k,j,2_ni)+8_ni*dat(k,j,1_ni)-8_ni*dat(k,j,nx-1_ni)+dat(k,j,nx-2_ni))/(6_ni*dx(j,1_ni))
          res(k,j,1_ni) = (-dat(k,j,3_ni)+8_ni*dat(k,j,2_ni)-8_ni*dat(k,j,nx)+dat(k,j,nx-1_ni))/(6_ni*dx(j,1_ni))
          res(k,j,2_ni) = (-dat(k,j,4_ni)+8_ni*dat(k,j,3_ni)-8_ni*dat(k,j,1_ni)+dat(k,j,nx))/(6_ni*dx(j,1_ni))
       end forall
    else 
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,1_ni) = nan
          res(k,j,2_ni) = nan
          res(k,j,nx-1_ni) = nan
          res(k,j,nx) = nan
       end forall
    end if
  end subroutine
  !
  !@ Calculates the second partial x derivative: b = partial^2(dat)/partial(x)^2
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     Second x-derivative of ``dat``.
  subroutine ddx2(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 1_ni:ny, i = 2_ni:nx-1_ni)
       res(k,j,i) = (dat(k,j,i+1_ni)-2_ni*dat(k,j,i)+dat(k,j,i-1_ni))/dx(j,i)**2_ni
    end forall
    if (grid_cyclic_ew) then
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,1_ni) = (dat(k,j,2_ni)-2_ni*dat(k,j,1_ni)+dat(k,j,nx     ))/dx(j,1_ni)**2_ni
          res(k,j,nx  ) = (dat(k,j,1_ni)-2_ni*dat(k,j,nx  )+dat(k,j,nx-1_ni))/dx(j,nx)**2_ni
       end forall
    else 
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,1_ni) = nan
          res(k,j,nx  ) = nan
       end forall
    end if
  end subroutine
  !
  !@ Calculates partial y derivative: ddatdy = partial(dat)/partial(y)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     y-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddy_o4`, :meth:`ddy_on_q`
  subroutine ddy(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 2_ni:ny-1_ni, i = 1_ni:nx)
       res(k,j,i) = (dat(k,j+1_ni,i)-dat(k,j-1_ni,i))/dy(j,i)
    end forall
    forall(k = 1_ni:nz, i = 1_ni:nx)
       res(k,1_ni,i)=nan
       res(k,ny,i)=nan
    end forall
  end subroutine
  !
  !@ Calculates partial x derivative, keeping datq constant:
  !@
  !@    ddatdy = partial(dat)/partial(y) - partial(dat)/partial(z)*partial(z)/partial(datq)*partial(datq)/dat(y)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ datq : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array for the vertical coordinate to be kept constant
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
  !@     x-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddy`, :meth:`ddy_o4`
  subroutine ddy_on_q(res,nx,ny,nz,dat,datq,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), datq(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) datq, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! Main part of the field
    do k = 2_ni,nz-1_ni
       do j = 2_ni,ny-1_ni
          do i = 1_ni,nx
             res(k,j,i) = (dat(k,j+1_ni,i) - dat(k,j-1_ni,i))/dx(j,i) &
                  &     - (dat(k+1_ni,j,i) - dat(k-1_ni,j,i))/(datq(k+1_ni,j,i)-datq(k-1_ni,j,i)) &
                  &      *(datq(k,j+1_ni,i) - datq(k,j-1_ni,i))
          enddo
       enddo
    enddo
    ! Northern/Southern boundaries
    do k = 1_ni,nz
       do i = 1_ni,ny
          res(k,1_ni,i) = nan
          res(k,ny  ,i) = nan
       end do
    end do
    ! top/bottom boundaries
    do j = 1_ni,ny
       do i = 1_ni,nx
          res(1_ni,j,i) = nan
          res(nz  ,j,i) = nan
       end do
    end do
    !
  end subroutine
  !
  !@ Calculates partial y derivative: ddatdy = partial(dat)/partial(y)
  !@
  !@ The routine uses 4nd-order centered differences. Returns NaN on first and last 
  !@ two latitudes.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     y-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddy`, :meth:`ddy_on_q`
  subroutine ddy_o4(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 3_ni:ny-2_ni, i = 1_ni:nx)
       res(k,j,i) = (-dat(k,j+2_ni,i)+8_ni*dat(k,j+1_ni,i)-8_ni*dat(k,j-1_ni,i)+dat(k,j-2_ni,i))/(6_ni*dy(j,i))
    end forall
    forall(k = 1_ni:nz, i = 1_ni:nx)
       res(k,1_ni,i)=nan
       res(k,2_ni,i)=nan
       res(k,ny-1_ni,i)=nan
       res(k,ny,i)=nan
    end forall
  end subroutine
  !
  !@ Calculates second partial y derivative: b = partial^2(dat)/partial(y)^2
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     Second y-derivative of ``dat``.
  subroutine ddy2(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 2_ni:ny-1_ni, i = 1_ni:nx)
       res(k,j,i) = (dat(k,j+1_ni,i)-2_ni*dat(k,j,i)+dat(k,j-1_ni,i))/dy(j,i)**2_ni
    end forall
    forall(k = 1_ni:nz, i = 1_ni:nx)
       res(k,1_ni,i)=nan
       res(k,ny,i)=nan
    end forall
  end subroutine
  !
  !@ Calculates second partial derivative in x and y directions::
  !@
  !@     b = partial^2(dat)/(partial(y)*partial(x)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude and longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     x- and y-derivative of ``dat``.
  subroutine ddxy(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 2_ni:ny-1_ni, i = 2_ni:nx-1_ni)
       res(k,j,i) = ( dat(k,j+1_ni,i+1_ni)-dat(k,j+1_ni,i-1_ni) - &
            &        (dat(k,j-1_ni,i+1_ni)-dat(k,j-1_ni,i-1_ni)) )/(dx(j,i)*dy(j,i))
    end forall
    if (grid_cyclic_ew) then
       forall(k = 1_ni:nz, j = 2_ni:ny-1_ni)
          res(k,j,1_ni) = ( dat(k,j+1_ni,2_ni)-dat(k,j+1_ni,nx) - &
               &           (dat(k,j-1_ni,2_ni)-dat(k,j-1_ni,nx)) )/(dx(j,1_ni)*dy(j,1_ni))
          res(k,j,nx)   = ( dat(k,j+1_ni,1_ni)-dat(k,j+1_ni,nx-1_ni) - &
               &           (dat(k,j-1_ni,1_ni)-dat(k,j-1_ni,nx-1_ni)) )/(dx(j,nx)*dy(j,nx))
       end forall
    else
       forall(k = 1_ni:nz, j = 2_ni:ny-1_ni)
          res(k,j,1_ni) = nan
          res(k,j,nx) = nan
       end forall
    end if
    forall(k = 1_ni:nz, i = 1_ni:nx)
       res(k,1_ni,i) = nan
       res(k,ny  ,i) = nan
    end forall
  end subroutine
  !
  !@ Calculates partial derivative in z direction ddatdz = partial(dat)/partial(z)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on top and
  !@ bottom level.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ dz : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(k,j,i)`` is expected to contain the x-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
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
  !@     z-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddz_on_q`
  subroutine ddz(res,nx,ny,nz,dat,dz)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dz(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 2_ni:nz-1_ni, j = 1_ni:ny, i = 1_ni:nx)
       res(k,j,i) = (dat(k+1_ni,j,i)-dat(k-1_ni,j,i))/dz(k,j,i)
    end forall
    forall(j = 1_ni:ny, i = 1_ni:nx)
       res(1_ni,j,i)=nan
       res(nz  ,j,i)=nan
    end forall
  end subroutine
  !
  !@ Calculates partial vertical derivative using datq as a vertical coordinate:
  !@
  !@    ddatdz = partial(dat)/partial(datq)
  !@
  !@ The routine uses 2nd-order centered differences with a rather curious discretisation.
  !@ If the vertical grid spacing is uneven, the grid emphasises the estimate from the closer
  !@ set of points over the estimate from the more distant set of points. 
  !@ 
  !@ The function returns nan on the upper and lowermost level.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ datq : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array for the vertical coordinate
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
  !@     x-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddz`
  subroutine ddz_on_q(res,nx,ny,nz,dat,datq,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), datq(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dq_upper, dq_lower
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) datq, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! Main part of the field
    do k = 2_ni,nz-1_ni
       do j = 1_ni,ny
          do i = 1_ni,nx
             dq_upper = datq(k,j,i) - datq(k+1_ni,j,i)
             dq_lower = datq(k-1_ni,j,i) - datq(k,j,i)
             res(k,j,i) = 1.0_nr/(dq_upper + dq_lower) &
                    &    * ( (dq_upper/dq_lower) * (dat(k-1_ni,j,i) - dat(k,j,i))  &
                    &      + (dq_lower/dq_upper) * (dat(k,j,i) - dat(k+1_ni,j,i)) )
          enddo
       enddo
    enddo
    ! top/bottom boundaries
    do j = 1_ni,ny
       do i = 1_ni,nx
          res(1_ni,j,i) = nan
          res(nz  ,j,i) = nan
       end do
    end do
    !
  end subroutine
  !
  !@ Calculates second partial derivative in x and y directions::
  !@
  !@     bx, by = grad(dat)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude and longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     x-derivative of ``dat``.
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     y-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`grad_3d`
  subroutine grad(resx,resy,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(resx,nx,ny,nz,dat,dx,dy)
    call ddy(resy,nx,ny,nz,dat,dx,dy)
    ! set to 0 where grad result not valid: j=1 and ny
    resx(:,1_ni,:) = nan
    resy(:,1_ni,:) = nan
    resx(:,ny,:  ) = nan
    resy(:,ny,:  ) = nan
    if (.not.(grid_cyclic_ew)) then
       ! set to 0 where grad result not valid: i=1 and nx
       resx(:,:,1_ni) = nan 
       resy(:,:,1_ni) = nan
       resx(:,:,nx  ) = nan
       resy(:,:,nx  ) = nan
    end if
  end subroutine
  !
  !@ Calculates second partial derivative in x and y directions::
  !@
  !@     bx, by, bz = grad3(dat)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude, longitude for non-cyclic grids and bottom and top levels.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Data array
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dz : np.ndarray with shape (nt,nz-2,ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(k,j,i)`` is expected to contain the x-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
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
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     x-derivative of ``dat``.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     y-derivative of ``dat``.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     z-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`grad`
  subroutine grad_3d(resx,resy,resz,nx,ny,nz,nt,dat,dx,dy,dz)
    real(kind=nr), intent(in)  :: dat(nt,nz,ny,nx), dx(ny,nx), dy(ny,nx), dz(nt,nz,ny,nx)
    real(kind=nr), intent(out) :: resx(nt,nz,ny,nx), resy(nt,nz,ny,nx), resz(nt,nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz,nt
    !f2py depend(nx,ny,nz,nt) resx, resy, resz
    !f2py depend(nx,ny) dx, dy
    !
    integer(kind=ni) :: t
    ! -----------------------------------------------------------------
    !
    do t = 1_ni,nt
       call ddx(resx(t,:,:,:),nx,ny,nz,dat(t,:,:,:),dx,dy)
       call ddy(resy(t,:,:,:),nx,ny,nz,dat(t,:,:,:),dx,dy)
       call ddz(resz(t,:,:,:),nx,ny,nz,dat(t,:,:,:),dz(t,:,:,:))
    end do
  end subroutine
  !
  !@ Calculates 2d Laplacian ldat = lap2(dat)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude and longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     Laplacian of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`antilap2`
  subroutine lap2(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! TODO: use ddx2 + ddy2 instead
    forall(k = 1_ni:nz, j = 2_ni:ny-1_ni, i = 2_ni:nx-1_ni)
       res(k,j,i) = (dat(k,j,i+1_ni)+dat(k,j,i-1_ni)-2_ni*dat(k,j,i))/dx(j,i)**2_ni + &
                    (dat(k,j+1_ni,i)+dat(k,j-1_ni,i)-2_ni*dat(k,j,i))/dy(j,i)**2_ni
    end forall
    if (grid_cyclic_ew) then
       forall(k = 1_ni:nz, j = 2_ni:ny-1_ni)
          res(k,j,1_ni) = (dat(k,j,2_ni)+dat(k,j,nx)-2_ni*dat(k,j,1_ni))/dx(j,1_ni)**2_ni + &
              (dat(k,j+1_ni,1_ni)+dat(k,j-1_ni,1_ni)-2_ni*dat(k,j,1_ni))/dy(j,1_ni)**2_ni
          res(k,j,nx) = (dat(k,j,1_ni)+dat(k,j,nx-1_ni)-2_ni*dat(k,j,nx))/dx(j,nx)**2_ni + &
                     (dat(k,j+1_ni,nx)+dat(k,j-1_ni,nx)-2_ni*dat(k,j,nx))/dy(j,nx)**2_ni
       end forall
    else
       forall(k = 1_ni:nz, j = 2_ni:ny-1_ni)
          res(k,j,1_ni) = nan
          res(k,j,nx) = nan
       end forall
    end if
    forall(k = 1_ni:nz, i = 1_ni:nx)
       res(k,1_ni,i) = nan
       res(k,ny  ,i) = nan
    end forall
  end subroutine
  !
  !@ Calculates 2d antisymmetric analogue of the Laplacian::
  !@ 
  !@    ldat = partial^2(dat)/partial(x)^2 - partial^2(dat)/partial(y)^2
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude and longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
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
  !@     Antisymmetric Laplacian of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`lap2`
  subroutine antilap2(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! TODO: use ddx2 - ddy2 instead
    forall(k = 1_ni:nz, j = 2_ni:ny-1_ni, i = 2_ni:nx-1_ni)
       res(k,j,i) = (dat(k,j,i+1_ni)+dat(k,j,i-1_ni)-2_ni*dat(k,j,i))/dx(j,i)**2_ni - &
                    (dat(k,j+1_ni,i)+dat(k,j-1_ni,i)-2_ni*dat(k,j,i))/dy(j,i)**2_ni
    end forall
    if (grid_cyclic_ew) then
       forall(k = 1_ni:nz, j = 2_ni:ny-1_ni)
          res(k,j,1_ni) = (dat(k,j,2_ni)+dat(k,j,nx)-2_ni*dat(k,j,1_ni))/dx(j,1_ni)**2_ni - &
              (dat(k,j+1_ni,1_ni)+dat(k,j-1_ni,1_ni)-2_ni*dat(k,j,1_ni))/dy(j,1_ni)**2_ni
          res(k,j,nx) = (dat(k,j,1_ni)+dat(k,j,nx-1_ni)-2_ni*dat(k,j,nx))/dx(j,nx)**2_ni - &
                     (dat(k,j+1_ni,nx)+dat(k,j-1_ni,nx)-2_ni*dat(k,j,nx))/dy(j,nx)**2_ni
       end forall
    else
       forall(k = 1_ni:nz, j = 2_ni:ny-1_ni)
          res(k,j,1_ni) = nan
          res(k,j,nx) = nan
       end forall
    end if
    forall(k = 1_ni:nz, i = 1_ni:nx)
       res(k,1_ni,i) = nan
       res(k,ny  ,i) = nan
    end forall
  end subroutine
  !
end module

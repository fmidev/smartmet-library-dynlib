! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- utility functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module utils
  use kind
  use config
  use consts !, only: nan
  !
  implicit none
contains
  !
  !@ Mirror global data in meridional direction to prepare for a Fourier transform 
  !@
  !@ Returns the data extended along complementary meridians (for fft)
  !@ For each lon, the reflected (lon+180) is attached below
  !@ so that data is periodic in x and y.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The data to be mirrorred. 
  !@     *Note*, the input data is expected to have a y-axis running from 90 degS to 90 degN,
  !@     *not* following the ERA-Interim convention. Furthermore, the grid size in 
  !@     x-direction (nx) must be even.
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
  !@ np.ndarray with shape (nz,2*ny-2,nx) and dtype float64
  !@     The mirrored data.
  subroutine mirror_y_domain(res, nx,ny,nz, dat)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,2_ni*ny-2_ni,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz,iextra
    !f2py depend(nx,ny,nz) res
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do i=1_ni,nx
          iextra=i+nx/2_ni
          if (iextra.GT.nx) then
             iextra=iextra-nx          
          end if
          do j=2_ni,ny-1_ni
             res(k,j,i) = dat(k,j,i)
             res(k,2*ny-j,i)=dat(k,j,iextra)
          end do
          res(k,1_ni,i)=dat(k,1_ni,i)
          res(k,ny,i)=dat(k,ny,i) 
       end do
    end do
  end subroutine
  !
  !@ Conservative (diffusive) smoothing in 2D, minimising the local Laplacian
  !@
  !@ The routine is taken from NCL 6.1.2, where it is called [d]filter2d and 
  !@ used in wrf_smooth_2d. It applies the same filter for all (x,y)-slices in 
  !@ a 3D (Z/t,y,x) field.
  !@
  !@ The smoothing coefficient ``smooth_coeff`` can be configured using the 
  !@ ``config`` module. Following NCL, the default is ``0.25``.
  !@
  !@ Licensed as-is, details at http://www.ncl.ucar.edu/Download/NCL_source_license.shtml .
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input data to be smoothed.
  !@ niter : int
  !@     Number of passes of the filter.
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
  !@     Smoothed data.
  subroutine smooth_xy(res,nx,ny,nz,dat,niter)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: niter
    integer(kind=ni) :: i,j,k, n, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !
    real(kind=nr) :: tmp(nz,ny,nx) ! temporary helper array to avoid changing dat
    ! -----------------------------------------------------------------
    !
    res(:,:,:) = dat(:,:,:)
    !
    do n = 1_ni,niter
       tmp(:,:,:) = res(:,:,:)
       do i = 1_ni,nx
          do j = 2_ni,ny-1_ni
             do k = 1_ni,nz
                res(k,j,i) = res(k,j,i) + smooth_coef * (tmp(k,j-1_ni,i)-2_ni*tmp(k,j,i)+tmp(k,j+1_ni,i))
             end do
          end do
       end do
       do i = 2_ni,nx-1_ni
          do j = 1_ni,ny
             do k = 1_ni,nz
                res(k,j,i) = res(k,j,i) + smooth_coef * (tmp(k,j,i-1_ni)-2_ni*tmp(k,j,i)+tmp(k,j,i+1_ni))
             end do
          end do
       end do
       if ( grid_cyclic_ew ) then
          do j = 1_ni,ny
             do k = 1_ni,nz
                res(k,j,1_ni) = res(k,j,1_ni) + smooth_coef * (tmp(k,j,nx)-2_ni*tmp(k,j,1_ni)+tmp(k,j,2_ni))
                res(k,j,nx) = res(k,j,nx) + smooth_coef * (tmp(k,j,nx-1_ni)-2_ni*tmp(k,j,nx)+tmp(k,j,1_ni))
             end do
          end do
       end if
    end do
    !
    return
  end subroutine
  !
  !@ Generic filtering in 2D
  !@
  !@ The filter is applied for all (x,y)-slices in a 3D (Z/t,y,x) field.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input data to be filtered.
  !@ filtr : np.ndarray with shape (nf,) and dtype float64
  !@     Odd length array of filter weights. The filter array should sum up to 1
  !@     if the filter is to be conservative.
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
  !@ nf : int
  !@     Length of the filter, must be odd.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Filtered data.
  subroutine filter_xy(res,nx,ny,nz,nf,dat,filtr)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), filtr(nf)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz, nf
    integer(kind=ni) :: i,j,k, n,m
    !f2py depend(nx,ny,nz) res
    !
    real(kind=nr) :: tmp(nz,ny,nx) ! temporary helper array to avoid changing dat
    ! -----------------------------------------------------------------
    !
    n = (nf - 1_ni)/2_ni
    !
    !res(:,:,:) = dat(:,:,:)
    !
    do i = 1_ni,nx
       do j = 1_ni,n
          do k = 1_ni,nz
             res(k,j,i) = nan
          end do
       end do
       do j = 1_ni+n,ny-n
          do k = 1_ni,nz
             res(k,j,i) = sum(filtr(:) * dat(k,j-n:j+n,i))
          end do
       end do
       do j = ny-n+1_ni,ny
          do k = 1_ni,nz
             res(k,j,i) = nan
          end do
       end do
    end do
    !
    tmp(:,:,:) = res(:,:,:)
    do i = 1_ni+n,nx-n
       do j = 1_ni+n,ny-n
          do k = 1_ni,nz
             res(k,j,i) = sum(filtr(:) * tmp(k,j,i-n:i+n))
          end do
       end do
    end do
    if ( grid_cyclic_ew ) then
       do j = 1_ni+n,ny-n
          do k = 1_ni,nz
             do i = 1_ni,n
                m = n-i+1_ni
                res(k,j,i) = sum(filtr(1_ni:m)*tmp(k,j,nx-m+1_ni:nx)) &
                      &    + sum(filtr(m+1_ni:)*dat(k,j,:nf-m))
             end do
             do i = nx-n+1_ni,nx
                m = nf - (i-nx+n + 1_ni)
                res(k,j,i) = sum(filtr(1_ni:m)*tmp(k,j,nx-m+1_ni:nx)) &
                      &    + sum(filtr(m+1_ni:)*tmp(k,j,:nf-m))
             end do
          end do
       end do
    else
       do j = 1_ni+n,ny-n
          do k = 1_ni,nz
             do i = 1_ni,n
                res(k,j,i) = nan
             end do
             do i = nx-n+1_ni,nx
                res(k,j,i) = nan
             end do
          end do
       end do
    end if
    !
    return
  end subroutine
  !
  !@ Generic filtering in the first dimension (typically: time)
  !@
  !@ The filter is in the time dimension of a 3D (t,y,x) field.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nt,ny,nx) and dtype float64
  !@     Input data to be filtered.
  !@ filtr : np.ndarray with shape (nf,) and dtype float64
  !@     Odd length array of filter weights. The filter array should sum up to 1
  !@     if the filter is to be conservative.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in t-direction.
  !@ nf : int
  !@     Length of the filter, must be odd.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Filtered data.
  subroutine filter_t(res,nx,ny,nz,nf,dat,filtr)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), filtr(nf)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz, nf
    integer(kind=ni) :: i,j,k, n
    !f2py depend(nx,ny,nz) res
    ! -----------------------------------------------------------------
    !
    n = (nf - 1_ni)/2_ni
    !
    !res(:,:,:) = dat(:,:,:)
    !
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k=1,n
             res(k,j,i) = nan
          end do
          do k = 1_ni+n,nz-n
             res(k,j,i) = sum(filtr(:) * dat(k-n:k+n,j,i))
          end do
          do k=nz-n+1_ni,nz
             res(k,j,i) = nan
          end do
       end do
    end do
    !
    return
  end subroutine
  !
  !@ Mark line locations on the map
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ nx : int
  !@     Grid size of the output array in y-direction.
  !@ ny : int
  !@     Grid size of the output array in y-direction.
  !@ line : np.ndarray with shape (nz,np,3) and dtype float64
  !@     Position and additional information about each point along the lines.
  !@ lineoff : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Point index for the first point of each line.
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
  !@ np : int
  !@     Maximum number of points in all lines.
  !@ nl : int
  !@     Maximum number of lines.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Mask map, indicating by 1 where lines are on the map.
  subroutine mask_lines(res,nx,ny,nz,np,nl,line,lineoff)
    real(kind=nr), intent(in) :: line(nz,np,3_ni)
    integer(kind=ni), intent(in) :: lineoff(nz,nl)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,np,nl
    !f2py depend(nz) res, lineoff
    !
    integer(kind=ni) :: i,j,k, n
    ! -----------------------------------------------------------------
    !
    res(:,:,:) = 0.0_nr
    !
    do k = 1_ni,nz
       do n = 1_ni,maxval(lineoff(k,:))
          i = nint(line(k,n,1_ni), ni)
          if (i > nx) i = i - nx
          j = nint(line(k,n,2_ni), ni)
          res(k,j,i) = 1.0_nr
       end do
    end do
  end subroutine
  !
  !@ Normalized line frequencies: Line length per area
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ line : np.ndarray with shape (nz,np,3) and dtype float64
  !@     Position and additional information about each point along the lines.
  !@ lineoff : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Point index for the first point of each line.
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
  !@ np : int
  !@     Maximum number of points in all lines.
  !@ nl : int
  !@     Maximum number of lines.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Length map, indicating the length of the line segments through each
  !@     grid point
  subroutine normalize_lines(res,nx,ny,nz,np,nl,line,lineoff,dx,dy)
    real(kind=nr), intent(in) :: line(nz,np,3_ni), dx(ny,nx), dy(ny,nx)
    integer(kind=ni), intent(in) :: lineoff(nz,nl)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,np,nl
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dy
    !f2py depend(nz) lineoff
    !
    real(kind=nr) :: toi, toj, fromi, fromj, distx, disty, dist
    real(kind=nr) :: gridsize(ny,nx)
    integer(kind=ni) :: i,j,k, m,n
    ! -----------------------------------------------------------------
    !
    res(:,:,:) = 0.0_nr
    gridsize(:,:) = abs(dx(:,:)*dy(:,:))/4.0_nr
    !
    do k = 1_ni,nz
       ! Loop over lines
       do m = 2_ni,nl
          if ( lineoff(k,m) == 0_ni ) then
             exit
          end if
          !
          ! Loop over points in line: +1 for Fortran indices, +1 for starting with second point
          do n = lineoff(k,m-1_ni)+2_ni,lineoff(k,m)
             ! Line from where to where?
             fromi = line(k,n-1_ni,1_ni)
             toi = line(k,n,1_ni)
             fromj = line(k,n-1_ni,2_ni)
             toj = line(k,n,2_ni)
             !
             ! Indices of the originating point
             i = nint(fromi, ni)
             if (i > nx) i = i - nx
             j = nint(fromj, ni)
             !
             ! Total length of the line segment
             distx = abs(toi-fromi)
             if ( distx > nx/2 ) then
                distx = nx - distx
             end if
             distx = distx*dx(j,i)/2.0
             disty = (toj-fromj)*dy(j,i)/2.0
             dist = sqrt(distx**2 + disty**2)
             !
             ! Half of the length to be counted for the originating point
             res(k,j,i) = res(k,j,i) + dist/2.0_nr
             !
             ! Indices of the terminating point
             i = nint(toi, ni)
             if (i > nx) i = i - nx
             j = nint(toj, ni)
             !
             ! Half of the length to be counted for the terminating point
             res(k,j,i) = res(k,j,i) + dist/2.0_nr
          end do
       end do
       !
       ! Normalize by size of the grid cells
       res(k,:,:) = res(k,:,:) / gridsize
    end do
  end subroutine
  !
  !@ Linear regression, for example of EOF patterns on instantaneous data
  !@
  !@ Constructs ``nn`` time series by calculating the projection of a varying
  !@ field ``datproj`` onto the given patterns ``pat``. The given data ``dat`` 
  !@ is regressed onto the different time series.
  !@
  !@ Neither the time series nor the regressed patterns are normalised.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input data to be regressed.
  !@ datproj : np.ndarray with shape (nz,nyp,nxp) and dtype float64
  !@     Input data to be projected onto the given patterns to construct the time series.
  !@ pat : np.ndarray with shape (nn,nyp,nxp) and dtype float64
  !@     Patterns to be projected onto.
  !@ mask : np.ndarray with shape (nz) and dtype bool
  !@     Which indexes along the third (typically time) dimension to take along
  !@     in the calculations.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ nxp : int
  !@     Grid size in x-direction of the patterns.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nyp : int
  !@     Grid size in y-direction of the patterns.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@ nn : int
  !@     Number of patterns.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nn,ny,nx) and dtype float64
  !@     Regressed data onto the constructed time series.
  !@ np.ndarray with shape (nn,nz) and dtype float64
  !@     Time series of the projection onto the different patterns.
  subroutine regress(res,ts,nx,nxp,nyp,ny,nz,nn,dat,datproj,pat,mask)
    real(kind=nr),intent(in) :: dat(nz,ny,nx), datproj(nz,nyp,nxp), pat(nn,nyp,nxp)
    real(kind=nr),intent(out) :: res(nn,ny,nx), ts(nn,nz)
    logical, intent(in) :: mask(nz)
    integer(kind=ni), intent(in) :: nx,nxp,ny,nyp,nz,nn
    !f2py depend(nyp,nxp) pat
    !f2py depend(nn,nz) ts
    !f2py depend(nn,ny,nx) res
    !f2py depend(nz) mask
    !
    integer(kind=ni) :: k,n
    ! -----------------------------------------------------------------
    !
    do k = 1_ni,nz
       if ( mask(k) ) then
          do n = 1_ni,nn
             ts(n,k) = sum(datproj(k,:,:)*pat(n,:,:))
             res(n,:,:) = res(n,:,:) + ts(n,k)*dat(k,:,:)
          end do
       end if
    end do
  end subroutine
  !
  !@ Make a 2d-histogram based on given linear ranges
  !@
  !@ The two data sets are connected using the same array indexes.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat1 : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data describing the first dimension in the histogram.
  !@ scale1 : np.ndarray with shape (ns1) and dtype float64
  !@     Linear scale to be used for the histogram along the first dimension.
  !@ dat2 : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data describing the second dimension in the histogram.
  !@ scale2 : np.ndarray with shape (ns2) and dtype float64
  !@     Linear scale to be used for the histogram along the second dimension.
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
  !@ ns1 : int
  !@     Length of the scale in the first dimension. Minimum 2.
  !@ ns2 : int
  !@     Length of the scale in the first dimension. Minimum 2.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ns1,ns2) and dtype float64
  !@     2d histogram of the given data set.
  !@ int
  !@     Number of data points within the given scales and hence represented in
  !@     the histogram.
  !@ int 
  !@     Number of data points that fell outside the histogram ranges.
  subroutine hist2d(res, incnt,outcnt, nx,ny,nz, dat1,scale1,ns1, dat2,scale2,ns2)
    real(kind=nr), intent(in) :: dat1(nz,ny,nx), dat2(nz,ny,nx), &
            &                    scale1(ns1), scale2(ns2)
    integer(kind=ni), intent(in) :: nx,ny,nz, ns1, ns2
    integer(kind=ni), intent(out) :: res(ns2,ns1), incnt, outcnt
    !f2py depend(nx,ny,nz) dat2
    !f2py depend(ns1,ns2) res
    !
    real(kind=nr) :: lower1, lower2, step1, step2
    integer(kind=ni) :: i,j,k, idx1, idx2
    ! -----------------------------------------------------------------
    !
    res(:,:) = 0_ni
    !
    lower1 = scale1(1_ni)
    lower2 = scale2(1_ni)
    step1  = scale1(2_ni) - scale1(1_ni)
    step2  = scale2(2_ni) - scale2(1_ni)
    !
    do k = 1_ni,nz
       do j = 1_ni,ny
          do i = 1_ni,nx
             idx1 = int((dat1(k,j,i)-lower1)/step1,ni) + 1_ni
             idx2 = int((dat2(k,j,i)-lower2)/step2,ni) + 1_ni
             !
             if ( idx1 < 1_ni .or. idx1 > ns1 .or. idx2 < 1_ni .or. idx2 > ns2 ) then
                outcnt = outcnt + 1_ni
             else 
                incnt = incnt + 1_ni
                res(idx2,idx1) = res(idx2,idx1) + 1_ni
             end if
          end do
       end do
    end do
    !
    return
  end subroutine
  !
  !@ Mask areas above given threshold and conntected to list of points
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (ny,nx) and dtype float64
  !@     Data to which the threshold is applied.
  !@ seeds : np.ndarray with shape (ns,2) and dtype int32
  !@     Array of j and i indexes of seed points.
  !@ thres : float
  !@     Minimum threshold.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ ns : int
  !@     Number of seed points.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype bool
  !@     Mask of areas connected to any of the seeds and above the given
  !@     threshold.
  subroutine mask_minimum_connect(res, nx,ny, ns, dat, seeds,thres)
    real(kind=nr), intent(in) :: dat(ny,nx), thres
    integer(kind=ni), intent(in) :: seeds(ns, 2_ni), nx,ny,ns
    logical, intent(out) :: res(ny,nx)
    !f2py depend(nx,ny) res
    !
    integer(kind=ni) :: n, i,j
    ! -----------------------------------------------------------------
    !
    res(:,:) = .false.
    do n = 1_ni,ns
       j = seeds(n,1_ni)
       i = seeds(n,2_ni)
       if ( dat(j,i) < thres ) then
          write(*,'(A18,I4,A24)') 'FATAL Error: seed ', n, ' itself below threshold.'
          stop 1
       end if
       if ( .not. res(j,i) ) then
          res(j,i) = .true.
          call minimum_connect(res, nx,ny, dat, i,j, thres)
       end if
    end do
  end subroutine
  !
  !@ Recursively find points above a threshold connected to a given point
  !@
  !@ This routine is used internally by :meth:`mask_minimum_connect` and is not
  !@ intended to be called directly.
  recursive subroutine minimum_connect(res, nx,ny, dat, i,j, thres)
    real(kind=nr), intent(in) :: dat(ny,nx), thres
    integer(kind=ni), intent(in) :: nx,ny, i,j
    logical, intent(inout) :: res(ny,nx)
    !f2py depend(nx,ny) res
    ! -----------------------------------------------------------------
    !
    ! North
    if ( j > 1_ni ) then
       if ( dat(j-1_ni,i) >= thres .and. .not. res(j-1_ni,i) ) then
          res(j-1_ni,i) = .true.
          call minimum_connect(res, nx,ny, dat, i,j-1_ni, thres)
       end if
    end if
    ! South
    if ( j < ny ) then
       if ( dat(j+1_ni,i) >= thres .and. .not. res(j+1_ni,i) ) then
          res(j+1_ni,i) = .true.
          call minimum_connect(res, nx,ny, dat, i,j+1_ni, thres)
       end if
    end if
    ! West
    if ( i > 1_ni ) then
       if ( dat(j,i-1_ni) >= thres .and. .not. res(j,i-1_ni) ) then
          res(j,i-1_ni) = .true.
          call minimum_connect(res, nx,ny, dat, i-1_ni,j, thres)
       end if
    end if
    ! East
    if ( i < nx ) then
       if ( dat(j,i+1_ni) >= thres .and. .not. res(j,i+1_ni) ) then
          res(j,i+1_ni) = .true.
          call minimum_connect(res, nx,ny, dat, i+1_ni,j, thres)
       end if
    end if
    ! Periodic boundary West
    if ( grid_cyclic_ew .and. i == 1_ni ) then
       if ( dat(j,nx) >= thres .and. .not. res(j,nx) ) then
          res(j,nx) = .true.
          call minimum_connect(res, nx,ny, dat, nx,j, thres)
       end if
    end if
    ! Periodic boundary East
    if ( grid_cyclic_ew .and. i == nx ) then
       if ( dat(j,1_ni) >= thres .and. .not. res(j,1_ni) ) then
          res(j,1_ni) = .true.
          call minimum_connect(res, nx,ny, dat, 1_ni,j, thres)
       end if
    end if
  end subroutine
  !
  !@ Connect 3D features in a fourth dimension by overlap
  !@
  !@ The routine returnes list-like arrays of connections between features, as 
  !@ well as features which are not connected. In direction of the fourth 
  !@ dimension, if feature come into being they are "born", if they cease to exist
  !@ they "die". 
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ mask : np.ndarray with shape (nt,nz,ny,nx) and dtype integer
  !@     Input data field containing object IDs
  !@ nn : int
  !@     Maximum number of features, births, deaths or connections.
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
  !@ np.ndarray with shape (nt-1,nn) and dtype int
  !@     List of features being born. If the fourth dimension is time: Features appearing in later 
  !@     time steps but not in earlier ones.
  !@ np.ndarray with shape (nt-1,nn,2) and dtype int
  !@     The size of the given area, measured in grid points.
  !@ np.ndarray with shape (nt-1,nn) and dtype int
  !@     List of dying features. E.g. features appearing in earlier time steps but not in later ones.
  subroutine connect_4d_from_3d(birth,life,death, nx,ny,nz,nt, nn, dat)
    integer(kind=ni), intent(in) :: dat(nt,nz,ny,nx), nx,ny,nz,nt, nn
    integer(kind=ni), intent(out) :: birth(nt-1_ni,nn), life(nt-1_ni,nn,2_ni), &
       &                             death(nt-1_ni,nn)
    !f2py depend(nt,nn) birth, life, death
    !
    integer(kind=ni) :: exists(nt,nn), i,j,k, m,n, t,t1, nbirth,ndeath
    logical :: found, found1
    ! -----------------------------------------------------------------
    !
    exists(:,:) = 0_ni
    !
    birth(:,:) = 0_ni
    life(:,:,:) = 0_ni
    death(:,:) = 0_ni
    !
    ! Find connections between time steps
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k = 1_ni,nz
             do t = 1_ni,nt
                ! Record which feature IDs exist
                if ( dat(t,k,j,i) > 0_ni ) then
                   found = .false.
                   do n = 1_ni,nn
                      if ( exists(t,n) == 0_ni ) then
                         found = .true.
                         exists(t,n) = dat(t,k,j,i)
                         exit
                      else if ( exists(t,n) == dat(t,k,j,i) ) then
                         found = .true.
                         exit
                      end if
                   end do
                   if ( .not. found ) then
                      write(*,*) 'Found more feature IDs then parameter nn allows.'
                      stop 1
                   end if
                end if
                !
                ! Record connection between features
                if ( t < nt ) then
                   t1 = t + 1_ni
                   if ( dat(t,k,j,i) > 0_ni .and. dat(t1,k,j,i) > 0_ni ) then
                      found = .false.
                      do m = 1_ni,nn
                         if ( life(t,m,1_ni) == 0_ni ) then
                            found = .true.
                            life(t,m,1_ni) = dat(t,k,j,i)
                            life(t,m,2_ni) = dat(t1,k,j,i)
                            exit
                         else if ( life(t,m,1_ni) == dat(t,k,j,i) .and. life(t,m,2_ni) == dat(t1,k,j,i) ) then
                            found = .true.
                            exit
                         end if
                      end do
                      if ( .not. found ) then
                         write(*,*) 'Found more feature connections then parameter nn allows.'
                         stop 1
                      end if
                   end if
                end if
             end do
          end do
       end do
    end do
    !
    ! Register what's not connected -> births and deaths
    found = .false.
    found1 = .false.
    do t = 1_ni,nt-1_ni
       ndeath = 0_ni
       nbirth = 0_ni
       !
       t1 = t + 1_ni
       do n = 1_ni,nn
          if ( exists(t,n) == 0_ni .and. exists(t1,n) == 0_ni ) exit
          !
          found = .false.
          found1 = .false.
          do m = 1_ni,nn
             if ( life(t,m,1_ni) == 0_ni ) then
                exit
             else if ( exists(t,n) > 0_ni .and. life(t,m,1_ni) == exists(t,n) ) then
                !write(*,*) 'Continuation of', exists(t,n), 'by', life(t,m,2_ni)
                found = .true.
             end if
             if ( exists(t1,n) > 0_ni .and. life(t,m,2_ni) == exists(t1,n) ) then
                found1 = .true.
             end if
          end do
          !
          !write(*,*) t, n, exists(t, n), found, found1
          !
          if ( exists(t,n) > 0_ni .and. .not. found ) then
             ndeath = ndeath + 1_ni
             if ( ndeath > nn ) then
                write(*,*) 'Found more deaths then parameter nn allows.'
                stop 1
             end if
             death(t,ndeath) = exists(t,n)
          end if
          !
          if ( exists(t1,n) > 0_ni .and. .not. found1 ) then
             nbirth = nbirth + 1_ni
             if ( nbirth > nn ) then
                write(*,*) 'Found more births then parameter nn allows.'
                stop 1
             end if
             birth(t,nbirth) = exists(t1,n)
          end if
       end do
    end do
    !
  end subroutine
  !
  !@ Label connected areas by mask input field
  !@
  !@ For a given mask, all connected areas (taking the 26 surrounding grid cells as neighbours),
  !@ will be labeled by the same integer identifer number. All grid points not belonging to any 
  !@ masked area will be labeled zero.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ mask : np.ndarray with shape (nz,ny,nx) and dtype bool
  !@     Input mask field.
  !@ gsize : np.ndarray with shape (ny,nx)
  !@     Size of the grid cells.
  !@ nn : int
  !@     Maxmimum number of features returned. Will stop with an error if too small.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype int
  !@     Connected areas labeled with a unique integer.
  !@ np.ndarray with shape (nn) and dtype int
  !@     The size of the given area, measured in grid points.
  !@ np.ndarray with shape (nn,6) and dtype int
  !@     Minimum and maximum grid point indices for the connected area, in
  !@     order: kmin, kmax, jmin, jmax, imin, imax 
  subroutine label_connected_3d(label,sizes,extent, nx,ny,nz, mask, gsize, nn)
    logical, intent(in) :: mask(nz,ny,nx)
    real(kind=nr), intent(in) :: gsize(ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nn
    integer(kind=ni), intent(out) :: label(nz,ny,nx)
    integer(kind=ni), intent(out) :: extent(nn,6)
    real(kind=nr), intent(out) :: sizes(nn)
    !f2py depend(nx,ny,nz) label
    !f2py depend(nx,ny) gsize
    !
    integer(kind=ni) :: i,j,k, cnt
    ! -----------------------------------------------------------------
    !
    cnt = 0_ni
    sizes(:) = 0_ni
    extent(:,:) = 0_ni
    label(:,:,:) = 0_ni
    !
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k = 1_ni,nz
             ! Found a new "feature"
             if ( mask(k,j,i) .and. label(k,j,i) == 0_ni ) then
                cnt = cnt + 1_ni
                if ( cnt > nn ) then
                   write(*,*) 'Found more features than allowed by input array:', nn
                   stop 1
                end if
                call label_connected_3d_single(label,sizes(cnt),extent(cnt,:), nx,ny,nz, mask, gsize, cnt, i,j,k)
             end if
          end do
       end do
    end do
    !
  end subroutine
  !
  !@ Find points belonging to a masked area and label them, using a breadth-first search 
  !@ implemented via a stack
  !@
  !@ This routine is used internally by :meth:`label_connected_2d` and is not
  !@ intended to be called directly.
  subroutine label_connected_3d_single(clabel,csize,extent, nx,ny,nz, mask, gsize, cnt, i,j,k)
    logical, intent(in) :: mask(nz,ny,nx)
    real(kind=nr), intent(in) :: gsize(ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz, cnt, i,j,k
    integer(kind=ni), intent(inout) :: clabel(nz,ny,nx)
    integer(kind=ni), intent(inout) :: extent(6)
    real(kind=nr), intent(inout) :: csize
    !f2py depend(nx,ny,nz) clabel
    !f2py depend(nx,ny) gsize
    !
    integer(kind=ni) :: si(nx*ny*nz), sj(nx*ny*nz), sk(nx*ny*nz)
    integer(kind=ni) :: nstack, ii,jj,kk, ci,cj,ck, ax,px,zx,ay,zy,az,zz
    ! -----------------------------------------------------------------
    !
    ! Initialise stack
    nstack = 1_ni
    si(nstack) = i
    sj(nstack) = j
    sk(nstack) = k
    !
    ! Record starting point
    clabel(k,j,i) = cnt
    csize = csize + gsize(j,i)
    extent(1_ni) = k
    extent(2_ni) = k
    extent(3_ni) = j
    extent(4_ni) = j
    extent(5_ni) = i
    extent(6_ni) = i
    !
    do while ( nstack > 0_ni )
       ! Remove current point from stack
       ii = si(nstack)
       jj = sj(nstack)
       kk = sk(nstack)
       nstack = nstack - 1_ni
       !
       ! Set index ranges of points to check in the neighbourhood
       ax = max(ii-1_ni, 1_ni)
       zx = min(ii+1_ni, nx)
       ay = max(jj-1_ni, 1_ni)
       zy = min(jj+1_ni, ny)
       az = max(kk-1_ni, 1_ni)
       zz = min(kk+1_ni, nz)
       !
       ! If the domain is periodic, this x index should be taken into account as well
       if ( ii == 1_ni) then
          px = nx
       else if ( ii == nx ) then
          px = 1_ni
       else 
          px = 0_ni
       end if
       !
       ! Check neighbours, record them and put them on stack
       do ci = ax,zx
          do cj = ay,zy
             do ck = az,zz
                if ( mask(ck,cj,ci) .and. clabel(ck,cj,ci) == 0_ni ) then
                   nstack = nstack + 1_ni
                   si(nstack) = ci
                   sj(nstack) = cj
                   sk(nstack) = ck
                   csize = csize + gsize(cj,ci)
                   clabel(ck,cj,ci) = cnt
                   if ( ck .lt. extent(1_ni) ) then
                      extent(1_ni) = ck
                   end if
                   if ( ck .gt. extent(2_ni) ) then
                      extent(2_ni) = ck
                   end if
                   if ( cj .lt. extent(3_ni) ) then
                      extent(3_ni) = cj
                   end if
                   if ( cj .gt. extent(4_ni) ) then
                      extent(4_ni) = cj
                   end if
                   if ( ci .lt. extent(5_ni) ) then
                      extent(5_ni) = ci
                   end if
                   if ( ci .gt. extent(6_ni) ) then
                      extent(6_ni) = ci
                   end if
               end if
             end do
          end do
       end do
       ! Respect periodic boundaries in x
       if ( grid_cyclic_ew .and. px > 0_ni ) then
          do cj = ay,zy
             do ck = az,zz
                if ( mask(ck,cj,px) .and. clabel(ck,cj,px) == 0_ni ) then
                   nstack = nstack + 1_ni
                   si(nstack) = px
                   sj(nstack) = cj
                   sk(nstack) = ck
                   csize = csize + gsize(cj,px)
                   clabel(ck,cj,px) = cnt
                   if ( ck .lt. extent(1_ni) ) then
                      extent(1_ni) = ck
                   end if
                   if ( ck .gt. extent(2_ni) ) then
                      extent(2_ni) = ck
                   end if
                   if ( cj .lt. extent(3_ni) ) then
                      extent(3_ni) = cj
                   end if
                   if ( cj .gt. extent(4_ni) ) then
                      extent(4_ni) = cj
                   end if
                   if ( px .lt. extent(5_ni) ) then
                      extent(5_ni) = px
                   end if
                   if ( px .gt. extent(6_ni) ) then
                      extent(6_ni) = px
                   end if
                end if
             end do
          end do
       end if
    end do
    !
  end subroutine
  !
  !@ Label connected areas by mask input field
  !@
  !@ For a given mask, all connected areas (taking the 8 surrounding grid cells as neighbours),
  !@ will be labeled by the same integer identifer number. All grid points not belonging to any 
  !@ masked area will be labeled zero.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ mask : np.ndarray with shape (ny,nx) and dtype bool
  !@     Input mask field.
  !@ nn : int
  !@     Maxmimum number of features returned. Will stop with an error if too small.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype int
  !@     Connected areas labeled with a unique integer.
  !@ np.ndarray with shape (nn) and dtype int
  !@     The size of the given area, measured in grid points.
  subroutine label_connected_2d(label,sizes, nx,ny, mask, nn)
    logical, intent(in) :: mask(ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nn
    integer(kind=ni), intent(out) :: label(ny,nx), sizes(nn)
    !f2py depend(nx,ny) label
    !
    integer(kind=ni) :: i,j, cnt
    ! -----------------------------------------------------------------
    !
    cnt = 0_ni
    sizes(:) = 0_ni
    label(:,:) = 0_ni
    !
    do i = 1_ni,nx
       do j = 1_ni,ny
          ! Found a new "feature"
          if ( mask(j,i) .and. label(j,i) == 0_ni ) then
             cnt = cnt + 1_ni
             if ( cnt > nn ) then
                write(*,*) 'Found more features than allowed by input array:', nn
                stop 1
             end if
             call label_connected_2d_single(label,sizes(cnt), nx,ny, mask, cnt, i,j)
          end if
       end do
    end do
    !
  end subroutine
  !
  !@ Recursively find points belonging to a masked area and label them
  !@
  !@ This routine is used internally by :meth:`label_connected_2d` and is not
  !@ intended to be called directly.
  recursive subroutine label_connected_2d_single(label,csize, nx,ny, mask, cnt, i,j)
    logical, intent(in) :: mask(ny,nx)
    integer(kind=ni), intent(in) :: nx,ny, cnt, i,j
    integer(kind=ni), intent(inout) :: label(ny,nx), csize
    !f2py depend(nx,ny) label
    ! -----------------------------------------------------------------
    !
    ! Record current point
    label(j,i) = cnt
    csize = csize + 1_ni
    !
    ! Look for neighbours
    ! -> North
    if ( j > 1_ni ) then
       if ( mask(j-1_ni,i) .and. label(j-1_ni,i) == 0_ni ) then
          call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i,j-1_ni)
       end if
    end if
    ! -> South
    if ( j < ny ) then
       if ( mask(j+1_ni,i) .and. label(j+1_ni,i) == 0_ni ) then
          call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i,j+1_ni)
       end if
    end if
    ! -> West
    if ( i > 1_ni ) then
       if ( mask(j,i-1_ni) .and. label(j,i-1_ni) == 0_ni ) then
          call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i-1_ni,j)
       end if
       ! North-West
       if ( j > 1_ni ) then
          if ( mask(j-1_ni,i-1_ni) .and. label(j-1_ni,i-1_ni) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i-1_ni,j-1_ni)
          end if
       end if
       ! South-West
       if ( j < ny ) then
          if ( mask(j+1_ni,i-1_ni) .and. label(j+1_ni,i-1_ni) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i-1_ni,j+1_ni)
          end if
       end if
    end if
    ! -> East
    if ( i < nx ) then
       if ( mask(j,i+1_ni) .and. label(j,i+1_ni) == 0_ni ) then
          call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i+1_ni,j)
       end if
       ! North-East
       if ( j > 1_ni ) then
          if ( mask(j-1_ni,i+1_ni) .and. label(j-1_ni,i+1_ni) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i+1_ni,j-1_ni)
          end if
       end if
       ! South-East
       if ( j < ny ) then
          if ( mask(j+1_ni,i+1_ni) .and. label(j+1_ni,i+1_ni) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, i+1_ni,j+1_ni)
          end if
       end if
    end if
    ! -> Periodic boundary West
    if ( grid_cyclic_ew .and. i == 1_ni ) then
       if ( mask(j,nx) .and. label(j,nx) == 0_ni ) then
          call label_connected_2d_single(label,csize, nx,ny, mask, cnt, nx,j)
       end if
       ! North-East
       if ( j > 1_ni ) then
          if ( mask(j-1_ni,nx) .and. label(j-1_ni,nx) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, nx,j-1_ni)
          end if
       end if
       ! South-East
       if ( j < ny ) then
          if ( mask(j+1_ni,nx) .and. label(j+1_ni,nx) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, nx,j+1_ni)
          end if
       end if
    end if
    ! -> Periodic boundary East
    if ( grid_cyclic_ew .and. i == nx ) then
       if ( mask(j,1_ni) .and. label(j,1_ni) == 0_ni ) then
          call label_connected_2d_single(label,csize, nx,ny, mask, cnt, 1_ni,j)
       end if
       ! North-East
       if ( j > 1_ni ) then
          if ( mask(j-1_ni,1_ni) .and. label(j-1_ni,1_ni) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, 1_ni,j-1_ni)
          end if
       end if
       ! South-East
       if ( j < ny ) then
          if ( mask(j+1_ni,1_ni) .and. label(j+1_ni,1_ni) == 0_ni ) then
             call label_connected_2d_single(label,csize, nx,ny, mask, cnt, 1_ni,j+1_ni)
          end if
       end if
    end if
  end subroutine
  !
  !@ Grow a given mask by including adjacent points N times
  !@
  !@ The masks are grown in 2D, independently of the first dimension.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ mask : np.ndarray with shape (nz,ny,nx) and dtype bool
  !@     Input mask.
  !@ nn : int
  !@     Number of interations, each growing the mask by one row of adjacent
  !@     grid points.
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
  !@ np.ndarray with shape (ny,nx) and dtype bool
  !@     Extended mask including N rows of adjacent points.
  subroutine mask_grow(res, nx,ny,nz, mask, nn)
    logical, intent(in) :: mask(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz, nn
    logical, intent(out) :: res(nz,ny,nx)
    !f2py depend(nx,ny,nz) res
    !
    logical :: tmp(nz,ny,nx)
    integer(kind=ni) :: n, i,j,k
    ! -----------------------------------------------------------------
    !
    ! Temporary copy to avoid overwriting mask
    tmp(:,:,:) = mask(:,:,:)
    res(:,:,:) = mask(:,:,:)
    !
    ! Grow the mask, one layer at a time
    do n = 1_ni,nn
       do k = 1_ni,nz
          do j = 1_ni,ny
             do i = 1_ni,nx
                if ( .not. res(k,j,i) ) then
                   res(k,j,i) = tmp(k,min(j+1_ni,ny),i) .or. tmp(k,max(j-1_ni,1_ni),i) .or. &
                              & tmp(k,j,min(i+1_ni,nx)) .or. tmp(k,j,max(i-1_ni,1_ni))
                end if
             end do
             ! Respect periodic boundaries in West and East
             if ( grid_cyclic_ew ) then
                res(k,j,1_ni) = res(k,j,1_ni) .or. tmp(k,j,nx)
                res(k,j,nx) = res(k,j,nx) .or. tmp(k,j,1_ni)
             end if
          end do
       end do
       tmp(:,:,:) = res(:,:,:)
    end do
  end subroutine
  !
  !@ Fill missing values (=NaNs) by the Poisson equation
  !@ 
  !@ The Poisson equation is iteratively solved for all missing values using 
  !@ the Successive-Over-Relaxation (SOR) algorithm to minimise the the Laplacian.
  !@
  !@ The subroutine is based on code provided by Sebastian Schemm.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float
  !@     Input data containing NaN values
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
  !@ np.ndarray with shape (nz,ny,nx) and dtype float
  !@     A copy of the input data with NaN values replaced.
  !@ float
  !@     The maximum correction during the last SOR iteration as a measure of convergence.
  subroutine fill_nan(res, conv, nx,ny,nz, dat)
    real(kind=nr), intent(in) :: dat(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx), conv
    integer(kind=ni), intent(in) :: nz,ny,nx
    !f2py depend(nx,ny,nz) res
    !
    logical :: lfill(nz,ny,nx)
    real(kind=nr) :: smooth, corr
    real(kind=nr), parameter :: omega = 1.5_nr
    integer(kind=ni) :: i,j,k, il,ir,jd,ju, iter
    integer(kind=ni), parameter :: maxiter = 300_ni
    ! -----------------------------------------------------------------
    !
    ! Create a mask array for where NaN's occur and initialise output,
    ! using zero as first guess to fill the field
    do i = 1_ni,nx
       do j = 1_ni,ny
          do k = 1_ni,nz
             lfill(k,j,i) = isnan(dat(k,j,i))
             if ( lfill(k,j,i) ) then
                res(k,j,i) = 0.0_nr
             else
                res(k,j,i) = dat(k,j,i)
             endif
          end do
       end do
    end do
    !
    ! Apply the Poisson filling - successive over-relaxation (SOR)
    iter = 0_ni
    do while ( iter < maxiter )
       conv = 0.0_nr
       !
       ! Loop over all points
       do i = 1_ni,nx
          do j = 1_ni,ny
             do k = 1_ni,nz
                ! Apply the updating only for specified points
                if ( lfill(k,j,i) ) then
                   ! Get neighbouring points - respecting periodicity in x
                   if ( grid_cyclic_ew ) then
                      il = i-1_ni
                      if ( il < 1_ni ) il = nx
                      ir = i+1_ni
                      if ( ir > nx ) ir = 1_ni
                   else
                      il = max(i-1_ni, 1_ni)
                      ir = min(i+1_ni, nx)
                   end if
                   jd = max(j-1_ni, 1_ni)
                   ju = min(j+1_ni, ny)
                   !
                   ! Update field
                   smooth = 0.25_nr * (res(k,j,il) + res(k,j,ir) + res(k,ju,i) + res(k,jd,i))
                   res(k,j,i) = omega * smooth + (1.0_nr - omega) * res(k,j,i)
                   !
                   ! Remember maximum change
                   if ( res(k,j,i) /= 0.0_nr ) then
                      corr = abs( omega * (smooth - res(k,j,i)) )
                      if ( corr > conv ) then
                         conv = corr
                      endif
                   endif
                endif
                !
             enddo
          enddo
       enddo
       iter = iter + 1_ni
    enddo
    !
  end subroutine
  !
  !@ Kernel density estimator in 2d for histograms and averages
  !@ 
  !@ For a given set of input data points, this codes estimates densities and
  !@ averages using a Gaussian kernel as smoothing.
  !@
  !@ This code is an adaptation of python code provided by Andrea Marcheggiani.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ xi : np.ndarray with shape (nn) and dtype float
  !@     The x-coordinate of the input
  !@ yi : np.ndarray with shape (nn) and dtype float
  !@     The y-coordinate of the input
  !@ zi : np.ndarray with shape (nq,nn) and dtype float
  !@     Data values to be kernel-averaged. The nq-axis is to accomodate several
  !@     indepdendent variables.
  !@ xo : np.ndarray with shape (nx) and dtype float
  !@     The x-coordinate of the output grid
  !@ yo : np.ndarray with shape (ny) and dtype float
  !@     The y-coordinate of the output grid
  !@ Lx : float
  !@     Averaging length scale in x. If unsure use 0.5 standard deviations of xi for a start.
  !@ Ly : float
  !@     Averaging length scale in y. If unsure use 0.5 standard deviations of yi for a start.
  !@ maskbelow : float
  !@     Mask results where the kernel density average is below this value. If unsure start with 1.0e-6.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nn : int
  !@     Number of data points in the dataset
  !@ nq : int
  !@     Number of independent values for which the kernel average is estimated.
  !@ nx : int
  !@     Grid size of the output grid in x-direction.
  !@ ny : int
  !@     Grid size of the output grid in y-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype float
  !@     Kernel density estimate (corresponds to histogram)
  !@ np.ndarray with shape (nq,ny,nx)
  !@     Kernel-averaged input values
  subroutine kernel_average_2d_fortran(kdens, kavg, nx,ny,nn,nq, xi,yi,zi, xo,yo, Lx, Ly, maskbelow)
    real(kind=nr), intent(in) :: xi(nn), yi(nn), zi(nq,nn), xo(nx), yo(ny), Lx, Ly, maskbelow
    real(kind=nr), intent(out) :: kdens(ny,nx), kavg(nq,ny,nx)
    integer(kind=ni), intent(in) :: nq,nn,ny,nx
    !f2py depend(nq,ny,nz) kavg
    !f2py depend(ny,nz) kdens
    !
    real(kind=nr) :: wgt, total_wgt, xc, yc
    integer(kind=ni) :: m, n, j, i
    ! -----------------------------------------------------------------
    !
    kdens(:,:) = 0.0_nr
    kavg(:,:,:) = 0.0_nr
    !
    ! For all points in the output grid
    do i = 1_ni,nx
       xc = xo(i)
       do j = 1_ni,ny
          yc = yo(j)
          total_wgt = 0.0_nr
          !
          ! For all records in the input
          do n = 1_ni,nn
             wgt = exp( -0.5_nr * ( ((xi(n)-xc)/Lx)**2.0_nr + ((yi(n)-yc)/Ly)**2.0_nr) )
             kdens(j,i) = kdens(j,i) + wgt
             total_wgt = total_wgt + wgt
             do m = 1_ni,nq
                kavg(m,j,i) = kavg(m,j,i) + wgt * zi(m,n)
             end do
          end do
          !
          if ( total_wgt .ge. maskbelow ) then
             do m = 1_ni,nq
                kavg(m,j,i) = kavg(m,j,i) / total_wgt
             end do
          else
             kdens(j,i) = nan
             do m = 1_ni,nq
                kavg(m,j,i) = nan
             end do
          end if
       end do
    end do
    !
  end subroutine
  !
end module

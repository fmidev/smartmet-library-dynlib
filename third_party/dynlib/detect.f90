! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- Feature detection algorithms
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Detection routines for various atmospheric features.
!@ 
!@ Instead of just calculation a diagnostic, these routines detect features
!@ in the atmosphere as objects.
module detect
  use kind
  use config
  use derivatives
  use diag
  !
  implicit none
contains
  !
  !@ Identify Rossby wave breaking by local gradient reversals
  !@
  !@ At each (i,j,k) grid point, finds the reversals of pv y-gradient and classes them as
  !@ either cyclonic (c) or anticyclonic (a).
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ pv : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     PV (or any other suitable field) to be used to detect gradient reversals.
  !@ mask : np.ndarray with shape (nz,ny,nx) and dtype int8
  !@     Whether or not to test for a gradient reveral at a given point. Only points 
  !@     where mask == 1 are tested.
  !@ latitudes : np.ndarray with shape (ny) and dtype float64
  !@     Latitudes describing the ``pv`` array.
  !@ ddythres : float
  !@     Cutoff y-gradient for pv. The magnitude of (negative) d(pv)/dy must be above 
  !@     ddythres for reversal to be detected.
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
  !@ np.ndarray with shape (nz,ny,nx) and dtype int8
  !@     Flag for anticyclonic reversals (threshold test applied).
  !@ np.ndarray with shape (nz,ny,nx) and dtype int8
  !@     Flag for cyclonic reversals (threshold test applied).
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Absolute gradient anticyclonic reversals (threshold test applied).
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Absolute gradient cyclonic reversals (threshold test applied).
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Absolute y-gradient anticyclonic reversals (no threshold test applied).
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Absolute y-gradient cyclonic reversals (no threshold test applied).
  !@ np.ndarray with shape (nz,ny,nx) and dtype int8
  !@     Flag for points testes for reversals.
  subroutine rwb_by_grad_rev(resa,resc,resai,resci,resaiy,resciy,tested,nx,ny,nz, &
       pv,mask,latitudes,ddythres,dx,dy)
    real(kind=nr), intent(in)  :: pv(nz,ny,nx), latitudes(ny), ddythres, dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resai(nz,ny,nx), resci(nz,ny,nx), & 
       resaiy(nz,ny,nx), resciy(nz,ny,nx)
    integer(kind=1), intent(out) :: resa(nz,ny,nx), resc(nz,ny,nx), tested(nz,ny,nx)
    integer(kind=1), intent(in) :: mask(nz,ny,nx)
    real(kind=nr) :: gradx, grady, gradmag, gradang, hemi, sig, neg, tempval,edge
    real(kind=nr) :: xgradient(nz,ny,nx), ygradient(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) resa,resc,resai,resci,resaiy,resciy,tested, mask,xgradient,ygradient
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) latitudes
    ! -----------------------------------------------------------------
    !
    call grad(xgradient,ygradient,nx,ny,nz,pv,dx,dy) ! returns grad at all testable points, 0,0 elsewhere
    !
    do i=1_ni,nx
       do j=1_ni,ny
          if (latitudes(j) > 0._nr) then !NH
             hemi = 1._nr
          else if (latitudes(j) < 0._nr) then !SH
             hemi = -1._nr
          else !equator
             hemi = 0._nr
          end if
          if (abs(latitudes(j))==90._nr) then 
             edge = 1._nr !at edge of grid, must skip
          elseif ((.NOT.(grid_cyclic_ew)).AND.((i==1).OR.(i==nx))) then
             edge=1._nr !at edge of grid, must skip
          else
             edge=0._nr
          end if
          do k=1_ni,nz
             !initialise all to 0, change later if reversal
             tested(k,j,i) = 0 
             resa(k,j,i)   = 0
             resai(k,j,i)  = 0._nr
             resaiy(k,j,i) = 0._nr
             resc(k,j,i)   = 0
             resci(k,j,i)  = 0._nr
             resciy(k,j,i) = 0._nr
             if ((edge==0._nr).AND.(mask(k,j,i)==1)) then ! not on grid edge, and high enough: test it
                tested(k,j,i) = 1
                gradx=xgradient(k,j,i)
                grady=ygradient(k,j,i)
                if (grady<0._nr) then 
                   neg=1._nr !y-gradient is negative
                else
                   neg=0._nr
                end if
                if (grady<(-ddythres)) then  
                   sig=1._nr ! ygradient is significantly negative 
                else
                   sig=0._nr
                end if
                gradmag=sqrt(gradx**2.+grady**2.)
                if (gradmag .NE. 0._nr) then
                   ! angle exists if gradmag not 0
                   gradang=atan2(grady,gradx) !NB slow! faster to do it by signs
                   if ((gradang+3.141592/2.)*hemi>0._nr) then ! CW NH or CCW SH -> ACYC
                      resa(k,j,i)=int(sig)                 ! flag significant ygradient reversal
                      resai(k,j,i)=gradmag*sig        ! |gradPV| of significant reversal
                      resaiy(k,j,i)=-grady*neg        ! sum |dPV/dy| of all reversals
                   else ! CCW NH or CCW SH -> CYC
                      resc(k,j,i)=int(sig)  
                      resci(k,j,i)=gradmag*sig
                      resciy(k,j,i)=-grady*neg
                   end if
                end if
             end if  ! end test
          end do !k
       end do !j
    end do !i
  end subroutine
  !
  !@ RWB detection by contour tracking, adapted from algorithm by Gwendal Riviere
  !@
  !@ Detects the occurrence of anticyclonic and cyclonic wave-breaking 
  !@ events from a PV field on isentropic coordinates.
  !@
  !@ The implementation assumes to be given a grid with longitudes starting at
  !@ -180, and with equal longitudial as latitudinal resolution (dlon == dlat). 
  !@
  !@ Reference: Riviere (2009, hereafter R09): Effect of latitudinal 
  !@ variations in low-level baroclinicity on eddy life cycles and upper-
  !@ tropospheric wave-breaking processes. J. Atmos. Sci., 66, 1569-1592.
  !@ See the appendix C.
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ pv_in : np.ndarray with shape (nz,ny,nx) and float64
  !@     Isentropic PV. Should be on a regular lat-lon grid and 180W must be the 
  !@     first longitude (If 180W is not the first longitude, the outputs will have 
  !@     180W as the first, so must be rearranged).
  !@ variable : string of length 2
  !@     Abbreviated variable, 'PT' for potential temperature on PV2, 'PV' for
  !@     potential vorticity on an isentropic level and 'Z2' for 250 hPa geopotential
  !@ lonvalues : np.ndarray with shape (nx) and dtype float64
  !@     Longitudes describing ``pv_in``.
  !@ latvalues : np.ndarray with shape (ny) and dtype float64
  !@     Latitudes describing ``pv_in``.
  !@ ncon : int 
  !@     Number of contours to test, normally 41 or 21.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype int32 
  !@     Flag array, the value 1 indicates anticyclonic wave breaking 
  !@ np.ndarray with shape (nz,ny,nx) and dtype int32 
  !@     Flag array, the value 1 indicates cyclonic wave breaking 
  subroutine rwb_by_contour(beta_a_out,beta_c_out,nx,ny,nz,pv_in,variable,lonvalues,latvalues,ncon,dx,dy)
    use detect_rwb_contour
    !
    implicit none
    !
    real(kind=nr), intent(in) :: pv_in(nz,ny,nx),lonvalues(nx),latvalues(ny),dx(ny,nx),dy(ny,nx)
    integer(kind=ni), intent(in) :: ncon
    integer(kind=ni), intent(out) :: beta_a_out(nz,ny,nx), beta_c_out(nz,ny,nx)
    integer(kind=ni) :: beta_a_outint(nz,ny,nx), beta_c_outint(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) beta_a_out, beta_c_out,beta_a_outint, beta_c_outint
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nx) lonvalues
    !f2py depend(ny) latvalues
    ! -----------------------------------------------------------------
    !
    integer(kind=ni) :: nlon, nlat, ntime, nplev
    integer(kind=ni) :: i,it,sd,j,l,nn,idecal
    !
    real(kind=nr), dimension(:,:,:,:), allocatable :: pv,pvint
    real(kind=nr), dimension(:), allocatable :: rdate
    !
    character :: variable*2
    character :: cx*18,cy*18
    integer(kind=ni), dimension(:,:,:,:), allocatable :: beta_a,beta_c,beta_cint
    real(kind=nr), dimension(:,:,:), allocatable :: xvalf,yvalf
    real(kind=nr) :: incr,incr2,incrlat,incrlat2
    integer(kind=ni) :: nlon2,nlat2,choicegrid
    real(kind=nr), dimension(:), allocatable :: lonval2,latval2,lonval,latval
    real(kind=nr), dimension(:,:,:,:), allocatable :: pv2
    real(kind=nr), dimension(:,:,:), allocatable :: gamma_a,gamma_c
    real(kind=nr), dimension(:,:), allocatable :: moy_ga,moy_gc
    real(kind=nr), dimension(:,:,:), allocatable :: gamma_a1,gamma_c1
    real(kind=nr), dimension(:,:), allocatable :: moy_ga1,moy_gc1
    real(kind=nr),dimension(:,:),allocatable :: zon_clim_u
    real(kind=nr),dimension(:,:,:),allocatable :: clim_u
    !
    real(kind=nr) :: lonleft,latleft,lonright,latright,xlength,ylength
    !
    nlon = nx
    nlat = ny
    ntime = nz
    nplev = 1
    incr=lonvalues(2)-lonvalues(1)
    incrlat=latvalues(2)-latvalues(1)
    !
    ! CHANGING THE ORDER OF THE VALUES TO GET AN INCREASE IN LATITUDES 
    ! AND TO GET -180._nr0 AS THE FIRST LONGITUDE    
    !
    ! In the case incrlat<0 that is if the latitudes (latvalues) decrease with the
    ! corresponding index, we have to reverse the order to orientate from the South
    ! to the North
    !
    allocate(lonval(nlon))
    allocate(latval(nlat))
    allocate(pv(nlon,nlat,nplev,ntime))
    allocate(pvint(nlon,nlat,nplev,ntime))
    !
    ! Arrange pv into format required for this routine
    do it=1,ntime
       nn=1 ! only 1 level
       do j=1,nlat
          do i=1,nlon
             pv(i,j,nn,it) = pv_in(it,j,i)
          end do
       end do
    end do
    !
    pvint=pv
    lonval=lonvalues
    latval=latvalues
    !
    if (incrlat.lt.0._nr) then
       do j=1,nlat
          pvint(:,j,:,:)=pv(:,nlat-j+1,:,:)
          latval(j)=latvalues(nlat-j+1)
       end do
    end if
    !
    ! In the case, lonvalues(1) different from -180.0, we translate the field
    ! to have the first longitude equal to -180.0
    pv=pvint
    idecal=int((180._nr+lonvalues(1))/incr)
    !
    if (lonvalues(1).ne.(-180._nr)) then
       do i=1,nlon
          if (i.le.idecal) then
             pv(i,:,:,:)=pvint(i+nlon-idecal,:,:,:)
             if (lonvalues(i+nlon-idecal).lt.180._nr) then
                lonval(i)=lonvalues(i+nlon-idecal)
             else
                lonval(i)=lonvalues(i+nlon-idecal)-360._nr
             end if
          else
             pv(i,:,:,:)=pvint(i-idecal,:,:,:)
             if (lonvalues(i-idecal).lt.180._nr) then
                lonval(i)=lonvalues(i-idecal)
             else
                lonval(i)=lonvalues(i-idecal)-360._nr
             end if
          end if
       end do
    end if
    !
    ! DEFINING THE SPATIAL GRID TO WHICH THE WAVE-BREAKING
    ! DETECTION ALGORITHM WILL BE APPLIED
    !
    ! If choicegrid=1 then INCR2=INCR (we keep all the points of the original grid)
    ! If choicegrid=2 then INCR2=2*INCR (we take one grid point over two)
    ! If choicegrid=3 then INCR2=3*INCR (we take one grid point over three) this is the 
    ! choice made in Michel and Riviere (2012). Starting from a 1.5x1.5 grid (ERA40),
    ! we applied our algorithm to a 4.5x4.5 grid
    choicegrid=1
    !
    nlon2=1  ! One more point to get the periodicity in longitudes 
    do i=1,nlon,choicegrid
       nlon2=nlon2+1
    end do
    nlat2=0
    do j=1,nlat,choicegrid
       nlat2=nlat2+1
    end do
    incr2=incr*choicegrid
    incrlat2=abs(incrlat)*choicegrid
    !
    allocate(lonval2(nlon2))
    lonval2=0._nr
    do i=1,nlon2
       lonval2(i)=lonval(1)+incr2*real(i-1)
    end do
    !
    allocate(latval2(nlat2))
    latval2=0._nr
    do j=1,nlat2
       latval2(j)=latval(1)+incrlat2*real(j-1)
    end do
    ! 
    allocate(pv2(nlon2,nlat2,nplev,ntime))
    do it=1,ntime
       do nn=1,nplev
          do i=1,nlon,choicegrid
             do j=1,nlat,choicegrid     
                if (choicegrid.eq.1) then
                   pv2(i,j,nn,it)=pv(i,j,nn,it)
                else
                   pv2(int(i/choicegrid)+1,int(j/choicegrid)+1,nn,it)=pv(i,j,nn,it)
                end if
             end do
          end do
       end do
    end do
    !
    ! periodicity in longitudes (the first and last longitudes are the same)
    do it=1,ntime
       do nn=1,nplev
          do j=1,nlat,choicegrid
             if (choicegrid.eq.1) then
                pv2(nlon2,j,nn,it)=pv(1,j,nn,it)
             else
                pv2(nlon2,int(j/choicegrid)+1,nn,it)=pv(1,j,nn,it)
             end if
          end do
       end do
    end do
    deallocate(pv)
    !
    ! Computation of the beta_a and beta_c as defined in Riviere (2009) and Riviere et al (2010)
    ! beta_a=1 if the grid point belongs to an anticyclonic breaking (AWB) and 0 otherwise
    ! beta_c=1 if the grid point belongs to a cyclonic breaking (CWB) and 0 otherwise
    !
    allocate(beta_a(nlon2,nlat2,nplev,ntime))
    allocate(beta_c(nlon2,nlat2,nplev,ntime))
    beta_a=0
    beta_c=0
    !
    allocate(xvalf(8000,nlat2,ncon))
    allocate(yvalf(8000,nlat2,ncon))
    xvalf=0._nr
    yvalf=0._nr
    !print *,'before detection'
    !print *,'variable=',variable
    !print *,'nlon2=',nlon2,'nlat2=',nlat2,'ncon=',ncon,' incr2=',incr2, ' incrlat2=', incrlat2
    !print *,'lonval2(1)=',lonval2(1),'latval2(1)=',latval2(1)
    !
    do it=1,ntime
       !print *,'time=', it
       do nn=1,nplev
          xvalf=0._nr
          yvalf=0._nr
          call wb_detection(variable,pv2(:,:,nn,it),nlon2,nlat2,ncon,&
                       &    incr2,incrlat2,lonval2,latval2,&
                       &    xvalf,yvalf,beta_a(:,:,nn,it),beta_c(:,:,nn,it))
       end do
    end do
    !
    ! In the Northern Hemisphere, the anticyclonic wave breaking is part
    ! of the contour which is locally overturned and oriented NE-SW while 
    ! in the Southern Hemisphere it is oriented SE-NW. 
    ! The output fields of the above routine WB_DETECTION beta_a and beta_c 
    ! correspond respectively to NE-SW and SE-NW oriented segments over the whole earth.
    ! Therefore, we have to reverse beta_a and beta_c in the Southern Hemisphere
    ! in order to have beta_a equivalent to an anticyclonic breaking and
    ! beta_c to a cyclonic wave-breaking over the whole earth.
    !
    allocate(beta_cint(nlon2,nlat2,nplev,ntime))
    beta_cint=beta_c
    do j=1,nlat2
       if (latval2(j).lt.0._nr) then
          beta_c(:,j,:,:)=beta_a(:,j,:,:)
          beta_a(:,j,:,:)=beta_cint(:,j,:,:)
       end if
    end do
    !
    ! allocate dynlib output arrays (NOTE!! ASSUMES CHOICEGRID=1!!)
    do it=1,ntime
       nn=1 ! only 1 level
       do j=1,nlat
          do i=1,nlon
             beta_a_out(it,j,i)=beta_a(i,j,nn,it)
             beta_c_out(it,j,i)=beta_c(i,j,nn,it)
          end do
       end do
    end do
    !
    ! If latitudes were reversed for calculation, reverse them back for outputs
    if (incrlat.lt.0._nr) then
       do j=1,nlat
          beta_a_outint(:,j,:)=beta_a_out(:,nlat-j+1,:)
          beta_c_outint(:,j,:)=beta_c_out(:,nlat-j+1,:)
       end do
    end if
    beta_a_out=beta_a_outint
    beta_c_out=beta_c_outint
    !
    ! Note: Need to rearrange lons if starting lon not -180 - not done yet!
    !
    deallocate(lonval) 
    deallocate(latval) 
    deallocate(pvint) 
    deallocate(lonval2) 
    deallocate(latval2) 
    deallocate(pv2) 
    deallocate(beta_a) 
    deallocate(beta_c) 
    deallocate(xvalf) 
    deallocate(yvalf) 
    deallocate(beta_cint) 
  end subroutine ! end of subroutine contour_rwb
  !
  !@ Blocking indicator using large-scale gradient reversals
  !@
  !@ Base data can for example be potential temperature on the PV2-surface (following 
  !@ e.g. Masato et al. 2012).
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     PV (or any other suitable field) to be used to detect gradient reversals.
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
  !@     Local block intensity indicator field. Maxima in this field are blocking centers 
  !@     that can be tracked in time and space to verify the minimum persistence and 
  !@     stationarity requirements.
  subroutine block_indicator_grad_rev(res, nx,ny,nz, dat, dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: tmp(nz,ny,nx)
    integer(kind=ni) :: i,j,k
    ! -----------------------------------------------------------------
    !
    res(:,:,:) = 0.0_nr
    !
    do k = 1_ni,nz
       ! Local blocking index
       do j = ny-block_dj,block_dj+1_ni,-1_ni
          do i = 1_ni,nx
             tmp(k,j,i) = (sum(dat(k,j-block_dj:j-1_ni,i)) - sum(dat(k,j+1_ni:j+block_dj,i))) / block_dj
          end do
       end do
       ! Longitudinal smoothing by running mean
       do j = ny-block_dj,block_dj+1_ni,-1_ni
          do i = 1_ni+block_di,nx-block_di
             res(k,j,i) = sum(tmp(k,j,i-block_di:i+block_di))/(2_ni*block_di+1_ni)
          end do
       end do
    end do
    !
    ! Taking periodic grid into account for longitudinal smoothing by running mean
    if (grid_cyclic_ew) then
       do k = 1_ni,nz
          do j = ny-block_dj,block_dj+1_ni,-1_ni
             do i = 1_ni,block_di
                res(k,j,i) = (sum(tmp(k,j,nx-(block_di-i):nx)) + sum(tmp(k,j,1_ni:i+block_di)) )/(2_ni*block_di+1_ni)
             end do
             do i = nx-block_di+1_ni,nx
                res(k,j,i) = (sum(tmp(k,j,i-block_di:nx)) + sum(tmp(k,j,1_ni:i-(nx-block_di))) )/(2_ni*block_di+1_ni)
             end do
          end do
       end do
    end if
    !
  end subroutine
  !
  !@ Cyclone masks by finding the outermost closed contour
  !@
  !@ Fortran part of a reimplementation of the Wernli and Schwierz (2006) 
  !@ algorithm, including the modifications described in Sprenger et al. (2017). 
  !@ To avoid the technically difficult contour tracing, we base the detection on 
  !@ sorted sea-level pressure values. 
  !@
  !@ For more info refer to the python function. This routine is not intended to
  !@ be called directly.
  !@
  !@ The algorithm uses the following configration parameters
  !@  - `dynfor.config.cyc_minsize`: Minimum size in km^2
  !@  - `dynfor.config.cyc_maxsize`: Maximum size in km^2
  !@  - `dynfor.config.cyc_maxoro`: Maximum orographic height in the units of
  !@    the `oro` parameter to this function.
  !@  - `dynfor.config.cyc_mindist`: Minimum distance between two cyclone centres in km
  !@  - `dynfor.config.cyc_minprominence`: Minimum prominence of the sea-level
  !@    pressure minimum in the unit of the `msl`/`msls` parameters to this function.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ nn : int
  !@     Maximum number of cyclones to be detected
  !@ msl : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Sea-level pressure (or any other suitable field).
  !@ msls : np.ndarray with shape (nz,ny*nx) and dtype float64
  !@     Sorted sea-level pressure values for each time step.
  !@ iis : np.ndarray with shape (nz,ny*nx) and dtype float64
  !@     Longitudinal grid index of the respective sea-level pressure value.
  !@ jjs : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Latitudinal grid index of the respective sea-level pressure value.
  !@ oro : np.ndarray with shape (ny,nx) and dtype float64
  !@     Orographic height to mask sea-level pressure minima over high terrain.
  !@ lon : np.ndarray with shape(nx)
  !@     Longitudes in degrees
  !@ lat : np.ndarray with shape(ny)
  !@     Latitudes in degrees
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
  !@ np.ndarray with shape (nz,ny,nx) and dtype int32
  !@     Cyclone mask with different integers designating different cyclones
  !@ np.ndarray with shape (nz,nn,5)
  !@     Meta data about each cyclone: (1) latitutde and (2) longitude of its centre,
  !@     (3) minimum SLP, (4) SLP at the outermost contour, and (5) cyclone size.
  subroutine cyclone_by_contour_fortran(mask,meta, nx,ny,nz,nn, msl, msls,iis,jjs, oro, lon,lat, dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: msl(nz,ny,nx), msls(nz,ny*nx), oro(ny,nx),  &
       &                          lon(nx), lat(ny), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: meta(nz,nn,5_ni)
    integer(kind=ni), intent(in) :: iis(nz,ny*nx), jjs(nz,ny*nx), nn,nx,ny,nz
    integer(kind=ni), intent(out) :: mask(nz,ny,nx)
    !f2py depend(nx,ny,nz) mask, msls, iis, jjs
    !f2py depend(nn,nz) meta
    !f2py depend(nx,ny) oro, dx, dy
    !f2py depend(nx) lon
    !f2py depend(ny) lat
    !
    real(kind=nr) :: cellsize(ny,nx), min_around, dist, mindist, dlon, lat1r, lat2r, val, prominence
    integer(kind=ni) :: loc_mask(8_ni), cycidx, ocycidx, minidx, ndiffloc, nmaxloc
    logical :: done(nn)
    integer(kind=ni) :: i,j,k,n,m, ip1,im1,jp1,jm1
    ! -----------------------------------------------------------------
    !
    ! Initialize stuff
    meta(:,:,:) = 0_ni
    !
    do i = 1_ni,nx
       do j = 1_ni,ny
          cellsize(j,i) = abs(dx(j,i) * dy(j,i)) / 4.0e6 ! double grid spacing in dx/dy and conversion from m^2 to km^2
          if ( isnan(cellsize(j,i)) ) then
             cellsize(j,i) = 0.0_nr
          end if
       end do
    end do
    !
    ! Loops over time and sorted arrays
    do k = 1,nz
       cycidx = 0_ni
       done(:) = .false.
       !
       do n = 1,ny*nx
          i = iis(k,n) + 1_ni ! conversion to Fortran indexes
          j = jjs(k,n) + 1_ni
          val = msls(k,n)
          !
          call find_min_around(val, k, j, i, nx, ny, msl, mask, min_around, loc_mask)
          ! How many different mask values in the vicinity?  0, 1, 2-or-more?
          ndiffloc = 0_ni
          nmaxloc = maxval(loc_mask)
          do m = 1_ni,8_ni
             if ( loc_mask(m) == 0_ni ) then
                cycle
             end if

             if ( loc_mask(m) == nmaxloc ) then
                ndiffloc = max(ndiffloc,1_ni)
             else
                ndiffloc = max(ndiffloc,2_ni)
             end if
          end do
          !
          ! Nothing around -> potential new cyclone center
          if ( minval(loc_mask) == 0_ni .and. maxval(loc_mask) == 0_ni ) then
             ! If SLP minimum and not above too high orography
             if ( val < min_around .and. oro(j,i) < cyc_maxoro ) then
                ! Check proximity to previously detected cyclone centres
                mindist = cyc_mindist
                minidx = -1_ni
                do m = 1_ni,cycidx
                   dlon = pi/180.0_nr * ( lon(i) - meta(k,m,2_ni) )
                   lat1r = pi/180.0_nr * lat(j)
                   lat2r = pi/180.0_nr * meta(k,m,1_ni)
                   dist = 6370.0_nr * acos(sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r)*cos(dlon))
                   if ( dist < mindist ) then
                      mindist = dist
                      minidx = m
                   end if
                end do

                ! Minimum belongs to existing cyclone
                if ( minidx > 0 ) then
                   ! ... but existing cyclone already complete
                   if ( done(minidx) ) then
                      mask(k,j,i) = -1_ni
                   ! ... but existing cyclone would become too large
                   else if ( meta(k,minidx,5_ni) + cellsize(j,i) > cyc_maxsize ) then
                      mask(k,j,i) = -1_ni
                      done(minidx) = .true.
                  ! Point really belongs to existing cyclone
                   else
                      mask(k,j,i) = minidx
                      meta(k,minidx,5_ni) = meta(k,minidx,5_ni) + cellsize(j,i)
                      meta(k,minidx,4_ni) = val
                   endif

                ! Actually new cyclone
                else
                   cycidx = cycidx + 1_ni
                   mask(k,j,i) = cycidx
                   meta(k,cycidx,1_ni) = lat(j)
                   meta(k,cycidx,2_ni) = lon(i)
                   meta(k,cycidx,3_ni) = val
                   meta(k,cycidx,4_ni) = val
                   meta(k,cycidx,5_ni) = cellsize(j,i)
                end if

             else
                mask(k,j,i) = -1_ni
             end if
          !
          ! One cyclone around -> potential extension of existing cyclone
          else if ( ndiffloc == 1_ni ) then
             ocycidx = nmaxloc
             
             ! Only close to previously discarded grid cell
             if ( ocycidx < 0_ni ) then
                mask(k,j,i) = -1_ni
             ! Other cyclone already complete
             else if ( done(ocycidx) ) then
                mask(k,j,i) = -1_ni
             ! Other cyclone would become too large
             else if ( meta(k,ocycidx,5_ni) + cellsize(j,i) > cyc_maxsize ) then
                mask(k,j,i) = -1_ni
                done(ocycidx) = .true.
             ! Point extends other cyclone
             else
                mask(k,j,i) = ocycidx
                meta(k,ocycidx,5_ni) = meta(k,ocycidx,5_ni) + cellsize(j,i)
                meta(k,ocycidx,4_ni) = val
             end if

          ! More than one cyclone around -> mark all of them as complete
          else
             mask(k,j,i) = -1_ni
             do m = 1_ni,8_ni
                if ( loc_mask(m) > 0_ni ) then
                   done(loc_mask(m)) = .true.
                end if
             end do
          end if
       end do ! Loop over sorted grid cells
       !
       !
       ! Clean up negative masks
       do j = 1_ni,ny
          do i = 1_ni,nx
             if ( mask(k,j,i) < 0_ni ) then
                mask(k,j,i) = 0_ni
             end if
          end do
       end do
       !
       ! Enforce minimum size and minimum prominence 
       do n = 1_ni,cycidx
          prominence = meta(k,n,4_ni) - meta(k,n,3_ni)
          ! Too weak or too small -> not a cyclone after all
          if ( prominence < cyc_minprominence .or. meta(k,n,5_ni) < cyc_minsize ) then
             meta(k,n,:) = 0.0_nr
             !
             do j = 1_ni,ny
                do i = 1_ni,nx
                   if ( mask(k,j,i) == n ) then
                      mask(k,j,i) = 0_ni
                   end if
                end do
             end do
          end if
       end do  
       ! 
    end do ! Loop over time
    !
  end subroutine
  !
  !@ Find jetaxes by zero-shear condition and a mask for well-defined wind maxima
  !@
  !@ The mask for well-defined maxima is defined by d/dn(U * dU/dn) < K. The
  !@ threshold K is to be defined in config module of dynfor.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ no : int
  !@     Maximum number of points along jet axes
  !@ nf : int
  !@     Maximum number of jet axis lines
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     u wind velocity component, typically filtered to T42 resolution.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     v wind velocity component, typically filtered to T42 resolution.
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
  !@ 
  !@ np.ndarray with shape (nz,no,3) and dtype float64
  !@     List of points points belonging to jet axes for each time step. For each point, 
  !@     (0) the j-index, (1) the i-index and (2) the wind speed is saved.
  !@ np.ndarray with shape (nz,nf) and dtype float64
  !@     List of point indexes marking the beginning of jet axes within the point array.
  subroutine jetaxis(ja,jaoff,nx,ny,nz,no,nf,u,v,dx,dy)
    use detect_lines
    use utils, only: smooth_xy
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: ja(nz,no,3_ni), jaoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) ja, jaoff
    !
    !real(kind=nr) :: us(nz,ny,nx), vs(nz,ny,nx)
    real(kind=nr) :: dsheardx(nz,ny,nx), dsheardy(nz,ny,nx), &
                 &   jetint(nz,ny,nx), shear(nz,ny,nx), ff(nz,ny,nx)
    integer(kind=ni) :: i,j,k, ip1,im1
    ! -----------------------------------------------------------------
    !
    write(*,*) 'preparing'
    !
    ff = sqrt( u(:,:,:)**2_ni + v(:,:,:)**2_ni )
    !
    ! Based on shear in natural coordinates, and the absolute wind speed u calculate d/dn(u * du/dn)
    call shear_nat(shear, nx,ny,nz, u,v, dx,dy)
    call ddx(dsheardx, nx,ny,nz, shear*ff, dx,dy)
    call ddy(dsheardy, nx,ny,nz, shear*ff, dx,dy)
    jetint(:,:,:) = (u(:,:,:)*dsheardy(:,:,:) - v(:,:,:)*dsheardx(:,:,:))/ff(:,:,:)
    !
    where(jetint > -jetint_thres)
       shear = nan
    end where
    !
    ! Save the wind speed along the jet axis in the third field along with the
    ! jet axis position. 
    !  (the original jetint is not needed anymore after the masking `where` block)
    jetint(:,:,:) = ff(:,:,:)
    !
    call line_locate(ja,jaoff, nx,ny,nz,no,nf, shear,jetint, dx,dy)
    !
    return
  end subroutine
  !
  !@ Find jetaxes by zero-shear condition and masks selecting wind maxima above a 
  !@ certain wind speed threshold.
  !@
  !@ The mask for wind maxima is defined by d^2U/dn < 0. The wind speed threshold
  !@ is currently fixed at 30 m/s.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ no : int
  !@     Maximum number of points along jet axes
  !@ nf : int
  !@     Maximum number of jet axis lines
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     u wind velocity component, typically filtered to T42 resolution.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     v wind velocity component, typically filtered to T42 resolution.
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
  !@ no : int
  !@     Maximum number of points in all lines.
  !@ nf : int
  !@     Maximum number of lines.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,no,3) and dtype float64
  !@     List of points points belonging to jet axes for each time step. For each point, 
  !@     (0) the j-index, (1) the i-index and (2) the wind speed is saved.
  !@ np.ndarray with shape (nz,nf) and dtype float64
  !@     List of point indexes marking the beginning of jet axes within the ``ja`` array.
  subroutine jetaxis_ff_thres(ja,jaoff,nx,ny,nz,no,nf,u,v,dx,dy)
    use detect_lines
    use utils, only: smooth_xy
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: ja(nz,no,3_ni), jaoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) ja, jaoff
    !
    real(kind=nr) :: us(nz,ny,nx), vs(nz,ny,nx), ffs(nz,ny,nx), & 
                 &   ddangdx(nz,ny,nx), ddangdy(nz,ny,nx), &
                 &   jetint(nz,ny,nx), jaloc(nz,ny,nx)
    integer(kind=ni) :: i,j,k, ip1,im1
    ! -----------------------------------------------------------------
    !
    write(*,*) 'preparing'
    !
    call smooth_xy(us, nx,ny,nz, u, nsmooth)
    call smooth_xy(vs, nx,ny,nz, v, nsmooth)
    !
    ffs(:,:,:) = sqrt(us(:,:,:)**2_ni + vs(:,:,:)**2_ni)
    !
    ! Based on shear in natural  coordinates
    call shear_nat(jaloc, nx,ny,nz, us,vs, dx,dy)
    call ddx(ddangdx, nx,ny,nz, jaloc, dx,dy)
    call ddy(ddangdy, nx,ny,nz, jaloc, dx,dy)
    !
    jetint(:,:,:) = us(:,:,:)*ddangdy(:,:,:) - vs(:,:,:)*ddangdx(:,:,:)
    !
    ! Mask wind minima
    where(jetint >= 0.0_nr)
       jaloc = nan
    end where
    ! Mask low wind zones
    where(ffs < 30.0_nr)
       jaloc = nan
    end where
    !
    ! Save the wind speed along the jet axis in the third field along with the
    ! jet axis position. 
    !  (the original jetint is not needed anymore after the masking `where` block)
    jetint(:,:,:) = sqrt(us(:,:,:)**2.0_nr + vs(:,:,:)**2.0_nr)
    !
    call line_locate(ja,jaoff, nx,ny,nz,no,nf, jaloc,jetint, dx,dy)
    !
    return
  end subroutine
  !
  !@ Front intensity function after G. Berry et al, based on the frontal detection
  !@ algorithm by Hewson (1998) which uses the Laplacian of the equivalent potential
  !@ temperature as dat
  !@
  !@ Used only internally in dynlib to detect front lines at the maximum curvature
  subroutine front_intensity_speed_maxcurv(frint,frspd,frloc,nx,ny,nz,dat,u,v,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: frint(nz,ny,nx), frspd(nz,ny,nx), frloc(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) frint, frspd, u, v
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: datx (nz,ny,nx), daty (nz,ny,nx), absgrad(nz,ny,nx), &
                 &   absx (nz,ny,nx), absy (nz,ny,nx), abslap (nz,ny,nx), &
                 &   absxx(nz,ny,nx), absyy(nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    call grad(datx,daty, nx,ny,nz, dat, dx,dy)
    absgrad(:,:,:) = sqrt(datx**2.0_nr + daty**2.0_nr)
    call grad(absx,absy, nx,ny,nz, absgrad, dx,dy)
    abslap (:,:,:) = sqrt(absx**2.0_nr + absy**2.0_nr)
    !
    frint(:,:,:) = (datx(:,:,:)*absx(:,:,:) + daty(:,:,:)*absy(:,:,:)) / absgrad(:,:,:)
    frspd(:,:,:) = (u(:,:,:)*absx(:,:,:) + v(:,:,:)*absy(:,:,:)) / abslap(:,:,:)
    !
    ! determine front line location type after Hewson 1998, eq. 5
    call ddx(absxx, nx,ny,nz, absx, dx,dy)
    call ddy(absyy, nx,ny,nz, absy, dx,dy)
    frloc(:,:,:) = absxx(:,:,:) + absyy(:,:,:)
    !
    return
  end subroutine
  !
  !@ Front intensity function after Jenkner et al. (2010), loosely based on the frontal 
  !@ detection algorithm by Hewson (1998) which uses the gradient of the equivalent 
  !@ potential temperature as dat
  !@
  !@ Used only internally in dynlib to detect front lines at the maximum gradient
  subroutine front_intensity_speed_maxgrad(frint,frspd,frloc,nx,ny,nz,dat,u,v,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: frint(nz,ny,nx), frspd(nz,ny,nx), frloc(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) frint, frspd, u, v
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: datx (nz,ny,nx), daty (nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    call grad(datx,daty, nx,ny,nz, dat, dx,dy)
    call lap2(frloc, nx,ny,nz, dat, dx,dy)
    !
    frint(:,:,:) = sqrt(datx**2.0_nr + daty**2.0_nr)
    frspd(:,:,:) = (u(:,:,:)*datx(:,:,:) + v(:,:,:)*daty(:,:,:)) / frint(:,:,:)
    !
    return
  end subroutine
  !
  !@ Front line detection after Berry et al. (2007)
  !@ 
  !@ The scheme is based on the frontal detection algorithm by Hewson (1998) 
  !@ which uses the Laplacian of the equivalent potential temperature as dat
  !@ 
  !@ The algorithm was developed with ERA40 data of 2.5x2.5 degree resolution, 
  !@ so some smoothing might be necessary to obtain useful results for other
  !@ data sets. Smoothing can be done within this function if configured by  
  !@ the ``dynfor.config.nsmooth`` parameter.
  !@ 
  !@ The front intensity threshold can be overridden via 
  !@ ``dynfor.config.frint_thres``.
  !@ The front movement threshold to separate between warm, stationary and cold
  !@ fronts can be set via ``dynfor.config.frspd_thres``.
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ no : integer
  !@     Maximum number of points along all front lines combined, 
  !@     determinining the output array size
  !@ nf : integer
  !@     Maximum number of front lines, determining the output array size
  !@ dat : np.ndarray with shape (nz,ny,nx)
  !@     Thermal frontal parameter (TFP), often potential or equivalent 
  !@     potential temperature
  !@ u : np.ndarray with shape (nz,ny,nx)
  !@     Wind velocity component u
  !@ v : np.ndarray with shape (nz,ny,nx)
  !@     Wind velocity component v
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
  !@ np.ndarray with shape (nz,no,3) and dtype float64
  !@     List of points points belonging to front lines for each time step. For each point, 
  !@     (0) the j-index, (1) the i-index and (2) the TFP is saved.
  !@ np.ndarray with shape (nz,nf) and dtype float64
  !@     List of point indexes marking the beginning of front lines within the point array.
  !@
  !@ See also
  !@ --------
  !@ :meth:`frontalzone_smallscale`, :meth:`frontalzone_largescale`, :meth:`frontline_at_maxgrad`
  subroutine frontline_at_maxcurv(fr,froff,nx,ny,nz,no,nf,dat,u,v,dx,dy)
    use config
    use detect_lines
    use utils, only: smooth_xy
    !
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: fr(nz,3_ni,no,3_ni), froff(nz,3_ni,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) u, v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) fr, froff
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:)
    integer(kind=ni), allocatable :: linelen(:)
    !
    real(kind=nr) :: us     (nz,ny,nx), vs   (nz,ny,nx), dats (nz,ny,nx), &
                 &   frint  (nz,ny,nx), frspd(nz,ny,nx), frloc(nz,ny,nx), &
                 &   frloc_cws(nz,ny,nx)
    integer(kind=ni) :: typ
    ! -----------------------------------------------------------------
    !
    write(*,*) 'preparing'
    !
    call smooth_xy(dats, nx,ny,nz, dat, nsmooth)
    call smooth_xy(us, nx,ny,nz, u, nsmooth)
    call smooth_xy(vs, nx,ny,nz, v, nsmooth)
    !
    call front_intensity_speed_maxcurv(frint,frspd,frloc,nx,ny,nz,dats,us,vs,dx,dy)
    !
    ! frint must be negative and below a configurable threshold for front
    where(frint > frint_thres)
       frloc = nan
    end where
    ! 
    do typ = 1_ni,3_ni
       ! cold fronts
       if (typ == 1_ni) then
          where(frspd(:,:,:) < -1_ni*frspd_thres)
             frloc_cws = frloc(:,:,:)
          elsewhere
             frloc_cws = nan
          end where
       ! warm fronts
       elseif(typ == 2_ni) then
          where(frspd(:,:,:) > frspd_thres)
             frloc_cws = frloc(:,:,:)
          elsewhere
             frloc_cws = nan
          end where
       ! stationary fronts
       else
          where(frspd(:,:,:) > -1_ni*frspd_thres .and. frspd(:,:,:) < frspd_thres)
             frloc_cws = frloc(:,:,:)
          elsewhere
             frloc_cws = nan
          end where
       end if
       !
       call line_locate(fr(:,typ,:,:),froff(:,typ,:), nx,ny,nz,no,nf, frloc_cws,frint, dx,dy)
       !
    end do ! loop over front type
    !
    return
  end subroutine
  !
  !@ Front line detection after Jenkner et al. (2010)
  !@ 
  !@ The scheme is loosely based on the frontal detection algorithm by Hewson (1998),
  !@ but uses the maximum in gradient rather than the maximum in curvature to locate
  !@ the front. Hence fronts will not be placed at the (leading) edge of an area of
  !@ strong gradients, but rather in it's centre. 
  !@
  !@ This algorithm is more robust for higher-resolution datasets that the one by 
  !@ Berry et al. (2007). Nevertheless, some smoothing might be necessary to obtain 
  !@ useful results. Smoothing can be done within this function if configured by  
  !@ the ``dynfor.config.nsmooth`` parameter.
  !@ 
  !@ The front intensity threshold can be overridden via 
  !@ ``dynfor.config.frint_thres``.
  !@ The front movement threshold to separate between warm, stationary and cold
  !@ fronts can be set via ``dynfor.config.frspd_thres``.
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ no : integer
  !@     Maximum number of points along all front lines combined, 
  !@     determinining the output array size
  !@ nf : integer
  !@     Maximum number of front lines, determining the output array size
  !@ dat : np.ndarray with shape (nz,ny,nx)
  !@     Thermal frontal parameter (TFP), often potential or equivalent 
  !@     potential temperature
  !@ u : np.ndarray with shape (nz,ny,nx)
  !@     Wind velocity component u
  !@ v : np.ndarray with shape (nz,ny,nx)
  !@     Wind velocity component v
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
  !@ np.ndarray with shape (nz,no,3) and dtype float64
  !@     List of points points belonging to front lines for each time step. For each point, 
  !@     (0) the j-index, (1) the i-index and (2) the TFP is saved.
  !@ np.ndarray with shape (nz,nf) and dtype float64
  !@     List of point indexes marking the beginning of front lines within the point array.
  !@
  !@ See also
  !@ --------
  !@ :meth:`frontalzone_smallscale`, :meth:`frontalzone_largescale`, :meth:`frontline_at_maxcurv`
  subroutine frontline_at_maxgrad(fr,froff,nx,ny,nz,no,nf,dat,u,v,dx,dy)
    use config
    use detect_lines
    use utils, only: smooth_xy
    !
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: fr(nz,3_ni,no,3_ni), froff(nz,3_ni,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) u, v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) fr, froff
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:)
    integer(kind=ni), allocatable :: linelen(:)
    !
    real(kind=nr) :: us     (nz,ny,nx), vs   (nz,ny,nx), dats (nz,ny,nx), &
                 &   frint  (nz,ny,nx), frspd(nz,ny,nx), frloc(nz,ny,nx), &
                 &   frloc_cws(nz,ny,nx)
    integer(kind=ni) :: typ
    ! -----------------------------------------------------------------
    !
    write(*,*) 'preparing'
    !
    call smooth_xy(dats, nx,ny,nz, dat, nsmooth)
    call smooth_xy(us, nx,ny,nz, u, nsmooth)
    call smooth_xy(vs, nx,ny,nz, v, nsmooth)
    !
    call front_intensity_speed_maxgrad(frint,frspd,frloc,nx,ny,nz,dats,us,vs,dx,dy)
    !
    ! frint must be positive and above a configurable threshold for front
    where(frint < frint_thres)
       frloc = nan
    end where
    ! 
    do typ = 1_ni,3_ni
       ! cold fronts
       if (typ == 1_ni) then
          where(frspd(:,:,:) > frspd_thres)
             frloc_cws = frloc(:,:,:)
          elsewhere
             frloc_cws = nan
          end where
       ! warm fronts
       elseif(typ == 2_ni) then
          where(frspd(:,:,:) < -1_ni*frspd_thres)
             frloc_cws = frloc(:,:,:)
          elsewhere
             frloc_cws = nan
          end where
       ! stationary fronts
       else
          where(frspd(:,:,:) > -1_ni*frspd_thres .and. frspd(:,:,:) < frspd_thres)
             frloc_cws = frloc(:,:,:)
          elsewhere
             frloc_cws = nan
          end where
       end if
       !
       call line_locate(fr(:,typ,:,:),froff(:,typ,:), nx,ny,nz,no,nf, frloc_cws,frint, dx,dy)
       !
    end do ! loop over front type
    !
    return
  end subroutine
  !
  !@ Convergence line detection based on the front detection algorithm, using
  !@ convergence instead of grad(dat)
  !@ 
  !@ *Note*: Routine not sufficiently tested for general applicability! More thorough 
  !@ documentation will be added once properly tested.
  subroutine convline(fr,froff,nx,ny,nz,no,nf,u,v,dx,dy)
    use config
    use detect_lines
    use utils, only: smooth_xy
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: fr(nz,no,3_ni), froff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) fr, froff
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:)
    integer(kind=ni), allocatable :: linelen(:)
    !
    real(kind=nr) :: us(nz,ny,nx), vs(ny,ny,nx), &
                 &   absgrad(nz,ny,nx), absx (nz,ny,nx), absy (nz,ny,nx), &
                 &   absxx(nz,ny,nx), absyy(nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    write(*,*) 'preparing'
    !
    ! In case of segfaults "double free" on return
    ! comment out the calls to smooth_xy below and use u and v instead of us and vs 
    !call smooth_xy(us, nx,ny,nz, u, nsmooth)
    !call smooth_xy(vs, nx,ny,nz, v, nsmooth)
    !
    call div(absgrad, nx,ny,nz, u,v, dx,dy)
    !call grad(absx,absy, nx,ny,nz, absgrad, dx,dy)
    !call ddx2(absxx, nx,ny,nz, absgrad, dx,dy)
    !call ddy2(absyy, nx,ny,nz, absgrad, dx,dy)
    !
    ! frint must be negative and below a configurable threshold for front
    !where(absgrad > div_thres)
    !   absx = nan
    !   absy = nan
    !end where
    ! Masking which criterion to use
    !where(absxx > absyy)
    !   absx = nan
    !else where
    !   absy = nan
    !end where
    !
    absgrad = -absgrad
    !
    call maxline_locate(fr,froff, nx,ny,nz,no,nf, absgrad,-div_thres, dx,dy)
    !
    return
  end subroutine
  !
  !@ Vorticity line detection based on the front detection algorithm, using
  !@ vorticity instead of grad(dat)
  !@ 
  !@ *Note*: Routine not sufficiently tested for general applicability! More thorough 
  !@ documentation will be added once properly tested.
  subroutine vorline(fr,froff,nx,ny,nz,no,nf,u,v,dx,dy)
    use config
    use detect_lines
    use utils, only: smooth_xy
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: fr(nz,no,3_ni), froff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) fr, froff
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:)
    integer(kind=ni), allocatable :: linelen(:)
    !
    real(kind=nr) :: us(nz,ny,nx), vs(ny,ny,nx), &
                 &   absgrad(nz,ny,nx), absx(nz,ny,nx), absy(nz,ny,nx), &
                 &   absxx(nz,ny,nx), absyy(nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    write(*,*) 'preparing'
    !
    ! In case of segfaults "double free" on return
    ! comment out the calls to smooth_xy below and use u and v instead of us and vs 
    !call smooth_xy(us, nx,ny,nz, u, nsmooth)
    !call smooth_xy(vs, nx,ny,nz, v, nsmooth)
    !
    call vor(absgrad, nx,ny,nz, u,v, dx,dy)
    !call grad(absx,absy, nx,ny,nz, absgrad, dx,dy)
    !call ddx2(absxx, nx,ny,nz, absgrad, dx,dy)
    !call ddy2(absyy, nx,ny,nz, absgrad, dx,dy)
    !
    ! frint must be negative and below a configurable threshold for front
    !where(absgrad < vor_thres)
    !   absx = nan
    !   absy = nan
    !end where
    ! Masking which criterion to use
    !where(absxx > absyy)
    !   absx = nan
    !else where
    !   absy = nan
    !end where
    ! 
    call maxline_locate(fr,froff, nx,ny,nz,no,nf, absgrad,vor_thres, dx,dy)
    !
    return
  end subroutine
  !
  !@ Deformation line detection based on the front detection algorithm, using
  !@ deformation instead of grad(dat)
  !@ 
  !@ *Note*: Routine not sufficiently tested for general applicability! More thorough 
  !@ documentation will be added once properly tested.
  subroutine defline(fr,froff,nx,ny,nz,no,nf,u,v,dx,dy)
    use config
    use detect_lines
    use utils, only: smooth_xy
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: fr(nz,no,3_ni), froff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) fr, froff
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:)
    integer(kind=ni), allocatable :: linelen(:)
    !
    real(kind=nr) :: us(nz,ny,nx), vs(ny,ny,nx),  &
                 &   absgrad(nz,ny,nx), absx(nz,ny,nx), absy(nz,ny,nx), &
                 &   absxx(nz,ny,nx), absyy(nz,ny,nx)
    ! -----------------------------------------------------------------
    !
    ! In case of segfaults "double free" on return
    ! comment out the calls to smooth_xy below and use u and v instead of us and vs 
    !call smooth_xy(us, nx,ny,nz, u, nsmooth)
    !call smooth_xy(vs, nx,ny,nz, v, nsmooth)
    !
    call def_total(absgrad, nx,ny,nz, u,v, dx,dy)
    !call grad(absx,absy, nx,ny,nz, absgrad, dx,dy)
    !call ddx2(absxx, nx,ny,nz, absgrad, dx,dy)
    !call ddy2(absyy, nx,ny,nz, absgrad, dx,dy)
    !
    ! frint must be negative and below a configurable threshold for front
    !where(absgrad < def_thres)
    !   absx = nan
    !   absy = nan
    !end where
    ! Masking which criterion to use
    !where(absxx > absyy)
    !   absx = nan
    !else where
    !   absy = nan
    !end where
    ! 
    call maxline_locate(fr,froff, nx,ny,nz,no,nf, absgrad,def_thres, dx,dy)
    !
    return
  end subroutine
  !
  ! TODO: Reduce duplication with line_locate
  subroutine maxline_locate(lines,lnoff,nx,ny,nz,no,nf,dat,thres,dx,dy)
    use consts
    use detect_lines, only: find_maxloc, linejoin
    !
    real(kind=nr), intent(in) :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx), thres
    real(kind=nr), intent(out) :: lines(nz,no,3_ni), lnoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) lines, lnoff
    !
    integer(kind=ni), parameter :: nn = 30000_ni
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:), linelen(:)
    integer(kind=ni), allocatable :: lineptcnt(:)
    !
    real(kind=nr) :: maxloc_(nz,nn,2_ni)
    integer(kind=ni) :: i,j,k, di, m, n, maxcnt(nz), ptcnt, linecnt, off
    ! -----------------------------------------------------------------
    !
    ! find lines by zero-criterion
    call find_maxloc(dat(:,:,:),thres, nx,ny,nz, nn,  maxloc_(:,:,:),maxcnt(:), dx,dy)
    !
    do k = 1_ni,nz
       write(*,'(I5,A4,I5,A)', advance='no') k, 'of', nz, cr
       !
       ! searchrad is in grid point indexes, as zero locations are found at grid
       ! resolution, hence the number of neighbours for a given searchrad does not
       ! depend on location within the grid
       !
       allocate(recj(maxcnt(k),maxcnt(k)), reci(maxcnt(k),maxcnt(k)), lineptcnt(maxcnt(k)), linelen(maxcnt(k)) )
       reci(:,:) = nan
       recj(:,:) = nan
       linecnt    = 0_ni ! number of lines
       ptcnt      = 0_ni ! total numer of points
       lineptcnt(:) = 0_ni ! number of points per line
       call linejoin(maxcnt(k), maxcnt(k), nx,ny, maxloc_(k,1_ni:maxcnt(k),2_ni), maxloc_(k,1_ni:maxcnt(k),1_ni), &
          &          searchrad, recj, reci, linelen, lineptcnt, dx,dy) 
       !
       off = 0_ni
       do n = 1_ni,maxcnt(k)
          if ( isnan(recj(n,1_ni)) ) then
             exit
          end if
          !
          ! filter fronts by length
          if (linelen(n) >= minlen) then
          !if (.true.) then
             linecnt = linecnt + 1_ni
             ptcnt = ptcnt + lineptcnt(n)
             !
             ! check if results larger than output array
             if (ptcnt > no) then
                write(*,*) 'Found more points than output array allows: ', no
                stop 1
             end if
             if (linecnt > nf) then
                write(*,*) 'Found more fronts than output array allows: ', nf
                stop 1
             end if
             !
             ! write into output arrays fr and froff
             do m = 1_ni,lineptcnt(n)
                lines(k,off+m,1_ni) = reci(n,m)
                lines(k,off+m,2_ni) = recj(n,m)
                lines(k,off+m,3_ni) = dat(k,int(recj(n,m),ni),int(reci(n,m),ni))
             end do
             lnoff(k,linecnt) = off
             off = off + lineptcnt(n)
          end if
       end do
       ! Save the ending of the last front by saving the beginning of the
       ! first non-existant
       lnoff(k,linecnt+1_ni) = off
       !
       deallocate(reci, recj, lineptcnt, linelen)
       !
    end do ! loop over k
    !
    return
  end subroutine
   !@ ---------------------------------------------------------------------------
   !@ Subroutine: find_min_around
   !@ Description: This subroutine determines the minimum value in
   !@              the vicinity of a grid cell and updates the corresponding
   !@              mask and local mask arrays.
   !@
   !@ Arguments:
   !@   val       - Input, the value of the current grid cell.
   !@   k         - Input, the current time index.
   !@   j         - Input, the row index of the current grid cell.
   !@   i         - Input, the column index of the current grid cell.
   !@   nx        - Input, the total number of columns in the grid.
   !@   ny        - Input, the total number of rows in the grid.
   !@   msl       - Input, the precipitation or pressure field array.
   !@   mask      - Input/output, the mask array indicating connected components.
   !@   loc_mask  - Output, an array indicating the presence of masks in the vicinity.
   !@   extremum  - Output, the determined maximum or minimum value in the vicinity.
   !@---------------------------------------------------------------------------
  subroutine find_min_around(val, k, j, i, nx, ny, msl, mask, extremum, loc_mask)
   real(kind=nr), intent(in) :: val, msl(:,:,:)  ! Assuming msl is a 3D array
   integer(kind=ni), intent(in) :: k, j, i, nx, ny
   integer(kind=ni), intent(inout) :: mask(:,:,:), loc_mask(8_ni)
   real(kind=nr), intent(out) :: extremum
   integer(kind=ni) :: ip1, im1, jp1, jm1

   ip1 = modulo(i, nx) + 1_ni
   im1 = modulo(i - 2_ni, nx) + 1_ni
   jp1 = min(j + 1_ni, ny)
   jm1 = max(j - 1_ni, 1_ni)

   extremum = min(msl(k,j,ip1), msl(k,jp1,ip1), msl(k,jp1,i), msl(k,jp1,im1), &
                  msl(k,j,im1), msl(k,jm1,im1), msl(k,jm1,i), msl(k,jm1,ip1) )

   loc_mask(:) = (/ mask(k,j,ip1), mask(k,jp1,ip1), mask(k,jp1,i), mask(k,jp1,im1), &
                    mask(k,j,im1), mask(k,jm1,im1), mask(k,jm1,i), mask(k,jm1,ip1) /)
  end subroutine
  !@ Precipitation blob masks by finding the outermost contour
  !@
  !@ ADAPTED FROM THE CYCLONE CONTOUR DETECTION
  !@ 
  !@ 
  !@ To avoid the technically difficult contour tracing, we base the detection on 
  !@ sorted precipitation values. 
  !@
  !@ For more info refer to the python function. This routine is not intended to
  !@ be called directly.
  !@
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ nn : int
  !@     Maximum number of precipitation blobs (objects) to be detected
  !@ precip : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Precipitation field.
  !@ precips : np.ndarray with shape (nz,ny*nx) and dtype float64
  !@     Sorted precipitation values for each time step.
  !@ iis : np.ndarray with shape (nz,ny*nx) and dtype float64
  !@     Longitudinal grid index of the respective precipitation  value.
  !@ jjs : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Latitudinal grid index of the respective precipitation value.
  !@ lon : np.ndarray with shape(nx)
  !@     Longitudes in degrees
  !@ lat : np.ndarray with shape(ny)
  !@     Latitudes in degrees
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ blob_mindist : float 
  !@     Minimum number of the distance between two maxima. 
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
  !@ np.ndarray with shape (nz,ny,nx) and dtype int32
  !@     Precipitation object mask with different integers designating different blobs (objects)
  !@ np.ndarray with shape (nz,nn,5)
  !@     Meta data about each blob: (1) latitutde and (2) longitude of its centre,
  !@     (3) maximum P, (4) P at the outermost contour, and (5) blob size.
  subroutine blobs_fortran(mask,meta, nx,ny,nz, nn, precip, precips,iis,jjs,  lon,lat, dx,dy, blob_mindist)
    use consts
    !
    real(kind=8), intent(in)  :: precip(nz,ny,nx), precips(nz,ny*nx),  &
       &                          lon(nx), lat(ny), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: meta(nz,nn,5_ni)
    integer(kind=ni), intent(in) :: iis(nz,ny*nx), jjs(nz,ny*nx), nn,nx,ny,nz
    integer(kind=ni), intent(out) :: mask(nz,ny,nx)
    !f2py depend(nx,ny,nz) mask, msls, iis, jjs
    !f2py depend(nn,nz) meta
    !f2py depend(nx,ny) oro, dx, dy
    !f2py depend(nx) lon
    !f2py depend(ny) lat
    !
    real(kind=nr) :: cellsize(ny,nx), max_around, dist, mindist, dlon, lat1r, lat2r, val, blob_mindist
    integer(kind=ni) :: loc_mask(8_ni), blobidx, oblobidx, maxidx, ndiffloc, nmaxloc
    integer(kind=ni) :: i,j,k,n,m, ip1,im1,jp1,jm1
    ! -----------------------------------------------------------------
    !
    ! Initialize stuff
    meta(:,:,:) = 0_ni
    !
    do i = 1_ni,nx
       do j = 1_ni,ny
          cellsize(j,i) = abs(dx(j,i) * dy(j,i)) / 4.0e6 ! double grid spacing in dx/dy and conversion from m^2 to km^2
          if ( isnan(cellsize(j,i)) ) then
             cellsize(j,i) = 0.0_nr
          end if
       end do
    end do
    !
    ! Loops over time and sorted arrays
    do k = 1,nz
       blobidx = 0_ni
       !
       do n = 1,ny*nx
          i = iis(k,n) + 1_ni ! conversion to Fortran indexes
          j = jjs(k,n) + 1_ni
          val = precips(k,n)
          !
          ! We only want to sort through precip with a little bit of magnitude 
          if ( val < 0) then  ! But set to negative because we pass a negative field. 
              call find_min_around(val, k, j, i, nx, ny, precip, mask, max_around, loc_mask)
              
              ndiffloc = 0_ni
              nmaxloc = maxval(loc_mask)
              do m = 0_ni,8_ni
                 if ( loc_mask(m) == 0_ni ) then
                    cycle
                 end if

                 if ( loc_mask(m) == nmaxloc ) then
                    ndiffloc = max(ndiffloc,1_ni)
                 else
                    ndiffloc = max(ndiffloc,2_ni)
                 end if
              end do

              ! Nothing around -> potential new blob center
              if ( minval(loc_mask) == 0_ni .and. maxval(loc_mask) == 0_ni ) then
               if ( val < max_around ) then
                  ! Check proximity to previously detected blob centres
                  mindist = blob_mindist
                  maxidx = -1_ni
                  do m = 1_ni,blobidx
                     dlon = pi/180.0_nr * ( lon(i) - meta(k,m,2_ni) )
                     lat1r = pi/180.0_nr * lat(j)
                     lat2r = pi/180.0_nr * meta(k,m,1_ni)
                     dist = 6370.0_nr * acos(sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r)*cos(dlon))
                     if ( dist < mindist ) then
                        mindist = dist
                        maxidx = m
                     end if
                  end do

                  ! Minimum belongs to existing blob
                  if ( maxidx > 0 ) then
                    ! Point really belongs to existing blob
                        mask(k,j,i) = maxidx
                        meta(k,maxidx,5_ni) = meta(k,maxidx,5_ni) + cellsize(j,i)
                        meta(k,maxidx,4_ni) = val

                  ! Actually new blob
                  else
                     blobidx = blobidx + 1_ni
                     mask(k,j,i) = blobidx
                     meta(k,blobidx,1_ni) = lat(j)
                     meta(k,blobidx,2_ni) = lon(i)
                     meta(k,blobidx,3_ni) = val
                     meta(k,blobidx,4_ni) = val
                     meta(k,blobidx,5_ni) = cellsize(j,i)
                  end if

               else
                  mask(k,j,i) = -1_ni
               end if
            !
            ! One blob around -> extension of existing blob
            else 
               oblobidx = nmaxloc
               ! Point extends to other blob
               mask(k,j,i) = oblobidx
               meta(k,oblobidx,5_ni) = meta(k,oblobidx,5_ni) + cellsize(j,i)
               meta(k,oblobidx,4_ni) = val
            end if
        end if 
     end do ! Loop over sorted grid cells
  end do ! Loop over time

  end subroutine

end module

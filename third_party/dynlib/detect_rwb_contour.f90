! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- contour tracking and RWB detection
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module adopted from code by Gwendal Riviere
! Module maintained by <?>
module detect_rwb_contour
  use kind
  use config
  !
  implicit none
  ! The module is used (only) in the contour_rwb subroutine in diag.f95
contains
!======================================================================
! This section contains all subroutines needed to detect RWBs
!======================================================================
! This subroutine is that called in contour_rwb subroutine to detect RWBs
        SUBROUTINE WB_DETECTION(grd,PV2d,NX,NY,ncon,INCRX,INCRY,&
        lon,lat,XVALF,YVALF,beta_a2d,beta_c2d)
! The input variables of the subroutine (given at the subroutine)
! grd= PV (for potential vorticity) or VO (for absolute vorticity)
        CHARACTER(LEN=2),INTENT(IN) :: grd
! NX and NY are the numbers of longitude and latitude
! ncon is the number of contours we want to detect (generally =21)
        integer(kind=ni),INTENT(IN) :: NX,NY,ncon
! ndim is the highest number of points that can define a contour
        integer(kind=ni),PARAMETER :: ndim=8000
! incrx,incry is the grid step in degrees in x and y directions
        real(kind=nr),INTENT(IN) :: INCRX, INCRY
! lon is the table containing all longitudes (from 180 degW to 180 degE)
        real(kind=nr),DIMENSION(NX),INTENT(IN) :: lon
! lat is the table containing all latitudes
        real(kind=nr),DIMENSION(NY),INTENT(IN) :: lat
! PV2d is the initial field on which we want to detect RWBs (the first
! and last points must be equal because it is the same longitude 180W
! and 180E)
        real(kind=nr),DIMENSION(NX,NY),INTENT(IN) :: PV2d
! XVALF and YVALF are the longitudes and latitudes defining the contours
! It depends on the latitude NY and on the number of the contour
! detected ncon
        real(kind=nr),DIMENSION(ndim,NY,ncon),INTENT(INOUT) :: XVALF,YVALF
! beta_a2d is a table of integers which are equal to 1 when there is an
! AWB and 0 when there is not
! beta_c2d is the same table but for CWB
        integer(kind=ni),DIMENSION(NX,NY),INTENT(INOUT) :: beta_a2d,beta_c2d

        integer(kind=ni) :: i,j,pp,p,dd,iii,jjj
        integer(kind=ni) :: numc,icompt,ee,nbptx,nbc,icomptprecedent
        integer(kind=ni) :: il,jl,ilongl,ilongu,jlatl,jlatu
        integer(kind=ni) :: nbpta,nbptc,iv_err
        integer(kind=ni) :: nbpta_temp,nbptc_temp
        integer(kind=ni),PARAMETER :: npt=8
        real(kind=nr),DIMENSION(ncon) :: zz
        real(kind=nr),DIMENSION(ndim) :: XVAL,YVAL
        real(kind=nr),DIMENSION(npt) :: xseg,yseg
        integer(kind=ni),DIMENSION(NY) :: flag
        integer(kind=ni),DIMENSION(NX,NY) :: beta_a_temp,beta_c_temp
        LOGICAL :: LRead
        integer(kind=ni) :: sd1
        integer(kind=ni) :: nc1,sd

! Initialisations
        zz=0.0
        xseg=0.0
        yseg=0.0
        XVAL=0.0
        YVAL=0.0
        flag=0
        beta_a_temp=0.0
        beta_c_temp=0.0
        XVALF=0.0
        YVALF=0.0
        sd1=0
        sd=0
        
! zz is the table containing contours values (with ncon=21 or 41)
! For PV: from 0 pvu to 10 pvu with an interval of 0.5 pvu
! For VO (absolute vorticity): from -20e-5 to 20e-5 with an interval of 2e-5 s-1
! For Z2 (250 hPa geopotential): from 90e3 to 105e3 in steps of 750 m^2 s^-2 for ncon=21
! For Z5 (500 hPa geopotential): from 49e3 to 56e3 in steps of 350 m^2 s^-2 for ncon=21
! For PT (potential temperature): from 280K to 380K in steps of 5K with ncon=21
        nc1=int(ncon/2)+1
!        print '(A4,I2)','nc1=',nc1
!        print '(A3,I2)','ncon=',ncon
        DO pp=1,ncon
! if the initial field studied is PV on the northern hemisphere
          if (grd.eq.'PV') then
            zz(pp)=0.5e-6*(real(pp-nc1))
          endif
! if the initial field studied is 250 hPa geopotential on the sphere
          if (grd.eq.'Z2') then
            zz(pp)=97500 + 750*(real(pp-nc1))
          endif
! if the initial field studied is 500 hPa geopotential on the sphere
          if (grd.eq.'Z5') then
            zz(pp)=52500 + 350*(real(pp-nc1))
          endif
! if the initial field studied is the absolute vorticity on the sphere
          if (grd.eq.'VO') then
            zz(pp)=2.0e-5*(real(pp-nc1))
          endif
! if the initial field studied is potential temperature on the sphere
          if (grd.eq.'PT') then
            zz(pp)=5.0*(real(pp-nc1)) + 330.0
          endif
        ENDDO
!        print *,'MAXVAL(PV2d)=',MAXVAL(PV2d)

! Scheme of a box with the coordinates of the summits
!     j+1      j+1
!     i -------- i+1
!       |      |
!       |      |
!     i -------- i+1
!       j      j
!
! Initialisations
        nbc=0
        nbpta=0
        nbptc=0

!        print *,'NX=',NX
!        print *,'NY=',NY
!        print *,'lat(NY)=',lat(NY)
!        print *,'lat(1)=',lat(1)
        
! Loop on the contours' values
        DO p=1,ncon
!          print *,'zz(p)=',zz(p)
          i=1
          flag=1
          numc=1
! Loop on the latitudes
          DO j=1,NY-1
            icompt=0
            XVAL=0.0
            YVAL=0.0
            xseg=0.0
            yseg=0.0
            ee=0
            nbptx=0
            if (flag(j).ne.2) then
! conrec is the subroutine that detects segments of the contours in the
! box studied
              call conrec(PV2d(i:i+1,j:j+1),lon(i:i+1),&
              lat(j:j+1),zz(p),xseg,yseg,npt,nbc)
              if (nbc.ne.0) then
                if ((xseg(nbc).ne.(-180.0)).and.(xseg(nbc-1).ne.(-180.0))) CYCLE 
                flag(j)=2
! orgseg is the subroutine that orders the points found by conrec in the
! west-east direction
                call orgseg(xseg(1:nbc),yseg(1:nbc),XVAL,YVAL,nbc,ndim,icompt)
                if (YVAL(icompt).eq.lat(1)) then
                  ee=1
                  EXIT
                endif
                XVALF(1:icompt,numc,p)=XVAL(1:icompt)
                YVALF(1:icompt,numc,p)=YVAL(1:icompt)
                il=i
                jl=j
                LRead=.TRUE.
                if (sd1.eq.0) then
                DO WHILE (LRead)
                  icomptprecedent=icompt
! boitesuiv is the subroutine that defines the coordinate of the next
! box
                  call boitesuiv(il,jl,lon(il:il+1),lat(jl:jl+1),NX-1,NY,&
                  XVAL(icompt),YVAL(icompt),ilongl,ilongu,jlatl,jlatu)
                  xseg=0.0
                  yseg=0.0
! Again we use conrec to find segments in this new box
                  call conrec(PV2d(ilongl:ilongu,jlatl:jlatu),&
                  lon(ilongl:ilongl+1),lat(jlatl:jlatl+1),zz(p),xseg,yseg,npt,nbc)
                  if (nbc.eq.0) then
                    LRead=.FALSE.
                  endif
                  if (ilongl.eq.(il+1)) then
                    nbptx=nbptx+1
                  else
                    if (ilongl.eq.(il-1)) then
                      nbptx=nbptx-1
                    else
                      nbptx=nbptx
                    endif
                  endif
                  il=ilongl
                  jl=jlatl
                  if (nbc.ne.0) then
                    DO dd=1,nbc
                      if ((xseg(dd).eq.(-180.0)).or.(ilongl.eq.1)) then
                        flag(jl)=2
                      endif
                    ENDDO
                  endif
                  if (nbc.ne.0) then
! Again the orgseg subroutine order the pints detected by conrec
                    call orgseg(xseg(1:nbc),yseg(1:nbc),XVAL,YVAL,nbc,ndim,icompt)
                    if (icompt.eq.icomptprecedent) then
                      XVAL=0.0
                      LRead=.FALSE.
                      CYCLE
                    endif
                    XVALF(icomptprecedent+1:icompt,numc,p)=XVAL(icomptprecedent+1:icompt)
                    YVALF(icomptprecedent+1:icompt,numc,p)=YVAL(icomptprecedent+1:icompt)
                  else
                    LRead=.FALSE.
                  endif
                  if ((XVALF(1,numc,p).eq.(-XVALF(icompt,numc,p))).and.&
                  (YVALF(1,numc,p).eq.YVALF(icompt,numc,p))) then
                    LRead=.FALSE.
                  endif
                  DO dd=icomptprecedent+1,icompt
                    if (YVALF(dd,numc,p).ge.lat(NY)) then
                      LRead=.FALSE.
                      ee=1
                    endif
                    if (YVALF(dd,numc,p).le.lat(1)) then
                      LRead=.FALSE.
                      ee=1
                    endif
                  ENDDO
                  if (icompt.ge.6990) then
                    LRead=.FALSE.
                    ee=1
                  endif
                ENDDO
! End of the DO WHILE loop
                if (ee.eq.1) then
                  XVAL=0.0
                  YVAL=0.0
                  XVALF(:,numc,p)=0.0
                  YVALF(:,numc,p)=0.0
                  numc=numc-1
                else
                  if (nbptx.eq.0) then
                    XVAL=0.0
                    YVAL=0.0
                    XVALF(:,numc,p)=0.0
                    YVALF(:,numc,p)=0.0
                    numc=numc-1
                  endif
                endif
                numc=numc+1
                endif
              endif
            endif
            if (sd.eq.0) then
            if (XVAL(1).ne.0.0) then
! Once the whole contour is detected, the eventfunctions subroutine
! finds the WB areas and fills in the two beta tables
              call eventfunctions(icompt,NX,NY,lon(1),lat(1),INCRX,INCRY,&
              XVAL(1:icompt),YVAL(1:icompt),beta_a_temp(1:NX,1:NY),&
              beta_c_temp(1:NX,1:NY),nbpta_temp,nbptc_temp)
              nbpta=nbpta_temp
              nbptc=nbptc_temp
              DO iii=1,NX
                DO jjj=1,NY
                  if (beta_a_temp(iii,jjj).ne.0) then
                    beta_a2d(iii,jjj)=beta_a_temp(iii,jjj)
                  endif
                  if (beta_c_temp(iii,jjj).ne.0) then
                    beta_c2d(iii,jjj)=beta_c_temp(iii,jjj)
                  endif
                ENDDO
              ENDDO
            endif
            endif
! End of the latitudes loop
          ENDDO
          if (p.eq.ncon) then
!          print *,'fin contours'
          endif
! End of the loop on the contours we want to detect
        ENDDO
!        print *,'Fin WB_DETECTION'
!        If (MAXVAL(beta_a2d).ne.0) THEN
!        print *,'Max de beta_a=',MAXVAL(beta_a2d)
!        print *,'Max de beta_c=',MAXVAL(beta_c2d)
!        print *,'Max de XVALF=',MAXVAL(XVALF)
!        print *,'Max de YVALF=',MAXVAL(YVALF)
!        Endif

        END SUBROUTINE

!======================================================================
! This subroutine finds segments of the contour in a rectangle box 
! divided in four triangles
        SUBROUTINE conrec(d,x,y,z,XVAL,YVAL,ndim,nbc)

        real(kind=nr),DIMENSION(2,2),INTENT(IN) :: d
        real(kind=nr),DIMENSION(2),INTENT(IN) :: x,y
        integer(kind=ni),INTENT(IN) :: ndim
        integer(kind=ni) :: i,j,k,m,m1,m2,m3
        real(kind=nr),INTENT(IN) :: z
        real(kind=nr),DIMENSION(ndim),INTENT(OUT) :: XVAL,YVAL
        integer(kind=ni),INTENT(INOUT) :: nbc
        real(kind=nr) x1,y1,x2,y2
!
!     Local declarations
!
        real(kind=nr) h(0:4)
        integer(kind=ni) sh(0:4)
        real(kind=nr) xh(0:4),yh(0:4)
        integer(kind=ni) im(1:4),jm(1:4)
        integer(kind=ni) caseid
        integer(kind=ni) castab(-1:1,-1:1,-1:1)
        integer(kind=ni) p1,p2,ii
        real(kind=nr) l1,l3,pi
        real(kind=nr) l11,l31,l41,l12,l32,l42
        real(kind=nr) dmin,dmax
        data im/0,1,1,0/
        data jm/0,0,1,1/
        data castab/0,0,9,  0,1,5,  7,4,8, &
       &            0,3,6,  2,3,2,  6,3,0, &
       &            8,4,7,  5,1,0,  9,0,0/


        pi=atan(1.)*4
        ii=0
        i=1
        j=1

        dmin = min(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))
        dmax = max(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))

        if ((dmax.ge.z).and.(dmin.le.z)) then
          do 22 m=4,0,-1
            if (m.gt.0) then
              h(m)=d(i+im(m),j+jm(m))-z
              xh(m)=x(i+im(m))
              yh(m)=y(j+jm(m))
            else
              l1=sqrt(((cos(y(j)*pi/180.)*x(i+1)-cos(y(j)*pi/180.)*x(i))/2)**2+((y(j+1)-y(j))/2)**2)
              l3=sqrt(((cos(y(j+1)*pi/180.)*x(i+1)-cos(y(j+1)*pi/180.)*x(i))/2)**2+((y(j+1)-y(j))/2)**2)
              h(0)=(h(1)/l1+h(2)/l1+h(3)/l3+h(4)/l3)/(2/l1+2/l3)
              xh(0)=0.5*(x(i)+x(i+1))
              yh(0)=0.5*(y(j)+y(j+1))
            endif
            if (h(m).gt.0.0) then
              sh(m)=+1
            else if (h(m).lt.0.0) then
              sh(m)=-1
            else
              sh(m)=0
            endif
22        continue
          do 60 m=1,4
            m1=m
            m2=0
            if (m.ne.4) then
              m3=m+1
            else
              m3=1
            endif
            caseid = castab(sh(m1),sh(m2),sh(m3))
            if (caseid.ne.0) then
              goto (31,32,33,34,35,36,37,38,39),caseid
!
!     Case 1 - Line between vertices 1 and 2
!
31            x1=xh(m1)
              y1=yh(m1)
              x2=xh(m2)
              y2=yh(m2)
              goto 40
!
!     Case 2 - Line between vertices 2 and 3
!
32            x1=xh(m2)
              y1=yh(m2)
              x2=xh(m3)
              y2=yh(m3)
              goto 40
!
!     Case 3 - Line between vertices 3 and 1
!
33            x1=xh(m3)
              y1=yh(m3)
              x2=xh(m1)
              y2=yh(m1)
              goto 40
!
!     Case 4 - Line between vertex 1 and side 2-3
!
34            x1=xh(m1)
              y1=yh(m1)
              call xsectsub(x2,h(m2),h(m3),xh(m2),xh(m3))
              call ysectsub(y2,h(m2),h(m3),yh(m2),yh(m3))
              goto 40
!
!     Case 5 - Line between vertex 2 and side 3-1
!
35            x1=xh(m2)
              y1=yh(m2)
              call xsectsub(x2,h(m3),h(m1),xh(m3),xh(m1))
              call ysectsub(y2,h(m3),h(m1),yh(m3),yh(m1))
              goto 40
!
!     Case 6 - Line between vertex 3 and side 1-2
!
36            x1=xh(m3)
              y1=yh(m3)
              call xsectsub(x2,h(m1),h(m2),xh(m1),xh(m2))
              call ysectsub(y2,h(m1),h(m2),yh(m1),yh(m2))
              goto 40
!
!     Case 7 - Line between sides 1-2 and 2-3
!
37            call xsectsub(x1,h(m1),h(m2),xh(m1),xh(m2))
              call ysectsub(y1,h(m1),h(m2),yh(m1),yh(m2))
              call xsectsub(x2,h(m2),h(m3),xh(m2),xh(m3))
              call ysectsub(y2,h(m2),h(m3),yh(m2),yh(m3))
              goto 40
!
!     Case 8 - Line between sides 2-3 and 3-1
!
38            call xsectsub(x1,h(m2),h(m3),xh(m2),xh(m3))
              call ysectsub(y1,h(m2),h(m3),yh(m2),yh(m3))
              call xsectsub(x2,h(m3),h(m1),xh(m3),xh(m1))
              call ysectsub(y2,h(m3),h(m1),yh(m3),yh(m1))
              goto 40
!
!     Case 9 - Line between sides 3-1 and 1-2
!
39            call xsectsub(x1,h(m3),h(m1),xh(m3),xh(m1))
              call ysectsub(y1,h(m3),h(m1),yh(m3),yh(m1))
              call xsectsub(x2,h(m1),h(m2),xh(m1),xh(m2))
              call ysectsub(y2,h(m1),h(m2),yh(m1),yh(m2))
40            continue
              XVAL(ii+1)=x1
              YVAL(ii+1)=y1
              XVAL(ii+2)=x2
              YVAL(ii+2)=y2
              ii=ii+2
            endif
60        continue
        endif
        nbc=ii

        end subroutine

!=====================================================================
! Function used in the conrec subroutine
 subroutine xsectsub(xsect,h1,h2,xh1,xh2)
        !from FUNCTION xsect(h1,h2,xh1,xh2)
        real(kind=nr), intent(out) :: xsect
        real(kind=nr), intent(in) :: h1,h2,xh1,xh2

        if (xh1.eq.xh2) then
          xsect=xh1
        else
          xsect=(h2*xh1-h1*xh2)/(h2-h1)
        endif

        END SUBROUTINE

!=============================================================
! Function used in the conrec subroutine
 subroutine ysectsub(ysect,h1,h2,yh1,yh2)
        !from FUNCTION ysect(h1,h2,yh1,yh2)
        real(kind=nr), intent(out) :: ysect
        real(kind=nr), intent(in) :: h1,h2,yh1,yh2

        if (yh1.eq.yh2) then
          ysect=yh1
        else
          ysect=(h2*yh1-h1*yh2)/(h2-h1)
        endif

        END SUBROUTINE

!===========================================================================
! orgseg is subroutine that orders the points, defined by their longitude and 
! latitude, found by conrec in the west-east direction
        SUBROUTINE orgseg(xseg,yseg,xxf,yyf,nbc,ndim,icompt)

        IMPLICIT NONE
        integer(kind=ni),INTENT(IN) :: ndim,nbc
        integer(kind=ni),INTENT(INOUT) :: icompt
        real(kind=nr),DIMENSION(nbc),INTENT(IN) :: xseg
        real(kind=nr),DIMENSION(nbc),INTENT(IN) :: yseg
        real(kind=nr),DIMENSION(ndim),INTENT(INOUT) :: xxf,yyf
        integer(kind=ni) kk,ll,rr,ptsauv,reste,reste1,llsauv,ptf,rrsuiv,ptsauvb
        integer(kind=ni) faux,icnew,icomptd
 !f2py depend(nbc) xseg,yseg
 !f2py depend(ndim) xxf,yyf

        if (icompt.eq.0) then
          DO ll=1,nbc
            if (xseg(ll).eq.(-180.0)) then
              xxf(1)=xseg(ll)
              yyf(1)=yseg(ll)
              ptsauv=ll
            endif        
          ENDDO
          call parite(ptsauv,ptsauvb)
          xxf(2)=xseg(ptsauvb)
          yyf(2)=yseg(ptsauvb)
          ptf=nbc/2+1
          if (ptf.gt.2) then
            icnew=icompt
            boucleorgpt : &
            DO kk=3,ptf
              faux=0
              DO ll=1,nbc
                if ((xseg(ll).eq.xxf(kk-1)).and.&
               (yseg(ll).eq.YYF(kk-1)).and.(ll.ne.ptsauvb)) then
                  call parite(ll,llsauv)
                  xxf(kk)=xseg(llsauv)
                  yyf(kk)=yseg(llsauv)
                  icnew=kk
                else
                  faux=faux+1
                endif
              ENDDO
              ptsauvb=llsauv
              if (faux.ge.nbc) then
                EXIT boucleorgpt
              endif
            ENDDO boucleorgpt
            if (faux.ge.nbc) then
              icompt=icnew
            else
              icompt=icompt+ptf
            endif
          endif
        else
          ptsauv=0
          DO rr=1,nbc
            if (((xseg(rr).eq.xxf(icompt)).or.&
            (xseg(rr).eq.(-xxf(icompt)))).and.(yseg(rr).eq.yyf(icompt))) then
              call parite(rr,rrsuiv)
              xxf(icompt+1)=xseg(rrsuiv)
              yyf(icompt+1)=yseg(rrsuiv)
              ptsauv=rrsuiv
            endif
          ENDDO
          if (ptsauv.eq.0) then
            goto 111
          endif
          ptf=nbc/2+1
          if (ptf.gt.2) then
            icnew = icompt
            boucleorgpts : &
            DO kk=icompt+2,icompt+ptf-1
              faux=0
              DO ll=1,nbc
                if ((xseg(ll).eq.xxf(kk-1)).and.&
                (yseg(ll).eq.YYF(kk-1)).and.(ll.ne.ptsauv)) then
                  call parite(ll,llsauv)
                  xxf(kk)=xseg(llsauv)
                  yyf(kk)=yseg(llsauv)
                  icnew=kk
                else
                  faux=faux+1
                endif
              ENDDO
              ptsauv=llsauv
              if (faux.ge.nbc) then
                EXIT boucleorgpts
              endif
            ENDDO boucleorgpts
            if (faux.ge.nbc) then      
              icompt=icnew
            else
              icompt=icompt+ptf-1
            endif
          endif
        endif
   
111     CONTINUE

        
        end subroutine

! ===============================================================================================
! PARITE is a subroutine checking if the number is even or odd
! It is used in the orgseg subroutine
        SUBROUTINE parite(nombrea,nombreb)

        IMPLICIT NONE
        integer(kind=ni),INTENT(IN) :: nombrea
        integer(kind=ni),INTENT(OUT) :: nombreb
        integer(kind=ni) :: reste

        reste=0
        reste=mod(nombrea,2)
        if (reste.eq.1) then
          nombreb=nombrea+1
        else
          nombreb=nombrea-1
        endif

        
        end subroutine


! =============================================================================================
! boitesuiv is a subroutine finding the following box where the contour is.
! It uses the last point found of the contour.
        SUBROUTINE boitesuiv(i,j,x,y,NXX,NYY,XVALF,YVALF,ilongl,ilongu,jlatl,jlatu)

        IMPLICIT NONE
        integer(kind=ni),INTENT(IN) :: i,j
        integer(kind=ni) :: inew
        integer(kind=ni),INTENT(OUT) :: ilongl,ilongu,jlatl,jlatu
        integer(kind=ni),INTENT(IN) :: NXX,NYY
        real(kind=nr),DIMENSION(2),INTENT(IN) :: x,y
        real(kind=nr),INTENT(IN) :: XVALF,YVALF
        integer(kind=ni) :: a,b


! 8 cases are possible, 0 is the initial box whose vertices are (i,j),(i+1,j),(i,j+1),(i+1,j+1)
!   |2|8|4|
!   |3|0|6|
!   |1|7|5|
!
        a=1
        b=1
! Case 1
        if ((XVALF.eq.x(a)).and.(YVALF.eq.y(b))) then
          if (i.eq.1) then
            inew=NXX+1
          else
            inew=i
          endif
          ilongl=inew-1
          ilongu=ilongl+1
          jlatl=j-1
          jlatu=jlatl+1
        else
! Case 2
          if ((XVALF.eq.x(a)).and.(YVALF.eq.y(b+1))) then
            if (i.eq.1) then
              inew=NXX+1
            else
              inew=i
            endif
            ilongl=inew-1
            ilongu=ilongl+1
            jlatl=j+1
            jlatu=jlatl+1
          else
! Case 3
            if (XVALF.eq.x(a)) then      
              if (i.eq.1) then
                inew=NXX+1
              else
                inew=i
              endif
              ilongl=inew-1
              ilongu=ilongl+1
              jlatl=j
              jlatu=jlatl+1
            else
! Case 4
              if ((XVALF.eq.x(a+1)).and.(YVALF.eq.y(b+1))) then
                if ((i+1).eq.(NXX+1)) then
                  inew=0
                else
                  inew=i
                endif
                ilongl=inew+1
                ilongu=ilongl+1
                jlatl=j+1
                jlatu=jlatl+1
              else
! Case 5
                if ((XVALF.eq.x(a+1)).and.(YVALF.eq.y(b))) then
                  if ((i+1).eq.(NXX+1)) then
                    inew=0
                  else
                    inew=i
                  endif
                  ilongl=inew+1
                  ilongu=ilongl+1
                  jlatl=j-1
                  jlatu=jlatl+1
                else
! Case 6
                  if (XVALF.eq.x(a+1)) then
                    if ((i+1).eq.(NXX+1)) then
                      inew=0
                    else
                      inew=i
                    endif
                    ilongl=inew+1
                    ilongu=ilongl+1
                    jlatl=j
                    jlatu=jlatl+1
                  else
! Case 7
                    if (YVALF.eq.y(b)) then
                      ilongl=i
                      ilongu=ilongl+1
                      jlatl=j-1
                      jlatu=jlatl+1
                    else
! Case 8
                      if (YVALF.eq.y(b+1)) then
                        ilongl=i
                        ilongu=ilongl+1
                        jlatl=j+1
                        jlatu=jlatl+1
                      endif
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif

        
        end subroutine
! ______________________________________________________________________
! ______________________________________________________________________
! This subroutine determines the points of a contour that belong to an
! AWB or to a CWB
! beta_a_temp is 1 when an AWB is detected at a grid point and 0 otherwise
! beta_c_temp is 1 when a CWB is detected at a grid point and 0 otherwise
        SUBROUTINE eventfunctions(icompt,NX,NY,longmin,latimin,incrx,incry, &
        XVAL,YVAL,beta_a_temp,beta_c_temp,nbpta_temp,nbptc_temp)

        integer(kind=ni),INTENT(IN) :: icompt,NX,NY
        real(kind=nr),INTENT(IN) :: longmin,latimin,incrx,incry
        real(kind=nr),DIMENSION(icompt),INTENT(IN) :: XVAL
        real(kind=nr),DIMENSION(icompt),INTENT(IN) :: YVAL
        integer(kind=ni),DIMENSION(NX,NY),INTENT(OUT) ::beta_a_temp,beta_c_temp
        integer(kind=ni),DIMENSION(icompt) :: flag1
        integer(kind=ni),INTENT(OUT) :: nbpta_temp,nbptc_temp
        integer(kind=ni) :: aa,bb
        integer(kind=ni) :: kk,ww,ptdpremier,indi,indj
        real(kind=nr) :: diffm1,diffp1
 !f2py depend(icompt) XVAL, YVAL, flag1
 !f2py depend(NX,NY) beta_a_temp,beta_c_temp 

        flag1=0
        aa=0
        bb=0
        beta_a_temp=0
        beta_c_temp=0
        nbpta_temp=0
        nbptc_temp=0
        DO kk=2,icompt-1
          if (flag1(kk).ne.2) then
            diffm1=XVAL(kk)-XVAL(kk-1)
            diffp1=XVAL(kk+1)-XVAL(kk)
            if ((XVAL(kk).eq.-180.0).and.(diffp1.gt.180.0)) then
              diffp1=diffp1-360.0
            endif
            if ((XVAL(kk-1).eq.-180.0).and.(diffm1.gt.180.0)) then
              diffm1=diffm1-360.0
            endif
            if ((diffm1.lt.0.0).and.(diffp1.lt.0.0)) then
              if (XVAL(kk-2).lt.XVAL(kk-1)) then
                ptdpremier=kk
                flag1(kk)=2
                if (YVAL(ptdpremier-1).le.YVAL(ptdpremier)) then
                  indi=int((XVAL(ptdpremier)-longmin)/incrx)+1
                  indj=int((YVAL(ptdpremier)-latimin)/incry)+1
                  if ((indi.gt.NX).or.(indj.gt.NY)) then
                    print *,'indi ou indj depassent NX ou NY'
                    STOP
                  endif
                  beta_c_temp(indi,indj)=1
                  aa=aa+1
                  ww=ptdpremier+1
                  diffm1=XVAL(ww)-XVAL(ww-1)
                  diffp1=XVAL(ww+1)-XVAL(ww)
                  if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                    diffp1=diffp1-360.0
                  endif
                  if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                    diffm1=diffm1-360.0
                  endif
                  DO WHILE((diffm1.lt.0.0).and.(diffp1.lt.0.0))
                    aa=aa+1
                    indi=int((XVAL(ww)-longmin)/incrx)+1
                    indj=int((YVAL(ww)-latimin)/incry)+1
                    if ((indi.gt.NX).or.(indj.gt.NY)) then
                      STOP
                    endif
                    beta_c_temp(indi,indj)=1
                    flag1(ww)=2
                    ww=ww+1
                    diffm1=XVAL(ww)-XVAL(ww-1)
                    diffp1=XVAL(ww+1)-XVAL(ww)
                    if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                      diffp1=diffp1-360.0
                    endif
                    if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                      diffm1=diffm1-360.0
                    endif
                  ENDDO
                else
                  if (YVAL(ptdpremier-1).gt.YVAL(ptdpremier)) then
                    indi=int((XVAL(ptdpremier)-longmin)/incrx)+1
                    indj=int((YVAL(ptdpremier)-latimin)/incry)+1
                    if ((indi.gt.NX).or.(indj.gt.NY)) then
                      STOP
                    endif
                    beta_a_temp(indi,indj)=1
                    bb=bb+1
                    ww=ptdpremier+1
                    diffm1=XVAL(ww)-XVAL(ww-1)
                    diffp1=XVAL(ww+1)-XVAL(ww)
                    if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                      diffp1=diffp1-360.0
                    endif
                    if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                      diffm1=diffm1-360.0
                    endif
                    DO WHILE((diffm1.lt.0.0).and.(diffp1.lt.0.0))
                      bb=bb+1
                      indi=int((XVAL(ww)-longmin)/incrx)+1
                      indj=int((YVAL(ww)-latimin)/incry)+1
                      beta_a_temp(indi,indj)=1
                      flag1(ww)=2
                      ww=ww+1
                      diffm1=XVAL(ww)-XVAL(ww-1)
                      diffp1=XVAL(ww+1)-XVAL(ww)
                      if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                        diffp1=diffp1-360.0
                      endif
                      if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                        diffm1=diffm1-360.0
                      endif
                    ENDDO
                  endif
                endif
              endif
            endif
          endif
        ENDDO
        nbptc_temp=aa
        nbpta_temp=bb

        
        END SUBROUTINE
end module

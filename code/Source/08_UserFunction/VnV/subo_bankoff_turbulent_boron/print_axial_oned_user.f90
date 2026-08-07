!
      SUBROUTINE print_axial_oned_user(time)
!
!     Print radial directional DATA for origin or excel graph
!
      USE Zconst1  , ONLY: vv_prob
!
      IMPLICIT NONE 
!
!.....Input
      REAL(8) :: time
!
      IF(vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d')THEN      
         CALL normalized_axial_out(time)
      ELSEIF(vv_prob.eq.'boron_trans')THEN      
         CALL normalized_axial_out(time)
      ELSE   
         RETURN
      ENDIF
!      
      END SUBROUTINE print_axial_oned_user
!      
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE normalized_axial_out(time)
!
!     Save lateral output for Origin
!
      USE VOL_DATA   , ONLY: cell            
      USE Zzone      , ONLY: ncell_fluid
      USE Zparam     , ONLY: ndim
      USE Zcore      , ONLY: np,myrank      
      USE Zconst1    , ONLY: vv_prob
      USE Zcoord1    , ONLY: xloc
      USE Zcoord2    , ONLY: xloc_xfc_min,xloc_xfc_max
      USE Ziat       , ONLY: ia
      USE Zndforce   , ONLY: d_bfc  
      USE Zturb      , ONLY: turb_ke,turb_dp,turb_keg,turb_dpg,pro_ke,pro_keg,yplus,yplusg,utau
      USE Zvector    , ONLY: vl_n,vg_n
      USE Zio_unit   , ONLY: unit_log
!
      IMPLICIT NONE 
!
!.....Size for prnvar along second dimension
      INTEGER, PARAMETER :: isize=19
!.....Input
      REAL(8) :: time
!.....Local variables
      INTEGER :: i,j,k,level
      INTEGER :: ip,i0,j0,ix
      INTEGER :: loop,nloop,nloop_l
      INTEGER :: filenumb
      INTEGER,SAVE :: nradius
      LOGICAL,SAVE :: INITIAL=.TRUE.
      CHARACTER*20 :: filename,lastname
      CHARACTER*20,SAVE :: varname(0:17)
      REAL(8) :: hdistance
      REAL(8) :: highloc,lowloc
      REAL(8),SAVE :: h1,h2
!.....Local arrays
      INTEGER :: iaa(np+1)
      INTEGER :: nloop_all(np),nloop_dsp(np)
      INTEGER :: nloopz_all(np),nloopz_dsp(np)
      REAL(8),SAVE :: x(2),y(2)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: icell_index_l,index
      REAL(8),DIMENSION(:),ALLOCATABLE :: hdistance_l
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar_l
!                  
      IF(INITIAL) THEN
         h1=0.0d0
         h2=1.0d0
!         
         IF(vv_prob.eq.'subo')THEN
            h1=0.0d0
            h2=3.883d0
            nradius = 2
            x(1)=0.005d0                                  !center  
            x(2)=0.017d0                                  !out-wall
            y(1)=0.0d0
            y(2)=0.0d0
         ELSEIF(vv_prob.eq.'bankoff')THEN
         ELSEIF(vv_prob.eq.'turbulent')THEN               !(0.4)x0.4x4m
            h1=0.0d0
            h2=40.0d0
            nradius = 2
            x(1)=0.0001d0                                 !left-wall
            x(2)=0.2d0-0.0001d0                           !center
            IF(ndim.eq.3)THEN
               y(1)=0.005d0/2.0d0                         !center
               y(2)=0.005d0/2.0d0                         !center
            ENDIF  
         ELSEIF(vv_prob.eq.'turbulent3d')THEN             !(0.4)x0.4x4m
            h1=0.0d0
            h2=40.0d0
            nradius = 2
            x(1)=0.0001d0                                 !left-wall
            x(2)=0.2d0-0.0001d0                           !center
            IF(ndim.eq.3)THEN
               y(1)=0.4d0/2.0d0                           !center
               y(2)=0.4d0/2.0d0                           !center
            ENDIF  
         ELSEIF(vv_prob.eq.'boron_trans')THEN
            h1=0.0d0
            h2=5.0d0
            nradius = 2
            x(1)=0.0001d0                                 !left-wall
            x(2)=0.25d0-0.0001d0                          !center
         ELSEIF(vv_prob.eq.'vawl')THEN
         ELSEIF(vv_prob.eq.'Hibiki')THEN
         ELSE
            WRITE(*,*)'           normalized_radial_out is not available for',vv_prob 
            WRITE(unit_log,*)'           normalized_radial_out is not available for',vv_prob 
            STOP
         ENDIF            
!      
         IF(myrank.eq.0) THEN
!
            varname(0)='0:Time'
            varname(1)='1:DistanceNR'  
            varname(2)='2:FractionG'
            varname(3)='3:VelocityL'                          !bankoff
            varname(4)='4:VelocityG'  
            varname(5)='5:IAC' 
            varname(6)='6:DiameterB'  
            varname(7)='7:TempLiq'                            !subo
!      
            varname(8)='8:YplusLiq'                           !turbulent 
            varname(9)='9:KeLiq'
            varname(10)='10:DpLiq'
            varname(11)='11:TvisLiq'
            varname(12)='12:TproLiq'
!      
            varname(13)='13:YplusGas'                         !turbulent   
            varname(14)='14:KeGas'
            varname(15)='15:DpGas'
            varname(16)='16:TvisGas'
            varname(17)='17:TproGas'                          !turbulent 
!
            IF(vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d') THEN
            ELSEIF(vv_prob.eq.'boron_trans')then
               varname(2)='2:BoronCon'
               varname(3)='3:TempLiq'
               varname(4)='4:DensLiq'                         !bankoff
            ENDIF         
         ENDIF         
         INITIAL=.FALSE.
      ENDIF
!
      ALLOCATE(icell_index_l(ncell_fluid),hdistance_l(ncell_fluid))
      DO level = 1, nradius 
         ip=myrank+1
         nloop_l=0
         DO i=1, ncell_fluid
!
!...........check the x-coordinate in a 2-D geometry
!
            highloc=xloc_xfc_max(i,1)
            lowloc=xloc_xfc_min(i,1)
            IF(highloc.le.x(level).or.lowloc.gt.x(level)) CYCLE
            hdistance_l(nloop_l+1)=xloc(i,2)  
!
!...........check the x-coordinate in a 3-D geometry
!
            IF(ndim.eq.3)then
               highloc=xloc_xfc_max(i,2)
               lowloc=xloc_xfc_min(i,2)
               IF(highloc.le.y(level).or.lowloc.gt.y(level)) CYCLE
               hdistance_l(nloop_l+1)=xloc(i,3)
            ENDIF        
            nloop_l=nloop_l+1
            icell_index_l(nloop_l)=i
         ENDDO
         IF(ALLOCATED(prnvar_l)) DEALLOCATE(prnvar_l)
         ALLOCATE(prnvar_l(nloop_l,isize))
         IF(vv_prob.eq.'boron_trans')then
            DO loop=1,nloop_l
               i=icell_index_l(loop)
               hdistance=hdistance_l(loop)
               prnvar_l(loop, 1)=(hdistance-h1)/(h2-h1)
               prnvar_l(loop, 2)=cell%cboron(i)
               prnvar_l(loop, 3)=cell%tl(i)
               prnvar_l(loop, 4)=cell%rhol(i)                             
            ENDDO
         ELSE
            DO loop=1,nloop_l
               i=icell_index_l(loop)
               hdistance=hdistance_l(loop)
               prnvar_l(loop, 1)=(hdistance-h1)/(h2-h1)
               prnvar_l(loop, 2)=cell%alphag(i)
               prnvar_l(loop, 3)=vl_n(i,ndim)                 
               prnvar_l(loop, 4)=vg_n(i,ndim)
               prnvar_l(loop, 5)=ia(i)
               prnvar_l(loop, 6)=cell%d1(i)*1000.0d0          !meter to mili-meter
               prnvar_l(loop, 7)=cell%tl(i)                   !subo
   !            
               prnvar_l(loop, 8)=yplus(i)
               prnvar_l(loop, 9)=turb_ke(i)
               prnvar_l(loop,10)=turb_dp(i)
               prnvar_l(loop,11)=cell%tviscosl(i)
               prnvar_l(loop,12)=pro_ke(i)
   !                                                          
               prnvar_l(loop,13)=yplusg(i)
               prnvar_l(loop,14)=turb_keg(i)
               prnvar_l(loop,15)=turb_dpg(i)
               prnvar_l(loop,16)=cell%tviscosg(i)
               prnvar_l(loop,17)=pro_keg(i)
   !
               prnvar_l(loop,18)=utau(i)
               prnvar_l(loop,19)=d_bfc(i)*cell%rhol(i)/cell%lviscosl(i)
            ENDDO
         ENDIF
         CALL allgather_i(nloop_l,nloop_all)
         iaa(1)=1
         DO ip=1,np
            iaa(ip+1)=iaa(ip)+nloop_all(ip)
            nloopz_all(ip)=nloop_all(ip)*isize
         ENDDO
         nloop_dsp(1)=0
         nloopz_dsp(1)=0
         DO ip=2,np
            nloop_dsp(ip)=nloop_dsp(ip-1)+nloop_all(ip-1)
            nloopz_dsp(ip)=nloopz_dsp(ip-1)+nloopz_all(ip-1)
         ENDDO
         nloop=iaa(np+1)-1
!
         IF(myrank.eq.0)THEN
            ALLOCATE(prnvar_all(nloop*isize))
         ELSE
            ALLOCATE(prnvar_all(1))
         ENDIF
         CALL gather_vec_r(prnvar_l,nloop_l*isize,prnvar_all,nloop*isize,nloopz_all,nloopz_dsp)
!
         IF(myrank.eq.0)THEN
            ALLOCATE(prnvar(nloop,isize),prnvar1(nloop),index(nloop))
            DO ip=1,np
               j0=nloopz_dsp(ip)
               i0=1
               DO i=iaa(ip),iaa(ip+1)-1
                  DO ix=1,isize
                     prnvar(i,ix)=prnvar_all(i0+j0+(ix-1)*nloop_all(ip))
                  ENDDO
                  i0=i0+1
               ENDDO
            ENDDO
         ENDIF
         DEALLOCATE(prnvar_all)
!
         IF(myrank.eq.0) THEN
!
!........Sort the extracted DATA      
!
            CALL sortx_r(prnvar(1,1),index,nloop)
            DO k=2,isize
               DO i=1,nloop
                  j=index(i)
                  prnvar1(i)=prnvar(j,k)
               ENDDO
               DO i=1,nloop
                  prnvar(i,k)=prnvar1(i)
               ENDDO
            ENDDO
!
!...........Print steady data at each level to check version    
!
            filenumb=70+level
            WRITE(lastname,*)level
            lastname=ADJUSTL(lastname)
            IF(vv_prob.eq.'subo')THEN         
               filename='VD3_rn'//TRIM(lastname)//'_vv.dat' 
               OPEN(filenumb,file=filename) 
               WRITE(filenumb,5000)(varname(k),k=1,7)
            ELSEIF(vv_prob.eq.'bankoff')THEN
               filename='VD4_rn'//TRIM(lastname)//'_vv.dat' 
               OPEN(filenumb,file=filename) 
               WRITE(filenumb,5000)(varname(k),k=1,3)
            ELSEIF(vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d')THEN
               filename='VFS8_rn'//TRIM(lastname)//'_vv.dat' 
               OPEN(filenumb,file=filename) 
               WRITE(filenumb,5000)varname(1),varname(3),(varname(k),k=8,12)
               DO j=1,nloop
                  prnvar(j,8)=prnvar(1,18)*prnvar(j,19)    !ut(1)*y(j)*rhol(j)/mu(j)
                  WRITE(filenumb,5001)prnvar(j,1),prnvar(j,3),(prnvar(j,k),k=8,12)
               ENDDO
            ELSEIF(vv_prob.eq.'boron_trans')THEN
               filename='VFS12_rn'//TRIM(lastname)//'_vv.dat' 
               OPEN(filenumb,file=filename) 
               WRITE(filenumb,5000)(varname(k),k=1,4)
               DO j=1,nloop
                  WRITE(filenumb,5001)(prnvar(j,k),k=1,4)
               ENDDO
            ENDIF   
!         
            WRITE(filenumb,5002)level,time,np
            CLOSE(filenumb)
!            
            DEALLOCATE(prnvar1,index)
            DEALLOCATE(prnvar)
         ENDIF
      ENDDO !level
      DEALLOCATE(icell_index_l,hdistance_l)
!         
5000  FORMAT(1x,20(A15))
5001  FORMAT(1x,20(1pe15.5))      
5002  FORMAT(1x,'level=',1I3,1x,'time=',1pe15.5,1x,'np=',1I3)
!
      END SUBROUTINE normalized_axial_out

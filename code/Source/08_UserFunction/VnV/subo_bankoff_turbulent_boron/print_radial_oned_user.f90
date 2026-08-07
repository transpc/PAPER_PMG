!
      SUBROUTINE print_radial_oned_user(time)
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
      IF(vv_prob.eq.'suboex')THEN
         CALL local_out_subo(time)
      ELSEIF(vv_prob.eq.'becker')THEN
         CALL local_out_becker
      ELSEIF(vv_prob.eq.'subo'.or.vv_prob.eq.'bankoff'.or.vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d')THEN      
         CALL normalized_lateral_out(time)
      ELSE  
         RETURN 
      ENDIF
!      
      END SUBROUTINE print_radial_oned_user
!      
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE normalized_lateral_out(time)
!
!     Save lateral output for Origin
!
      USE VOL_DATA   , ONLY: cell            
      USE Zzone      , ONLY: ncell_fluid
      USE Zcore      , ONLY: np,myrank      
      USE Zparam     , ONLY: ndim
      USE Zconst1    , ONLY: vv_prob
      USE Zcoord1    , ONLY: xloc
      USE Zcoord2    , ONLY: xloc_xfc_min,xloc_xfc_max
      USE Ziat       , ONLY: ia
      USE Zndforce   , ONLY: d_bfc  
      USE Zturb      , ONLY: turb_ke,turb_dp,turb_keg,turb_dpg,pro_ke,pro_keg,yplus,yplusg,utau,tauw,utaug,tauwg
      USE Zvector    , ONLY: vl_n,vg_n
      USE Zio_unit   , ONLY: unit_log
!
      IMPLICIT NONE 
!
!.....Size for prnvar
      INTEGER, PARAMETER :: isize=23
!.....Input
      REAL(8) :: time
!.....Local variables
      INTEGER :: i,j,k,level
      INTEGER :: ip,i0,j0,ix
      INTEGER :: loop,nloop_l,nloop
      INTEGER :: filenumb
      INTEGER,SAVE :: Nheight
      LOGICAL,SAVE :: INITIAL=.TRUE.
      LOGICAL,SAVE :: INITIAL_p=.TRUE.
      CHARACTER*20 :: filename,lastname
      CHARACTER*20,SAVE :: varname(0:23)
      REAL(8) :: rdistance
      REAL(8) :: hzloc,lzloc
      REAL(8) :: ymax,ymin,temp1     
      REAL(8),SAVE :: ymean,xc,yc
      REAL(8),SAVE :: r1,r2
!.....Local arrays
      INTEGER :: iaa(np+1)
      INTEGER :: nloop_all(np),nloop_dsp(np)
      INTEGER :: nloopz_all(np),nloopz_dsp(np)
      REAL(8),SAVE :: height(6)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: icell_index_l,index
      REAL(8),DIMENSION(:),ALLOCATABLE :: rdistance_l
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar_l
!      
      IF(INITIAL)THEN
         xc=0.0d0
         yc=0.0d0
         IF(vv_prob.eq.'subo')THEN
            nheight = 6
            r1=0.00499d0
            r2=0.01775d0
            height(1)=0.47d0  
            height(2)=1.1053d0  
            height(3)=1.741d0  
            height(4)=2.4211d0 
            height(5)=2.9983d0
            height(6)=3.1486d0  
         ELSEIF(vv_prob.eq.'bankoff')THEN
            nheight = 2
            r1=0.0d0
            r2=0.019d0
            height(1)=0.1d0
            height(2)=1.368d0  
         ELSEIF(vv_prob.eq.'turbulent')THEN                           !(0.4)x0.4x4m
            nheight = 2
            r1=0.0d0
            r2=0.4d0                                                  !40.0d0
            height(1)=26.7d0                                          !2/3 
            height(2)=30.1d0                                          !3/4
            IF(ndim.eq.3)THEN
               xc=0.0d0
               yc=0.005d0/2.0d0
            ENDIF                 
         ELSEIF(vv_prob.eq.'turbulent3d')THEN                         !(0.4)x0.4x4m
            nheight = 2
            r1=0.0d0
            r2=0.4d0                                                  !40.0d0
            height(1)=2.67d1                                          !2/3 
            height(2)=3.00d1                                          !3/4
            IF(ndim.eq.3)THEN
               xc=0.0d0
               yc=0.2d0
            ENDIF 
         ELSEIF(vv_prob.eq.'vawl')THEN
            nheight = 3
            r1=0.0d0
            r2=0.04d0
            height(1)=0.01d0
            height(2)=3.279d0
            height(3)=7.959d0
         ELSEIF(vv_prob.eq.'Hibiki')THEN
            nheight = 2
            r1=0.0d0
            r2=0.0254d0
            height(1)=0.305d0
            height(2)=2.718d0
         ELSE
            WRITE(*,*) '           normalized_lateral_out is not available for',vv_prob 
            WRITE(unit_log,*)'           normalized_lateral_out is not available for',vv_prob 
            CALL finalize_mpi
            STOP
         ENDIF            
!      
         IF(myrank.eq.0) THEN
            varname( 0)='0:Time'
            varname( 1)='1:DistanceNR'  
            varname( 2)='2:FractionG'
            varname( 3)='3:VelocityL'                            !bankoff
            varname( 4)='4:VelocityG'  
            varname( 5)='5:IAC' 
            varname( 6)='6:DiameterB'  
            varname( 7)='7:TempLiq'                              !subo
!      
            varname( 8)='8:YplusLiq'                             !turbulent 
            varname( 9)='9:KeLiq'
            varname(10)='10:DpLiq'
            varname(11)='11:TvisLiq'
            varname(12)='12:TproLiq'
            varname(13)='13:U_tau'                               !turbulent   
            varname(14)='14:Y_plus'
            varname(15)='15:Tau_wall' 
!
            varname(16)='16:YplusGas'                            !turbulent 
            varname(17)='17:KeGas'
            varname(18)='18:DpGas'
            varname(19)='19:TvisGas'
            varname(20)='20:TproGas'
            varname(21)='21:U_taug'                              !turbulent   
            varname(22)='22:Y_plusg'
            varname(23)='23:Tau_wallg' 
         ENDIF
!
         IF(vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d')varname(3)='3:U+'
!
         IF(vv_prob.eq.'turbulent3d')THEN ! The number of y-directional mesh should be odd!
            ymin=huge(0.d0)
            ymax=-ymin
            DO i=1,ncell_fluid
               ymin=MIN(ymin,xloc(i,2))
               ymax=MAX(ymax,xloc(i,2))
            ENDDO
            CALL allreducei_min_r1(ymin)
            CALL allreducei_max_r1(ymax)
            ymean=(ymin+ymax)/2.0d0
         ENDIF
         INITIAL=.FALSE.
      ENDIF
!      
      ALLOCATE(icell_index_l(ncell_fluid),rdistance_l(ncell_fluid))
      DO level=1,Nheight 
         ip=myrank+1
         nloop_l=0
         DO i=1,ncell_fluid
!
!..............Bypass if not a certain XY 2-D plane in a 3-D geometry
!
               IF(vv_prob.eq.'turbulent3d')THEN            
                  IF(xloc(i,2).gt.ymean+1.0d-3.or.xloc(i,2).lt.ymean-1.0d-3) CYCLE
               ENDIF
             hzloc=xloc_xfc_max(i,ndim)
             lzloc=xloc_xfc_min(i,ndim)
!
!..............Check the height
!
               IF(hzloc.le.height(level).or.lzloc.gt.height(level)) CYCLE
!
!..............Check the radius
!
               rdistance=ABS(xloc(i,1)-xc)
               IF(ndim.eq.3) rdistance=SQRT((xloc(i,1)-xc)**2+(xloc(i,2)-yc)**2)    
               IF(rdistance.lt.r1.or.rdistance.gt.r2) CYCLE
               nloop_l=nloop_l+1 
               icell_index_l(nloop_l)=i
               rdistance_l(nloop_l)=rdistance
         ENDDO
         IF(ALLOCATED(prnvar_l)) DEALLOCATE(prnvar_l)
         ALLOCATE(prnvar_l(nloop_l,isize))
         DO loop=1,nloop_l
            i=icell_index_l(loop)
            rdistance=rdistance_l(loop)
!
            prnvar_l(loop, 1)=(rdistance-r1)/(r2-r1)
            prnvar_l(loop, 2)=cell%alphag(i)
            prnvar_l(loop, 3)=vl_n(i,ndim) !bankoff
            prnvar_l(loop, 4)=vg_n(i,ndim)
            prnvar_l(loop, 5)=ia(i)
            prnvar_l(loop, 6)=cell%d1(i)*1000.0d0 !meter to mili-meter
            prnvar_l(loop, 7)=cell%tl(i) !subo
!
            prnvar_l(loop, 8)=yplus(i)
            prnvar_l(loop, 9)=turb_ke(i)
            prnvar_l(loop,10)=turb_dp(i)
            prnvar_l(loop,11)=cell%tviscosl(i)
            prnvar_l(loop,12)=pro_ke(i)
            prnvar_l(loop,13)=utau(i)
            prnvar_l(loop,14)=d_bfc(i)*cell%rhol(i)/cell%lviscosl(i)
            prnvar_l(loop,15)=tauw(i)
!
            prnvar_l(loop,16)=yplusg(i)
            prnvar_l(loop,17)=turb_keg(i)
            prnvar_l(loop,18)=turb_dpg(i)
            prnvar_l(loop,19)=cell%tviscosg(i)
            prnvar_l(loop,20)=pro_keg(i)
            prnvar_l(loop,21)=utaug(i)
            prnvar_l(loop,22)=d_bfc(i)*cell%rhog(i)/cell%lviscosg(i)
            prnvar_l(loop,23)=tauwg(i)
!
         ENDDO
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
!...........Sort the exctracted DATA      
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
!         
            IF(vv_prob.eq.'subo')THEN         
               filename='VD3_hn'//TRIM(lastname)//'_vv.dat' 
               OPEN(filenumb,file=filename) 
               WRITE(filenumb,5000)(varname(k),k=1,7)
               DO j=1,nloop
                  WRITE(filenumb,5001)(prnvar(j,k),k=1,7)
               ENDDO
            ELSEIF(vv_prob.eq.'bankoff')THEN
               filename='VD4_hn'//TRIM(lastname)//'_vv.dat' 
               OPEN(filenumb,file=filename) 
               WRITE(filenumb,5000)(varname(k),k=1,3)
               DO j=1,nloop
                  WRITE(filenumb,5001)(prnvar(j,k),k=1,3)
               ENDDO
            ELSEIF(vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d')THEN
               filename='VFS8_hn'//TRIM(lastname)//'_vv.dat' 
               OPEN(filenumb,file=filename) 
               WRITE(filenumb,5000)varname(1),varname(8),'U_plus    ',varname(3),varname(13),varname(15),(varname(k),k=9,12)   !Liquid
!               WRITE(filenumb,5000)varname(1),varname(16),'U_plus    ',varname(4),varname(21),varname(23),(varname(k),k=17,20)   !Gas
               DO j=1,nloop
!Liquid            
                  prnvar(j,8)=prnvar(1,13)*prnvar(j,14)                        !ut(1)*y(j)*rhol(j)/mu(j)
                  !prnvar(j,3)=prnvar(j,3)/prnvar(1,18)                         !U(j)/ut(1)
                  temp1=prnvar(j,3)/prnvar(1,13) 
                  !                         r        Y+        U+     vel            U_tau        tau_w          ke   dp   k_vis   k_pro
                  WRITE(filenumb,5001)prnvar(j,1),prnvar(j,8),temp1,prnvar(j,3),prnvar(j,13),prnvar(j,15),(prnvar(j,k),k=9,12)
!GAS
                  prnvar(j,16)=prnvar(1,21)*prnvar(j,22)                        !y+=ut(1)*y(j)*rhol(j)/mu(j)
!                  !prnvar(j,4)=prnvar(j,4)/prnvar(1,21)                         !U+=U(j)/ut(1)
!                  temp1=prnvar(j,4)/prnvar(1,21) 
!                  WRITE(filenumb,5001)prnvar(j,1),prnvar(j,16),temp1,prnvar(j,4),prnvar(j,21),prnvar(j,23),(prnvar(j,k),k=17,20)
               ENDDO
            ELSE
            ENDIF   
!         
            WRITE(filenumb,5002)level,time,np
            CLOSE(filenumb)
!            
            DEALLOCATE(prnvar1,index)
            IF(level.lt.Nheight) DEALLOCATE(prnvar)
         ENDIF 
      ENDDO !level
      DEALLOCATE(icell_index_l,rdistance_l)
!
      IF(myrank.eq.0) THEN
!
!........Print transient at top and center to check parallel computing
!
         IF(INITIAL_p)THEN
            INITIAL_p=.FALSE.
!         
            IF(vv_prob.eq.'subo')THEN            
               OPEN(69,file='VD3_ref.dat') 
               WRITE(69,5000)(varname(k),k=0,7)
            ELSEIF(vv_prob.eq.'bankoff')THEN               
               OPEN(69,file='VD4_ref.dat') 
               WRITE(69,5000)(varname(k),k=0,3)
            ELSEIF(vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d')THEN
               OPEN(69,file='VFS8_ref.dat') 
               WRITE(69,5000)(varname(k),k=0,1),varname(3),(varname(k),k=8,12)
            ENDIF
         ENDIF 
!      
         IF(vv_prob.eq.'subo')THEN            
            WRITE(69,5001)time,(prnvar(1,k),k=1,7) 
         ELSEIF(vv_prob.eq.'bankoff')THEN               
            WRITE(69,5001)time,(prnvar(1,k),k=1,3) 
         ELSEIF(vv_prob.eq.'turbulent'.or.vv_prob.eq.'turbulent3d')THEN    
            prnvar(1,8)=prnvar(1,18)*prnvar(1,19)                           !ut(1)*y(j)*rhol(j)/mu(j)
           !prnvar(1,3)=prnvar(1,3)/prnvar(1,18)                          !U(1)/ut(1)
! bug j ???
!          temp1=prnvar(j,3)/prnvar(1,13)                                 !U/ut
            temp1=prnvar(1,3)/prnvar(1,13)                                 !U/ut
            WRITE(69,5001)time,prnvar(1,1),temp1,(prnvar(1,k),k=8,12) 
         ENDIF
!
         DEALLOCATE(prnvar)
!
      ENDIF
!         
5000  FORMAT(1x,20(A15))
5001  FORMAT(1x,20(1pe15.5))      
5002  FORMAT(1x,'level=',1I3,1x,'time=',1pe15.5,1x,'np=',1I3)
!
      END SUBROUTINE normalized_lateral_out
!      
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE local_out_subo(time)
!
!     Save lateral output for Origin
!
      USE VOL_DATA , ONLY: cell
      USE Zmpi     , ONLY: ncell_fp
      USE Zparam   , ONLY: ndim
      USE Zconst1  , ONLY: vv_prob
      USE Zcoord1  , ONLY: xloc_tmp
      USE Zcore    , ONLY: myrank      
      USE Zcoord2  , ONLY: xloc_xfc_min,xloc_xfc_max
      USE Ziat     , ONLY: ia,ia_conv,iat_break,iat_coal,iat_nucl,iat_size
      USE Zqvol    , ONLY: h_ig,h_il       
      USE Zturb    , ONLY: turb_ke,turb_dp,turb_keg,turb_dpg,pro_ke,pro_keg
      USE Zvector  , ONLY: vl_n,vg_n
      USE Zzone    , ONLY: ncell_fluid,ncell_fluid_all   
!
      IMPLICIT NONE 
!
!     size for prnvar
      INTEGER, PARAMETER :: isize1=17
!.....Input
      REAL(8) :: time
!.....Local variables
      INTEGER :: i, j, k, na, wopt, level
      INTEGER :: iradius
      INTEGER :: Nheight
      INTEGER :: loop,nloop
      LOGICAL,SAVE :: INITIAL_p=.TRUE.
      CHARACTER*16 sradius
      REAL(8) alphag_sum, ia_sum, vg_sum, dsm_sum
      REAL(8) rdistance
      REAL(8) height(6), r1, r2
      REAL(8) hzloc,lzloc      
!.....Local allocatable arrays            
      CHARACTER*16,DIMENSION(:,:),ALLOCATABLE ::  varname
      INTEGER,DIMENSION(:),ALLOCATABLE :: icell_index,index
      REAL(8),DIMENSION(:),ALLOCATABLE :: rdistance_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: pr,ag,diab,cboron
      REAL(8),DIMENSION(:),ALLOCATABLE :: hil,hig,D1,eviscosl,lviscosl
      REAL(8),DIMENSION(:),ALLOCATABLE :: iae,iae_conv,iae_break,iae_coal,iae_nucl,iae_size
      REAL(8),DIMENSION(:),ALLOCATABLE :: tk,td,tp,tlq
      REAL(8),DIMENSION(:),ALLOCATABLE :: hzloc1_all,lzloc1_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: vlp,vgp
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar0
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
!
      na=ncell_fluid_all
      IF(myrank.eq.0) THEN
         ALLOCATE(pr(na),ag(na),diab(na),cboron(na))
         ALLOCATE(hil(na),hig(na),D1(na))
         ALLOCATE(iae(na),iae_conv(na),iae_break(na),     &
                  iae_coal(na),iae_nucl(na),iae_size(na))
!
         ALLOCATE(eviscosl(na),lviscosl(na))
         ALLOCATE(tk(na),td(na),tp(na),tlq(na))
         ALLOCATE(hzloc1_all(ncell_fluid_all),lzloc1_all(ncell_fluid_all))
         ALLOCATE(vlp(ncell_fluid_all,ndim),vgp(ncell_fluid_all,ndim))
      ELSE
         ALLOCATE(pr(1),ag(1),diab(1),cboron(1))
         ALLOCATE(hil(1),hig(1),D1(1))
         ALLOCATE(iae(1),iae_conv(1),iae_break(1),     &
                  iae_coal(1),iae_nucl(1),iae_size(1))
!
         ALLOCATE(eviscosl(1),lviscosl(1))
         ALLOCATE(tk(1),td(1),tp(1),tlq(1))
         ALLOCATE(hzloc1_all(1),lzloc1_all(1))
         ALLOCATE(vlp(1,1),vgp(1,1))
      ENDIF
!          
      CALL gatherv_r_2d(vl_n,ncell_fp,vlp,ncell_fluid,na,0)
      CALL gatherv_r_2d(vg_n,ncell_fp,vgp,ncell_fluid,na,0)
!
      CALL gatherv_r(cell%p     ,ncell_fluid,pr       ,na,0)
      CALL gatherv_r(cell%alphag,ncell_fluid,ag       ,na,0)
      CALL gatherv_r(cell%d1    ,ncell_fluid,diab     ,na,0)
      CALL gatherv_r(cell%cboron,ncell_fluid,cboron   ,na,0)
!
      CALL gatherv_r(H_il       ,ncell_fluid,hil      ,na,0)
      CALL gatherv_r(H_ig       ,ncell_fluid,hig      ,na,0)
      CALL gatherv_r(cell%D1    ,ncell_fluid,D1       ,na,0)
!
      CALL gatherv_r(ia         ,ncell_fluid,iae      ,na,0)
      CALL gatherv_r(ia_conv    ,ncell_fluid,iae_conv ,na,0)
      CALL gatherv_r(iat_break  ,ncell_fluid,iae_break,na,0)
      CALL gatherv_r(iat_coal   ,ncell_fluid,iae_coal ,na,0)
      CALL gatherv_r(iat_nucl   ,ncell_fluid,iae_nucl ,na,0)
      CALL gatherv_r(iat_size   ,ncell_fluid,iae_size ,na,0)
!
      IF(.true.)THEN  ! print liquid phase data
         CALL gatherv_r(cell%eviscosl,ncell_fluid,eviscosl,na,0)
         CALL gatherv_r(cell%lviscosl,ncell_fluid,lviscosl,na,0)
         CALL gatherv_r(turb_ke      ,ncell_fluid,tk      ,na,0)
         CALL gatherv_r(turb_dp      ,ncell_fluid,td      ,na,0)
         CALL gatherv_r(pro_ke       ,ncell_fluid,tp      ,na,0)
         CALL gatherv_r(cell%tl      ,ncell_fluid,tlq     ,na,0)
      ELSE       ! print gas phase data
         CALL gatherv_r(cell%eviscosg,ncell_fluid,eviscosl,na,0)
         CALL gatherv_r(cell%lviscosg,ncell_fluid,lviscosl,na,0)
         CALL gatherv_r(turb_keg     ,ncell_fluid,tk      ,na,0)
         CALL gatherv_r(turb_dpg     ,ncell_fluid,td      ,na,0)
         CALL gatherv_r(pro_keg      ,ncell_fluid,tp      ,na,0)
      ENDIF
!
      IF(ndim.eq.2) THEN
         CALL gatherv_r(xloc_xfc_max(1,2),ncell_fluid,hzloc1_all,ncell_fluid_all,0)
         CALL gatherv_r(xloc_xfc_min(1,2),ncell_fluid,lzloc1_all,ncell_fluid_all,0)
      ELSE
         CALL gatherv_r(xloc_xfc_max(1,3),ncell_fluid,hzloc1_all,ncell_fluid_all,0)
         CALL gatherv_r(xloc_xfc_min(1,3),ncell_fluid,lzloc1_all,ncell_fluid_all,0)
      ENDIF
!
      IF(myrank.eq.0) THEN
!      
         IF(vv_prob.eq.'suboex')THEN
            nheight = 6
            r1=0.00499d0
            r2=0.01775d0
            height(1)=0.47d0
            height(2)=1.1053d0
            height(3)=1.741d0
            height(4)=2.4211d0 
            height(5)=2.9983d0
            height(6)=3.1486d0
         ELSE
            nheight = 0
         ENDIF
!
         ALLOCATE(icell_index(na),rdistance_all(na))
         DO level=1,Nheight 
            loop=0
            alphag_sum=0.0d0
            ia_sum=0.0d0
            vg_sum=0.0d0
            dsm_sum=0.0d0
!        
            nloop=0
            DO i=1,ncell_fluid_all
               hzloc=hzloc1_all(i)
               lzloc=lzloc1_all(i)
!
!..............Check the height
!
               IF(hzloc.le.height(level).or.lzloc.gt.height(level)) CYCLE
!
!..............Check the radius
!
               wopt = 0
               rdistance=dabs(xloc_tmp(i,1))
               IF(ndim.eq.3) rdistance=SQRT(xloc_tmp(i,1)**2+xloc_tmp(i,2)**2)
               IF(rdistance.ge.r1.and.rdistance.le.r2) wopt=1
!
!..............Extract the DATA to be printed if wopt ==1
!
               IF(wopt.ne.1) CYCLE                                        !when height and radius is in proper value
               nloop=nloop+1
               icell_index(nloop)=i
               rdistance_all(nloop)=rdistance
            ENDDO
            ALLOCATE(prnvar(nloop,isize1),prnvar0(nloop))
            prnvar(:,:)=0.d0
            DO loop=1,nloop
               i=icell_index(loop)
               rdistance=rdistance_all(loop)
               prnvar(loop,1)=(rdistance-r1)/(r2-r1)
               prnvar(loop,2)=ag(i)
               prnvar(loop,3)=iae(i)
               prnvar(loop,4)=diab(i)*1000.0d0
               prnvar(loop,5)=vlp(i,ndim)  
               prnvar(loop,6)=vgp(i,ndim)
               prnvar(loop,7)=pr(i)
               prnvar(loop,8)=cboron(i)  
               prnvar(loop,9)=eviscosl(i)/lviscosl(i)
               prnvar(loop,6)=tk(i)
               prnvar(loop,7)=td(i)
               prnvar(loop,8)=tp(i)
               prnvar(loop,9)=tlq(i)
               prnvar(loop,1)=hil(i)
               prnvar(loop,11)=D1(i)
               prnvar(loop,12)=hig(i)   
               prnvar(loop,13)=iae_conv(i)
               prnvar(loop,14)=iae_break(i)
               prnvar(loop,15)=iae_coal(i)
               prnvar(loop,16)=iae_nucl(i)
               prnvar(loop,17)=iae_size(i)    
            ENDDO
!
!...........Sort the exctracted DATA      
!
            CALL sortx_r(prnvar(1,1),index,nloop)
            DO k=2,isize1
               DO i=1,nloop
                  j=index(i)
                  prnvar0(i)=prnvar(j,1)
               ENDDO
               DO i=1,nloop
                  prnvar(i,k)=prnvar0(i)
               ENDDO
            ENDDO
!         
!...........PRINT the sorted DATA
!
            IF(INITIAL_p)THEN
               INITIAL_p=.FALSE.
               IF(vv_prob.eq.'suboex')THEN
                  OPEN(71,file='h1subo_tl.dat')    
                  OPEN(72,file='h1subo_ag.dat')    
                  OPEN(73,file='h1subo_ia.dat') 
                  OPEN(74,file='h1subo_hil.dat') 
                  OPEN(75,file='h1subo_vel.dat')  
                  OPEN(76,file='h1subo_D1.dat') 
                  OPEN(77,file='h1subo_hig.dat')        
                  OPEN(81,file='h2subo_tl.dat')    
                  OPEN(82,file='h2subo_ag.dat')    
                  OPEN(83,file='h2subo_ia.dat') 
                  OPEN(84,file='h2subo_hil.dat') 
                  OPEN(85,file='h2subo_vel.dat')  
                  OPEN(86,file='h2subo_D1.dat') 
                  OPEN(87,file='h2subo_hig.dat')   
                  OPEN(91,file='h3subo_tl.dat')    
                  OPEN(92,file='h3subo_ag.dat')    
                  OPEN(93,file='h3subo_ia.dat') 
                  OPEN(94,file='h3subo_hil.dat') 
                  OPEN(95,file='h3subo_vel.dat')  
                  OPEN(96,file='h3subo_D1.dat') 
                  OPEN(97,file='h3subo_hig.dat')   
                  OPEN(101,file='h4subo_tl.dat')    
                  OPEN(102,file='h4subo_ag.dat')    
                  OPEN(103,file='h4subo_ia.dat') 
                  OPEN(104,file='h4subo_hil.dat') 
                  OPEN(105,file='h4subo_vel.dat')  
                  OPEN(106,file='h4subo_D1.dat') 
                  OPEN(107,file='h4subo_hig.dat')   
                  OPEN(111,file='h5subo_tl.dat')    
                  OPEN(112,file='h5subo_ag.dat')    
                  OPEN(113,file='h5subo_ia.dat') 
                  OPEN(114,file='h5subo_hil.dat') 
                  OPEN(115,file='h5subo_vel.dat')  
                  OPEN(116,file='h5subo_D1.dat') 
                  OPEN(117,file='h5subo_hig.dat')   
                  OPEN(121,file='h6subo_tl.dat')    
                  OPEN(122,file='h6subo_ag.dat')    
                  OPEN(123,file='h6subo_ia.dat') 
                  OPEN(124,file='h6subo_hil.dat') 
                  OPEN(125,file='h6subo_vel.dat')  
                  OPEN(126,file='h6subo_D1.dat') 
                  OPEN(127,file='h6subo_hig.dat')               
                  OPEN(78,file='h1subo_IAT.dat') 
                  OPEN(88,file='h2subo_IAT.dat')
                  OPEN(98,file='h3subo_IAT.dat')
                  OPEN(108,file='h4subo_IAT.dat')
                  OPEN(118,file='h5subo_IAT.dat')
                  OPEN(128,file='h6subo_IAT.dat')  
               ELSE
                  OPEN(69,file='etc.dat') 
               ENDIF
! 
               ALLOCATE(varname(8,nloop))
               DO j=1,nloop
                  iradius=prnvar(j,1)*1000
                  WRITE(sradius,*)iradius
                  sradius=adjustl(sradius)
                  IF(vv_prob.eq.'suboex')THEN
                     varname(1,j)='tl'//trim(sradius) 
                     varname(2,j)='ag_'//trim(sradius)   
                     varname(3,j)='ia_'//trim(sradius)   
                     varname(4,j)='Hil_'//trim(sradius)  
                     varname(5,j)='Vl_'//trim(sradius) 
                     varname(6,j)='D1_'//trim(sradius)  
                     varname(7,j)='Hig_'//trim(sradius)
                     varname(8,j)='ia_conv'//trim(sradius)
                  ELSE                 
                     varname(1,j)='alphag_'//trim(sradius)
                  ENDIF
                  IF(vv_prob.ne.'suboex')varname(2,j)='vlz_'//trim(sradius)
               ENDDO                
!
               IF(vv_prob.eq.'suboex')THEN
                  DO k=1,6
                     WRITE(71+(k-1)*10,5000)(varname(1,j), j=1,nloop)
                     WRITE(72+(k-1)*10,5000)(varname(2,j), j=1,nloop) 
                     WRITE(73+(k-1)*10,5000)(varname(3,j), j=1,nloop) 
                     WRITE(74+(k-1)*10,5000)(varname(4,j), j=1,nloop)  
                     WRITE(75+(k-1)*10,5000)(varname(5,j), j=1,nloop)  
                     WRITE(76+(k-1)*10,5000)(varname(6,j), j=1,nloop)
                     WRITE(77+(k-1)*10,5000)(varname(7,j), j=1,nloop)
                     WRITE(78+(k-1)*10,5000)(varname(8,j), j=1,nloop)
                  ENDDO
               ENDIF 
            ENDIF 
            DEALLOCATE(varname)
!
!...........Print DATA        
!
!           IF(level.eq.Nheight)THEN      ! use when the data at specific height is required
            IF(vv_prob.eq.'suboex')THEN
               WRITE(71+(level-1)*10,5001)time,(prnvar(j, 9),j=1,nloop)
               WRITE(72+(level-1)*10,5001)time,(prnvar(j, 2),j=1,nloop)
               WRITE(73+(level-1)*10,5001)time,(prnvar(j, 3),j=1,nloop)
               WRITE(74+(level-1)*10,5001)time,(prnvar(j,10),j=1,nloop)
               WRITE(75+(level-1)*10,5001)time,(prnvar(j, 5),j=1,nloop)
               WRITE(76+(level-1)*10,5001)time,(prnvar(j,11),j=1,nloop)
               WRITE(77+(level-1)*10,5001)time,(prnvar(j,12),j=1,nloop)
               WRITE(78+(level-1)*10,5001)time,(prnvar(j,13),j=1,nloop),(prnvar(j,14),j=1,nloop),                         &
                                               (prnvar(j,15),j=1,nloop),(prnvar(j,16),j=1,nloop),(prnvar(j,17),j=1,nloop)
            ELSE
               WRITE(69,5001)time,(prnvar(j,2),j=1,nloop),(prnvar(j,3),j=1,nloop)
            ENDIF        
!           ENDIF
            DEALLOCATE(prnvar0,index)
            DEALLOCATE(prnvar)
         ENDDO
         DEALLOCATE(icell_index,rdistance_all)
!
      ENDIF
!      
5000  FORMAT(9x,'Time            ',1x,2(100A16))
5001  FORMAT(1x,1e16.8,1x,5(1000e16.8))      
!
      DEALLOCATE(pr,ag,diab,cboron)
      DEALLOCATE(hil,hig,D1,eviscosl,lviscosl)
      DEALLOCATE(iae,iae_conv,iae_break,iae_coal,iae_nucl,iae_size)
      DEALLOCATE(tk,td,tp,tlq)
      DEALLOCATE(vlp,vgp)
!      
      END SUBROUTINE local_out_subo
!      
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE local_out_becker
!
!     Save lateral output for Origin
!
      USE VOL_DATA , ONLY: cell
      USE Wall_DATA, ONLY: face                    
      USE Zconst1  , ONLY: vv_prob
      USE Zcoord1  , ONLY: xloc_tmp
      USE Zcoord2  , ONLY: xloc_xfc_min,xloc_xfc_max
      USE Zcore    , ONLY: myrank  
      USE Zface    , ONLY: qqcell,qecell,qclcell,qcgcell,ndensitycell           
      USE Zqvol    , ONLY: dry_weight       
      USE Zzone    , ONLY: ncell_fluid,ncell_fluid_all   
!
      IMPLICIT NONE 
!
!     size for prnvar
      INTEGER, PARAMETER :: isize1=22
!.....Local variables
      INTEGER :: i, j, k, na, wopt
      INTEGER :: loop,nloop
      LOGICAL,SAVE:: INITIAL_p=.TRUE.
      REAL(8) :: height, h1, h2
      REAL(8) :: hzloc,lzloc      
!
      CHARACTER*16 varname(10)
!            
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: icell_index,index
      REAL(8),DIMENSION(:),ALLOCATABLE :: height_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: ag,twall,ddepart,ndensity,dry_area
      REAL(8),DIMENSION(:),ALLOCATABLE :: qq,qe,qcl,qcg,tl
      REAL(8),DIMENSION(:),ALLOCATABLE :: hzloc1_all,lzloc1_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar0
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
!      
      na=ncell_fluid_all
      IF(myrank.eq.0) THEN
         ALLOCATE(ag(na),twall(na),ddepart(na),ndensity(na),dry_area(na))
         ALLOCATE(qq(na),qe(na),qcl(na),qcg(na),tl(na))
         ALLOCATE(hzloc1_all(ncell_fluid_all),lzloc1_all(ncell_fluid_all))
      ELSE
         ALLOCATE(ag(1),twall(1),ddepart(1),ndensity(1),dry_area(1))
         ALLOCATE(qq(1),qe(1),qcl(1),qcg(1),tl(1))
         ALLOCATE(hzloc1_all(1),lzloc1_all(1))
      ENDIF
!          
      CALL gatherv_r(cell%alphag         ,ncell_fluid,ag      ,na,0)
      CALL gatherv_r(face%twall_partition,ncell_fluid,twall   ,na,0)
      CALL gatherv_r(cell%ddepart        ,ncell_fluid,ddepart ,na,0)
      CALL gatherv_r(ndensitycell        ,ncell_fluid,ndensity,na,0)
      CALL gatherv_r(dry_weight          ,ncell_fluid,dry_area,na,0)
      CALL gatherv_r(qqcell              ,ncell_fluid,qq      ,na,0)
      CALL gatherv_r(qecell              ,ncell_fluid,qe      ,na,0)
      CALL gatherv_r(qclcell             ,ncell_fluid,qcl     ,na,0)
      CALL gatherv_r(qcgcell             ,ncell_fluid,qcg     ,na,0)
      CALL gatherv_r(cell%tl             ,ncell_fluid,tl      ,na,0)
!
!...........Check the x-location (for Becker)
!
      CALL gatherv_r(xloc_xfc_max(1,1),ncell_fluid,hzloc1_all,ncell_fluid_all,0)
      CALL gatherv_r(xloc_xfc_min(1,1),ncell_fluid,lzloc1_all,ncell_fluid_all,0)
!
      IF(myrank.eq.0) THEN
!      
         IF(vv_prob.eq.'becker')THEN
            h1=0.0d0
            h2=7.0d0
         ENDIF
!     
         ALLOCATE(icell_index(na),height_all(na))
         nloop=0   
         DO i=1,ncell_fluid_all
!      
!...........Check the x-location (for Becker)
!
            hzloc=hzloc1_all(i)
            lzloc=lzloc1_all(i)
!
!...........Check the height
!
            IF(hzloc.le.0.00475d0.or.lzloc.gt.0.00475d0) CYCLE   !0.00475d0=x location of near wall cell
!
!...........Check the height
!
            wopt = 0
            height=dabs(xloc_tmp(i,3)) 
            IF(height.ge.h1.and.height.le.h2) wopt=1
!
!...........Extract the DATA to be printed if wopt ==1
!
            IF(wopt.ne.1) CYCLE                                        !when height and radius is in proper value
            nloop=nloop+1
            icell_index(nloop)=i
            height_all(nloop)=height
         ENDDO
         ALLOCATE(prnvar(nloop,isize1),prnvar0(nloop))
         DO loop=1,nloop
            i=icell_index(loop)
            height=height_all(loop)
            prnvar(loop,1)=height
            prnvar(loop,2)=ag(i)
!            prnvar(loop,3)=iae(i)
!            prnvar(loop,4)=diab(i)*1000.0d0
!            prnvar(loop,5)=vlp(ndim,i)  
!            prnvar(loop,6)=vgp(ndim,i)
!            prnvar(loop,7)=pr(i)
!            prnvar(loop,8)=cboron(i)  
!            prnvar(loop,9)=eviscosl(i)/lviscosl(i)
!            prnvar(loop,6)=tk(i)
!            prnvar(loop,7)=td(i)
!            prnvar(loop,8)=tp(i)
!            prnvar(loop,9)=tlq(i)
!            prnvar(loop,10)=hil(i)
!            prnvar(loop,11)=D1(i)
!            prnvar(loop,12)=hig(i)   
!            prnvar(loop,13)=iae_conv(i)
            prnvar(loop,14)=qcl(i)
            prnvar(loop,15)=qq(i)
            prnvar(loop,16)=qe(i)
            prnvar(loop,17)=qcg(i)
            prnvar(loop,18)=dry_area(i)
            prnvar(loop,19)=ndensity(i)
            prnvar(loop,20)=ddepart(i)
            prnvar(loop,21)=twall(i)
            prnvar(loop,22)=tl(i)
         ENDDO
!
!........Sort the exctracted DATA      
!
            CALL sortx_r(prnvar(1,1),index,loop)
            DO k=2,isize1
               DO i=1,nloop
                  j=index(i)
                  prnvar0(i)=prnvar(j,k)
               ENDDO
               DO i=1,nloop
                  prnvar(i,k)=prnvar0(i)
               ENDDO
            ENDDO
!         
!........PRINT the sorted DATA
!
         IF(INITIAL_p)THEN
            INITIAL_p=.FALSE.
            IF(vv_prob.eq.'becker')THEN
               OPEN(71,file='becker_ss.dat')    
            ELSE
               OPEN(69,file='etc.dat') 
            ENDIF
! 
            IF(vv_prob.eq.'becker')THEN
                 varname( 1)='ag'
                 varname( 2)='qcl'
                 varname( 3)='qq'  
                 varname( 4)='qe' 
                 varname( 5)='qcg'
                 varname( 6)='dry_area'
                 varname( 7)='ndensity'
                 varname( 8)='ddepart'  
                 varname( 9)='twall'  
                 varname(10)='tl'                              
                 WRITE(71,5000)(varname(j), j=1, 10)
            ENDIF 
         ENDIF 
!
!........Print DATA        
!
!         IF(level.eq.Nheight)THEN      ! use when the data at specific height is required
            IF(vv_prob.eq.'becker')THEN
               DO i=1,nloop               
                  WRITE(71,5001)prnvar(i,1),prnvar(i,2),(prnvar(i,j),j=14,22)
               ENDDO
            ENDIF        
!         ENDIF
         DEALLOCATE(prnvar0,index)
         DEALLOCATE(prnvar)
         DEALLOCATE(icell_index,height_all)
      ENDIF
      DEALLOCATE(ag,twall,ddepart,ndensity,dry_area)
      DEALLOCATE(qq,qe,qcl,qcg,tl)
      DEALLOCATE(hzloc1_all,lzloc1_all)
!      
5000  FORMAT(9x,'Height           ',1x,2(100A16))
5001  FORMAT(1x,1e16.8,1x,5(1000e16.8))      
!
!      
      END SUBROUTINE local_out_becker            
      

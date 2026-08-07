!------------------------------------------------------------------------------------------------------------------            
      SUBROUTINE siphon_print_ref     
      USE VOL_DATA        , ONLY: cell
      USE Zsiphon         , ONLY: topen
      USE Ztimecon        , ONLY: time
      USE Zzone           , ONLY: ncell_fluid
      USE Zcoord1         , ONLY: xloc
      USE Zcoord3         , ONLY: vol
      USE Zcore           , ONLY: myrank
      
      IMPLICIT NONE  
      INTEGER :: i 
      LOGICAL,SAVE :: initial=.true.
      REAL(8),SAVE :: xx1,xx2,yy1,yy2,tank_height !siphon
      REAL(8) :: av_sum,v_sum,water_level !siphon      
      REAL(8),SAVE :: print_time

          IF (initial) THEN
               initial=.false.
               IF(myrank.eq.0)OPEN(333, file='VD14_siphon_ref.dat')
               print_time=time-0.001d0
               tank_height=4.5d0
               xx1=-0.44d0
               xx2=xx1+0.1d0
               yy1=-1.77d0
               yy2=yy1+0.1d0
            ENDIF
            IF(time.ge.print_time)THEN
                print_time=print_time+0.01d0
                av_sum=0.0d0
                v_sum=0.0d0
                DO i=1,ncell_fluid
                   IF(xloc(i,1).ge.xx1.and.xloc(i,1).le.xx2.and.xloc(i,2).ge.yy1.and.xloc(i,2).le.yy2)THEN
                      av_sum=av_sum+cell%alphal(i)*vol(i)
                      v_sum=v_sum+vol(i)  
                   ENDIF
                ENDDO                
                CALL allreducei_r1(av_sum)
                CALL allreducei_r1(v_sum)
                IF(myrank.eq.0)THEN
                   water_level=av_sum/v_sum*tank_height
                   WRITE(333,"(2(1x,1pe12.5))")time-topen,water_level
                ENDIF
          ENDIF 
      RETURN           
      ENDSUBROUTINE siphon_print_ref      
!----------------------------------------------------------------------------------
      SUBROUTINE siphon_initialization     
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid,num_max_zone
      USE Zconst1      , ONLY: vv_prob
      USE Zcoord1      , ONLY: xloc
      USE Zcore        , ONLY: myrank
      USE Zio_unit   , ONLY: unit_log
            
      IMPLICIT NONE
      INTEGER i
      REAL(8) :: rmax0,xx,yy,rr,xx1,xx2,yy1,yy2,zz1,zz2,zsurf,xcen,ycen    
      
!
      IF(vv_prob.eq.'siphon')THEN    
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'Initilization for siphon in initialize_specific_variables!'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Initilization for siphon in initialize_specific_variables!'
         xx1=-0.445d0
         xx2=3.555d0
         yy1=-1.775d0
         yy2=1.775d0
         zz1=0.0d0
         zz2=4.5d0
         zsurf=3.855d0
         xcen=2.155d0
         ycen=0.0d0
         rmax0=69.0d-3/2.0d0 !2.5*2.54e-2/2.0d0
         xcen=xcen+0.8d0
         DO i=1,ncell_fluid
             xx=xloc(i,1)-xcen
             yy=xloc(i,2)-ycen
             rr=(xx*xx+yy*yy)**0.5d0
             IF(rr.le.rmax0)THEN
             ELSE
                IF((xloc(i,1).ge.xx1.and.xloc(i,1).le.xx2).and.(xloc(i,2).ge.yy1.and.xloc(i,2).le.yy2))THEN
                   IF(xloc(i,3).ge.zsurf.and.xloc(i,3).le.zz2)THEN
                      cell%alphal(i)=0.001d0
                      cell%alphag(i)=1.0-cell%alphal(i)
                      cell%quala(i)=1.0d0  
                   ENDIF
                ELSE
                ENDIF   
             ENDIF   
         ENDDO 
      ENDIF
!      
      ENDSUBROUTINE siphon_initialization     
!------------------------------------------------------------------------------------------------------------------
      SUBROUTINE siphon_control_pbnd
!
!     Calculate the pressure at outlets. Only for "i_horizontal_outlet" is used.
!
      USE VOL_DATA     , ONLY: cell          
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np,myrank
      USE Zparam       , ONLY: ndim,nin_max,nb_max
      USE Znum_cell    , ONLY: i_neigh,neigh
      USE Ztimecon     , ONLY: time,itim
      USE Zbc_index    , ONLY: nbcon,npb
      USE Zconst1      , ONLY: vv_prob
      USE Zconst2      , ONLY: grav
      USE Zcoord1      , ONLY: xloc
      USE Zpress       , ONLY: p
      USE Zb_condition , ONLY: pbnd
      USE Zsiphon      , ONLY: topen 
      USE Zbc_index    , ONLY: l_horizontal_outlet_init
      USE Zvec_geo     , ONLY: sv_nf,fac_non,fac1_non
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,k,i0,j0,k1
      INTEGER :: i1,ii,kk
      INTEGER :: kill
      INTEGER :: loc_top(1)
      INTEGER :: index1,index2,index3,nout,ipb    
      INTEGER :: out_opt_max 
      INTEGER,SAVE :: valveopt=1
      INTEGER,SAVE :: npb_max
      LOGICAL,SAVE :: INITIAL=.TRUE.
      LOGICAL,SAVE :: head_LOGICAL
      REAL(8) :: a,ffac,ffac1,svx,svy,svz
      REAL(8) :: head,projection
      REAL(8) :: dtopen   
!.....Local arrays
      INTEGER :: out_opt(nb_max-nin_max)
      INTEGER,SAVE,ALLOCATABLE :: outneigh(:),outneigh_tot(:,:)     
      INTEGER,SAVE,ALLOCATABLE :: outcell(:,:),outcell_tot(:,:,:)   
      REAL(8),SAVE,ALLOCATABLE::elev_out(:),elevout_tot(:,:)
      INTEGER,SAVE,ALLOCATABLE::nout_tot(:),pbnd_last(:),pbnd_first(:) 
      INTEGER,ALLOCATABLE::indx(:),itemp(:),itemp1(:),itemp2(:)
!         
!.....Horizontal_flow
!
      IF(INITIAL.or.l_horizontal_outlet_init)THEN
!      
         npb_max=nb_max-nin_max      
!
         IF(.not.ALLOCATED(nout_tot))   ALLOCATE(nout_tot(npb_max))
         IF(.not.ALLOCATED(pbnd_last))  ALLOCATE(pbnd_last(nb_max))
         IF(.not.ALLOCATED(pbnd_first)) ALLOCATE(pbnd_first(nb_max))
!     
         IF(myrank.eq.0)then
            WRITE(*,*)'          Gravity effect on P.B. will be controlled!'
            WRITE(*,*)'          Cells included in a P.B. should be in a subdomain!' 
            WRITE(unit_log,*)'          Gravity effect on P.B. will be controlled!'
            WRITE(unit_log,*)'          Cells included in a P.B. should be in a subdomain!'                
         ENDIF
!            
         nout_tot(:)=0
         DO i=1,ncell_fluid
            k=npb(i)
            IF(k.gt.0.and.k.le.npb_max)THEN
               nout_tot(k)=nout_tot(k)+1
            ENDIF
         ENDDO
         DO i=1,npb_max
            IF(nout_tot(i).gt.0) WRITE(*,"(11x,a,1i3,1x,1i8,1x,1i3)")'npb,#cells,rank=',i,nout_tot(i),myrank
         ENDDO           
         nout=MAXVAL(nout_tot(:))
!
!........Pressure boundary must be in a subdomain.
!         
         out_opt(:)=0
         out_opt_max=0
         DO ipb=1,npb_max
            IF(nout_tot(ipb).ge.1.and.npb(ipb).le.npb_max) out_opt(ipb)=1
         ENDDO
         IF(np.gt.1) CALL allreducei_i(out_opt,npb_max)
         kill=0
         DO ipb=1,npb_max
             out_opt_max=MAX(out_opt(ipb),out_opt_max)
             IF(out_opt(ipb).ge.2)THEN
                kill=1
                IF(myrank.eq.0)THEN
                   WRITE(*,"(11x,a,1i2,1x,a)")  'P.B. #',ipb+4,'is splitted into more than two subdomains!'
                   WRITE(unit_log,"(11x,a,1i2,1x,a)") 'P.B. #',ipb+4,'is splitted into more than two subdomains!'
                ENDIF   
             ENDIF 
         ENDDO
         IF(vv_prob.ne.'Horizontal_flow') THEN
            IF(np.gt.1) CALL allreducei_i1(kill)
            IF(kill.gt.0) THEN
               CALL barrier_mpi
               CALL finalize_mpi
               STOP
            ENDIF
         ENDIF
!         
         IF(ALLOCATED(outcell))      DEALLOCATE(outcell)
         IF(ALLOCATED(elev_out))     DEALLOCATE(elev_out)
         IF(ALLOCATED(outcell_tot))  DEALLOCATE(outcell_tot)
         IF(ALLOCATED(elevout_tot))  DEALLOCATE(elevout_tot)
         IF(ALLOCATED(outneigh))     DEALLOCATE(outneigh)
         IF(ALLOCATED(outneigh_tot)) DEALLOCATE(outneigh_tot)
!         
         ALLOCATE(outcell(nout,10),elev_out(nout))
         ALLOCATE(outcell_tot(nout,10,npb_max),elevout_tot(nout,npb_max))
         ALLOCATE(outneigh(nout),outneigh_tot(nout,npb_max))               !siphon
!
         ALLOCATE(indx(nout),itemp(nout),itemp1(nout),itemp2(nout))
!         
         DO ipb=1,npb_max         
             outcell(:,:)=0
             elev_out(:)=0.d0
             IF(nout_tot(ipb).le.0)CYCLE
             nout=0
             DO i=1, ncell_fluid
                IF(npb(i).eq.ipb)THEN
                   nout=nout+1
                   outcell(nout,1)=i
                   elev_out(nout)=DOT_PRODUCT(xloc(i,:),grav(:))
                   j0=i_neigh(i)-1
                   DO j=i_neigh(i),i_neigh(i+1)-1
                      IF(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max) outneigh(nout)=j-j0 !siphon
                   ENDDO                   
                ENDIF
             ENDDO
             IF(nout.le.0)CYCLE
!
!............Sort based on elev_out =>outcell(:,1),outneigh(:)
!
             CALL sortx_r(elev_out,indx,nout)
             do i=1,nout
                itemp(i)=outcell(i,1)
             enddo
             do i=1,nout
                k=indx(i)
                outcell(i,1)=itemp(k)
             enddo
!
             do i=1,nout
                itemp(i)=outneigh(i)
             enddo
             do i=1,nout
                k=indx(i)
                outneigh(i)=itemp(k)
             enddo
!
             outcell(:,2)=0
             outcell(:,3:10)=0
!
!             DO i=1,nout
!                WRITE(*,"(11x,a,1i3,1x,1e12.4,1x,1i8,1x,1i3,1x,1i3)")'index,elevation,cell,rank=',i,elev_out(i),outcell(i,1),npb(outcell(i,1)),myrank
!             ENDDO
             DO i0=1,nout
                i=outcell(i0,1)
                DO j=i_neigh(i),i_neigh(i+1)-1
                   k=neigh(j)                                      ! neighbor of outcell(i0,1)
! bug neigh array has 0 entry when nbcon # 0
                   IF(k.eq.0) CYCLE
                   IF(npb(k).ne.ipb)CYCLE
                   DO ii=1,nout
                      IF(k.eq.outcell(ii,1))THEN                     ! neighbor of outcell(i0,1)
                         IF(elev_out(i0).ge.elev_out(ii)+1.0e-8)THEN 
                            outcell(i0,2)=outcell(i0,2)+1
                            outcell(i0,outcell(i0,2)+2)=k
                            IF(outcell(i0,2).gt.8)then
                               WRITE(*,*) 'Error1 in udfn_horizontal_outlet!!!'
                               WRITE(unit_log,*)'Error1 in udfn_horizontal_outlet!!!'
                               STOP
                            ENDIF   
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
!             outlet_top=MINVAL(elev_out(:))
!             loc_top=MINLOC(elev_out(:))
             outcell_tot(1:nout,1:10,ipb)=outcell(1:nout,1:10)
             elevout_tot(1:nout,ipb)=elev_out(1:nout)
             outneigh_tot(1:nout,ipb)=outneigh(1:nout) !siphon
          ENDDO
          DEALLOCATE(indx,itemp,itemp1,itemp2)
          INITIAL=.FALSE.
          pbnd_last(:)=pbnd(:) !siphon
          pbnd_first(:)=pbnd(:) !siphon
      ENDIF 
!------------------------------------------------------------------------------------------------------------------            
!
!.....Control pressure boundary like a valve 
!  
      IF(vv_prob.eq.'siphon')THEN 
         dtopen=0.1d0
         nout_tot(2)=0 !vertical pressure boundary(nbcon=6) is not controlled.          
         IF(valveopt.eq.1)THEN !opened
            IF(time.le.topen)THEN !closing
               valveopt=0
               IF(myrank.eq.0)THEN
                  WRITE(*,"(11x,a,1pe11.3,a,1pe11.3)")'Valve is closed at ',time,',& will be opened at ',topen                   
                  WRITE(unit_log,"(11x,a,1pe11.3,a,1pe11.3)")'Valve is closed at ',time,',& will be opened at ',topen                   
               ENDIF 
               CALL nbcon_change_start    
               DO ipb=1,1 !npb_max
                   nout=nout_tot(ipb)
                   IF(nout.le.1)CYCLE            
                   DO ii=1,nout
                      i=outcell_tot(ii,1,ipb)
                      j0=i_neigh(i)-1
                      j=outneigh_tot(ii,ipb)
                      nbcon(j+j0)=-1
!                     npb(i)=0
                   ENDDO
               ENDDO              
               CALL nbcon_change_end    
               RETURN !closed
            ELSE
            !opened
            ENDIF
         ELSE !closed
            IF(time.gt.topen)THEN !opening
               valveopt=1
               IF(myrank.eq.0)WRITE(*,"(11x,a,1pe12.5)")'Valve is opened at time=',time                  
               IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1pe12.5)")'Valve is opened at time=',time                  
               CALL nbcon_change_start    
               DO ipb=1,1 !npb_max
                   nout=nout_tot(ipb)
                   pbnd_last(ipb)=0.0d0
                   IF(nout.le.0)CYCLE            
                   DO ii=1,nout
                      i=outcell_tot(ii,1,ipb)
                      j=outneigh_tot(ii,ipb)
                      j0=i_neigh(i)-1
                      nbcon(j+j0)=ipb+nin_max
!                     npb(i)=ipb
                      pbnd_last(ipb)=pbnd_last(ipb)+p(i)
                      WRITE(*,"(11x,a,1i2,1x,1i8,1x,1i3)")'myrank,i,j=',myrank,i,j
                   ENDDO
                   pbnd_last(ipb)=pbnd_last(ipb)/nout
               ENDDO    
               CALL nbcon_change_end    
            ELSE
               RETURN !closed
            ENDIF   
         ENDIF
!          
!.........set the pressure of the pressure boundary
          IF(time.gt.topen.and.time.lt.(topen+dtopen))THEN
             pbnd(:)=(pbnd_first(:)-pbnd_last(:))*(time-topen)/dtopen+pbnd_last(:)
             IF(mod(itim,10).eq.0)THEN
               DO ipb=1,npb_max
                   IF(nout_tot(ipb).gt.0)WRITE(*,"(11x,a,1x,1i2,1x,2(1pe10.3))")'myrank&pbnd(1,2)=',myrank,pbnd(1)/1.d6,pbnd(2)/1.d6      
               ENDDO    
             ENDIF  
          ELSE
             pbnd(:)=pbnd_first(:)
          ENDIF   
      ENDIF 
!------------------------------------------------------------------------------------------------------------------            
!
!.....Set a pressure at pressure boundary cells
!
          DO i=1,ncell_fluid
            IF(npb(i).gt.0) p(i)=pbnd(npb(i))
          ENDDO    
!
!....every time step
!
      DO ipb=1,npb_max
!         
         nout=nout_tot(ipb)
         IF(nout.le.1)CYCLE
         outcell(1:nout,1:10)=outcell_tot(1:nout,1:10,ipb)
         elev_out(1:nout)=elevout_tot(1:nout,ipb)
         DO ii=1,nout
            head=0.0d0
            projection=0.0d0
            i=outcell(ii,1)
            index1=outcell(ii,2)
!
!...........no neighbor in different level with i-cell 
!
            head_LOGICAL=.FALSE.
            IF(index1.eq.0)THEN
               IF(elev_out(ii).gt.MINVAL(elev_out(1:ii))) THEN
                  loc_top = MINLOC(elev_out(1:ii))
                  head=DOT_PRODUCT((xloc(outcell(loc_top(1),1),:)-xloc(i,:)),grav(:))*cell%rhom(outcell(loc_top(1),1))
                  p(i)=p(outcell(loc_top(1),1))-head
               ENDIF
            ENDIF
!
!...........has neighbor in different level with i-cell
!
            IF(index1.eq.0) CYCLE
            j0=i_neigh(i)-1 
            DO j=i_neigh(i),i_neigh(i+1)-1
               k1=0
               k=neigh(j)
               IF(k.eq.0) CYCLE
               IF(npb(k).ne.ipb) CYCLE
               DO kk=1,index1
                  index2=kk+2
                  index3=outcell(ii,index2)
                  IF(k.eq.index3)THEN
                     k1=index3
                     EXIT 
                  ENDIF
               ENDDO
               IF(k1.eq.0) CYCLE
!
!..............Get offset i1 in vector space of (j-j0,i)
!
               CALL get_vector_disp(j-j0,i,i1)
!
               IF(i1.gt.0) THEN
                  ffac=fac_non(i1)
                  ffac1=fac1_non(i1)
                  svx=sv_nf(i1,1)
                  svy=sv_nf(i1,2)
                  IF(ndim.eq.3) svz=sv_nf(i1,3)
               ELSE
                  i1=-i1
                  ffac=fac1_non(i1)
                  ffac1=fac_non(i1)
                  svx=-sv_nf(i1,1)
                  svy=-sv_nf(i1,2)
                  IF(ndim.eq.3) svz=-sv_nf(i1,3)
               ENDIF
               IF(index1.ge.2)THEN
                  a=svx*grav(1)+svy*grav(2)
                  IF(ndim.eq.3) THEN
                    a=a+svz*grav(3)
                  ENDIF 
                  head=head+(-(cell%rhom(i)*ffac1+cell%rhom(k1)*ffac)*DOT_PRODUCT(grav(:),xloc(k1,:)-xloc(i,:))+p(k1)) &
                             *a/DSQRT(DOT_PRODUCT(grav(:),grav(:)))
                   projection=projection+a/DSQRT(DOT_PRODUCT(grav(:),grav(:)))               
               ELSE
                  head=head+(-(cell%rhom(i)*ffac1+cell%rhom(k1)*ffac)*DOT_PRODUCT(grav(:),xloc(k1,:)-xloc(i,:))+p(k1)) 
                  projection=1.0d0                
               ENDIF                        
            ENDDO  
            head=head/projection
            p(i)=head
        ENDDO
      ENDDO
!     
      RETURN 
      END SUBROUTINE siphon_control_pbnd
!-------------------------------------------------------------------------------

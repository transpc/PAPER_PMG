!
      SUBROUTINE rocom_control_bc_user
!
!     Define a boundary conditoin for the ROCOM problem
!
      USE Zzone           , ONLY: ncell_fluid
      USE Zparam          , ONLY: nb_max
      USE Zcore           , ONLY: myrank      
      USE Znum_cell       , ONLY: i_neigh
      USE Zbc_index       , ONLY: nbcon,npb,ngrad,icell_type,iface_wall
      USE Ztimecon        , ONLY: time,t_end_ctrl 
      USE Zrocom_specific , ONLY: closeloop,injection
      USE Zio_unit        , ONLY: unit_log
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER i,j,j0
      INTEGER, SAVE :: idex
      INTEGER, SAVE :: changeopt=0
      LOGICAL, SAVE :: INITIAL=.TRUE.
!.....Local arrays
      INTEGER, SAVE :: istore(1000),jstore(1000),nbstore(1000),npbstore(1000),ngradstore(1000)
      INTEGER, SAVE :: icell_typestore(1000),iface_wallstore(1000)
!
      IF(INITIAL)THEN     
         WRITE(*,*)'          Control_bc_rocom are used now!' 
         IF(myrank.eq.0)WRITE(unit_log,*)'          Control_bc_rocom are used now!'               
         closeloop=0
         injection=0
      ENDIF        
!
!.....closeloop(9),injection=0-->injection start
!.....Injection start: injection=1(2=sustain),changeopt=1(0)
!.....Injection stop: injection=3(4-sustain),changeopt=1(0)  
!
      IF(closeloop.eq.0 .and. Time.ge.0.0d0) closeloop=1         
      IF(injection.eq.0 .and. Time.ge.t_end_ctrl(1)) injection=1 
      IF(injection.eq.2 .and. Time.ge.t_end_ctrl(2)) injection=3
!
!.....Just once, store nbcon and its index for debugging
!
      IF(INITIAL)THEN
         idex=0
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).ge.1 .and. nbcon(j).le.nb_max)THEN
                  idex=idex+1
                  istore(idex)=i
                  jstore(idex)=j-j0
                  nbstore(idex)=nbcon(j)
                  npbstore(idex)=npb(i)
                  ngradstore(idex)=ngrad(i)
                  icell_typestore(idex)=icell_type(i)
                  iface_wallstore(idex)=iface_wall(i)
               ENDIF
               IF(idex.gt.1000)THEN
                  WRITE(*,*)'Too many boundary cell in SUBROUTINE control_bc_rocom!'
                  IF(myrank.eq.0)WRITE(unit_log,*)'Too many boundary cell in SUBROUTINE control_bc_rocom!'
                  STOP
               ENDIF       
            ENDDO
         ENDDO
         idex=0
      ENDIF      
!
!.....Just once, close the loop by temperary p.b.c.( 9 --> 1)
!
      IF(closeloop.eq.1)THEN
         CALL nbcon_change_start
         closeloop=2                                                                !no more closeloop change
         changeopt=1
         WRITE(*,*)'          nbcon 9 => -1 at time=',time
         IF(myrank.eq.0)WRITE(unit_log,*)'          nbcon 9 => -1 at time=',time
         DO i=1,ncell_fluid
!           j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).eq.9)THEN
                  nbcon(j)=-1
!                 npb(i)=0
!                 ngrad(i)=1 
!                 icell_type(i)=1
!                 iface_wall(i)=j-j0
                ENDIF
            ENDDO
         ENDDO
         CALL nbcon_change_end
      ENDIF         
!
!.....Just once, stop injection by changing and storing nbcon of 5,1 into -1,-1
!
      IF(INITIAL)THEN
         CALL nbcon_change_start
         WRITE(*,*)'          Loop-I are closed: 5,1 => -1at time=',time
         IF(myrank.eq.0)WRITE(unit_log,*)'          Loop-I are closed: 5,1 => -1at time=',time
         INITIAL=.FALSE.                                                           !INITIAL is true for just one time.
         idex=0
         changeopt=1
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).eq.5)THEN
                  idex=idex+1
                  istore(idex)=i
                  jstore(idex)=j-j0
                  nbstore(idex)=nbcon(j)
                  npbstore(idex)=npb(i)
                  ngradstore(idex)=ngrad(i)
                  icell_typestore(idex)=icell_type(i)
                  iface_wallstore(idex)=iface_wall(i)                   
                  nbcon(j)=-1
!                 npb(i)=0
!                 ngrad(i)=1
!                 icell_type(i)=1
!                 IF(iface_wall(i).le.0) iface_wall(i)=j-j0
               ELSEIF(nbcon(j).eq.1)THEN
                  idex=idex+1
                  istore(idex)=i
                  jstore(idex)=j-j0
                  nbstore(idex)=nbcon(j)                 
                  npbstore(idex)=npb(i)
                  ngradstore(idex)=ngrad(i)
                  icell_typestore(idex)=icell_type(i)
                  iface_wallstore(idex)=iface_wall(i)  
                  nbcon(j)=-1
!                 npb(i)=0
!                 ngrad(i)=1
!                 icell_type(i)=1
!                 IF(iface_wall(i).le.0) iface_wall(i)=j-j0
               ENDIF
               IF(idex.gt.1000)THEN
                  WRITE(*,*)'Too many boundary cell in SUBROUTINE control_bc_rocom!'
                  IF(myrank.eq.0)WRITE(unit_log,*)'Too many boundary cell in SUBROUTINE control_bc_rocom!'
                  STOP
               ENDIF       
            ENDDO
         ENDDO
         CALL nbcon_change_end
      ENDIF 
!      
!.....Start injection by restoring nbcon which is changed into -1, -1 from 5, 1        
!
      IF(injection.eq.1)THEN
         CALL nbcon_change_start
         WRITE(*,*)'          Loop-I are opened: -1 => 5,1 at time=',time
         IF(myrank.eq.0)WRITE(unit_log,*)'          Loop-I are opened: -1 => 5,1 at time=',time
         injection=2                                                             !no more injection change       
         changeopt=1
         DO i=1,idex
            j0=i_neigh(istore(i))-1
            nbcon(jstore(i)+j0)=nbstore(i)
            npb(istore(i))=npbstore(i)
            ngrad(istore(i))=ngradstore(i)
            icell_type(istore(i))=icell_typestore(i)  
            iface_wall(istore(i))=iface_wallstore(i)
         ENDDO
         CALL nbcon_change_end
       ELSEIF(injection.eq.3)THEN
         CALL nbcon_change_start
         WRITE(*,*)'          Loop-I are closed:5,1 => -1 at time=',time 
         IF(myrank.eq.0)WRITE(unit_log,*)'          Loop-I are closed:5,1 => -1 at time=',time      
         injection=4                                                             !no more injection change       
         changeopt=1
         DO i=1,idex
            j0=i_neigh(istore(i))-1
            nbcon(jstore(i)+j0)=-1
            npb(istore(i))=0
            ngrad(istore(i))=1
            icell_type(istore(i))=1
            IF(iface_wall(istore(i)).le.0) iface_wall(istore(i))=iface_wallstore(i)
         ENDDO      
         CALL nbcon_change_end
      ENDIF      
! 
      END SUBROUTINE rocom_control_bc_user

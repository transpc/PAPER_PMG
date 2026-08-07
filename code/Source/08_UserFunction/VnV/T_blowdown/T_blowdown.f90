!
      SUBROUTINE T_blowdown_user
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: np,myrank
      USE Ztimecon     , ONLY: time
      USE Zcoord3      , ONLY: volp
      USE Zmass_conv   , ONLY: tot_mass
      USE Zpress       , ONLY: p
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,na,nout
      INTEGER,SAVE :: iout=0
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: tot_eng
!.....Local arrays
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: dat
!
      IF(initial.and.myrank.eq.0)THEN
         OPEN(990,file='VFT11_T_blowdown_ref.dat')
         initial=.false.
      ENDIF
!
      tot_eng=0.d0
      DO i=1,ncell_fluid
         tot_eng=tot_eng+( cell%alphag(i)*cell%rhog(i)*cell%eg(i)          &
                          +cell%alphal(i)*cell%rhol(i)*cell%el(i))*volp(i)
      ENDDO
!
      IF(time.le.0.01)THEN
         nout=1
      ELSE
         nout=10
      ENDIF
!
      iout=iout+1
      IF(iout.ge.nout)THEN
         na=ncell_fluid_all
         IF(myrank.eq.0) THEN
            ALLOCATE(dat(na,2))
         ELSE
            ALLOCATE(dat(1,2))
         ENDIF
         CALL gatherv_r(cell%alphag,ncell_fluid,dat(1,1),na,0)
         CALL gatherv_r(p          ,ncell_fluid,dat(1,2),na,0)
         IF(np.gt.1) CALL allreducei_r1(tot_eng)
         IF(myrank.eq.0) THEN
            WRITE(990,10) time,tot_mass,tot_eng,dat(25,1),dat(50,1),dat(610,1),dat(25,2),dat(50,2),dat(610,2)
         ENDIF
         iout=0
         DEALLOCATE(dat)
      ENDIF
   10 FORMAT(10(4e15.7,1x))
!
      END SUBROUTINE T_blowdown_user

!
      SUBROUTINE udfn_mat_prop
!      
!     User-defined material properties (only when 'udfl_mat_prop' is used.)
!     This is different from Initial_Variables_solid_user.f90. This is a kind of bc for solids after calculating T_solid
!  
      USE SOLID_DATA , ONLY: solid
      USE Zconst1    , ONLY: vv_prob      
      USE Zmpi       , ONLY: ncell_ps
      USE Zzone      , ONLY: nmaterial_c,icore,ncell_cond
      USE Zcore      , ONLY: myrank
      USE Znum_cell  , ONLY: n_fluid
      USE Zio_unit   , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER i
      INTEGER iOKr,iOKk
      LOGICAL,SAVE :: initial
      REAL(8) CpVol,Condu 
      REAL(8),ALLOCATABLE,SAVE :: multiplier(:)
      DATA initial/.true./
!      
      IF(vv_prob.eq.'PAFS-POOL') then
         DO i=1,ncell_ps !ncell_cond
            CALL mat_prop(ABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iOKr,iOKk)
            solid%rhocps(i)=CpVol/100.0d0
            solid%conds(i)=Condu*1.0d0
            solid%matnum(i)=ABS(nmaterial_c(i))
         ENDDO
      ELSEIF(vv_prob.eq.'2D_conduction'.or.vv_prob.eq.'Horizontal_flow')THEN
         DO i=1,ncell_ps !ncell_cond
            IF(ABS(nmaterial_c(i)).lt.50) THEN
               solid%rhocps(i)=4100.0d0
               solid%conds(i)=100.d0
            ENDIF
         ENDDO
      ELSEIF(vv_prob.eq.'apr1400_lbloca')THEN
         IF(initial)THEN  
            initial=.false.
            ALLOCATE(multiplier(ncell_ps))
            multiplier(:)=1.0d0
            DO i=1,ncell_cond 
               IF(icore(n_fluid(i)).lt.1)multiplier(i)=1.0d0 !sensitivity
            ENDDO
            IF(myrank.eq.0)WRITE(*,"(8x,a,1f5.1,a)")'*** In udfn_mat_prop, solid%rhocps(i) = CpVol * ', MAXVAL(multiplier(:)), '.'
            IF(myrank.eq.0)WRITE(unit_log,"(8x,a,1f5.1,a)")'*** In udfn_mat_prop, solid%rhocps(i) = CpVol * ', MAXVAL(multiplier(:)), '.'
      ENDIF   
         DO i=1,ncell_cond
            CALL mat_prop(ABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iOKr,iOKk)
            solid%rhocps(i)=CpVol*multiplier(i)
            solid%conds(i)=Condu
         ENDDO   
      ELSEIF(vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')then
         DO i=1,ncell_cond !ncell_ps-error in nmaterial_c
             IF(abs(nmaterial_c(i)).le.5.and.abs(nmaterial_c(i).gt.0)) THEN
                CALL mat_prop(abs(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iokr,iokk)
                solid%rhocps(i)=CpVol 
                solid%conds(i)=Condu 
             ENDIF
             solid%matnum(i)=abs(nmaterial_c(i))
             solid%rhocps(i)=0.4227193d7 !*5.0d0 !atlas_mc_porous sensitivity for 16 ton of ATLAS RPV except for heater
             solid%conds(i)=18.0d0 !
         ENDDO         
      ELSE
         IF(myrank.eq.0)WRITE(*,*)'udfl_mat_prop is true, but vv_prob is not set in udfn_mat_prop!!!'
         IF(myrank.eq.0)WRITE(unit_log,*)'udfl_mat_prop is true, but vv_prob is not set in udfn_mat_prop!!!'
         STOP         
      ENDIF   
!           
      END SUBROUTINE udfn_mat_prop

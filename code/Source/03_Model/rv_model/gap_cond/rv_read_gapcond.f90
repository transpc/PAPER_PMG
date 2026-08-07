!
      SUBROUTINE rv_read_gapcond(initial)
!
      USE Zrv_gap_cond      , ONLY:nr_gapi,nr_gapo,iplas,gap_P0,GapDisp_fission,GapRough,ngas, &
                                    indgas,molgas,rhocp_gap
      use Zio_unit          , only: unit_gap_cond
      use unitManager       , only: createUnit
!
      IMPLICIT NONE
!
      INTEGER i
      INTEGER err                                    
      REAL(8) GapRough_fuel,GapRough_clad
      LOGICAL initial
!
      IF(.not.initial) RETURN          
      initial=.false.
!
      ! OPEN(5,file='rv_gapcond.in',status='old',iostat=err)
      unit_gap_cond=createUnit("rv_gap_conductance")
      unit_gap_cond=5
      OPEN(unit_gap_cond,file='rv_gapcond.in',status='old',iostat=err)
!      
      IF(err.eq.0)THEN
! 
!........gap inside node index, gap outside node index
!      
         READ(unit_gap_cond,*) nr_gapi,nr_gapo  
!
!........initial gap pressure [Pa]
!         
         READ(unit_gap_cond,*) gap_P0
!
!........cladding plastic deformation (0:not included, 1:included)
!        
         READ(unit_gap_cond,*) iplas
!
!........Gap radial displacement due to fission gas-induced swelling and densification (must be >= 0) [m]
!      
         READ(unit_gap_cond,*) GapDisp_fission    
!
!........Gap roughness for fuel and clad [m]
!      
         READ(unit_gap_cond,*) GapRough_fuel,GapRough_clad 
         GapRough=GapRough_fuel+GapRough_clad 
!
!........gas component in the gap 
!        ngas=number of gases) must be >= 1
!      
         READ(unit_gap_cond,*) ngas 
         IF(ngas.gt.7.or.ngas.lt.0) THEN
            PRINT*,'***************************************'
            PRINT*,'              ERROR!                   '
            PRINT*,'ngas: 0<ngas<=7'
            PRINT*,'***************************************'
            PAUSE
            STOP
         ENDIF 
!
!........gas component in the gap 
!        indgas(i)=1:he, 2:h2, 3:n2, 4:kr, 5:xe, 6:o2, 7:ar
!        molgas(i)=!Mole fraction of He: mole fraction of ith gas/sum of mole fractions of all gases (mole fraction=weight fraction of ith gas/molecular weight of ith gas  
!
         indgas(:)=0
         molgas(:)=0.d0
         DO i=1,ngas
            READ(unit_gap_cond,*) indgas(i),molgas(i)
         ENDDO
!
!........rho*cp of gap material
!
         READ(unit_gap_cond,*) rhocp_gap
!     
      ELSE
         PRINT*,'***************************************'
         PRINT*,'              ERROR!                   '
         PRINT*,'If rv_gapcond options is turned on, rv_gapcond.in is necessary'
         PRINT*,''
         PRINT*,'***************************************'
         PAUSE
         STOP
! 
      ENDIF


      RETURN
      END SUBROUTINE rv_read_gapcond

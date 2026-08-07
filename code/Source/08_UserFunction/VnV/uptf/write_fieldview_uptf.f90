!
      SUBROUTINE write_fieldview_uptf
!
!     SAVE field view DATA for Paraview: 3 vectors, 48 scalars = 51 variables
!
      USE VOL_DATA    
      USE Zparam           , ONLY: nn,ndim
      USE Zcore            , ONLY: myrank,np
      USE Ztimecon         , ONLY: time
      USE Zvector          , ONLY: vg_n,vl_n
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all
      USE viewData_common  , ONLY: nframe,viwUnit,viwname,viewField
      USE Zio_unit         , ONLY: unit_log
!
      IMPLICIT NONE
! 
      INTEGER i,j,ivar,ix
      INTEGER nncup
      INTEGER nvector,nscalar   
!     
      REAL(8) crit_zero 
      REAL(8),ALLOCATABLE :: vg_all(:,:),vl_all(:,:)
      REAL(8),ALLOCATABLE :: alphag_all(:),alphal_all(:),alphad_all(:)         
      REAL(8),ALLOCATABLE :: vfgl_all(:),vfwl_all(:),vfwg_all(:)
      INTEGER,ALLOCATABLE :: regime_all(:)
!
      nncup=nn
      !IF(cupid_mars)nncup=ncell_old(1) 
      nvector=0
      nscalar=0          
!        
!.....OPEN file and WRITE header        
!
      IF (nframe == 0) THEN 
         IF(myrank.eq.0) THEN
            viwname ='somaPlot.viw'
            OPEN (unit=viwUnit, file=trim(viwname),status='replace', form='unformatted')
            WRITE(viwUnit) nncup, ncell_fluid
            WRITE(viwUnit) viewField%nVectors
            WRITE(viwUnit) (viewField%vectorVar(i), i=1,viewField%nVectors)
            WRITE(viwUnit) viewField%nScalars
            WRITE(viwUnit) (viewField%scalarVar(i), i=1,viewField%nScalars)
            nframe = nframe + 1
         ENDIF   
      ENDIF
!      
      IF(myrank.eq.0)WRITE(viwUnit) nframe,time
!
!.....Save alphag&alphal to check the zero fraction
!
      crit_zero=1.0d-4    ! criterion to check whether the phase fraction is zero or not
      ALLOCATE(alphag_all(nncup))
      alphag_all(:)=0.d0
      IF(np.gt.1) THEN               
         CALL allgatherv_r(cell%alphag,alphag_all,ncell_fluid,ncell_fluid_all,0)        
      ELSEIF(np.eq.1) THEN    
         alphag_all(1:ncell_fluid_all)=cell%alphag(1:ncell_fluid_all)
      ENDIF   
      ALLOCATE(alphal_all(nncup))
      alphal_all(:)=0.d0
      IF(np.gt.1) THEN               
         CALL allgatherv_r(cell%alphal,alphal_all,ncell_fluid,ncell_fluid_all,0)       
      ELSEIF(np.eq.1) THEN    
         alphal_all(1:ncell_fluid_all)=cell%alphal(1:ncell_fluid_all)
      ENDIF       
!
      DO ivar = 1, viewField%nVectors
!
!----------------------------------------------------------------------
!-----------------------Basic vector variables-------------------------
!----------------------------------------------------------------------               
!     
!........1. Gas-phase velocity
!
         IF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "vg" ) THEN
            nvector=nvector+1
            ALLOCATE(vg_all(nncup,ndim))
            vg_all(:,:)=0.d0
            IF(np.gt.1) THEN 
                DO ix=1,ndim
                   CALL allgatherv_r(vg_n(1,ix),vg_all(1,ix),ncell_fluid,ncell_fluid_all,0)
                ENDDO  
             ELSEIF(np.eq.1) THEN 
                DO ix=1,ndim
                   vg_all(1:ncell_fluid_all,ix)=vg_n(1:ncell_fluid_all,ix)
                ENDDO 
             ENDIF                              
!
!...........Make zero if the pahse fraction is almost zero
!            
            DO j=1, ncell_fluid_all
               IF(alphag_all(j).le.crit_zero)THEN
                  vg_all(j,:)=0.0d0
               ENDIF
            ENDDO
!
            IF(myrank.eq.0) THEN                   
               WRITE(viwUnit) (vg_all(i,1), i = 1, nncup)
               WRITE(viwUnit) (vg_all(i,2), i = 1, nncup)
               IF (ndim == 3) WRITE(viwUnit) (vg_all(i,3),  i = 1, nncup)
            ENDIF   
!
            IF(myrank.eq.0) WRITE(154,150) time,vg_all(536,3),vg_all(537,3),vg_all(538,3),vg_all(539,3),  &
                                                vg_all(320,3),vg_all(321,3),vg_all(322,3),vg_all(323,3),vg_all(324,3),vg_all(325,3),vg_all(326,3),vg_all(327,3),   &
                                                vg_all(662,3),vg_all(655,3)
150         FORMAT(20(e14.7,1x))                                                  
!
            DEALLOCATE(vg_all)
         ENDIF 
!
!........2. Liquid-phase velocity
!
         IF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "vl" ) THEN
            nvector=nvector+1
            ALLOCATE(vl_all(nncup,ndim))
            vl_all(:,:)=0.d0
            IF(np.gt.1) THEN 
               DO ix=1,ndim
                  CALL allgatherv_r(vl_n(1,ix),vl_all(1,ix),ncell_fluid,ncell_fluid_all,0)
               ENDDO         
            ELSEIF(np.eq.1) THEN    
               DO ix=1,ndim
                  vl_all(1:ncell_fluid_all,ix)=vl_n(1:ncell_fluid_all,ix)
               ENDDO   
            ENDIF            
!
!...........Make zero if the pahse fraction is almost zero
!            
            DO j=1, ncell_fluid_all
               IF(alphal_all(j).le.crit_zero)THEN
                  vl_all(j,:)=0.0d0
               ENDIF
            ENDDO
!            
            IF(myrank.eq.0) THEN   
               WRITE(viwUnit) (vl_all(i,1), i = 1, nncup)
               WRITE(viwUnit) (vl_all(i,2), i = 1, nncup)
               IF (ndim == 3)THEN
                  WRITE(viwUnit) (vl_all(i,3), i = 1, nncup)
               ENDIF   
             ENDIF   
!
            IF(myrank.eq.0) WRITE(155,150) time,vl_all(536,3),vl_all(537,3),vl_all(538,3),vl_all(539,3),  &
                                                vl_all(320,3),vl_all(321,3),vl_all(322,3),vl_all(323,3),vl_all(324,3),vl_all(325,3),vl_all(326,3),vl_all(327,3),   &
                                                vl_all(662,3),vl_all(655,3)
!
             
            DEALLOCATE(vl_all)
         ENDIF 
!
      ENDDO
!
      DO ivar = 1, viewField%nScalars
!      
!----------------------------------------------------------------------
!-----------------------Basic scalar variables-------------------------
!----------------------------------------------------------------------               
!         
!........4. Gas-phase volume fraction
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphag" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (alphag_all(i), i = 1, nncup)     
         ENDIF 
!
!.......5. Liquid-phase volume fraction
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphal" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (alphal_all(i), i = 1, nncup)   
         ENDIF 
!
!........6. Gas-phase volume fraction
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphad" ) THEN
            nscalar=nscalar+1
            ALLOCATE(alphad_all(nncup))
            alphad_all(:)=0.d0
            IF(np.gt.1) THEN               
               CALL allgatherv_r(cell%alphad,alphad_all,ncell_fluid,ncell_fluid_all,0)       
            ELSEIF(np.eq.1) THEN    
               alphad_all(1:ncell_fluid_all)=cell%alphad(1:ncell_fluid_all)               
            ENDIF 
            IF(myrank.eq.0) WRITE(viwUnit) (alphad_all(i), i = 1, nncup)      
            DEALLOCATE(alphad_all)
         ENDIF 
!        
!........24. Interfacial drag 
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "vfgl" ) THEN
            nscalar=nscalar+1
            ALLOCATE(vfgl_all(nncup))
            vfgl_all(:)=0.d0
            IF(np.gt.1) THEN
               CALL allgatherv_r(cell%vfgl,vfgl_all,ncell_fluid,ncell_fluid_all,0)       
            ELSEIF(np.eq.1) THEN    
               vfgl_all(1:ncell_fluid_all)=cell%vfgl(1:ncell_fluid_all)
            ENDIF 
            IF(myrank.eq.0) WRITE(viwUnit) (vfgl_all(i), i = 1, nncup)    
!
            IF(myrank.eq.0) WRITE(153,150) time,vfgl_all(536),vfgl_all(537),vfgl_all(538),vfgl_all(539),  &
                                                 vfgl_all(320),vfgl_all(321),vfgl_all(322),vfgl_all(323),vfgl_all(324),vfgl_all(325),vfgl_all(326),vfgl_all(327),   &
                                                 vfgl_all(662),vfgl_all(655)
!                   
            DEALLOCATE(vfgl_all)
         ENDIF
!        
!........52. Flow regime
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "regime" ) THEN
            nscalar=nscalar+1
            ALLOCATE(regime_all(nncup))
            regime_all(:)=0
            IF(np.gt.1) THEN             
               CALL allgatherv_i(cell%regime,regime_all,ncell_fluid,ncell_fluid_all,0)       
            ELSEIF(np.eq.1) THEN    
              regime_all(1:ncell_fluid_all)= cell%regime(1:ncell_fluid_all)
            ENDIF 
            IF(myrank.eq.0) WRITE(viwUnit) (REAL(regime_all(i)), i = 1, nncup)    
!
            IF(myrank.eq.0) WRITE(156,160) time,regime_all(536),regime_all(537),regime_all(538),regime_all(539),  &
                                                 regime_all(320),regime_all(321),regime_all(322),regime_all(323),regime_all(324),regime_all(325),regime_all(326),regime_all(327),   &
                                                 regime_all(662),regime_all(655)
!             
160         FORMAT(e14.7,1x,20(i4,1x)) 
            DEALLOCATE(regime_all)
         ENDIF
!         
!----------------------------------------------------------------------
!-----------------------Addition---------------------------------------
!----------------------------------------------------------------------
! 
!........58. Wall drag of liquid phase
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "vfwl" ) THEN
            nscalar=nscalar+1
            ALLOCATE(vfwl_all(nncup))
            vfwl_all(:)=0.d0
            IF(np.gt.1) THEN
               CALL allgatherv_r(cell%vfwl,vfwl_all,ncell_fluid,ncell_fluid_all,0)       
            ELSEIF(np.eq.1) THEN    
               vfwl_all(1:ncell_fluid_all)=cell%vfwl(1:ncell_fluid_all)
            ENDIF 
            IF(myrank.eq.0) WRITE(viwUnit) (vfwl_all(i), i = 1, nncup)     
!
            IF(myrank.eq.0) WRITE(151,150) time,vfwl_all(536),vfwl_all(537),vfwl_all(538),vfwl_all(539),  &
                                                 vfwl_all(320),vfwl_all(321),vfwl_all(322),vfwl_all(323),vfwl_all(324),vfwl_all(325),vfwl_all(326),vfwl_all(327),   &
                                                 vfwl_all(662),vfwl_all(655)
!            
            DEALLOCATE(vfwl_all)
         ENDIF   
!
!........80. Wall drag of gas phase
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "vfwg" ) THEN
            nscalar=nscalar+1
            ALLOCATE(vfwg_all(nncup))
            vfwg_all(:)=0.d0
            IF(np.gt.1) THEN
               CALL allgatherv_r(cell%vfwg,vfwg_all,ncell_fluid,ncell_fluid_all,0)       
            ELSEIF(np.eq.1) THEN    
               vfwg_all(1:ncell_fluid_all)=cell%vfwg(1:ncell_fluid_all)
            ENDIF 
            IF(myrank.eq.0) WRITE(viwUnit) (vfwg_all(i), i = 1, nncup)     
!
!
            IF(myrank.eq.0) WRITE(152,150) time,vfwg_all(536),vfwg_all(537),vfwg_all(538),vfwg_all(539),  &
                                                 vfwg_all(320),vfwg_all(321),vfwg_all(322),vfwg_all(323),vfwg_all(324),vfwg_all(325),vfwg_all(326),vfwg_all(327),   &
                                                 vfwg_all(662),vfwg_all(655)
!      
!            
            DEALLOCATE(vfwg_all)
         ENDIF 
!            
      ENDDO
!      
      DEALLOCATE(alphag_all)
      DEALLOCATE(alphal_all)
!           
      IF(nvector.ne.viewField%nVectors.or.nscalar.ne.viewField%nScalars)THEN
         IF(myrank.eq.0)THEN
            WRITE(*,*)'          Check variable names using Variable Info and write_fieldview.f90 !!!'
            WRITE(*,*)nvector,viewField%nVectors,nscalar,viewField%nScalars
            WRITE(unit_log,*)'          Check variable names using Variable Info and write_fieldview.f90 !!!'
            WRITE(unit_log,*)nvector,viewField%nVectors,nscalar,viewField%nScalars
         ENDIF
         STOP
      ENDIF
!   
      END SUBROUTINE write_fieldview_uptf

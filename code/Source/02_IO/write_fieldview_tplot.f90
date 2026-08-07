!
      SUBROUTINE write_fieldview_tplot
!
!     SAVE field view DATA for Paraview: 3 vectors, 48 scalars = 51 variables
!
      USE VOL_DATA         , ONLY: cell               
      USE Zmpi             , ONLY: ncell_fp
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nn,ndim
      USE Zconst1          , ONLY: restart
      USE Zio_unit         , ONLY: unit_tplotv,unit_tplots
      USE Ziat             , ONLY: ia
      USE Zpress           , ONLY: p
      USE Ztimecon         , ONLY: time
      USE Ztplot           , ONLY: time2view,tplot_cell,tplot_num       
      USE Zvector          , ONLY: vg_n,vl_n
      USE viewData_common  , ONLY: nvector,nscalar,crit_zero, &
                                   alphal_all,alphag_all, &
                                   viewField
!
      IMPLICIT NONE
! 
!.....Local variables
      INTEGER :: indx,indx_namelist,kk,ix !only for real-time-plot
      INTEGER :: ivar
      INTEGER :: ii,j 
      INTEGER :: nncup,na
      LOGICAL,SAVE :: initial_v=.true.
      LOGICAL,SAVE :: initial_s=.true.
!.....Local arrays
      CHARACTER*20  vector_list(100),scalar_list(100)
!.....Local allocatable arrays
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: vg_all,vl_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: tplot_all        !Only for real-time plot        
!      
      IF(tplot_num.eq.0) RETURN
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
      nvector=0
      nscalar=0       
!      
!.....Save alphag&alphal to check the zero fraction
!
      crit_zero=5.0d-3    ! criterion to check whether the phase fraction is zero or not
      IF(myrank.eq.0) THEN                   
         ALLOCATE(alphag_all(na),alphal_all(na))
      ELSE
         ALLOCATE(alphag_all(1),alphal_all(na))
      ENDIF
      CALL gatherv_r(cell%alphag,ncell_fluid,alphag_all,na,0)
      CALL gatherv_r(cell%alphal,ncell_fluid,alphal_all,na,0)
!
!
!.....Only for the real-time plot option (tplot_num is greater than 0)
!
      indx=0
      indx_namelist=0
      IF(tplot_num.gt.0) ALLOCATE(tplot_all(tplot_num,100))
      IF(initial_v) THEN
         IF(myrank.eq.0) THEN
            indx_namelist=indx_namelist+1
            vector_list(indx_namelist)="time"
         ENDIF      
      ENDIF      
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
            IF(myrank.eq.0) THEN
               ALLOCATE(vg_all(na,ndim))
            ELSE
               ALLOCATE(vg_all(1,ndim))
            ENDIF
            CALL gatherv_r_2d(vg_n,ncell_fp,vg_all,ncell_fluid,na,0)
            IF(myrank.eq.0) THEN                   
!
!..............Make zero if the pahse fraction is almost zero
!            
               DO j=1, na
                  IF(alphag_all(j).le.crit_zero) vg_all(j,:)=0.0d0
               ENDDO
!
!..............Only for real-time plot
!
               IF(initial_v) THEN
                  DO ix=1,ndim
                     indx_namelist=indx_namelist+1
                     IF(ix.eq.1) vector_list(indx_namelist)="vgx"
                     IF(ix.eq.2) vector_list(indx_namelist)="vgy"
                     IF(ix.eq.3) vector_list(indx_namelist)="vgz"
                  ENDDO
               ENDIF   
               DO ix=1,ndim
                  indx=indx+1                  
                  DO kk=1,tplot_num
                     tplot_all(kk,indx)=vg_all(tplot_cell(kk),ix)
                  ENDDO  
               ENDDO
            ENDIF   
            DEALLOCATE(vg_all)
         ENDIF 
!
!........2. Liquid-phase velocity
!
         IF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "vl" ) THEN
            nvector=nvector+1
            IF(myrank.eq.0) THEN
               ALLOCATE(vl_all(na,ndim))
            ELSE
               ALLOCATE(vl_all(1,ndim))
            ENDIF
            CALL gatherv_r_2d(vl_n,ncell_fp,vl_all,ncell_fluid,na,0)
            IF(myrank.eq.0) THEN                               
!
!..............Make zero if the pahse fraction is almost zero
!            
               DO j=1, na
                  IF(alphal_all(j).le.crit_zero) vl_all(j,:)=0.0d0
               ENDDO
!            
!..............Only for real-time plot
!
               IF(initial_v) THEN
                  DO ix=1,ndim     
                     indx_namelist=indx_namelist+1
                     IF(ix.eq.1) vector_list(indx_namelist)="vlx"
                     IF(ix.eq.2) vector_list(indx_namelist)="vly"
                     IF(ix.eq.3) vector_list(indx_namelist)="vlz"
                  ENDDO
               ENDIF  
               DO ix=1,ndim
                  indx=indx+1                  
                  DO kk=1,tplot_num
                     tplot_all(kk,indx)=vl_all(tplot_cell(kk),ix)
                  ENDDO  
               ENDDO
            ENDIF   
            DEALLOCATE(vl_all)
         ENDIF 
!         
      ENDDO
!
!.....only for real-time plot
!
      IF(myrank.eq.0.and.(time >= time2view.or.initial_v))THEN
         IF(initial_v.and.restart.eq.0) THEN
            DO kk=1,tplot_num
               ! write(550+kk,329) (vector_list(ii),ii=1,indx_namelist)
               write(unit_tplotv(kk),329) (vector_list(ii),ii=1,indx_namelist)
            ENDDO
         ENDIF      
         initial_v=.false. 
         DO kk=1,tplot_num
            ! write(550+kk,330) time,(tplot_all(kk,ii),ii=1,indx)
            write(unit_tplotv(kk),330) time,(tplot_all(kk,ii),ii=1,indx)
         ENDDO  
      ENDIF
 329  FORMAT(200(A14,1x))     
 330  FORMAT(200(e13.6,1x))      
      DEALLOCATE(tplot_all)
!
      indx=0
      indx_namelist=0
      IF(tplot_num.gt.0) ALLOCATE(tplot_all(tplot_num,200))   
      IF(initial_s.and.myrank.eq.0) THEN
         indx_namelist=indx_namelist+1
         scalar_list(indx_namelist)="time"
      ENDIF
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
!
            IF(myrank.eq.0) THEN              
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="alphag"
               ENDIF
               indx=indx+1                  
               DO kk=1,tplot_num
                  tplot_all(kk,indx)=alphag_all(tplot_cell(kk))
               ENDDO  
!               
            ENDIF 
         ENDIF 
!
!.......5. Liquid-phase volume fraction
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphal" ) THEN
            nscalar=nscalar+1
!
            IF(myrank.eq.0) THEN              
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="alphal"               
               ENDIF               
               indx=indx+1                  
               DO kk=1,tplot_num
                  tplot_all(kk,indx)=alphal_all(tplot_cell(kk))
               ENDDO  
!               
            ENDIF                        
         ENDIF 
!
!........6. Gas-phase volume fraction
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphad" ) THEN
            CALL wr1_1d(cell%alphad,tplot_all,indx)
            IF(myrank.eq.0) THEN              
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="alphad"               
               ENDIF               
            ENDIF                          
         ENDIF 
!
!........7. Gas-phase density
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "rhog" ) THEN
            CALL wr1_1d(cell%rhog,tplot_all,indx)
            IF(myrank.eq.0) THEN              
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="rhog"               
               ENDIF               
            ENDIF                          
         ENDIF 
!
!........8. Liquid-phase density
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "rhol" ) THEN
            CALL wr1_1d(cell%rhog,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="rhol"               
               ENDIF               
            ENDIF                                       
         ENDIF 
!
!........9. Gas-phase density
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "rhod" ) THEN
            CALL wr1_1d(cell%rhod,tplot_all,indx)
            IF(myrank.eq.0) THEN    
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="rhod"               
               ENDIF               
            ENDIF                                       
         ENDIF 
!
!........10. Pressure
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "p" ) THEN
            CALL wr1_1d(p,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="p"
               ENDIF               
            ENDIF                                     
         ENDIF  
!
!........11. Gas-phase temperature
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tg" ) THEN
            CALL wr1_1d_z(cell%tg,tplot_all,alphag_all,indx)
            IF(myrank.eq.0) THEN                
!               
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="tg"               
               ENDIF               
            ENDIF                                     
         ENDIF
!
!........12. Liquid-phase temperature
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tl" ) THEN
            CALL wr1_1d_z(cell%tl,tplot_all,alphal_all,indx)
            IF(myrank.eq.0) THEN                
!            
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="tl"               
               ENDIF               
            ENDIF               
         ENDIF 
!
!........13. Gas-phase energy
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "eg" ) THEN
            CALL wr1_1d(cell%eg,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="eg"               
               ENDIF               
            ENDIF                        
         ENDIF 
!
!........14. Liquid-phase energy
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "el" ) THEN
            CALL wr1_1d(cell%el,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="el"               
               ENDIF               
            ENDIF                                    
         ENDIF 
!        
!........15. Saturation temperature
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ts" ) THEN
            CALL wr1_1d(cell%ts,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="ts"               
               ENDIF               
            ENDIF                                    
         ENDIF
!        
!........16. Boron concentration
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "cboron" ) THEN
            CALL wr1_1d(cell%cboron,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="cboron"               
               ENDIF               
            ENDIF                                      
         ENDIF
!
!........17. Drop-phase energy
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ed" ) THEN
            CALL wr1_1d(cell%ed,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="ed"
               ENDIF               
            ENDIF                                     
         ENDIF          
!
!........18. Saturation gas enthalpy
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hg" ) THEN
            CALL wr1_1d(cell%hg,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="hg"
               ENDIF               
            ENDIF                                    
         ENDIF
!
!........19. Saturation liquid enthalpy
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hl" ) THEN
            CALL wr1_1d(cell%hl,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="hl"
               ENDIF               
            ENDIF                                      
         ENDIF  
!                
!----------------------------------------------------------------------
!-----------------------Interaface-Transfer-related variables----------
!----------------------------------------------------------------------
!        
!........20. Bubble diameter
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "dbubble" ) THEN
            CALL wr1_1d(cell%D1,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="dbubble"
               ENDIF               
            ENDIF                                      
         ENDIF
!        
!........21. Interfacial area concentration
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "iac" ) THEN
            CALL wr1_1d(ia,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............Only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="iac"
               ENDIF               
            ENDIF                                      
         ENDIF
!        
!........23. Quality 
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "quality" ) THEN
            CALL wr1_1d(cell%quala,tplot_all,indx)
            IF(myrank.eq.0) THEN                
!
!..............only for real-time plot
!
               IF(initial_s) THEN
                  indx_namelist=indx_namelist+1
                  scalar_list(indx_namelist)="quality"
               ENDIF               
            ENDIF                                      
         ENDIF
!
      ENDDO
!         
!........only for real-time plot
!
         IF(myrank.eq.0.and.(time >= time2view.or.initial_s))THEN
            IF(initial_s) THEN
               DO kk=1,tplot_num
                  ! write(660+kk,329) (scalar_list(ii),ii=1,indx_namelist)
                  write(unit_tplots(kk),329) (scalar_list(ii),ii=1,indx_namelist)
               ENDDO             
            ENDIF
            initial_s=.false.         
            DO kk=1,tplot_num
               write(unit_tplots(kk),330) time,(tplot_all(kk,ii),ii=1,indx)
            ENDDO  
         ENDIF   

!      
      DEALLOCATE(alphag_all)
      DEALLOCATE(alphal_all)
      DEALLOCATE(tplot_all)      
!           
      END SUBROUTINE write_fieldview_tplot
!
      SUBROUTINE wr1_1d(x,tplot_all,indx)
!
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: myrank
      USE Ztplot           , ONLY: tplot_cell,tplot_num       
      USE viewData_common  , ONLY: nscalar
!
      IMPLICIT NONE
! 
!.....Input
      INTEGER :: indx
      REAL(8) :: x(ncell_fluid)
      REAL(8) :: tplot_all(tplot_num,*)
!.....Local variables
      INTEGER :: i,k
      INTEGER :: na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp
!
      na=ncell_fluid_all
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN              
         ALLOCATE(tmp(na))
      ELSE
         ALLOCATE(tmp(1))
      ENDIF
      CALL gatherv_r(x,ncell_fluid,tmp,na,0)
!
      IF(myrank.eq.0) THEN              
         indx=indx+1                  
         DO i=1,tplot_num
            k=tplot_cell(i)
            tplot_all(i,indx)=tmp(k)
         ENDDO  
      ENDIF                          
      DEALLOCATE(tmp)
!           
      END SUBROUTINE wr1_1d
!
      SUBROUTINE wr1_1d_z(x,tplot_all,alpha_all,indx)
!
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: myrank
      USE Ztplot           , ONLY: tplot_cell,tplot_num       
      USE viewData_common  , ONLY: nscalar,crit_zero
!
      IMPLICIT NONE
! 
!.....Input
      INTEGER :: indx
      REAL(8) :: x(ncell_fluid)
      REAL(8) alpha_all(ncell_fluid_all)
      REAL(8) :: tplot_all(tplot_num,*)
!.....Local variables
      INTEGER :: i,k
      INTEGER :: na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp
!
      na=ncell_fluid_all
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN              
         ALLOCATE(tmp(na))
      ELSE
         ALLOCATE(tmp(1))
      ENDIF
      CALL gatherv_r(x,ncell_fluid,tmp,na,0)
!
      IF(myrank.eq.0) THEN              
!
!........Make zero if the pahse fraction is almost zero
!           
         DO i=1, na
            IF(alpha_all(i).le.crit_zero) tmp(i)=0.d0
          ENDDO
         indx=indx+1                  
         DO i=1,tplot_num
            k=tplot_cell(i)
            tplot_all(i,indx)=tmp(k)
         ENDDO  
      ENDIF                          
      DEALLOCATE(tmp)
!           
      END SUBROUTINE wr1_1d_z

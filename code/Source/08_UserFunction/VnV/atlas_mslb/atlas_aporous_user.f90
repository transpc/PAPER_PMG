!
      SUBROUTINE atlas_aporous_user
!
!     Define Aporous (only when "udfl_heat_partition_porous" is used)
!    
      USE SOLID_DATA   , ONLY: solid
      USE Zparam       , ONLY: pi
      USE Zconst1      , ONLY: vv_prob
      USE Zcoord3      , ONLY: aporous
      USE Zcoord3      , ONLY: aporous,vol,porosity  
      USE Zcore        , ONLY: np,myrank
      USE Znum_cell    , ONLY: n_fluid
      USE Zzone        , ONLY: ncell_fluid,ncell_cond
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER i,ii
!      
      REAL(8) strucvol,aporoussum,strucmas
      REAL(8) porous_s  !mcc-pik
!      
      LOGICAL,SAVE ::initial
!      
      DATA initial/.TRUE./
!      
      IF(initial)then
!      
         initial=.FALSE.
         strucvol=0.0d0
         strucmas=0.0d0 
         aporoussum=0.0d0
!         
         IF(vv_prob.eq.'atlas_mc')then
             DO i=1,ncell_fluid
                porous_s=1.0d0-porosity(ii)
                strucvol=strucvol+vol(i)*porous_s
             ENDDO  
         ELSEIF(vv_prob.eq.'atlas_mc_porous')then
             DO i=1,ncell_cond
                ii=n_fluid(i)
                porous_s=1.0d0-porosity(ii)
                strucvol=strucvol+vol(ii)*porous_s
                strucmas=strucmas+vol(ii)*porous_s*solid%rhocps(i)/500.d0
             ENDDO
         ENDIF
         IF(np.gt.1)then
            CALL allreducei_r1(strucvol)
            CALL allreducei_r1(strucmas)
         ENDIF   
!
        IF(vv_prob.eq.'atlas_mc'.or.vv_prob.eq.'atlas_mc_porous')then
             DO i=1,ncell_fluid
                porous_s=1.0d0-porosity(i)
                aporous(i)=2.0d0*pi*(0.095d0/2.d0)*1.9d0*390.0d0*vol(i)*porous_s/strucvol !must be saved
                aporoussum=aporoussum+aporous(i)
             ENDDO  
         ENDIF
         IF(np.gt.1)then
            CALL allreducei_r1(aporoussum)
         ENDIF          
!                    
         IF(myrank.eq.0)then
            WRITE(*,"(11x,a,1pe17.5,a)")'ATLAS RPV solid structure  volume=',strucvol,'m3' 
            WRITE(*,"(11x,a,1pe17.5,a)")'ATLAS RPV solid structure  mass=',strucmas/1000.d0,'ton' 
            WRITE(*,"(11x,a,1pe17.5,a)")'ATLAS RPV solid structure  area=',aporoussum,'m2'             
            WRITE(unit_log,"(11x,a,1pe17.5,a)")'ATLAS RPV solid structure volume=',strucvol,'m3' 
            WRITE(unit_log,"(11x,a,1pe17.5,a)")'ATLAS RPV solid structure mass=',strucmas/1000.d0,'ton' 
            WRITE(unit_log,"(11x,a,1pe17.5,a)")'ATLAS RPV solid structure area=',aporoussum,'m2' 
         ENDIF
!            
      ELSE
!      
      ENDIF 
!      
      RETURN
      END SUBROUTINE atlas_aporous_user         
        

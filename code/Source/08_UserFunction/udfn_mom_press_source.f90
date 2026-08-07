!
      SUBROUTINE udfn_mom_press_source
!
!     Modifies the momentum source terms at the free surface cells
!
      USE VOL_DATA    , ONLY: cell         
      USE Zzone       , ONLY: ncell_fluid
      USE Zparam      , ONLY: ndim
      USE Zconst1     , ONLY: vv_prob,topolsurface
      USE Zconst2     , ONLY: grav
      USE Zpress      , ONLY: dpdx
      USE Zvoid       , ONLY: gamma_void,gradient_void      
      USE Zm_src      , ONLY: src_gas,src_liq,src_drp
!
      IMPLICIT NONE
!            
!.....Local variables
      INTEGER i,ix,igrav
!   
      LOGICAL :: freesurf=.false.
!      
!.....Local arrays
!     
!........PAFS-POOL, SMALL-POOL, fluidic_device, PAFS-3D
!
!         IF(vv_prob.eq.'PAFS-POOL'.or.vv_prob.eq.'fluidic_device'.or.  &
!            vv_prob.eq.'DIVA-NEW'.or.vv_prob.eq.'UPTF-RV') THEN
         IF(vv_prob.eq.'PAFS-POOL'.or. vv_prob.eq.'fluidic_device'.or. &
            vv_prob.eq.'DIVA-NEW' .or. vv_prob.eq.'UPTF-RV'       .or. &
            vv_prob.eq.'ST2-CT-01'.or. vv_prob.eq.'ST2-CT-02'     .or. &
            vv_prob.eq.'ST2-CT-03'                                     ) THEN
!             
            gradient_void=0.4d0    
            
! 
!...........Find gravity direction
!
            DO igrav=1,ndim
               IF(grav(igrav)*grav(igrav).gt.0.5d0)EXIT
            ENDDO  
!                
            DO i=1,ncell_fluid
! 
!..............Free surface check for both topology map or not case
!
               IF(topolsurface.eq.1.and.cell%regime(i).eq.3) freesurf=.true. 
               IF(topolsurface.eq.0.and.gamma_void(i).gt.gradient_void) freesurf=.true.          
!            
               DO ix=1,ndim
                  IF(freesurf.and.ix.eq.igrav)THEN
                     src_gas(i,ix)=src_gas(i,ix)-cell%alphag(i)*(cell%rhog(i)/  &
                                   (cell%rhog(i)*cell%alphag(i)+cell%rhol(i)*cell%alphal(i)))*dpdx(i,ix)
                     src_liq(i,ix)=src_liq(i,ix)-cell%alphal(i)*(cell%rhol(i)/  &
                                   (cell%rhog(i)*cell%alphag(i)+cell%rhol(i)*cell%alphal(i)))*dpdx(i,ix)
                     src_drp(i,ix)=src_drp(i,ix)-cell%alphad(i)*dpdx(i,ix)
                  ELSE
                     src_gas(i,ix)=src_gas(i,ix)-cell%alphag(i)*dpdx(i,ix)
                     src_liq(i,ix)=src_liq(i,ix)-cell%alphal(i)*dpdx(i,ix)
                     src_drp(i,ix)=src_drp(i,ix)-cell%alphad(i)*dpdx(i,ix)
                  ENDIF
               ENDDO
               freesurf=.false.
            ENDDO
         ENDIF
!
      RETURN
      END SUBROUTINE udfn_mom_press_source

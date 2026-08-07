!
      SUBROUTINE vectorize_deallocate_face
!      
      USE Zcoord1       , ONLY: xloc_m
      USE Zcoord2       , ONLY: fac,fac1,xfc, &
                                fac_c,fac1_c
      USE Zcoord3       , ONLY: sv,svp,permeability
      USE Zcoord4       , ONLY: sa,saa,sap,dji,dji_x,dji_a,sad,sap_c
      USE Znormal       , ONLY: xn      
!      
      IMPLICIT NONE
!      
      DEALLOCATE(xfc   )
      DEALLOCATE(xn    )
      DEALLOCATE(xloc_m)
      DEALLOCATE(sv    ) 
      DEALLOCATE(svp   )
      DEALLOCATE(sa    )  
      DEALLOCATE(saa   )
      DEALLOCATE(sad   )
      DEALLOCATE(sap   )
!      
      DEALLOCATE(dji_x)
      DEALLOCATE(dji  )
      DEALLOCATE(dji_a)
      DEALLOCATE(fac  )
      DEALLOCATE(fac1 )
      DEALLOCATE(permeability)
!
      IF(ALLOCATED(sap_c )) DEALLOCATE(sap_c )
      IF(ALLOCATED(fac_c )) DEALLOCATE(fac_c )
      IF(ALLOCATED(fac1_c)) DEALLOCATE(fac1_c)
!
!     DEALLOCATE(ecnvc_l   )   !removed
!     DEALLOCATE(ecnvc_g   )   !removed
!     DEALLOCATE(ecnvc_d   )   !removed
!
!     DEALLOCATE(vap_conv  )   !removed   
!     DEALLOCATE(void_conv )   !removed
!     DEALLOCATE(liq_conv  )   !removed
!     DEALLOCATE(drp_conv  )   !removed
!     DEALLOCATE(quala_conv)   !removed
!     DEALLOCATE(al_conv   )   !removed
!     DEALLOCATE(ad_conv   )   !removed
!
      RETURN
!      
      ENDSUBROUTINE vectorize_deallocate_face 

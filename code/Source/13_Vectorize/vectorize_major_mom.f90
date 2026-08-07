      SUBROUTINE vectorize_major_mom
!
!     This routine vectorize volume fluxes.
!
      USE Zcore         ,ONLY:np
      USE Zvec_param
      USE Zvec_major    ,ONLY: fluxl_cx_non,fluxg_cx_non,fluxd_cx_non,&
                               fluxl_cy_non,fluxg_cy_non,fluxd_cy_non,&
                               fluxl_cz_non,fluxg_cz_non,fluxd_cz_non,&
                               fluxl_cx_inl,fluxg_cx_inl,fluxd_cx_inl,&
                               fluxl_cy_inl,fluxg_cy_inl,fluxd_cy_inl,&
                               fluxl_cz_inl,fluxg_cz_inl,fluxd_cz_inl,&
                               fluxl_cx_out,fluxg_cx_out,fluxd_cx_out,&
                               fluxl_cy_out,fluxg_cy_out,fluxd_cy_out,&
                               fluxl_cz_out,fluxg_cz_out,fluxd_cz_out,&
                               fluxl_cx_mcc,fluxg_cx_mcc,fluxd_cx_mcc,&
                               fluxl_cy_mcc,fluxg_cy_mcc,fluxd_cy_mcc,&
                               fluxl_cz_mcc,fluxg_cz_mcc,fluxd_cz_mcc 
      USE Zvec_major    ,ONLY: mfluxl_cup_non,mfluxg_cup_non,mfluxd_cup_non,&
                                mfluxl_cup_inl,mfluxg_cup_inl,mfluxd_cup_inl,&
                                mfluxl_cup_out,mfluxg_cup_out,mfluxd_cup_out,&
                                mfluxl_cup_mcc,mfluxg_cup_mcc,mfluxd_cup_mcc                          
!
      IMPLICIT NONE
!
      INTEGER i
!
      LOGICAL,SAVE:: initial
!
      DATA initial/.true./ 
!
      IF(initial)THEN
         IF(ALLOCATED(fluxl_cx_non)) DEALLOCATE(fluxl_cx_non)
         IF(ALLOCATED(fluxg_cx_non)) DEALLOCATE(fluxg_cx_non)
         IF(ALLOCATED(fluxd_cx_non)) DEALLOCATE(fluxd_cx_non)
         IF(ALLOCATED(fluxl_cy_non)) DEALLOCATE(fluxl_cy_non)
         IF(ALLOCATED(fluxg_cy_non)) DEALLOCATE(fluxg_cy_non)
         IF(ALLOCATED(fluxd_cy_non)) DEALLOCATE(fluxd_cy_non)
         IF(ALLOCATED(fluxl_cz_non)) DEALLOCATE(fluxl_cz_non)
         IF(ALLOCATED(fluxg_cz_non)) DEALLOCATE(fluxg_cz_non)
         IF(ALLOCATED(fluxd_cz_non)) DEALLOCATE(fluxd_cz_non)
!
         IF(ALLOCATED(fluxl_cx_inl)) DEALLOCATE(fluxl_cx_inl)
         IF(ALLOCATED(fluxg_cx_inl)) DEALLOCATE(fluxg_cx_inl)
         IF(ALLOCATED(fluxd_cx_inl)) DEALLOCATE(fluxd_cx_inl)
         IF(ALLOCATED(fluxl_cy_inl)) DEALLOCATE(fluxl_cy_inl)
         IF(ALLOCATED(fluxg_cy_inl)) DEALLOCATE(fluxg_cy_inl)
         IF(ALLOCATED(fluxd_cy_inl)) DEALLOCATE(fluxd_cy_inl)
         IF(ALLOCATED(fluxl_cz_inl)) DEALLOCATE(fluxl_cz_inl)
         IF(ALLOCATED(fluxg_cz_inl)) DEALLOCATE(fluxg_cz_inl)
         IF(ALLOCATED(fluxd_cz_inl)) DEALLOCATE(fluxd_cz_inl)
!
         IF(ALLOCATED(fluxl_cx_out)) DEALLOCATE(fluxl_cx_out)
         IF(ALLOCATED(fluxg_cx_out)) DEALLOCATE(fluxg_cx_out)
         IF(ALLOCATED(fluxd_cx_out)) DEALLOCATE(fluxd_cx_out)
         IF(ALLOCATED(fluxl_cy_out)) DEALLOCATE(fluxl_cy_out)
         IF(ALLOCATED(fluxg_cy_out)) DEALLOCATE(fluxg_cy_out)
         IF(ALLOCATED(fluxd_cy_out)) DEALLOCATE(fluxd_cy_out)
         IF(ALLOCATED(fluxl_cz_out)) DEALLOCATE(fluxl_cz_out)
         IF(ALLOCATED(fluxg_cz_out)) DEALLOCATE(fluxg_cz_out)
         IF(ALLOCATED(fluxd_cz_out)) DEALLOCATE(fluxd_cz_out)
!
         IF(ALLOCATED(fluxl_cx_mcc)) DEALLOCATE(fluxl_cx_mcc)
         IF(ALLOCATED(fluxg_cx_mcc)) DEALLOCATE(fluxg_cx_mcc)
         IF(ALLOCATED(fluxd_cx_mcc)) DEALLOCATE(fluxd_cx_mcc)
         IF(ALLOCATED(fluxl_cy_mcc)) DEALLOCATE(fluxl_cy_mcc)
         IF(ALLOCATED(fluxg_cy_mcc)) DEALLOCATE(fluxg_cy_mcc)
         IF(ALLOCATED(fluxd_cy_mcc)) DEALLOCATE(fluxd_cy_mcc)
         IF(ALLOCATED(fluxl_cz_mcc)) DEALLOCATE(fluxl_cz_mcc)
         IF(ALLOCATED(fluxg_cz_mcc)) DEALLOCATE(fluxg_cz_mcc)
         IF(ALLOCATED(fluxd_cz_mcc)) DEALLOCATE(fluxd_cz_mcc)
!
         IF(ALLOCATED(mfluxl_cup_non)) DEALLOCATE(mfluxl_cup_non)
         IF(ALLOCATED(mfluxg_cup_non)) DEALLOCATE(mfluxg_cup_non)
         IF(ALLOCATED(mfluxd_cup_non)) DEALLOCATE(mfluxd_cup_non)
!
         IF(ALLOCATED(mfluxl_cup_inl)) DEALLOCATE(mfluxl_cup_inl)
         IF(ALLOCATED(mfluxg_cup_inl)) DEALLOCATE(mfluxg_cup_inl)
         IF(ALLOCATED(mfluxd_cup_inl)) DEALLOCATE(mfluxd_cup_inl)
!
         IF(ALLOCATED(mfluxl_cup_out)) DEALLOCATE(mfluxl_cup_out)
         IF(ALLOCATED(mfluxg_cup_out)) DEALLOCATE(mfluxg_cup_out)
         IF(ALLOCATED(mfluxd_cup_out)) DEALLOCATE(mfluxd_cup_out)
!
         IF(ALLOCATED(mfluxl_cup_mcc)) DEALLOCATE(mfluxl_cup_mcc)
         IF(ALLOCATED(mfluxg_cup_mcc)) DEALLOCATE(mfluxg_cup_mcc)
         IF(ALLOCATED(mfluxd_cup_mcc)) DEALLOCATE(mfluxd_cup_mcc)
!
         ALLOCATE(fluxl_cx_non(nf_non),fluxg_cx_non(nf_non),fluxd_cx_non(nf_non)) 
         ALLOCATE(fluxl_cy_non(nf_non),fluxg_cy_non(nf_non),fluxd_cy_non(nf_non)) 
         ALLOCATE(fluxl_cz_non(nf_non),fluxg_cz_non(nf_non),fluxd_cz_non(nf_non)) 
         ALLOCATE(fluxl_cx_inl(nf_inl),fluxg_cx_inl(nf_inl),fluxd_cx_inl(nf_inl)) 
         ALLOCATE(fluxl_cy_inl(nf_inl),fluxg_cy_inl(nf_inl),fluxd_cy_inl(nf_inl)) 
         ALLOCATE(fluxl_cz_inl(nf_inl),fluxg_cz_inl(nf_inl),fluxd_cz_inl(nf_inl)) 
         ALLOCATE(fluxl_cx_out(nf_out),fluxg_cx_out(nf_out),fluxd_cx_out(nf_out)) 
         ALLOCATE(fluxl_cy_out(nf_out),fluxg_cy_out(nf_out),fluxd_cy_out(nf_out)) 
         ALLOCATE(fluxl_cz_out(nf_out),fluxg_cz_out(nf_out),fluxd_cz_out(nf_out)) 
         ALLOCATE(fluxl_cx_mcc(nf_mcc),fluxg_cx_mcc(nf_mcc),fluxd_cx_mcc(nf_mcc)) 
         ALLOCATE(fluxl_cy_mcc(nf_mcc),fluxg_cy_mcc(nf_mcc),fluxd_cy_mcc(nf_mcc)) 
         ALLOCATE(fluxl_cz_mcc(nf_mcc),fluxg_cz_mcc(nf_mcc),fluxd_cz_mcc(nf_mcc)) 
         ALLOCATE(mfluxl_cup_non(nf_non),mfluxg_cup_non(nf_non),mfluxd_cup_non(nf_non)) 
         ALLOCATE(mfluxl_cup_inl(nf_inl),mfluxg_cup_inl(nf_inl),mfluxd_cup_inl(nf_inl)) 
         ALLOCATE(mfluxl_cup_out(nf_out),mfluxg_cup_out(nf_out),mfluxd_cup_out(nf_out)) 
         ALLOCATE(mfluxl_cup_mcc(nf_mcc),mfluxg_cup_mcc(nf_mcc),mfluxd_cup_mcc(nf_mcc))          
      ENDIF
      fluxl_cx_non=0.0d0  
      fluxl_cy_non=0.0d0  
      fluxl_cz_non=0.0d0  
      fluxl_cx_inl=0.0d0  
      fluxl_cy_inl=0.0d0  
      fluxl_cz_inl=0.0d0  
      fluxl_cx_out=0.0d0  
      fluxl_cy_out=0.0d0  
      fluxl_cz_out=0.0d0  
      fluxl_cx_mcc=0.0d0  
      fluxl_cy_mcc=0.0d0  
      fluxl_cz_mcc=0.0d0  
      fluxg_cx_non=0.0d0  
      fluxg_cy_non=0.0d0  
      fluxg_cz_non=0.0d0  
      fluxg_cx_inl=0.0d0  
      fluxg_cy_inl=0.0d0  
      fluxg_cz_inl=0.0d0  
      fluxg_cx_out=0.0d0  
      fluxg_cy_out=0.0d0  
      fluxg_cz_out=0.0d0  
      fluxg_cx_mcc=0.0d0  
      fluxg_cy_mcc=0.0d0  
      fluxg_cz_mcc=0.0d0  
      fluxd_cx_non=0.0d0  
      fluxd_cy_non=0.0d0  
      fluxd_cz_non=0.0d0  
      fluxd_cx_inl=0.0d0  
      fluxd_cy_inl=0.0d0  
      fluxd_cz_inl=0.0d0  
      fluxd_cx_out=0.0d0  
      fluxd_cy_out=0.0d0  
      fluxd_cz_out=0.0d0  
      fluxd_cx_mcc=0.0d0  
      fluxd_cy_mcc=0.0d0  
      fluxd_cz_mcc=0.0d0  
      mfluxl_cup_non=0.0d0
      mfluxl_cup_inl=0.0d0
      mfluxl_cup_out=0.0d0
      mfluxg_cup_non=0.0d0
      mfluxg_cup_inl=0.0d0
      mfluxg_cup_out=0.0d0
      mfluxd_cup_non=0.0d0
      mfluxd_cup_inl=0.0d0
      mfluxd_cup_out=0.0d0
!
!      initial=.FALSE.
!      
      RETURN
      END SUBROUTINE vectorize_major_mom


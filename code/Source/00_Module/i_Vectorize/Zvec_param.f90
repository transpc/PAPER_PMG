      MODULE Zvec_param
! 
      IMPLICIT NONE
      SAVE
!
!DEC$IF defined (AVX512)
      INTEGER,PARAMETER :: vl_f=8
!DEC$ELSEIF defined (AVX2)
      INTEGER,PARAMETER :: vl_f=4
!DEC$ELSEIF defined (AVX)
      INTEGER,PARAMETER :: vl_f=4
!DEC$ELSEIF defined (SSE)
      INTEGER,PARAMETER :: vl_f=2
!DEC$ELSE
      INTEGER,PARAMETER :: vl_f=1
!DEC$ENDIF
!
      INTEGER :: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym, &
                 nf_hconv,nf_hcond,nf_hvad
      INTEGER :: nf_inl_old,nf_out_old,nf_adw_old,nf_flux_old ! for nbcon_change
      INTEGER :: nf_totk,nf_tot,nf_nbcon_change,nf_nbcon_change_flux,nf_tot_nbcon,nf_tot_svp
      INTEGER :: nf_flux,nf_fluxk,nf_flux1,nf_fluxk1,nf_flux2,nf_fluxk2
      INTEGER :: nf_tot_nb1
!
!.....Solid
!
      INTEGER :: nfc_nonk
      INTEGER :: nfc_non,nfc_fsw,nfc_ctw,nfc_chw,nfc_chtcw                
      INTEGER :: nfc_tot,nfc_tot1
      INTEGER :: nfc_tot_nb1
!
      END MODULE Zvec_param

!   
      SUBROUTINE convective_heat_solid
!      
      USE SOLID_DATA       , ONLY: solid
      USE Zvec_param       , ONLY: nfc_nonk,nfc_non,nfc_fsw,nfc_ctw,nfc_chw,nfc_tot
      USE Znum_cell        , ONLY: istartc_nf,                  &
                                   lens,nf_number_id,istart_nfs
      USE Zvec_index_solid , ONLY: left_solid_nf,    &
                                   nbcon_solid_chtcw
      USE Zqvol            , ONLY: qconv_sol,htc_convw,tb_convw,ha_convw
      USE Zvec_geo         , ONLY: sap_c_nf      
!
      IMPLICIT NONE
!
      INTEGER :: i,k,ii
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i1,i2
      LOGICAL,SAVE::initial=.true.
      REAL(8),SAVE,ALLOCATABLE::hap_c_nf(:)
!
!.....Build summation info for non(0),fsw(1),ctw(2),chw(3),chtcw(4)
!
      nf_number_id=2 !???
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3                     
      nf_number_id(4)=4
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nfc_nonk
      istart_nfs(1)=istart_nfs(0) +nfc_non
      istart_nfs(2)=istart_nfs(1) +nfc_fsw
      istart_nfs(3)=istart_nfs(2) +nfc_ctw 
      istart_nfs(4)=istart_nfs(3) +nfc_chw 
      lens         =istart_nfs(4)    
!     
!.....Set heat transfer area      
!.....Skip calculating qconv_sol because solid%tsol is not uninitialized
!      
      IF(initial)THEN
         initial=.FALSE.
         ALLOCATE(hap_c_nf(nfc_tot))
         hap_c_nf(:)=0.0d0
         nv=4                            
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart2=istart0-istart_nfs(1)
         istart=istartc_nf(1,nf_number)
         len   =istartc_nf(2,nf_number)
         DO i=1,len
            i1=istart+i  !cell 
            i2=istart2+i !source
            ii=left_solid_nf(i1)
            k=-(nbcon_solid_chtcw(i)+30)
            IF(ha_convw(k).le.0.0d0)THEN
              hap_c_nf(i1)=sap_c_nf(i1)
            ELSE
              hap_c_nf(i1)=ha_convw(k)
            ENDIF   
         ENDDO          
         RETURN
      ENDIF
!      
!.....Set convetive heat loss from solid wall
!      
      qconv_sol(:)=0.0d0
      nv=4                            
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart2=istart0-istart_nfs(1)
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      DO i=1,len
         i1=istart+i          !a face 
         i2=istart2+i         !source for a face
         ii=left_solid_nf(i1) !a cell
         k=-(nbcon_solid_chtcw(i)+30)
         qconv_sol(ii)=qconv_sol(ii)-htc_convw(k)*(solid%tsol(ii)-tb_convw(k))*hap_c_nf(i1) 
      ENDDO   
!      
      RETURN
      END SUBROUTINE convective_heat_solid
!

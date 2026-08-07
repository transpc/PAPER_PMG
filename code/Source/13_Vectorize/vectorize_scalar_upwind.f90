      SUBROUTINE vectorize_scalar_upwind
!
!     This routine vectorizes basic scalar variables
!
      USE Znum_cell     ,ONLY: istart_nf
      USE Zvec_param    ,ONLY: nf_flux,nf_non
      USE Znum_cell     ,ONLY: i_neigh
      USE Zvec_index    ,ONLY: jneigh_nf,left_nf,right_non
      USE Zbc_index     ,ONLY: nbcon
      USE Zare          ,ONLY: ar_gas,ar_liq,ar_drp      
      USE Zb_condition  ,ONLY: alphab_gas,alphab_liq,alphab_drp,rhob_gas,rhob_liq,rhob_drp
      USE c3com_cupid   ,ONLY: i3invtbl
      USE Zvec_major    ,ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zvec_scalar   ,ONLY: arli_nf,argi_nf,ardi_nf,           &
                               arli_noi,argi_noi,arli_nok,argi_nok
      
!
      IMPLICIT NONE
!
      INTEGER i,k,j0
      INTEGER :: ii,jj,kk,idx
      INTEGER :: nf_number,istart,isize,i1
!
      LOGICAL,SAVE:: initial
!
      DATA initial/.true./ 
!
      IF(initial)THEN
         ALLOCATE(arli_nf(nf_flux),argi_nf(nf_flux),ardi_nf(nf_flux))
         ALLOCATE(arli_noi(nf_non),arli_nok(nf_non),argi_noi(nf_non),argi_nok(nf_non)) !turb_ke_convection.f90 lsj
         initial=.FALSE.     
!!!         DO i=1,ncell_fluid
!!!             ar_liq(i)=cell%alphal(i)*cell%rhol(i)
!!!             ar_gas(i)=cell%alphag(i)*cell%rhog(i)
!!!         ENDDO
!!!         CALL communicate_1d(ar_liq)
!!!         CALL communicate_1d(ar_gas)
      ENDIF
! bug
!     ar_liq(:)=cell%alphal(:)*cell%rhol(:)
!     ar_gas(:)=cell%alphag(:)*cell%rhog(:)
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)

           IF    (flux_l_nf(i1).eq.0.0d0)THEN
              arli_nf(i1)=.5d0*(ar_liq(ii)+ar_liq(kk))
           ELSEIF(flux_l_nf(i1).gt.0.0d0)THEN
              arli_nf(i1)=ar_liq(ii)
           ELSE
              arli_nf(i1)=ar_liq(kk)
           ENDIF
           IF    (flux_g_nf(i1).eq.0.0d0)THEN
              argi_nf(i1)=.5d0*(ar_gas(ii)+ar_gas(kk))
           ELSEIF(flux_g_nf(i1).gt.0.0d0)THEN
              argi_nf(i1)=ar_gas(ii)
           ELSE
              argi_nf(i1)=ar_gas(kk)
           ENDIF
           IF    (flux_d_nf(i1).eq.0.0d0)THEN
              ardi_nf(i1)=.5d0*(ar_drp(ii)+ar_drp(kk))
           ELSEIF(flux_d_nf(i1).gt.0.0d0)THEN
              ardi_nf(i1)=ar_drp(ii)
           ELSE
              ardi_nf(i1)=ar_drp(kk)
           ENDIF
!
!          turb_ke_convection.f90 lsj
           arli_noi(i)=ar_liq(ii)          
           argi_noi(i)=ar_gas(ii)          
           arli_nok(i)=ar_liq(kk)         
           argi_nok(i)=ar_gas(kk)  
!                               
      ENDDO
!
      nf_number=2
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         j0=i_neigh(ii)-1
         jj=jneigh_nf(i1)
           k=nbcon(jj+j0)
           arli_nf(i1)=rhob_liq(k)*alphab_liq(k)
           argi_nf(i1)=rhob_gas(k)*alphab_gas(k)
           ardi_nf(i1)=rhob_drp(k)*alphab_drp(k)
!
      ENDDO           
!
      nf_number=3
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
           arli_nf(i1)=ar_liq(ii)
           argi_nf(i1)=ar_gas(ii)
           ardi_nf(i1)=ar_drp(ii)
      ENDDO 
!
      nf_number=1
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         idx=i3invtbl(i)
         arli_nf(i1)=ar_liq(ii)
         argi_nf(i1)=ar_gas(ii)
         ardi_nf(i1)=ar_drp(ii)
      ENDDO         
!              
      RETURN
      END SUBROUTINE vectorize_scalar_upwind

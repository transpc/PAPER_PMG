!
      SUBROUTINE turb_ke_convection(ke,dp,ar,keb,dpb,ab,rb,flux_nf,conv_k,conv_e)
!
!     This routine calculates convection terms of the k-epsilon transport equation.
!   
      USE Zinterface
      USE Zmpi       , ONLY: ncell_fp
      USE Zzone      , ONLY: ncell_fluid
      USE Zparam     , ONLY: nin_max
      USE Zvec_param , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_flux,nf_fluxk
      USE Znum_cell  , ONLY: istart_nf,istart_nbcon_nf, &
                             nf_number_nb,lens,                           &
                             right_nb_k,istart_nfs,nf_number_id,istart_nf
      USE Zvec_index , ONLY: left_nf,right_non,nbcon_nf
      USE Z2nd_order , ONLY: turb_conv_2nd
!
      IMPLICIT NONE
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: ke,dp,ar
      REAL(8),DIMENSION(nin_max) :: keb,dpb,ab,rb
      REAL(8),DIMENSION(nf_flux) :: flux_nf
!.....Output
      REAL(8),DIMENSION(ncell_fluid) ::  conv_k,conv_e
!.....Local variables
      INTEGER :: nv,nf_number,istart,len,istart0,istart2,i0,i1,i2
      INTEGER :: i,k
      INTEGER :: ii,kk
      REAL(8) ark1,ark2,are1,are2,ar1,ar2
      REAL(8) conv_k1,conv_e1,conv_m1
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: conv_k_non,conv_e_non,conv_m_non
      REAL(8),DIMENSION(nf_fluxk) :: conv_k_nf,conv_e_nf
!     
!.....Cells non with 2nd order convection
!
      IF(turb_conv_2nd.gt.0) THEN
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!       
            ar1=ar(ii)
            ark1=ke(ii)*ar1
            are1=dp(ii)*ar1
!       
            ar2=ar(kk)
            ark2=ke(kk)*ar2
            are2=dp(kk)*ar2
!
            conv_k_non(i)=MIN(flux_nf(i1),0.d0)*ark2+MAX(flux_nf(i1),0.d0)*ark1
            conv_e_non(i)=MIN(flux_nf(i1),0.d0)*are2+MAX(flux_nf(i1),0.d0)*are1
            conv_m_non(i)=MIN(flux_nf(i1),0.d0)*ar2 +MAX(flux_nf(i1),0.d0)*ar1
         ENDDO
!
         CALL turb_2nd_conv(ke,dp,ar,flux_nf,conv_k_non,conv_e_non)
      ENDIF
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_mcc
      istart_nfs(3)=istart_nfs(2) +nf_inl
      lens         =istart_nfs(3) +nf_out
!
!
!.....Computing Cell
!
      IF(turb_conv_2nd.gt.0) THEN
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            conv_k_nf(i0)=conv_k_non(i)-ke(ii)*conv_m_non(i)
            conv_e_nf(i0)=conv_e_non(i)-dp(ii)*conv_m_non(i)
         ENDDO
!
         nv=-1
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(0)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            conv_k_nf(i)  =-conv_k_non(k)+ke(ii)*conv_m_non(k)
            conv_e_nf(i)  =-conv_e_non(k)+dp(ii)*conv_m_non(k)
         ENDDO
      ELSE
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!       
            ar1=ar(ii)
            ark1=ke(ii)*ar1
            are1=dp(ii)*ar1
!       
            ar2=ar(kk)
            ark2=ke(kk)*ar2
            are2=dp(kk)*ar2
!
            conv_k1=MIN(flux_nf(i1),0.d0)*ark2+MAX(flux_nf(i1),0.d0)*ark1
            conv_e1=MIN(flux_nf(i1),0.d0)*are2+MAX(flux_nf(i1),0.d0)*are1
            conv_m1=MIN(flux_nf(i1),0.d0)*ar2 +MAX(flux_nf(i1),0.d0)*ar1
            conv_k_nf(i0)=conv_k1-ke(ii)*conv_m1
            conv_e_nf(i0)=conv_e1-dp(ii)*conv_m1
         ENDDO
!
         nv=-1
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k)
!       
            ar1=ar(kk)
            ark1=ke(kk)*ar1
            are1=dp(kk)*ar1
!       
            ar2=ar(ii)
            ark2=ke(ii)*ar2
            are2=dp(ii)*ar2
!
            conv_k1=MIN(flux_nf(k),0.d0)*ark2+MAX(flux_nf(k),0.d0)*ark1
            conv_e1=MIN(flux_nf(k),0.d0)*are2+MAX(flux_nf(k),0.d0)*are1
            conv_m1=MIN(flux_nf(k),0.d0)*ar2 +MAX(flux_nf(k),0.d0)*ar1
            conv_k_nf(i)=-conv_k1+ke(ii)*conv_m1
            conv_e_nf(i)=-conv_e1+dp(ii)*conv_m1
         ENDDO
      ENDIF
!      
!.....Cells mcc
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
!       
         ar1=ar(ii)
         ark1=ke(ii)*ar1
         are1=dp(ii)*ar1
!
         conv_k1=MIN(flux_nf(i1),0.d0)*ark2+MAX(flux_nf(i1),0.d0)*ark1
         conv_e1=MIN(flux_nf(i1),0.d0)*are2+MAX(flux_nf(i1),0.d0)*are1
         conv_m1=MIN(flux_nf(i1),0.d0)*ar2 +MAX(flux_nf(i1),0.d0)*ar1
         conv_k_nf(i0)=conv_k1-ke(ii)*conv_m1
         conv_e_nf(i0)=conv_e1-dp(ii)*conv_m1
      ENDDO
!
!.....Cells inl
!                          
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
!       
         ar1=ar(ii)
         ark1=ke(ii)*ar1
         are1=dp(ii)*ar1
!
         ar2=ab(k)*rb(k)
         ark2=keb(k)*ar2
         are2=dpb(k)*ar2
!
         conv_k1=MIN(flux_nf(i1),0.d0)*ark2+MAX(flux_nf(i1),0.d0)*ark1
         conv_e1=MIN(flux_nf(i1),0.d0)*are2+MAX(flux_nf(i1),0.d0)*are1
         conv_m1=MIN(flux_nf(i1),0.d0)*ar2 +MAX(flux_nf(i1),0.d0)*ar1
         conv_k_nf(i0)=conv_k1-ke(ii)*conv_m1
         conv_e_nf(i0)=conv_e1-dp(ii)*conv_m1
      ENDDO
!      
!.....Cells out
!
      nv=3
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
!       
         ar1=ar(ii)
         ark1=ke(ii)*ar1
         are1=dp(ii)*ar1
!
         conv_k1=MIN(flux_nf(i1),0.d0)*ark2+MAX(flux_nf(i1),0.d0)*ark1
         conv_e1=MIN(flux_nf(i1),0.d0)*are2+MAX(flux_nf(i1),0.d0)*are1
         conv_m1=MIN(flux_nf(i1),0.d0)*ar2 +MAX(flux_nf(i1),0.d0)*ar1
         conv_k_nf(i0)=conv_k1-ke(ii)*conv_m1
         conv_e_nf(i0)=conv_e1-dp(ii)*conv_m1
      ENDDO
!
      CALL sum_nf(0,0,              &
                  conv_k_nf,conv_k, &
                  conv_e_nf,conv_e)
!
      END SUBROUTINE turb_ke_convection

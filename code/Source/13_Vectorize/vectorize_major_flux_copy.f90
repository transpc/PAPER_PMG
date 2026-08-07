      SUBROUTINE vectorize_major_flux_copy
!
!     This routine vectorizes volume fluxes.
!
      USE Zflux        ,ONLY: fluxvol_l,fluxvol_g,fluxvol_d
      USE Zvec_index   ,ONLY: left_nf,jneigh_nf, &
                              right_non,kneigh_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
!
      IMPLICIT NONE
!
      INTEGER :: nf_number,istart,isize,i1
      INTEGER :: i,jj,ii,kk,jk
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          flux_l_nf(i1)=fluxvol_l(jj,ii)
          flux_g_nf(i1)=fluxvol_g(jj,ii)
          flux_d_nf(i1)=fluxvol_g(jj,ii)
       ENDDO
      nf_number=2
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          flux_l_nf(i1)=fluxvol_l(jj,ii)
          flux_g_nf(i1)=fluxvol_g(jj,ii)
          flux_d_nf(i1)=fluxvol_d(jj,ii)
       ENDDO
      nf_number=3
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          flux_l_nf(i1)=fluxvol_l(jj,ii)
          flux_g_nf(i1)=fluxvol_g(jj,ii)
          flux_d_nf(i1)=fluxvol_d(jj,ii)
       ENDDO
      nf_number=2
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          flux_l_nf(i1)=fluxvol_l(jj,ii)
          flux_g_nf(i1)=fluxvol_g(jj,ii)
          flux_d_nf(i1)=fluxvol_d(jj,ii)
       ENDDO
!     make fluxvol symetric
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
          i1=istart+i
         kk=right_non(i)
         jk=kneigh_non(i)
         fluxvol_l(jk,kk)=-flux_l_nf(i1)
         fluxvol_g(jk,kk)=-flux_g_nf(i1)
         fluxvol_d(jk,kk)=-flux_d_nf(i1)
      ENDDO
!
      RETURN
      END SUBROUTINE vectorize_major_flux_copy
!-----------------------------------------------------------------------------
      SUBROUTINE vectorize_major_flux_store
!
!     This routine vectorizes volume fluxes.
!
      USE Zflux        ,ONLY: fluxvol_l,fluxvol_g,fluxvol_d
      USE Zvec_index   ,ONLY: left_nf,jneigh_nf, &
                              right_non,kneigh_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
!
      IMPLICIT NONE
!
      INTEGER :: nf_number,istart,isize,i1
      INTEGER :: i,jj,ii,kk,jk
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          fluxvol_l(jj,ii)=flux_l_nf(i1)
          fluxvol_g(jj,ii)=flux_g_nf(i1)
          fluxvol_g(jj,ii)=flux_d_nf(i1)
       ENDDO
      nf_number=2
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          fluxvol_l(jj,ii)=flux_l_nf(i1)
          fluxvol_g(jj,ii)=flux_g_nf(i1)
          fluxvol_g(jj,ii)=flux_d_nf(i1)
       ENDDO
      nf_number=3
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          fluxvol_l(jj,ii)=flux_l_nf(i1)
          fluxvol_g(jj,ii)=flux_g_nf(i1)
          fluxvol_g(jj,ii)=flux_d_nf(i1)
       ENDDO
      nf_number=2
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
       DO i=1,isize
          i1=istart+i
          ii=left_nf(i1)
          jj=jneigh_nf(i1)
          fluxvol_l(jj,ii)=flux_l_nf(i1)
          fluxvol_g(jj,ii)=flux_g_nf(i1)
          fluxvol_g(jj,ii)=flux_d_nf(i1)
       ENDDO
!     make fluxvol symetric
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
          i1=istart+i
         kk=right_non(i)
         jk=kneigh_non(i)
         fluxvol_l(jk,kk)=-flux_l_nf(i1)
         fluxvol_g(jk,kk)=-flux_g_nf(i1)
         fluxvol_d(jk,kk)=-flux_d_nf(i1)
      ENDDO
!
      RETURN
      END SUBROUTINE vectorize_major_flux_store

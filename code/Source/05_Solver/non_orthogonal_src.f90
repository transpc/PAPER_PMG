!
      SUBROUTINE non_orthogonal_src(source2,poiss_non_i,poiss_non_k)
!
!     This routine calculate additional source to account the non-orthogonal grid
!     This source will be used for the second pressure correction
!
      USE Zinterface
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_nonk,nf_non
      USE Znum_cell     , ONLY: istart_nf, &
                                nf_number_nb,lens,                           &
                                right_nb_k,istart_nfs,nf_number_id,istart_nf
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zbc_index     , ONLY: npb
      USE Zgradoption   , ONLY: non_orth
      USE Zpress        , ONLY: dpdx
      USE Zvec_geo      , ONLY: xn_nf,djia_nf,dji_x_nf,         &
                                xloc_m_non_i,xloc_m_non_k,      &
                                xn_non_k,djia_non_k,dji_x_non_k
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(nf_non) :: poiss_non_i
      REAL(8),DIMENSION(nf_nonk) :: poiss_non_k
!.....Output
      REAL(8) :: source2(ncell_fluid)
!.....Local variables
      INTEGER :: i,k,ii,kk
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      REAL(8) :: grdPx1,grdPx2,grdPx3
      REAL(8) :: dpi,dpj,dp,dpx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non) :: source_non
!
!.....Build summation info for non,mcc,inl
!
      nf_number_nb=0
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      lens         =istart_nfs(0) +nf_non
      IF(non_orth.eq.1)THEN
         IF(ndim.eq.2)THEN
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               IF(npb(ii).eq.0) THEN
                  kk=right_non(i)
!
!.................Arithmetic average on cell face between cell centers
!               
                  grdPx1=0.5d0*(dpdx(ii,1)+dpdx(kk,1))
                  grdPx2=0.5d0*(dpdx(ii,2)+dpdx(kk,2))
!
!.................Amplitude of gradient(P) in cell normal direction to get !'dp'
!
                  dpx=grdPx1*xn_nf(i1,1)+grdPx2*xn_nf(i1,2)
                  dp=dpx*djia_nf(i1)
!
!.................Amplitude of gradient(P) between cell centers to get !'dpx'
!
                  dpx=grdPx1*dji_x_nf(i1,1)+grdPx2*dji_x_nf(i1,2)
!
!.................Define deviation between 'dp' and 'dpx'
!
                  dp=dpx-dp
!
!.................Correction with the deviation
!              
                  source_non(i0)=poiss_non_i(i)*dp
               ELSE 
                  source_non(i0)=0.d0
               ENDIF
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               IF(npb(ii).eq.0) THEN
                  kk=left_nf(k)
                  grdPx1=0.5d0*(dpdx(kk,1)+dpdx(ii,1))
                  grdPx2=0.5d0*(dpdx(kk,2)+dpdx(ii,2))
                  dpx=grdPx1*xn_non_k(i,1)+grdPx2*xn_non_k(i,2)
                  dp=dpx*djia_non_k(i)
                  dpx=grdPx1*dji_x_non_k(i,1)+grdPx2*dji_x_non_k(i,2)
                  dp=dpx-dp
                  source_non(i)=-poiss_non_k(i)*dp
               ELSE 
                  source_non(i)=0.d0
               ENDIF
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
               IF(npb(ii).eq.0) THEN
                  kk=right_non(i)
                  grdPx1=0.5d0*(dpdx(ii,1)+dpdx(kk,1))
                  grdPx2=0.5d0*(dpdx(ii,2)+dpdx(kk,2))
                  grdPx3=0.5d0*(dpdx(ii,3)+dpdx(kk,3))
                  dpx=grdPx1*xn_nf(i1,1)+grdPx2*xn_nf(i1,2)+grdPx3*xn_nf(i1,3)
                  dp=dpx*djia_nf(i1)
                  dpx=grdPx1*dji_x_nf(i1,1)+grdPx2*dji_x_nf(i1,2)+grdPx3*dji_x_nf(i1,3)
                  dp=dpx-dp
                  source_non(i0)=poiss_non_i(i)*dp
               ELSE 
                  source_non(i0)=0.d0
               ENDIF
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               IF(npb(ii).eq.0) THEN
                  kk=left_nf(k)
                  grdPx1=0.5d0*(dpdx(kk,1)+dpdx(ii,1))
                  grdPx2=0.5d0*(dpdx(kk,2)+dpdx(ii,2))
                  grdPx3=0.5d0*(dpdx(kk,3)+dpdx(ii,3))
                  dpx=grdPx1*xn_non_k(i,1)+grdPx2*xn_non_k(i,2)+grdPx3*xn_non_k(i,3)
                  dp=dpx*djia_non_k(i)
                  dpx=grdPx1*dji_x_non_k(i,1)+grdPx2*dji_x_non_k(i,2)+grdPx3*dji_x_non_k(i,3)
                  dp=dpx-dp
                  source_non(i)=-poiss_non_k(i)*dp
               ELSE 
                  source_non(i)=0.d0
               ENDIF
            ENDDO
         ENDIF
      ELSEIF(non_orth.eq.2)THEN
         IF(ndim.eq.2)THEN
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               IF(npb(ii).eq.0) THEN
                  kk=right_non(i)
                  dpi=dpdx(ii,1)*xloc_m_non_i(i,1)+dpdx(ii,2)*xloc_m_non_i(i,2)
                  dpj=dpdx(kk,1)*xloc_m_non_k(i,1)+dpdx(kk,2)*xloc_m_non_k(i,2)
                  dp=dpi-dpj
                  source_non(i0)=poiss_non_i(i)*dp
               ELSE 
                  source_non(i0)=0.d0
               ENDIF
            ENDDO
!           
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
                  k=right_nb_k(i)
               ii=right_non(k)
               IF(npb(ii).eq.0) THEN
                  kk=left_nf(k)
                  dpi=dpdx(kk,1)*xloc_m_non_i(k,1)+dpdx(kk,2)*xloc_m_non_i(k,2)
                  dpj=dpdx(ii,1)*xloc_m_non_k(k,1)+dpdx(ii,2)*xloc_m_non_k(k,2)
                  dp=dpi-dpj
                  source_non(i)=-poiss_non_k(i)*dp
               ELSE 
                  source_non(i)=0.d0
               ENDIF
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
               IF(npb(ii).eq.0) THEN
                  kk=right_non(i)
                  dpi=dpdx(ii,1)*xloc_m_non_i(i,1)+dpdx(ii,2)*xloc_m_non_i(i,2)+dpdx(ii,3)*xloc_m_non_i(i,3)
                  dpj=dpdx(kk,1)*xloc_m_non_k(i,1)+dpdx(kk,2)*xloc_m_non_k(i,2)+dpdx(kk,3)*xloc_m_non_k(i,3)
                  dp=dpi-dpj
                  source_non(i0)=poiss_non_i(i)*dp
               ELSE 
                  source_non(i0)=0.d0
               ENDIF
            ENDDO
!           
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               IF(npb(ii).eq.0) THEN
                  kk=left_nf(k)
                  dpi=dpdx(kk,1)*xloc_m_non_i(k,1)+dpdx(kk,2)*xloc_m_non_i(k,2)+dpdx(kk,3)*xloc_m_non_i(k,3)
                  dpj=dpdx(ii,1)*xloc_m_non_k(k,1)+dpdx(ii,2)*xloc_m_non_k(k,2)+dpdx(ii,3)*xloc_m_non_k(k,3)
                  dp=dpi-dpj
                  source_non(i)=-poiss_non_k(i)*dp
               ELSE 
                  source_non(i)=0.d0
               ENDIF
            ENDDO
         ENDIF
      ELSE
         STOP '"non_orth" should be 1 or 2'
      ENDIF
!
      CALL sum_nf(0,0,                &
                  source_non,source2)
!
      END SUBROUTINE non_orthogonal_src

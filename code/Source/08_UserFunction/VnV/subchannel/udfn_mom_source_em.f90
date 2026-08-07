!
      SUBROUTINE udfn_mom_source_em
!
!     Modifies the momentum source terms at the free surface cells
!
      USE Zinterface
      USE VOL_DATA   , ONLY: cell
      USE Zzone      , ONLY: ncell_fluid
      USE Zparam     , ONLY: ndim
      USE Zvec_param , ONLY: nf_nonk,nf_non
      USE Znum_cell  , ONLY: istart_nf
      USE Zvec_index , ONLY: left_nf,right_non
      USE Znum_cell  , ONLY: istart_nf,istart_nb1,                 &
                             ia_nb,icell_nb,right_nb_k,            &
                             nf_number_nb,lens,                    &
                             right_nb_k,istart_nfs,nf_number_id
      USE Zvector    , ONLY: vl_o
      USE Zporous    , ONLY: s_gapij_non_i,s_gapij_non_k
      USE Zporous    , ONLY: ftm,beta
      USE Zm_src     , ONLY: src_liq
      USE Zvec_geo   , ONLY: xn_nf,saa_nf
!
      IMPLICIT NONE
!            
!.....Local variables
      INTEGER :: i,ix,k
      INTEGER :: ii,kk,nb
      INTEGER :: nv,nf_number,istart,len,istart0,istart1,i0,i1
      REAL(8) :: u_i,u_j,s_gap,g_gap,sum_turb_mixing,a_subchannel
      REAL(8) :: turb_mixing_non
      REAL(8) :: s1_s
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: turb_mixing,a_sub
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non) :: g_gap_non
!
!.....Build summation info for non
!
      nf_number_nb=1
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      lens         =istart_nfs(0)+nf_non
!
!.....Cells non
!
      ix=ndim
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
         IF(abs(xn_nf(i1,1)).eq.1.d0 .or.abs(xn_nf(i1,2)).eq.1.d0) THEN
            u_i=vl_o(ii,ix)
            u_j=vl_o(kk,ix)
            turb_mixing_non=(u_i-u_j)
            g_gap=(cell%rhol(ii)*u_i+cell%rhol(kk)*u_j)*0.5d0 
            s_gap=s_gapij_non_i(i)
            g_gap_non(i0)=ftm*beta*s_gap*g_gap*turb_mixing_non
         ELSE
            g_gap_non(i0)=0.d0
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
         IF(abs(xn_nf(k,1)).eq.1.d0 .or.abs(xn_nf(k,2)).eq.1.d0) THEN
            u_i=vl_o(kk,ix)
            u_j=vl_o(ii,ix)
            turb_mixing_non=(u_i-u_j)
            g_gap=(cell%rhol(kk)*u_i+cell%rhol(ii)*u_j)*0.5d0 
            s_gap=s_gapij_non_k(k)
            g_gap_non(i)=-ftm*beta*s_gap*g_gap*turb_mixing_non
         ELSE
            g_gap_non(i)=0.d0
         ENDIF
      ENDDO
!
      CALL sum_nf(0,0,                   &
                  g_gap_non,turb_mixing)
!
!.....Compute a_subchannel into a_sub(i)
!
      DO i=1,ncell_fluid
         a_sub(i)=1.d0
      ENDDO
!
!.....Cells non_k
!
      nf_number=-1
      len   =istart_nb1(2,nf_number)
      DO nb=1,len
         ii=icell_nb(nb)
!DIR$ NOVECTOR
         DO i=ia_nb(nb),ia_nb(nb+1)-1
            k=right_nb_k(i)
!           IF(abs(xn_nf(k,3)).eq.1.d0) s1_s=saa_nf(k)
!           a_sub(ii)=s1_s
            IF(abs(xn_nf(k,3)).eq.1.d0) a_sub(ii)=saa_nf(k)
         ENDDO
      ENDDO
!
!.....Cells non_k
!
      nf_number=0
      istart =istart_nf(1,nf_number)
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
         ii=icell_nb(i1)
         s1_s=a_sub(ii)
!DIR$ NOVECTOR
         DO i=ia_nb(i1),ia_nb(i1+1)-1
            i0=istart+i
!           IF(abs(xn_nf(i0,3)).eq.1.d0) s1_s=saa_nf(i0)
!           a_sub(ii)=s1_s
            IF(abs(xn_nf(i0,3)).eq.1.d0) a_sub(ii)=saa_nf(i0)
         ENDDO
      ENDDO
!
      DO i=1,ncell_fluid
         sum_turb_mixing=turb_mixing(i)
         a_subchannel=a_sub(i)
         sum_turb_mixing=sum_turb_mixing/a_subchannel
         src_liq(i,ix)=src_liq(i,ix)-sum_turb_mixing
      ENDDO
!
      END SUBROUTINE udfn_mom_source_em

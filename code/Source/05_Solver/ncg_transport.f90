!
      SUBROUTINE ncg_transport
!
      USE Zncg         , ONLY: n_ncg_sp,imp_ncg,qn_cell
!
      IMPLICIT NONE
!
      IF(n_ncg_sp.eq.1)THEN
         qn_cell(:,1)=1.0d0
      ELSE
!
         IF(imp_ncg.eq.0)THEN
!
!...........Explicit transport
!
            CALL ncg_transport_exp
!
         ELSE
!
!...........Implicit transport
!
            CALL ncg_transport_imp
!
         ENDIF
      ENDIF
!
!.....Normalize each ncg fractions
!
      CALL norm_ncg
!
!.....CAlculate ncg properties
!
      CALL ncg_cell
!      
      RETURN
      END SUBROUTINE ncg_transport
!
      SUBROUTINE ncg_transport_exp
!
!     This routine calculates ncg transport explicitly
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zconst2      , ONLY: dt
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf, &
                                nf_number_nb,lens,nf_number_id,istart_nfs,     &
                               right_nb_k
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Z2nd_order   , ONLY: ncg_conv_2nd
      USE Zare         , ONLY: ar_gas
      USE Zb_condition , ONLY: alphab_gas,rhob_gas,qualab
      USE Zcoord3      , ONLY: volp
      USE Zncg         , ONLY: n_ncg_sp,qn_cell,qn_cell_o,qn_nvin,ncg_diff
      USE Ztimecon     , ONLY: alpha_min
      USE Zvec_param   , ONLY: nf_nonk,nf_non,nf_inl
      USE Zvec_major   , ONLY: flux_g_nf
      USE Zvec_geo     , ONLY: fac1_non,fac_non,sap_nf
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i,k,nc
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: arx2
      REAL(8) :: ardi
      REAL(8) :: avt
      REAL(8) :: ncg_conv,ncg_diff1,ncg_diff2
      REAL(8) :: temp_i,temp_k
!.....Local arrays
      REAL(8),DIMENSION(ncell_fp) :: arx,qn
      REAL(8),DIMENSION(ncell_fluid) :: ncg_diff1v,ncg_diff2v
      REAL(8),DIMENSION(ncell_fluid) :: ncg_convv
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: temp_non,ardi_non1,ardi_non2
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_inl) :: temp_nf
!
!.....Save convective variables
!
      DO i=1,ncell_fluid
         arx(i)=ar_gas(i)*cell%quala(i)
      ENDDO
!
!.....Communicate convective variables
!
      IF(np.gt.1) CALL communicate_1d(arx,        &
                                      ar_gas,     &
                                      cell%quala, &
                                      cell%mdiff)
!
!.....Calculate convective flux through cell face
!
      DO nc=1,n_ncg_sp
!
         DO i=1,ncell_fluid
            qn(i)=qn_cell_o(i,nc)
         ENDDO
!
         IF(np.gt.1) CALL communicate_1d(qn)
!
         IF(ncg_conv_2nd.gt.0) THEN
            nf_number=0
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               IF(flux_g_nf(i1).lt.0) THEN
                  temp_non(i)=qn(kk)
               ELSEIF(flux_g_nf(i1).gt.0) THEN
                  temp_non(i)=qn(ii)
               END IF
            ENDDO
!
!...........Move the multiply by flux_g_nfaway from mult_ncg_2nd_conv to enhance simplify the code
!
            CALL mult_ncg_2nd_conv1(qn,temp_non)
!
!...........Build summation info for non,inl
!
            nf_number_nb=1
            nf_number_id(-1)=-1
            nf_number_id(0)=0
            nf_number_id(1)=2
            istart_nfs(-1)=0
            istart_nfs(0)=istart_nfs(-1)+nf_nonk
            istart_nfs(1)=istart_nfs(0) +nf_non
            lens         =istart_nfs(1) +nf_inl
!
!...........Computing Cell
!
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
               IF(flux_g_nf(i1).lt.0) THEN
                  temp_i=(temp_non(i)-qn(ii))*flux_g_nf(i1)
                  temp_nf(i0)=temp_i*arx(kk)
               ELSE IF(flux_g_nf(i1).gt.0) THEN
                  temp_i=(temp_non(i)-qn(ii))*flux_g_nf(i1)
                  temp_nf(i0)=temp_i*arx(ii)
               ELSE
                  temp_nf(i0)=0.d0
               END IF
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               kk=left_nf(k)
               IF(flux_g_nf(k).lt.0) THEN
                  temp_k=(temp_non(k)-qn(ii))*flux_g_nf(k)
                  temp_nf(i)=-temp_k*arx(ii)
               ELSEIF(flux_g_nf(k).gt.0) THEN
                  temp_k=(temp_non(k)-qn(ii))*flux_g_nf(k)
                  temp_nf(i)=-temp_k*arx(kk)
               ELSE
                  temp_nf(i)=0.d0
               END IF
            ENDDO
         ELSE
!
!...........Build summation info for non,inl
!
            nf_number_nb=1
            nf_number_id(-1)=-1
            nf_number_id(0)=0
            nf_number_id(1)=2
            istart_nfs(-1)=0
            istart_nfs(0)=istart_nfs(-1)+nf_nonk
            istart_nfs(1)=istart_nfs(0) +nf_non
            lens         =istart_nfs(1) +nf_inl
!
!...........Computing Cell
!
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
               IF(flux_g_nf(i1).lt.0) THEN
                  temp_nf(i0)=(qn(kk)-qn(ii))*flux_g_nf(i1)*arx(kk)
               ELSE
                  temp_nf(i0)=0.d0
               END IF
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               kk=left_nf(k)
               IF(flux_g_nf(k).le.0) THEN
                  temp_nf(i)=0.d0
               ELSEIF(flux_g_nf(k).gt.0) THEN
                  temp_nf(i)=-(qn(kk)-qn(ii))*flux_g_nf(k)*arx(kk)
               END IF
            ENDDO
         ENDIF
!
!........Inlet
!                          
         nv=1
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
            IF(flux_g_nf(i1).lt.0.d0)THEN
               temp_i=qn_nvin(k,nc)*flux_g_nf(i1)
               arx2=alphab_gas(k)*rhob_gas(k)*qualab(k)
               temp_nf(i0)=(temp_i-qn(ii)*flux_g_nf(i1))*arx2
            ELSE
               temp_nf(i0)=0.d0
            ENDIF
         ENDDO
!
      CALL sum_nf(0,0,               &
                  temp_nf,ncg_convv)
!
         IF(ncg_diff.gt.0) THEN
            nf_number_nb=0
            nf_number_id(0)=0
            istart_nfs(0)=0
            lens         =istart_nfs(0)+nf_non
!
            nv=0
            nf_number=nf_number_id(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               ardi=(fac1_non(i)*ar_gas(ii)*cell%mdiff(ii)+fac_non(i)*ar_gas(kk)*cell%mdiff(kk))*sap_nf(i1)
               ardi_non1(i)=ardi*(cell%quala(kk)*qn(kk)-cell%quala(ii)*qn(ii))
               ardi_non2(i)=ardi*(cell%quala(kk)-cell%quala(ii))
            ENDDO 
!
            CALL sum_nf(0,-1,                 &
                        ardi_non1,ncg_diff1v, &
                        ardi_non2,ncg_diff2v)
!
         ENDIF
!
         DO i=1,ncell_fluid
            ncg_conv=ncg_convv(i)
            ncg_diff1=ncg_diff1v(i)
            ncg_diff2=ncg_diff2v(i)
            IF(cell%alphag(i).le.alpha_min .or. cell%alphag_o(i).le.alpha_min .or. cell%quala(i).le.alpha_min)THEN
               qn_cell(i,nc)=qn(i)
            ELSE
               avt=dt/(volp(i)*arx(i))
               qn_cell(i,nc)=qn(i)-ncg_conv*avt
!
!..............Mass diffusion
!
               IF(ncg_diff.gt.0) qn_cell(i,nc)=qn_cell(i,nc)+(ncg_diff1-qn(i)*ncg_diff2)*avt
!
            ENDIF
            qn_cell(i,nc)=MAX(qn_cell(i,nc),0.0d0)
            qn_cell(i,nc)=MIN(qn_cell(i,nc),1.0d0)
         ENDDO
!
      ENDDO
!
      END SUBROUTINE ncg_transport_exp
!
      SUBROUTINE norm_ncg
!
      USE Zncg         , ONLY: n_ncg_sp,qn_cell
      USE Zzone        , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,nc
      REAL(8) :: sum_nc,eps
!
      DO i=1,ncell_fluid
         sum_nc=0.0d0
         DO nc=1,n_ncg_sp
            sum_nc=sum_nc+qn_cell(i,nc)
         ENDDO
         eps=sum_nc-1.0d0
         IF(ABS(eps).gt.0.0d0)THEN
            DO nc=1,n_ncg_sp
               qn_cell(i,nc)=qn_cell(i,nc)-eps*qn_cell(i,nc)/sum_nc
            ENDDO
         ENDIF
      ENDDO
!
      END SUBROUTINE norm_ncg
!

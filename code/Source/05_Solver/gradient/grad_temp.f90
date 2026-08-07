!
      SUBROUTINE grad_temp(sl,dtldx,tinl, &
                           sg,dtgdx,ting)
!
!     This routine calculate gradient of a scalar using the Green-Gauss theorem
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim,nb_max
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_nonk
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf,right_nb_k, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zcoord3      , ONLY: volr
      USE Zb_condition , ONLY: twall
      USE Zvec_geo     , ONLY: sv_nf,fac1_non,fac_non
!
      IMPLICIT NONE
      INCLUDE '../../10_LinkToMARS/c3com.h'
!
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: sl,sg
      REAL(8),DIMENSION(nb_max) :: tinl,ting
!.....Output
      REAL(8),DIMENSION(ncell_fp,ndim) :: dtldx,dtgdx
!.....Local variables
      INTEGER :: i,k,ii,kk
      INTEGER :: nv,nf_number,istart,len,istart2,i1,i2,istart0,i0
      REAL(8) :: c1,c2,fie
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_mcc+nf_inl+nf_out+nf_adw+nf_fsw+nf_ctw+nf_chw+nf_sym,ndim) :: fiel_nf,fieg_nf
!
!.....Build summation info for all nf
!
      nf_number_nb=8
      nf_number_id(-1)=-1      
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      nf_number_id(4)=4
      nf_number_id(5)=5
      nf_number_id(6)=6
      nf_number_id(7)=7
      nf_number_id(8)=8
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      istart_nfs(4)=istart_nfs(3)+nf_out
      istart_nfs(5)=istart_nfs(4)+nf_adw
      istart_nfs(6)=istart_nfs(5)+nf_fsw
      istart_nfs(7)=istart_nfs(6)+nf_ctw
      istart_nfs(8)=istart_nfs(7)+nf_chw
      lens         =istart_nfs(8)+nf_sym  
!
      IF(ndim.eq.2) THEN
!          
         nv=-1
         nf_number=nf_number_id(nv)
         len=istart_nf(2,nf_number)
         DO i=1,len  
            k=right_nb_k(i)
            !ii=right_non(k) !the result is different from symmetric sum.
            !kk=left_nf(k)
            kk=right_non(k)
            ii=left_nf(k)            
            c1=fac1_non(k)
            c2=fac_non(k)
            fie=c1*sl(ii)+c2*sl(kk)
            fiel_nf(i,1)=-fie*sv_nf(k,1)
            fiel_nf(i,2)=-fie*sv_nf(k,2)
            fie=c1*sg(ii)+c2*sg(kk)
            fieg_nf(i,1)=-fie*sv_nf(k,1)
            fieg_nf(i,2)=-fie*sv_nf(k,2)
         ENDDO                      
!          
         DO nv=0,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
!
!...........Cells non
!
            IF(nf_number.eq.0) THEN
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  c1=fac1_non(i)
                  c2=fac_non(i)
                  fie=c1*sl(ii)+c2*sl(kk)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fie=c1*sg(ii)+c2*sg(kk)
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
!
!...........Cells inl
!
            ELSEIF(nf_number.eq.2) THEN
               istart2=istart_nbcon_nf(nf_number)
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  i2=istart2+i
                  ii=left_nf(i1)
                  k=nbcon_nf(i2)
                  fie=tinl(k)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fie=ting(k)
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
!
!...........Cells ctw
!
            ELSEIF(nf_number.eq.6) THEN
               istart2=istart_nbcon_nf(nf_number)
               DO i=1,len 
                  i0=istart0+i
                  i1=istart+i
                  i2=istart2+i
                  ii=left_nf(i1)
                  k=-nbcon_nf(i2)
                  fie=twall(k)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
!
!...........The rest
!
            ELSE
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  fie=sl(ii)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fie=sg(ii)
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
            ENDIF
         ENDDO
      ELSE

         nv=-1
         nf_number=nf_number_id(nv)
         len=istart_nf(2,nf_number)
         DO i=1,len  
            k=right_nb_k(i)
            !ii=right_non(k) !the result is different from symmetric sum.
            !kk=left_nf(k)
            kk=right_non(k) 
            ii=left_nf(k)            
            c1=fac1_non(k)
            c2=fac_non(k)            
            fie=c1*sl(ii)+c2*sl(kk)
            fiel_nf(i,1)=-fie*sv_nf(k,1)
            fiel_nf(i,2)=-fie*sv_nf(k,2)
            fiel_nf(i,3)=-fie*sv_nf(k,3)
            fie=c1*sg(ii)+c2*sg(kk)
            fieg_nf(i,1)=-fie*sv_nf(k,1)
            fieg_nf(i,2)=-fie*sv_nf(k,2)
            fieg_nf(i,3)=-fie*sv_nf(k,3)
         ENDDO           
          
         DO nv=0,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
!
!...........Cells non
!
            IF(nf_number.eq.0) THEN
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  c1=fac1_non(i)
                  c2=fac_non(i)
                  fie=c1*sl(ii)+c2*sl(kk)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fiel_nf(i0,3)=fie*sv_nf(i1,3)
                  fie=c1*sg(ii)+c2*sg(kk)
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
                  fieg_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
!
!...........Cells inl
!
            ELSEIF(nf_number.eq.2) THEN
               istart2=istart_nbcon_nf(nf_number)
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  i2=istart2+i
                  ii=left_nf(i1)
                  k=nbcon_nf(i2)
                  fie=tinl(k)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fiel_nf(i0,3)=fie*sv_nf(i1,3)
                  fie=ting(k)
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
                  fieg_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
!
!........Cells ctw
!
            ELSEIF(nf_number.eq.6) THEN
               istart2=istart_nbcon_nf(nf_number)
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  i2=istart2+i
                  ii=left_nf(i1)
                  k=-nbcon_nf(i2)
                  fie=twall(k)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fiel_nf(i0,3)=fie*sv_nf(i1,3)              
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
                  fieg_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
!
!...........The rest
!
            ELSE
               DO i=1,len  
                  i0=istart0+i 
                  i1=istart+i
                  ii=left_nf(i1)
                  fie=sl(ii)
                  fiel_nf(i0,1)=fie*sv_nf(i1,1)
                  fiel_nf(i0,2)=fie*sv_nf(i1,2)
                  fiel_nf(i0,3)=fie*sv_nf(i1,3)
                  fie=sg(ii)
                  fieg_nf(i0,1)=fie*sv_nf(i1,1)
                  fieg_nf(i0,2)=fie*sv_nf(i1,2)
                  fieg_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!  
      CALL sum_nf_ndim(0,0,ncell_fp, &
                       fiel_nf,dtldx,   &
                       fieg_nf,dtgdx)  
!
      IF(ndim.eq.2)THEN
         DO i=1,ncell_fluid
            dtldx(i,1)=dtldx(i,1)*volr(i)
            dtldx(i,2)=dtldx(i,2)*volr(i)
            dtgdx(i,1)=dtgdx(i,1)*volr(i)
            dtgdx(i,2)=dtgdx(i,2)*volr(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            dtldx(i,1)=dtldx(i,1)*volr(i)
            dtldx(i,2)=dtldx(i,2)*volr(i)
            dtldx(i,3)=dtldx(i,3)*volr(i)
            dtgdx(i,1)=dtgdx(i,1)*volr(i)
            dtgdx(i,2)=dtgdx(i,2)*volr(i)
            dtgdx(i,3)=dtgdx(i,3)*volr(i)
         ENDDO
      ENDIF
!
      END SUBROUTINE grad_temp

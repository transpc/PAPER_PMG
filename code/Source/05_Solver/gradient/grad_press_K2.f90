!
      SUBROUTINE grad_pressK2(s,dsdx,idg)
!
!     double-loop
!
!     This routine calculates the components of the gradient
!     vector of pressure at the cell center, using conservative
!     scheme based on the gauss theorem.
!
!     idg (only for MARS interfaces): 1=p, 2=dp, 3=ELSE 
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_nonk
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Znum_cell    , ONLY: istart_nf,right_nb_k, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zcoord3      , ONLY: volr
      USE Zgradoption  , ONLY: iter_grad
      USE c3com_cupid  , ONLY: i3invtbl
      USE Zvec_geo     , ONLY: sv_nf,dxfc_nf,   &
                               dxfc_non_k,      &
                               fac1_non,fac_non
      USE Zrv_model    , ONLY: rv_choke,rv_mcp,rv_valve      
!
      IMPLICIT NONE
!      
      INCLUDE '../../10_LinkToMARS/c3com.h'
!
!.....Input
      INTEGER :: idg
      REAL(8),DIMENSION(ncell_fp) :: s
!.....Output
      REAL(8),DIMENSION(ncell_fp,ndim) :: dsdx
!.....Local variables
      INTEGER :: i,idx
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,isize,i1,len,k,i0,istart0
      INTEGER :: itergx
      REAL(8) :: fie
      REAL(8) :: sk
      REAL(8) :: dp,p_i,p_k      
!.....Local arrays
      REAL(8),DIMENSION(ncell_fp,ndim) :: dgdx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_mcc+nf_inl+nf_out+nf_adw+nf_fsw+nf_ctw+nf_chw+nf_sym,ndim) :: fie_nf   
!
!.....Build summation info for non,inl
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
      IF(ndim.eq.2)THEN
!
!........Cells non
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            fie=fac1_non(i)*s(ii)+fac_non(i)*s(kk)
            fie_nf(i0,1)=fie*sv_nf(i1,1)
            fie_nf(i0,2)=fie*sv_nf(i1,2)
         ENDDO
!         
         nv=-1
         nf_number=nf_number_id(nv)
         len=istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            kk=right_non(k)
            ii=left_nf(k)            
            fie=fac1_non(k)*s(ii)+fac_non(k)*s(kk)
            fie_nf(i,1)=-fie*sv_nf(k,1)
            fie_nf(i,2)=-fie*sv_nf(k,2)     
         ENDDO
!
!........Cells mcc
!
         nv=1
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         IF(idg.eq.0)THEN
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3pa(1,idx)
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ELSEIF(idg.eq.1)THEN
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3delp(1,idx)
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)     
            ENDDO
         ELSEIF(idg.eq.3)THEN
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3rtp(1,idx,1) !al
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ELSEIF(idg.eq.4)THEN
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3rtp(1,idx,2) !ag                  
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ELSE
            DO i=1,isize
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               sk=s(ii)
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ENDIF
!
!........The rest
!
         DO nv=2,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               fie=s(ii)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2) 
            ENDDO
         ENDDO
      ELSE
!
!........Cells non
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            fie=fac1_non(i)*s(ii)+fac_non(i)*s(kk)
            fie_nf(i0,1)=fie*sv_nf(i1,1)
            fie_nf(i0,2)=fie*sv_nf(i1,2)
            fie_nf(i0,3)=fie*sv_nf(i1,3)
         ENDDO
!         
         nv=-1
         nf_number=nf_number_id(nv)
         len=istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            kk=right_non(k)
            ii=left_nf(k)            
            fie=fac1_non(k)*s(ii)+fac_non(k)*s(kk)
            fie_nf(i,1)=-fie*sv_nf(k,1)
            fie_nf(i,2)=-fie*sv_nf(k,2)
            fie_nf(i,3)=-fie*sv_nf(k,3)
         ENDDO            
!
!........Cells mcc
!
         nv=1
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         IF(idg.eq.0)THEN
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3pa(1,idx)
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSEIF(idg.eq.1)THEN
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3delp(1,idx)
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSEIF(idg.eq.3)THEN
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3rtp(1,idx,1) !al
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1) 
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSEIF(idg.eq.4)THEN
            DO i=1,isize
               i0=istart0+i  
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               sk=c3rtp(1,idx,2) !ag                  
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSE
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               sk=s(ii)
               fie=0.5d0*(s(ii)+sk)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ENDIF
!
!........The rest
!
         DO nv=2,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               fie=s(ii)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ENDDO
      ENDIF
!
!.....fluxBC: choke model, mcp model, valve model
!  
      IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_gradp(s,fie_nf)      
!
      CALL sum_nf_ndim(0,0,ncell_fp,fie_nf,dsdx) 
!
!.....Calculate gradient components at cell centers
!
      IF(ndim.eq.2)THEN
         DO i=1,ncell_fluid
            dgdx(i,1)=dsdx(i,1)*volr(i)
            dgdx(i,2)=dsdx(i,2)*volr(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            dgdx(i,1)=dsdx(i,1)*volr(i)
            dgdx(i,2)=dsdx(i,2)*volr(i)
            dgdx(i,3)=dsdx(i,3)*volr(i)
         ENDDO
      ENDIF
!
!
!.....MPI communication
!
      IF(np.gt.1) CALL communicate_2d(dgdx)
!
!.....Additional correction of the gradient for a mesh non-orthogonality
!
      DO itergx=1,iter_grad
         IF(itergx.eq.iter_grad) THEN
            IF(ndim.eq.2)THEN
               DO i=1,ncell_fp 
                  dsdx(i,1)=dgdx(i,1)
                  dsdx(i,2)=dgdx(i,2)
               ENDDO
            ELSE
               DO i=1,ncell_fp 
                  dsdx(i,1)=dgdx(i,1)
                  dsdx(i,2)=dgdx(i,2)
                  dsdx(i,3)=dgdx(i,3)
               ENDDO
            ENDIF
            EXIT
         ENDIF
!
         IF(ndim.eq.2)THEN
!
!........Cells non
!
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            DO i=1,isize
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               dp= dgdx(ii,1)*dxfc_nf(i1,1) &
                  +dgdx(ii,2)*dxfc_nf(i1,2)
               p_i=s(ii)+dp
               dp= dgdx(kk,1)*dxfc_non_k(i,1) &
                  +dgdx(kk,2)*dxfc_non_k(i,2)
               p_k=s(kk)+dp
               fie=0.5d0*(p_i+p_k)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len=istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               kk=right_non(k)
               ii=left_nf(k)               
               dp= dgdx(ii,1)*dxfc_nf(k,1) &
                  +dgdx(ii,2)*dxfc_nf(k,2)
               p_i=s(ii)+dp
               dp= dgdx(kk,1)*dxfc_non_k(k,1) &
                  +dgdx(kk,2)*dxfc_non_k(k,2)
               p_k=s(kk)+dp
               fie=0.5d0*(p_i+p_k)
               fie_nf(i,1)=-fie*sv_nf(k,1)
               fie_nf(i,2)=-fie*sv_nf(k,2)  
            ENDDO            
!
!........Cells mcc
!
            nv=1
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            IF(idg.eq.0)THEN
               DO i=1,isize
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  idx=i3invtbl(i)
                  sk=c3pa(1,idx)
                  fie=0.5d0*(s(ii)+sk)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)    
               ENDDO
            ELSEIF(idg.eq.1)THEN
               DO i=1,isize
                  i0=istart0+i 
                  i1=istart+i
                  ii=left_nf(i1)
                  idx=i3invtbl(i)
                  sk=c3delp(1,idx)
                  fie=0.5d0*(s(ii)+sk)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
            ELSE
               DO i=1,isize
                  i0=istart0+i  
                  i1=istart+i
                  ii=left_nf(i1)
                  sk=s(ii)
                  fie=0.5d0*(s(ii)+sk)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
            ENDIF
!
!........Cells out
!
            nv=3
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               dp= dgdx(ii,1)*dxfc_nf(i1,1) &
                  +dgdx(ii,2)*dxfc_nf(i1,2)
               fie=s(ii)+dp
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
!
!...........The rest
!
            DO nv=2,8
               nf_number=nf_number_id(nv)
               IF(nf_number.eq.3) cycle
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               isize =istart_nf(2,nf_number)
               DO i=1,isize
                  i0=istart0+i   
                  i1=istart+i
                  ii=left_nf(i1)
                  fie=s(ii)
                  fie_nf(i0,1)=fie*sv_nf(i1,1) 
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
            ENDDO
         ELSE
!
!...........Cells non
!
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
!
!...........Computing cell
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               dp= dgdx(ii,1)*dxfc_nf(i1,1) &
                  +dgdx(ii,2)*dxfc_nf(i1,2) &
                  +dgdx(ii,3)*dxfc_nf(i1,3)
               p_i=s(ii)+dp
               dp= dgdx(kk,1)*dxfc_non_k(i,1) &
                  +dgdx(kk,2)*dxfc_non_k(i,2) &
                  +dgdx(kk,3)*dxfc_non_k(i,3)
               p_k=s(kk)+dp
               fie=0.5d0*(p_i+p_k)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3) 
            ENDDO
!            
            nv=-1
            nf_number=nf_number_id(nv)
            len=istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               kk=right_non(k)
               ii=left_nf(k)               
               dp= dgdx(ii,1)*dxfc_nf(k,1) &
                  +dgdx(ii,2)*dxfc_nf(k,2) &
                  +dgdx(ii,3)*dxfc_nf(k,3)
               p_i=s(ii)+dp
               dp= dgdx(kk,1)*dxfc_non_k(k,1) &
                  +dgdx(kk,2)*dxfc_non_k(k,2) &
                  +dgdx(kk,3)*dxfc_non_k(k,3)
               p_k=s(kk)+dp
               fie=0.5d0*(p_i+p_k)
               fie_nf(i,1)=-fie*sv_nf(k,1)
               fie_nf(i,2)=-fie*sv_nf(k,2)
               fie_nf(i,3)=-fie*sv_nf(k,3)
            ENDDO                 
!
!...........Cells mcc
!
            nv=1
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            IF(idg.eq.0)THEN
               DO i=1,isize
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  idx=i3invtbl(i)
                  sk=c3pa(1,idx)
                  fie=0.5d0*(s(ii)+sk)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3) 
               ENDDO
            ELSEIF(idg.eq.1)THEN
               DO i=1,isize
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  idx=i3invtbl(i)
                  sk=c3delp(1,idx)
                  fie=0.5d0*(s(ii)+sk)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3)   
               ENDDO
            ELSE
               DO i=1,isize
                  i0=istart0+i !!  
                  i1=istart+i
                  ii=left_nf(i1)
                  sk=s(ii)
                  fie=0.5d0*(s(ii)+sk)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
            ENDIF
!
!...........Cells out
!
            nv=3
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            DO i=1,isize
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               dp= dgdx(ii,1)*dxfc_nf(i1,1) &
                  +dgdx(ii,2)*dxfc_nf(i1,2) &
                  +dgdx(ii,3)*dxfc_nf(i1,3)
               fie=s(ii)+dp
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
!
!...........The rest
!
            DO nv=2,8
               nf_number=nf_number_id(nv)
               IF(nf_number.eq.3) cycle
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               isize =istart_nf(2,nf_number)
               DO i=1,isize
                  i0=istart0+i 
                  i1=istart+i
                  ii=left_nf(i1)
                  fie=s(ii)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
            ENDDO
!
         ENDIF  
!
!........fluxBC: choke model, mcp model, valve model
!  
         IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_gradpK2(s,dgdx,fie_nf)
!
         CALL sum_nf_ndim(0,0,ncell_fp,fie_nf,dsdx)  
!
!........Calculate gradient components at cell centers
!
         IF(ndim.eq.2)THEN
            DO i=1,ncell_fluid
               dgdx(i,1)=dsdx(i,1)*volr(i)
               dgdx(i,2)=dsdx(i,2)*volr(i)
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               dgdx(i,1)=dsdx(i,1)*volr(i)
               dgdx(i,2)=dsdx(i,2)*volr(i)
               dgdx(i,3)=dsdx(i,3)*volr(i)
            ENDDO
         ENDIF
!
!........MPI communication
!
         IF(np.gt.1) CALL communicate_2d(dgdx)
!
      ENDDO
!      
      END SUBROUTINE grad_pressK2

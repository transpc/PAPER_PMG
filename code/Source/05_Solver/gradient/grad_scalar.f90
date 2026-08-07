      SUBROUTINE grad_scalar(s,dsdx,ncell)
!
!     This routine calculates the components of the gradient based on the gauss theorem.
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_tot
      USE Znum_cell    , ONLY: istart_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zcoord3      , ONLY: volr
      USE Zvec_geo     , ONLY: sv_nf,fac1_non,fac_non
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: ncell
      REAL(8),DIMENSION(ncell_fp) :: s
!.....Output
      REAL(8),DIMENSION(ncell,ndim) :: dsdx ! ncell is ncell_fluid or ncell_fp if dsdx needs to be communicated
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,len,i1 
      REAL(8) :: fie
!.....Local vector arrays
      REAL(8),DIMENSION(nf_tot,ndim) :: fie_nf
!
!.....Build summation info for non,inl
!
      nf_number_nb=8
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      nf_number_id(4)=4
      nf_number_id(5)=5
      nf_number_id(6)=6
      nf_number_id(7)=7
      nf_number_id(8)=8
      istart_nfs(0)=0
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
!.........Computing cell
!
         nv=0
         nf_number=nf_number_id(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            fie=fac1_non(i)*s(ii)+fac_non(i)*s(kk)
            fie_nf(i1,1)=fie*sv_nf(i1,1)
            fie_nf(i1,2)=fie*sv_nf(i1,2)
         ENDDO
!
!.......The rest
!
         DO nv=1,8
            nf_number=nf_number_id(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               fie=s(ii)
               fie_nf(i1,1)=fie*sv_nf(i1,1)
               fie_nf(i1,2)=fie*sv_nf(i1,2)
            ENDDO
         ENDDO
      ELSE
!
!.........Computing cell
!
         nv=0
         nf_number=nf_number_id(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            fie=(fac1_non(i)*s(ii)+fac_non(i)*s(kk))
            fie_nf(i1,1)=fie*sv_nf(i1,1)
            fie_nf(i1,2)=fie*sv_nf(i1,2)
            fie_nf(i1,3)=fie*sv_nf(i1,3)
         ENDDO
!
!........The rest
!
         DO nv=1,8
            nf_number=nf_number_id(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               fie=s(ii)
               fie_nf(i1,1)=fie*sv_nf(i1,1)
               fie_nf(i1,2)=fie*sv_nf(i1,2)
               fie_nf(i1,3)=fie*sv_nf(i1,3)
            ENDDO
         ENDDO
      ENDIF
!
      CALL sum_nf_ndim(0,-1,ncell, &
                       fie_nf,dsdx)
!
      IF(ndim.eq.2)THEN
         DO i=1,ncell_fluid
            dsdx(i,1)=dsdx(i,1)*volr(i)
            dsdx(i,2)=dsdx(i,2)*volr(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            dsdx(i,1)=dsdx(i,1)*volr(i)
            dsdx(i,2)=dsdx(i,2)*volr(i)
            dsdx(i,3)=dsdx(i,3)*volr(i)
         ENDDO
      ENDIF
!
      END SUBROUTINE grad_scalar

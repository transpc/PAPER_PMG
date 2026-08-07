!
      SUBROUTINE press_work_face(arvg_nf,avg_nf,arvl_nf,avl_nf,arvg_p,avg_p,arvl_p,avl_p)
!
      USE Zinterface
      USE Zzone         , ONLY: ncell_fluid
      USE Zvec_param    , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_flux
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf, &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zpress        , ONLY: p
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf
      USE Zb_condition  , ONLY: p_fb
!
      IMPLICIT NONE
!
!     input
      REAL(8) arvg_nf(nf_flux),avg_nf(nf_flux),arvl_nf(nf_flux),avl_nf(nf_flux)
!     output
      REAL(8) arvg_p(ncell_fluid),avg_p(ncell_fluid),arvl_p(ncell_fluid),avl_p(ncell_fluid)
!     local variables
      INTEGER i,ii,kk,k
      INTEGER nv,nf_number,istart,len,i1,istart2,i2
!     local vector array
      REAL(8) arvg_p_nf(nf_flux),avg_p_nf(nf_flux),arvl_p_nf(nf_flux),avl_p_nf(nf_flux)
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      lens         =istart_nfs(3)+nf_out
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            avg_p_nf(i1)=avg_nf(i1)*p(ii)
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            avg_p_nf(i1)=avg_nf(i1)*p(kk)
         ELSE
            avg_p_nf(i1)=0.0d0
         ENDIF
         arvg_p_nf(i1)=arvg_nf(i1)
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            avl_p_nf(i1)=avl_nf(i1)*p(ii)
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            avl_p_nf(i1)=avl_nf(i1)*p(kk)
         ELSE
            avl_p_nf(i1)=0.0d0
         ENDIF
         arvl_p_nf(i1)=arvl_nf(i1)
      ENDDO
!      
      nv=1
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         avg_p_nf(i1)=avg_nf(i1)*p(ii)
         avl_p_nf(i1)=avl_nf(i1)*p(ii)           
      ENDDO
!
      nv=2
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,len
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         avg_p_nf (i1)=avg_nf (i1)*p_fb(k)
         avl_p_nf (i1)=avl_nf (i1)*p_fb(k)
      ENDDO
!      
      nv=3
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         avg_p_nf(i1)=avg_nf(i1)*p(ii)
         avl_p_nf(i1)=avl_nf(i1)*p(ii)           
      ENDDO
!
      CALL sum_nf(0,-1,            &
                  arvg_nf ,arvg_p, &
                  avg_p_nf,avg_p,  &
                  arvl_nf ,arvl_p, &
                  avl_p_nf,avl_p)
!
      DO i=1,ncell_fluid
         arvg_p(i)=arvg_p(i)*p(i)
         arvl_p(i)=arvl_p(i)*p(i)
      ENDDO
!
      END SUBROUTINE press_work_face

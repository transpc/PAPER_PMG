      SUBROUTINE face_alpha
 
      USE VOL_DATA
      USE Zzone         , ONLY: ncell_fluid
      USE Zb_condition  , ONLY: alphab_gas,alphab_liq,rhob_liq,rhob_gas
      USE Zvec_geo      , ONLY: f0,f1
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index    , ONLY: left_nf,right_non, &
                                nbcon_nf
      USE Zvec_param    , ONLY: nf_non,nf_inl,nf_out
      
      !
      IMPLICIT NONE
      !
!     INTEGER :: get_nf_number
      INTEGER :: i
      INTEGER :: ii,kk,k
      INTEGER :: nf_number,istart,len,istart2,i1,i2,isize
      !
      REAL(8) :: a,b
!
      REAL(8),ALLOCATABLE :: rhogfluxg_non(:),rhogfluxg_inl(:),rhogfluxg_out(:)
      REAL(8),ALLOCATABLE :: rholfluxl_non(:),rholfluxl_inl(:),rholfluxl_out(:)
      REAL(8),ALLOCATABLE :: alphag_up_non(:),alphag_up_inl(:),alphag_up_out(:)
      REAL(8),ALLOCATABLE :: alphal_up_non(:),alphal_up_inl(:),alphal_up_out(:)
      REAL(8),ALLOCATABLE :: rhogfluxg(:),rholfluxl(:)
!
!-------------------------------------------------------------------------------------------
!
!Case1
! yjm
! define face void fractions
!      IF(initial)THEN
!         initial=.false.
!         ALLOCATE(alphagj(nf_flux),alphalj(nf_flux))
!      ENDIF
!      !
!      nf_number=0
!      istart=istart_nf(1,nf_number)
!      len   =istart_nf(2,nf_number)
!      DO i=1,len  
!         i1=istart+i
!         ii=left_nf(i1)
!         kk=right_non(i)          
!         a=f1(i)*porosity(ii)
!         b=f0(i)*porosity(kk)
!         c=1.d0/perm_non(i)
!         IF(vfporous.eq.0) THEN
!            a=f1(i)
!            b=f0(i)
!            c=1.d0
!         ENDIF
!         alphagj(i1)=(a*cell%alphag(ii)+b*cell%alphag(kk))*c
!         alphalj(i1)=1.0d0-alphagj(i1)
!      ENDDO
!      !
!      nf_number=2  ! inlet  i1=20
!      istart=istart_nf(1,nf_number)
!      istart2=istart_nbcon_nf(nf_number)
!      len   =istart_nf(2,nf_number)
!      DO i=1,len  
!         i1=istart+i
!         i2=istart2+i
!         ii=left_nf(i1)
!         k=nbcon_nf(i2)
!         alphagj(i1)=alphab_gas(k)
!         alphalj(i1)=1.0d0-alphagj(i1)
!      ENDDO
!      !
!      nf_number=3   !  outlet i1=21
!      istart=istart_nf(1,nf_number)
!      len   =istart_nf(2,nf_number)
!      DO i=1,len  
!         i1=istart+i
!         ii=left_nf(i1)
!         a=porosity(ii)/perm_out(i)
!         IF(vfporous.eq.0) THEN
!            a=1.0d0
!         ENDIF
!         alphagj(i1)=cell%alphag(ii)*a
!         alphalj(i1)=1.0d0-alphagj(i1)
!      ENDDO
!!      
!      DO i=1,ncell_fluid
!         DO j=i,num_neigh(i)
!            nf_number=get_nf_number(nbcon(j,i))
!            IF(nf_number.eq.3) goto 300  ! outlet
!            IF(nf_number.eq.2) goto 200  ! inlet
!         ENDDO
!         IF(vg_o(i,ndim).ge.0.0d0)THEN
!            cell%alphagf(i)=alphagj(i)
!            !cell%alphalf(i)=1.0d0-cell%alphagf(i)
!         ELSE
!            cell%alphagf(i)=alphagj(i-1)
!            !cell%alphalf(i)=1.0d0-cell%alphagf(i)
!         ENDIF
!         IF(vl_o(i,ndim).ge.0.0d0)THEN
!            cell%alphalf(i)=alphalj(i)
!         ELSE
!            cell%alphalf(i)=alphalj(i-1)
!         ENDIF
!         goto 100
!300      CONTINUE
!         IF(vg_o(i,ndim).ge.0.0d0)THEN
!            cell%alphagf(i)=alphagj(i)
!         ELSE
!            cell%alphagf(i)=alphagj(21)
!         ENDIF
!         IF(vl_o(i,ndim).ge.0.0d0)THEN
!            cell%alphalf(i)=alphalj(i)
!         ELSE
!            cell%alphalf(i)=alphalj(21)
!         ENDIF
!         goto 100
!200      CONTINUE
!         IF(vg_o(i,ndim).ge.0.0d0)THEN
!            cell%alphagf(i)=alphagj(i)
!         ELSE
!            cell%alphagf(i)=alphagj(21)
!         ENDIF
!         IF(vl_o(i,ndim).ge.0.0d0)THEN
!            cell%alphalf(i)=alphalj(i)
!         ELSE
!            cell%alphalf(i)=alphalj(i-1)
!         ENDIF
!100      CONTINUE
!      ENDDO
!
!-------------------------------------------------------------------------------------------
!
!!Case2
!************************************
! find upwind void fraction: face mass flux
!************************************
      ALLOCATE(rhogfluxg_non(nf_non),rholfluxl_non(nf_non))
      ALLOCATE(rhogfluxg_inl(nf_inl),rholfluxl_inl(nf_inl))
      ALLOCATE(rhogfluxg_out(nf_out),rholfluxl_out(nf_out))
!
      ALLOCATE(alphag_up_non(nf_non),alphal_up_non(nf_non))
      ALLOCATE(alphag_up_inl(nf_inl),alphal_up_inl(nf_inl))
      ALLOCATE(alphag_up_out(nf_out),alphal_up_out(nf_out))    
!
      ALLOCATE(rhogfluxg(ncell_fluid),rholfluxl(ncell_fluid))
!
      rhogfluxg(:)=0.d0
      rholfluxl(:)=0.d0
!     
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)     
         a=f1(i)
         b=f0(i)
!         
         rhogfluxg_non(i)=a*cell%rhog(ii)*flux_g_nf(i1) + b*cell%rhog(kk)*flux_g_nf(i1)               
!         rhogfluxg_non(i)=a*flux_g_nf(i1) + b*flux_g_nf(i1)
         rholfluxl_non(i)=a*cell%rhol(ii)*flux_l_nf(i1) + b*cell%rhol(kk)*flux_l_nf(i1)
!         rholfluxl_non(i)=a*flux_l_nf(i1) + b*flux_l_nf(i1)         
!         
         IF(flux_g_nf(i1).lt.0.d0) THEN
            rhogfluxg(ii)=rhogfluxg(ii)+DABS(rhogfluxg_non(i))
         ELSE
            rhogfluxg(ii)=rhogfluxg(ii)
         ENDIF
         IF(flux_l_nf(i1).lt.0.d0) THEN
            rholfluxl(ii)=rholfluxl(ii)+DABS(rholfluxl_non(i))
         ELSE
            rholfluxl(ii)=rholfluxl(ii)     
         ENDIF         
      ENDDO
!    
      nf_number=2
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,isize
         i1=istart+i        !global face number 
         i2=istart2+i            
         ii=left_nf(i1)     !global cell number 
         k=nbcon_nf(i2)     !inlet number
!
         rhogfluxg_inl(i)=rhob_gas(k)*flux_g_nf(i1)
!         rhogfluxg_inl(i)=flux_g_nf(i1)
         rholfluxl_inl(i)=rhob_liq(k)*flux_l_nf(i1)
!         rholfluxl_inl(i)=flux_l_nf(i1)
!         
         IF(flux_g_nf(i1).lt.0.d0) THEN
            rhogfluxg(ii)=rhogfluxg(ii)+DABS(rhogfluxg_inl(i))
         ELSE
            rhogfluxg(ii)=rhogfluxg(ii)   
         ENDIF
         IF(flux_l_nf(i1).lt.0.d0) THEN
            rholfluxl(ii)=rholfluxl(ii)+DABS(rholfluxl_inl(i))
         ELSE
            rholfluxl(ii)=rholfluxl(ii)              
         ENDIF         
      ENDDO         
!
!************************************
!find upwind void fraction: face void fraction
!************************************
!
      cell%alphagf(:)=cell%alphag(:)
      cell%alphalf(:)=cell%alphal(:)
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)

      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)          
         a=f1(i)
         b=f0(i)
!         
         alphag_up_non(i)=a*cell%alphag(ii)*cell%rhog(ii)*flux_g_nf(i1)   &
                         +b*cell%alphag(kk)*cell%rhog(kk)*flux_g_nf(i1)
!         alphag_up_non(i)=a*cell%alphag(ii)*flux_g_nf(i1)   &
!                         +b*cell%alphag(kk)*flux_g_nf(i1)    
         alphal_up_non(i)=a*cell%alphal(ii)*cell%rhol(ii)*flux_l_nf(i1)   &
                         +b*cell%alphal(kk)*cell%rhol(kk)*flux_l_nf(i1)
!         alphal_up_non(i)=a*cell%alphal(ii)*flux_l_nf(i1)   &
!                         +b*cell%alphal(kk)*flux_l_nf(i1)                           
!         
         IF(flux_g_nf(i1).lt.0.d0) THEN
            cell%alphagf(ii)=cell%alphagf(ii)+DABS(alphag_up_non(i))/rhogfluxg(ii)
         ENDIF  
         IF(flux_l_nf(i1).lt.0.d0) THEN                 
            cell%alphalf(ii)=cell%alphalf(ii)+DABS(alphal_up_non(i))/rholfluxl(ii)
         ENDIF
      ENDDO
!
      nf_number=2
     
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,isize
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
!
         alphag_up_inl(i)=alphab_gas(k) !*rhob_gas(k)*flux_g_nf(i1)
         alphal_up_inl(i)=alphab_liq(k) !*rhob_liq(k)*flux_l_nf(i1)
!         alphag_up_inl(i)=alphab_gas(k)*flux_g_nf(i1)
!         alphal_up_inl(i)=alphab_liq(k)*flux_l_nf(i1)
!         
!         IF(flux_g_nf(i1).lt.0.d0) THEN
            cell%alphagf(ii)=cell%alphagf(ii)+DABS(alphag_up_inl(i)) !/rhogfluxg(ii)
!         ENDIF  
!         IF(flux_l_nf(i1).lt.0.d0) THEN                 
            cell%alphalf(ii)=cell%alphalf(ii)+DABS(alphal_up_inl(i)) !/rholfluxl(ii)
!         ENDIF
      ENDDO 
!      
      nf_number=3
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         alphag_up_out(i)=cell%alphag(ii) !*cell%rhog(ii)*flux_g_nf(i1)
         alphal_up_out(i)=cell%alphal(ii) !*cell%rhol(ii)*flux_l_nf(i1)
!         alphag_up_out(i)=cell%alphag(ii)*flux_g_nf(i1)
!         alphal_up_out(i)=cell%alphal(ii)*flux_l_nf(i1)
!
         IF(flux_g_nf(i1).lt.0.d0) THEN
            cell%alphagf(ii)=cell%alphagf(ii)+DABS(alphag_up_out(i)) !/rhogfluxg(ii)
         ENDIF  
         IF(flux_l_nf(i1).lt.0.d0) THEN                 
            cell%alphalf(ii)=cell%alphalf(ii)+DABS(alphal_up_out(i)) !/rholfluxl(ii)
         ENDIF          
      ENDDO 
!      
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         IF(cell%alphagf(ii).eq.0.d0) THEN
            cell%alphagf(ii)=cell%alphag(ii)
         ENDIF
         IF(cell%alphalf(ii).eq.0.d0) THEN
            cell%alphalf(ii)=cell%alphal(ii)
         ENDIF
      ENDDO 
!
      DEALLOCATE(rhogfluxg_non,rhogfluxg_inl,rhogfluxg_out)
      DEALLOCATE(rholfluxl_non,rholfluxl_inl,rholfluxl_out)
      DEALLOCATE(alphag_up_non,alphag_up_inl,alphag_up_out)
      DEALLOCATE(alphal_up_non,alphal_up_inl,alphal_up_out)
      DEALLOCATE(rhogfluxg,rholfluxl) 
!
!-------------------------------------------------------------------------------------------
!
! Case3
!      ALLOCATE(maxfluxg(ncell_fluid),maxfluxl(ncell_fluid))
!      ALLOCATE(max_neg_fluxg_i1_non(ncell_fluid),max_neg_fluxl_i1_non(ncell_fluid))
!      ALLOCATE(max_neg_fluxg_i1_inl(ncell_fluid),max_neg_fluxl_i1_inl(ncell_fluid))
!!
!      maxfluxg(:)=0.d0
!      maxfluxl(:)=0.d0
!      max_neg_fluxg_i1_non(:)=0
!      max_neg_fluxl_i1_non(:)=0      
!!     
!      nf_number=0
!      istart=istart_nf(1,nf_number)
!      len   =istart_nf(2,nf_number)
!      DO i=1,len  
!         i1=istart+i
!         ii=left_nf(i1)
!         kk=right_non(i)     
!!         
!         IF(flux_g_nf(i1).lt.0.d0) THEN
!            IF(DABS(flux_g_nf(i1)).gt.maxfluxg(ii)) THEN
!               maxfluxg(ii)=DABS(flux_g_nf(i1))
!               max_neg_fluxg_i1_non(ii)=i1
!            ENDIF   
!         ENDIF
!         IF(flux_l_nf(i1).lt.0.d0) THEN
!            IF(DABS(flux_l_nf(i1)).gt.maxfluxl(ii)) THEN
!               maxfluxl(ii)=DABS(flux_l_nf(i1))
!               max_neg_fluxl_i1_non(ii)=i1
!            ENDIF   
!         ENDIF   
!!              
!      ENDDO
!      
!
!!************************************
!!find upwind void fraction: face void fraction
!!************************************
!!
!      cell%alphagf(:)=cell%alphag(:)
!      cell%alphalf(:)=cell%alphal(:)
!!
!      nf_number=0
!      istart=istart_nf(1,nf_number)
!      len   =istart_nf(2,nf_number)
!      DO i=1,len  
!         i1=istart+i
!         ii=left_nf(i1)
!         kk=right_non(i)          
!         a=f1(i)
!         b=f0(i)
!!                     
!         IF(i1.eq.max_neg_fluxg_i1_non(ii)) THEN
!            cell%alphagf(ii)=a*cell%alphag(ii)+b*cell%alphag(kk)
!         ENDIF
!         IF(i1.eq.max_neg_fluxl_i1_non(ii)) THEN
!            cell%alphalf(ii)=a*cell%alphal(ii)+b*cell%alphal(kk)
!         ENDIF         
!!
!      ENDDO
!!
!      nf_number=2
!      istart=istart_nf(1,nf_number)
!      isize =istart_nf(2,nf_number)
!      istart2=istart_nbcon_nf(nf_number)
!      DO i=1,isize
!         i1=istart+i
!         i2=istart2+i
!         ii=left_nf(i1)
!         k=nbcon_nf(i2)
!!
!!         IF(i1.eq.max_neg_fluxg_i1_inl(ii)) THEN
!            cell%alphagf(ii)=alphab_gas(k)
!!         ENDIF
!!         IF(i1.eq.max_neg_fluxl_i1_inl(ii)) THEN
!            cell%alphalf(ii)=alphab_liq(k)
!!         ENDIF         
!!
!      ENDDO 
!!
!      nf_number=3
!      istart=istart_nf(1,nf_number)
!      len   =istart_nf(2,nf_number)
!      DO i=1,len
!         i1=istart+i
!         ii=left_nf(i1)
!         cell%alphagf(ii)=cell%alphag(ii)
!         cell%alphalf(ii)=cell%alphal(ii)
!      ENDDO       
!!      
!      DEALLOCATE(maxfluxg,maxfluxl) 
!      DEALLOCATE(max_neg_fluxg_i1_non,max_neg_fluxl_i1_non)
!      DEALLOCATE(max_neg_fluxg_i1_inl,max_neg_fluxl_i1_inl)
! 
      RETURN
      END SUBROUTINE face_alpha


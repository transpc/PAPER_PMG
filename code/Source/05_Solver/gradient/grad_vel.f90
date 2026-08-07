      SUBROUTINE grad_vel(iflag,s,dsdx,vb,vin)
!
!     This routine calculates the components of the gradient
!     vector of velocity at the cell center, using conservative
!     scheme based on the gauss theorem.
!
      USE Zzone       , ONLY: ncell_fluid
      USE Zmpi        , ONLY: ncell_fp
      USE Zparam      , ONLY: ndim,nb_max
      USE Zvec_param  , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_sym,nf_nonk
      USE Znum_cell   , ONLY: istart_nf,istart_nbcon_nf,right_nb_k, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index  , ONLY: left_nf,right_non,nbcon_nf
      USE Zcoord3     , ONLY: volr
      USE Zuserdefined, ONLY: vel_bc_profile_inl
      USE c3com_cupid , ONLY: i3invtbl
      USE Zbc_index   , ONLY: vin_norm
      USE Zvec_geo    , ONLY: xn_nf,sv_nf,     &
                              fac_non,fac1_non
      USE Zcore       , ONLY: myrank
!
      IMPLICIT NONE 
      INCLUDE '../../10_LinkToMARS/c3com.h' 
!
!.....Input
      INTEGER :: iflag 
      REAL(8),DIMENSION(ncell_fp,ndim) :: s
      REAL(8),DIMENSION(nb_max,ndim) :: vb
      REAL(8),DIMENSION(nb_max) :: vin
!.....Output
      REAL(8),DIMENSION(ncell_fp,ndim,ndim) :: dsdx
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk,idx
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: c1,c2,contra 
      REAL(8) :: f_profile
      REAL(8) :: ggs1,ggs2,ggs3
      REAL(8) :: fie1,fie2,fie3
!.....Local vector arrays
      !REAL(8),DIMENSION(nf_non+nf_mcc+nf_inl+nf_out+nf_sym,ndim,ndim) :: fie_nf
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_mcc+nf_inl+nf_out+nf_sym,ndim,ndim) :: fie_nf
!
!.....Build summation info for all nf
!
      nf_number_nb=4
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      nf_number_id(4)=8
!      istart_nfs(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk      
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      istart_nfs(4)=istart_nfs(3)+nf_out
      lens         =istart_nfs(4)+nf_sym
!
!.....Cells non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2) THEN
         DO i=1,len  
            i0=istart0+i 
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            c1=fac1_non(i)
            c2=fac_non(i)
            fie1=c1*s(ii,1)+c2*s(kk,1)
            fie2=c1*s(ii,2)+c2*s(kk,2)
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
         ENDDO
      ELSE
         DO i=1,len  
            i0=istart0+i 
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            c1=fac1_non(i)
            c2=fac_non(i)
            fie1=c1*s(ii,1)+c2*s(kk,1)
            fie2=c1*s(ii,2)+c2*s(kk,2)
            fie3=c1*s(ii,3)+c2*s(kk,3)
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,3,1)=fie1*sv_nf(i1,3)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
            fie_nf(i0,3,2)=fie2*sv_nf(i1,3)
            fie_nf(i0,1,3)=fie3*sv_nf(i1,1)
            fie_nf(i0,2,3)=fie3*sv_nf(i1,2)
            fie_nf(i0,3,3)=fie3*sv_nf(i1,3)
         ENDDO
      ENDIF
!      
      nv=-1
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2) THEN
         DO i=1,len  
            k=right_nb_k(i)
            kk=right_non(k)
            ii=left_nf(k) 
         !ii=right_non(k)
         !kk=left_nf(k)            
            c1=fac1_non(k)
            c2=fac_non(k)
            fie1=c1*s(ii,1)+c2*s(kk,1)
            fie2=c1*s(ii,2)+c2*s(kk,2)
            fie_nf(i,1,1)=-fie1*sv_nf(k,1)
            fie_nf(i,2,1)=-fie1*sv_nf(k,2)
            fie_nf(i,1,2)=-fie2*sv_nf(k,1)
            fie_nf(i,2,2)=-fie2*sv_nf(k,2)
         ENDDO
      ELSE
         DO i=1,len  
            k=right_nb_k(i) 
            kk=right_non(k)
            ii=left_nf(k) 
         !ii=right_non(k)
         !kk=left_nf(k)            
            c1=fac1_non(k)
            c2=fac_non(k)
            fie1=c1*s(ii,1)+c2*s(kk,1)
            fie2=c1*s(ii,2)+c2*s(kk,2)
            fie3=c1*s(ii,3)+c2*s(kk,3)
            fie_nf(i,1,1)=-fie1*sv_nf(k,1)
            fie_nf(i,2,1)=-fie1*sv_nf(k,2)
            fie_nf(i,3,1)=-fie1*sv_nf(k,3)
            fie_nf(i,1,2)=-fie2*sv_nf(k,1)
            fie_nf(i,2,2)=-fie2*sv_nf(k,2)
            fie_nf(i,3,2)=-fie2*sv_nf(k,3)
            fie_nf(i,1,3)=-fie3*sv_nf(k,1)
            fie_nf(i,2,3)=-fie3*sv_nf(k,2)
            fie_nf(i,3,3)=-fie3*sv_nf(k,3)
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
      IF(iflag.eq.1) THEN
         IF(ndim.eq.2) THEN
            DO i=1,len  
               i0=istart0+i  
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               fie1=c3vg(1,idx)*xn_nf(i1,1)
               fie2=c3vg(1,idx)*xn_nf(i1,2)
               fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
               fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
               fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
               fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
            ENDDO
         ELSE
            DO i=1,len  
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               fie1=c3vg(1,idx)*xn_nf(i1,1)
               fie2=c3vg(1,idx)*xn_nf(i1,2)
               fie3=c3vg(1,idx)*xn_nf(i1,3)
               fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
               fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
               fie_nf(i0,3,1)=fie1*sv_nf(i1,3)
               fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
               fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
               fie_nf(i0,3,2)=fie2*sv_nf(i1,3)
               fie_nf(i0,1,3)=fie3*sv_nf(i1,1)
               fie_nf(i0,2,3)=fie3*sv_nf(i1,2)
               fie_nf(i0,3,3)=fie3*sv_nf(i1,3)
            ENDDO
         ENDIF
      ELSEIF(iflag.eq.2)THEN
         IF(ndim.eq.2) THEN
            DO i=1,len  
               i0=istart0+i   
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               fie1=c3vl(1,idx)*xn_nf(i1,1)
               fie2=c3vl(1,idx)*xn_nf(i1,2)
               fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
               fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
               fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
               fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
            ENDDO
         ELSE
            DO i=1,len  
               i0=istart0+i   
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               fie1=c3vl(1,idx)*xn_nf(i1,1)
               fie2=c3vl(1,idx)*xn_nf(i1,2)
               fie3=c3vl(1,idx)*xn_nf(i1,3)
               fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
               fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
               fie_nf(i0,3,1)=fie1*sv_nf(i1,3)
               fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
               fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
               fie_nf(i0,3,2)=fie2*sv_nf(i1,3)
               fie_nf(i0,1,3)=fie3*sv_nf(i1,1)
               fie_nf(i0,2,3)=fie3*sv_nf(i1,2)
               fie_nf(i0,3,3)=fie3*sv_nf(i1,3)
            ENDDO
         ENDIF
      ELSE
         IF(myrank.eq.0) THEN
            WRITE(*,"(11x,a,1i3)")'iflag should be 1 or 2, but iflag=',iflag
         ENDIF
         CALL finalize_mpi
         STOP
      ENDIF
!
!.....Cells inl
!
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      IF(ndim.eq.2) THEN
         DO i=1,len  
            i0=istart0+i 
            i1=istart+i
            i2=istart2+i
            k=nbcon_nf(i2)
            f_profile=vel_bc_profile_inl(i)
            IF(vin_norm(k).eq.0)THEN
               fie1=vb(k,1)*f_profile
               fie2=vb(k,2)*f_profile
            ELSE
               fie1=vin(k)*xn_nf(i1,1)*f_profile
               fie2=vin(k)*xn_nf(i1,2)*f_profile
            ENDIF
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
         ENDDO
      ELSE
         DO i=1,len  
            i0=istart0+i 
            i1=istart+i
            i2=istart2+i
            k=nbcon_nf(i2)
            f_profile=vel_bc_profile_inl(i)
            IF(vin_norm(k).eq.0)THEN
               fie1=vb(k,1)*f_profile
               fie2=vb(k,2)*f_profile
               fie3=vb(k,3)*f_profile
            ELSE
               fie1=vin(k)*xn_nf(i1,1)*f_profile
               fie2=vin(k)*xn_nf(i1,2)*f_profile
               fie3=vin(k)*xn_nf(i1,3)*f_profile
            ENDIF
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,3,1)=fie1*sv_nf(i1,3)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
            fie_nf(i0,3,2)=fie2*sv_nf(i1,3)
            fie_nf(i0,1,3)=fie3*sv_nf(i1,1)
            fie_nf(i0,2,3)=fie3*sv_nf(i1,2)
            fie_nf(i0,3,3)=fie3*sv_nf(i1,3)
         ENDDO
      ENDIF
!
!.....Cells out
!
      nv=3
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2) THEN
         DO i=1,len 
            i0=istart0+i  
            i1=istart+i
            ii=left_nf(i1)
            fie1=s(ii,1)
            fie2=s(ii,2)
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
         ENDDO
      ELSE
         DO i=1,len  
            i0=istart0+i   
            i1=istart+i
            ii=left_nf(i1)
            fie1=s(ii,1)
            fie2=s(ii,2)
            fie3=s(ii,3)
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,3,1)=fie1*sv_nf(i1,3)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
            fie_nf(i0,3,2)=fie2*sv_nf(i1,3)
            fie_nf(i0,1,3)=fie3*sv_nf(i1,1)
            fie_nf(i0,2,3)=fie3*sv_nf(i1,2)
            fie_nf(i0,3,3)=fie3*sv_nf(i1,3)
         ENDDO
      ENDIF
!
!.....Cells sym
!
      nv=4
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2) THEN
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            contra= s(ii,1)*xn_nf(i1,1) &
                   +s(ii,2)*xn_nf(i1,2)
            ggs1=s(ii,1)-2.d0*contra*xn_nf(i1,1)
            ggs2=s(ii,2)-2.d0*contra*xn_nf(i1,2)
            fie1=0.5d0*(s(ii,1)+ggs1)
            fie2=0.5d0*(s(ii,2)+ggs2)
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
         ENDDO
      ELSE
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            contra= s(ii,1)*xn_nf(i1,1) &
                   +s(ii,2)*xn_nf(i1,2) &
                   +s(ii,3)*xn_nf(i1,3)
            ggs1=s(ii,1)-2.d0*contra*xn_nf(i1,1)
            ggs2=s(ii,2)-2.d0*contra*xn_nf(i1,2)
            ggs3=s(ii,3)-2.d0*contra*xn_nf(i1,3)
            fie1=0.5d0*(s(ii,1)+ggs1)
            fie2=0.5d0*(s(ii,2)+ggs2)
            fie3=0.5d0*(s(ii,3)+ggs3)
            fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
            fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
            fie_nf(i0,3,1)=fie1*sv_nf(i1,3)
            fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
            fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
            fie_nf(i0,3,2)=fie2*sv_nf(i1,3)
            fie_nf(i0,1,3)=fie3*sv_nf(i1,1)
            fie_nf(i0,2,3)=fie3*sv_nf(i1,2)
            fie_nf(i0,3,3)=fie3*sv_nf(i1,3)
         ENDDO
      ENDIF
!
      !CALL sum_nf_ndim2(0,-1,ncell_fp, &
      !                  fie_nf,dsdx)
      CALL sum_nf_ndim2(0,0,ncell_fp, &
                        fie_nf,dsdx)      
!
      IF(ndim.eq.2)THEN
         DO i=1,ncell_fluid
            dsdx(i,1,1)=dsdx(i,1,1)*volr(i)
            dsdx(i,1,2)=dsdx(i,1,2)*volr(i)
            dsdx(i,2,1)=dsdx(i,2,1)*volr(i)
            dsdx(i,2,2)=dsdx(i,2,2)*volr(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            dsdx(i,1,1)=dsdx(i,1,1)*volr(i)
            dsdx(i,1,2)=dsdx(i,1,2)*volr(i)
            dsdx(i,1,3)=dsdx(i,1,3)*volr(i)
            dsdx(i,2,1)=dsdx(i,2,1)*volr(i)
            dsdx(i,2,2)=dsdx(i,2,2)*volr(i)
            dsdx(i,2,3)=dsdx(i,2,3)*volr(i)
            dsdx(i,3,1)=dsdx(i,3,1)*volr(i)
            dsdx(i,3,2)=dsdx(i,3,2)*volr(i)
            dsdx(i,3,3)=dsdx(i,3,3)*volr(i)
         ENDDO
      ENDIF
!
      END SUBROUTINE grad_vel

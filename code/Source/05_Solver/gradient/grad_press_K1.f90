!
      SUBROUTINE grad_pressK1(s,dsdx,idg)
!
!     gravity-weighted, gradP only
!     
!     This routine calculates the components of the gradient
!     vector of pressure at the cell center, using conservative
!     scheme based on the gauss theorem.
!
!     idg (only for MARS interfaces): 1=p, 2=dp, 3=ELSE 
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell            
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_nonk
      USE Znum_cell    , ONLY: istart_nf,right_nb_k, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zconst2      , ONLY: grav,gfactor
      USE c3com_cupid  , ONLY: i3invtbl
      USE Zcoord3      , ONLY: volr
      USE Zvec_geo     , ONLY: sv_nf,dxfc_nf,        &
                               dxfc_non_k
      USE Zmcp 
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
      INTEGER :: nv,nf_number,istart,len,i1,istart0
      INTEGER :: ii,kk,i0,k
      REAL(8) :: c1,c2,fie,p_i,p_k,dp
      REAL(8) :: mars_rhom,sk
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
         nv=-1
         nf_number=nf_number_id(nv)
         len=istart_nf(2,nf_number)
         DO i=1,len  
            k=right_nb_k(i)
            kk=right_non(k) 
            ii=left_nf(k)            
            p_i= grav(1)*dxfc_nf(k,1)  *gfactor(ii) &
                +grav(2)*dxfc_nf(k,2)  *gfactor(ii)
            p_k= grav(1)*dxfc_non_k(k,1)*gfactor(ii) &
                +grav(2)*dxfc_non_k(k,2)*gfactor(ii)
            fie=0.5d0*((s(ii)+s(kk))+(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))
            fie_nf(i,1)=-fie*sv_nf(k,1)
            fie_nf(i,2)=-fie*sv_nf(k,2)
         ENDDO   
!         
         DO nv=0,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            IF(nf_number.eq.0) THEN
!
!..............Cells non
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  p_i= grav(1)*dxfc_nf(i1,1)  *gfactor(ii) &
                      +grav(2)*dxfc_nf(i1,2)  *gfactor(ii)
                  p_k= grav(1)*dxfc_non_k(i,1)*gfactor(ii) &
                      +grav(2)*dxfc_non_k(i,2)*gfactor(ii)
                  fie=0.5d0*((s(ii)+s(kk))+(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
!
            ELSEIF(nf_number.eq.1) THEN
!
!..............Cells mcc
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  idx=i3invtbl(i)
                  mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
                  c1=0.5d0/cell%rhom(ii)
                  c2=0.5d0/mars_rhom
                  IF(idg.eq.0)THEN
                     sk=c3pa(1,idx)
                  ELSEIF(idg.eq.1)THEN
                     sk=c3delp(1,idx)
                  ELSEIF(idg.eq.3)THEN
                     sk=c3rtp(1,idx,1) !al
                  ELSEIF(idg.eq.4)THEN
                     sk=c3rtp(1,idx,2) !ag                  
                  ELSE
                     sk=s(ii)
                  ENDIF
                  fie=(c1*s(ii)+c2*sk)/(c1+c2)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
!
            ELSEIF(nf_number.eq.8 .or. nf_number.eq.3) THEN
!
!..............Cells out,sym
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  fie=s(ii)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)  
               ENDDO
!
            ELSE
!
!..............The rest 
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  dp= grav(1)*dxfc_nf(i1,1)*gfactor(ii)  &
                     +grav(2)*dxfc_nf(i1,2)*gfactor(ii)
                  dp=cell%rhom(ii)*dp
                  fie=s(ii)+dp
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
               ENDDO
            ENDIF
         ENDDO
      ELSE
!          
         nv=-1
         nf_number=nf_number_id(nv)
         len=istart_nf(2,nf_number)
         DO i=1,len  
            k=right_nb_k(i)
            kk=right_non(k)  
            ii=left_nf(k)            
            p_i= grav(1)*dxfc_nf(k,1)  *gfactor(ii) &
                +grav(2)*dxfc_nf(k,2)  *gfactor(ii) &
                +grav(3)*dxfc_nf(k,3)  *gfactor(ii)
            p_k= grav(1)*dxfc_non_k(k,1)*gfactor(ii) &
                +grav(2)*dxfc_non_k(k,2)*gfactor(ii) &
                +grav(3)*dxfc_non_k(k,3)*gfactor(ii)
            fie=0.5d0*((s(ii)+s(kk))+(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))
            fie_nf(i,1)=-fie*sv_nf(k,1)
            fie_nf(i,2)=-fie*sv_nf(k,2)
            fie_nf(i,3)=-fie*sv_nf(k,3)
         ENDDO                     
!          
         DO nv=0,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            IF(nf_number.eq.0) THEN
!
!..............Cells non
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  p_i= grav(1)*dxfc_nf(i1,1)  *gfactor(ii) &
                      +grav(2)*dxfc_nf(i1,2)  *gfactor(ii) &
                      +grav(3)*dxfc_nf(i1,3)  *gfactor(ii)
                  p_k= grav(1)*dxfc_non_k(i,1)*gfactor(ii) &
                      +grav(2)*dxfc_non_k(i,2)*gfactor(ii) &
                      +grav(3)*dxfc_non_k(i,3)*gfactor(ii)
                  fie=0.5d0*((s(ii)+s(kk))+(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
!
            ELSEIF(nf_number.eq.1) THEN
!
!..............Cells mcc
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  idx=i3invtbl(i)
                  mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
                  c1=0.5d0/cell%rhom(ii)
                  c2=0.5d0/mars_rhom
                  IF(idg.eq.0)THEN
                     sk=c3pa(1,idx)
                  ELSEIF(idg.eq.1)THEN
                     sk=c3delp(1,idx)
                  ELSEIF(idg.eq.3)THEN
                     sk=c3rtp(1,idx,1) !al
                  ELSEIF(idg.eq.4)THEN
                     sk=c3rtp(1,idx,2) !ag                  
                  ELSE
                     sk=s(ii)
                  ENDIF
                  fie=(c1*s(ii)+c2*sk)/(c1+c2)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
!
            ELSEIF(nf_number.eq.8 .or. nf_number.eq.3) THEN
!
!..............Cells out,sym
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  fie=s(ii)
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
!
            ELSE 
!
!..............The rest 
!
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  dp= grav(1)*dxfc_nf(i1,1)*gfactor(ii) &
                     +grav(2)*dxfc_nf(i1,2)*gfactor(ii) &
                     +grav(3)*dxfc_nf(i1,3)*gfactor(ii)
                  dp=cell%rhom(ii)*dp
                  fie=s(ii)+dp
                  fie_nf(i0,1)=fie*sv_nf(i1,1)
                  fie_nf(i0,2)=fie*sv_nf(i1,2)
                  fie_nf(i0,3)=fie*sv_nf(i1,3)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!.....fluxBC: choke model, mcp model, valve model
!  
      IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_gradpK1(s,fie_nf)      
!
      CALL sum_nf_ndim(0,0,ncell_fp,fie_nf,dsdx)
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
      END SUBROUTINE grad_pressK1

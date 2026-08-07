      SUBROUTINE grad_press(s,dsdx,idg)
!
!     This routine calculates the components of the gradient vector of pressure at the cell center
!     using conservative scheme based on the gauss theorem.
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
      USE Zconst1      , ONLY: lsquareoff
      USE Zcoord3      , ONLY: volr
      USE c3com_cupid  , ONLY: i3invtbl      
      USE Zvec_geo     , ONLY: sv_nf,                 &
                               fac1_non,fac_non,xn_nf
      USE Zgrad_ls_c3d , ONLY: lsindex
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
      INTEGER :: i,idx,i0,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,len,i1,istart0  
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: c1,c2,fie
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
!........Cells non
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv) !!
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            c1=fac1_non(i)*cell%rhomr(ii)
            c2=fac_non(i) *cell%rhomr(kk)
            fie=(c1*s(ii)+c2*s(kk))/(c1+c2)
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
            c1=fac1_non(k)*cell%rhomr(ii)  !f1*left + f*right
            c2=fac_non(k) *cell%rhomr(kk)
            fie=(c1*s(ii)+c2*s(kk))/(c1+c2)
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
         len   =istart_nf(2,nf_number)
         IF(idg.eq.0)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3pa(1,idx)
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ELSEIF(idg.eq.1)THEN
            DO i=1,len  
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3delp(1,idx)
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ELSEIF(idg.eq.3)THEN
            DO i=1,len  
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3rtp(1,idx,1) !al
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ELSEIF(idg.eq.4)THEN
            DO i=1,len  
               i0=istart0+i  
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3rtp(1,idx,2) !ag                  
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDDO
         ELSE
            DO i=1,len  
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=s(ii)
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
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
            len   =istart_nf(2,nf_number)
            DO i=1,len  
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
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            c1=fac1_non(i)*cell%rhomr(ii)
            c2=fac_non(i) *cell%rhomr(kk)
            fie=(c1*s(ii)+c2*s(kk))/(c1+c2)
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
            c1=fac1_non(k)*cell%rhomr(ii)
            c2=fac_non(k) *cell%rhomr(kk)
            fie=(c1*s(ii)+c2*s(kk))/(c1+c2)
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
         len   =istart_nf(2,nf_number)
         IF(idg.eq.0)THEN
            DO i=1,len  
               i0=istart0+i  
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3pa(1,idx)
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSEIF(idg.eq.1)THEN
            DO i=1,len  
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3delp(1,idx)
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSEIF(idg.eq.3)THEN
            DO i=1,len  
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3rtp(1,idx,1) !al
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSEIF(idg.eq.4)THEN
            DO i=1,len
               i0=istart0+i 
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=c3rtp(1,idx,2) !ag                  
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDDO
         ELSE
            DO i=1,len  
               i0=istart0+i  
               i1=istart+i
               ii=left_nf(i1)
               idx=i3invtbl(i)
               mars_rhom=c3rtp(1,idx,2)*c3rtp(1,idx,5)+c3rtp(1,idx,1)*c3rtp(1,idx,4)
               c1=0.5d0*cell%rhomr(ii)
               c2=0.5d0/mars_rhom
               sk=s(ii)
               fie=(c1*s(ii)+c2*sk)/(c1+c2)
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
            len   =istart_nf(2,nf_number)
            DO i=1,len 
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
!......fluxBC: choke model, mcp model, valve model
!  
      IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_gradp(s,fie_nf)
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
!.....Apply least square & linear interpoltion for the pressure gradient at the boundary cells
!
      IF(initial.and.lsquareoff.eq.2) THEN
         initial=.false.
         nf_number=4
         istart=istart_nf(1,nf_number)
         len=istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            IF(dabs(xn_nf(i1,3)).lt.0.5d0) lsindex(ii)=0
         END DO
      END IF
!      
      IF(lsquareoff.eq.0 .or. lsquareoff.eq.2) CALL grad_press_ls(s,dsdx)
!
      END SUBROUTINE grad_press
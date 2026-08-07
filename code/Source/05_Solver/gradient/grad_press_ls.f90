!
      SUBROUTINE grad_press_ls(s,dsdx)
!
!     Calculates the components of the gradient of pressure at the cell center, using least square method and linear interpolation in 3-DIMENSION.
!     (only when "lsquareoff" option is 1)
!
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: myrank
      USE Zparam       , ONLY: ndim      
      USE Zio_unit     , ONLY: unit_screen
      USE Zcoord1      , ONLY: xloc
      USE Zgrad_ls_c2d , ONLY: a11_2,a12_2,a_2,a22_2,det_2
      USE Zgrad_ls_c3d , ONLY: lsindex,det_3,a11_3,a12_3,a13_3,a_3,a22_3,a23_3,a33_3
      USE Zmpi         , ONLY: jperm      
      USE Zbc_index    , ONLY: ngrad
!
      USE Zvec_param   , ONLY: nf_non,nf_sym
      USE Znum_cell    , ONLY: i_neigh,neigh,istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: xn_nf,dxfc_nf
      USE Zrv_model    , ONLY: rv_valve        
!
      IMPLICIT NONE
!
      LOGICAL,SAVE::INITIAL     
!         
!.....Input
      REAL(8) :: s(ncell_fp)
!.....Output
      REAL(8) :: dsdx(ncell_fp,ndim)
!.....Local variables
      INTEGER :: i,ix,j
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,isize,i1
      REAL(8) :: eps,contra,det_tmp      
      REAL(8) :: xxc1,xxc2,xxc3
      REAL(8) :: xloc21,xloc22,xloc23
      REAL(8) :: f_non1,f_non2,f_non3
      REAL(8) :: dsdx1,dsdx2,dsdx3
      REAL(8) :: ds_non
!.....Local arrays
      REAL(8) :: ff_non(ncell_fluid,ndim)
      REAL(8) :: ssdx(ncell_fluid,ndim)
!.....Local vector arrays
      REAL(8) :: sdx(ndim)
      REAL(8) :: dx_non(nf_non,ndim)
      REAL(8) :: dx_sym(nf_sym,ndim)
      REAL(8) :: f_non(nf_non,ndim)
!      
      DATA INITIAL/.TRUE./
!
      eps=1.e-5
!        
      IF(INITIAL)THEN
         IF(ndim.eq.2) CALL grad_press_ls_c_2d
         IF(ndim.eq.3) CALL grad_press_ls_c_3d
         INITIAL=.FALSE.
      ENDIF
!
      IF(ndim.eq.2) THEN
!
!........Cells non
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF((ngrad(ii).gt.0 .and. lsindex(ii).ne.0) .or.  &
               (ngrad(kk).gt.0 .and. lsindex(kk).ne.0)) THEN
               ds_non=s(kk)-s(ii)
               dx_non(i,1)=xloc(kk,1)-xloc(ii,1)
               dx_non(i,2)=xloc(kk,2)-xloc(ii,2)
               f_non(i,1)=dx_non(i,1)*ds_non
               f_non(i,2)=dx_non(i,2)*ds_non
            ENDIF
         ENDDO
!
!........Cells sym
!
         nf_number=8
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            IF(ngrad(ii).gt.0 .and. lsindex(ii).ne.0) THEN
               xxc1=dxfc_nf(i1,1)
               xxc2=dxfc_nf(i1,2)
               contra=xxc1*xn_nf(i1,1)+xxc2*xn_nf(i1,2)
! ==>warning no runvv case to go through this path not tested
!     xxc = xfcx-xloc  -xxc+xfc = xloc
!              xloc21=-xxc1-2.d0*contra*xn_nf(i1,1)+xfcx_nf(i1)
!              xloc22=-xxc2-2.d0*contra*xn_nf(i1,2)+xfcy_nf(i1)
!              dx_sym(i,1)=xloc21-xloc(ii,1)
!              dx_sym(i,2)=xloc22-xloc(ii,2)
               xloc21=-2.d0*contra*xn_nf(i1,1)
               xloc22=-2.d0*contra*xn_nf(i1,2)
               dx_sym(i,1)=xloc21
               dx_sym(i,2)=xloc22
            ENDIF
         ENDDO
      ELSE
!
!........Cells non
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF((ngrad(ii).gt.0 .and. lsindex(ii).ne.0) .or.  &
               (ngrad(kk).gt.0 .and. lsindex(kk).ne.0)) THEN
               ds_non=s(kk)-s(ii)
               dx_non(i,1)=xloc(kk,1)-xloc(ii,1)
               dx_non(i,2)=xloc(kk,2)-xloc(ii,2)
               dx_non(i,3)=xloc(kk,3)-xloc(ii,3)
               f_non(i,1)=dx_non(i,1)*ds_non
               f_non(i,2)=dx_non(i,2)*ds_non
               f_non(i,3)=dx_non(i,3)*ds_non
            ENDIF
         ENDDO
!
!........Cells sym
!
         nf_number=8
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            IF((ngrad(ii).gt.0 .and. lsindex(ii).ne.0)) THEN
               xxc1=dxfc_nf(i1,1)
               xxc2=dxfc_nf(i1,2)
               xxc3=dxfc_nf(i1,3)
               contra=xxc1*xn_nf(i1,1)+xxc2*xn_nf(i1,2)+xxc3*xn_nf(i1,3)
! ==>warning no runvv case to go through this path not tested
!     xxc = xfcx-xloc  -xxc+xfc = xloc
!              xloc21=-xxc1-2.d0*contra*xn_nf(i1,1)+xfcx_nf(i1)
!              xloc22=-xxc2-2.d0*contra*xn_nf(i1,2)+xfcy_nf(i1)
!              xloc23=-xxc3-2.d0*contra*xn_nf(i1,3)+xfcz_nf(i1)
!              dx_sym(i,1)=xloc21-xloc(ii,1)
!              dx_sym(i,2)=xloc22-xloc(ii,2)
!              dx_sym(i,3)=xloc23-xloc(ii,3)
               xloc21=-2.d0*contra*xn_nf(i1,1)
               xloc22=-2.d0*contra*xn_nf(i1,2)
               xloc23=-2.d0*contra*xn_nf(i1,3)
               dx_sym(i,1)=xloc21
               dx_sym(i,2)=xloc22
               dx_sym(i,3)=xloc23
            ENDIF
         ENDDO
      ENDIF
!         
!.....valve model         
!      
      IF(rv_valve.eq.1) CALL valve_model_grad_press_ls(dx_non,f_non)
!
      CALL grad_press_ls_sum(eps,f_non,dx_non,dx_sym,ff_non,ssdx)
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            IF((ngrad(i).gt.0 .and. lsindex(i).ne.0)) THEN
               IF(det_2(i).eq.1.0) THEN
                  f_non1=ff_non(i,1)
                  f_non2=ff_non(i,2)
                  dsdx(i,1)=(a11_2(i)*f_non1+a12_2(i)*f_non2)/a_2(i)
                  dsdx(i,2)=(a12_2(i)*f_non1+a22_2(i)*f_non2)/a_2(i)
               ELSE
                  dsdx1=ff_non(i,1)
                  dsdx2=ff_non(i,2)
                  dsdx(i,1)=dsdx1
                  dsdx(i,2)=dsdx2
               ENDIF
            ENDIF
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            IF((ngrad(i).gt.0 .and. lsindex(i).ne.0)) THEN
               IF(det_3(i).eq.1.0) THEN
                  f_non1=ff_non(i,1)
                  f_non2=ff_non(i,2)
                  f_non3=ff_non(i,3)
                  dsdx(i,1)=(a11_3(i)*f_non1+a12_3(i)*f_non2+a13_3(i)*f_non3)/a_3(i)
                  dsdx(i,2)=(a12_3(i)*f_non1+a22_3(i)*f_non2+a23_3(i)*f_non3)/a_3(i)
                  dsdx(i,3)=(a13_3(i)*f_non1+a23_3(i)*f_non2+a33_3(i)*f_non3)/a_3(i)
               ELSE
                  dsdx1=ff_non(i,1)
                  dsdx2=ff_non(i,2)
                  dsdx3=ff_non(i,3)
                  dsdx(i,1)=dsdx1
                  dsdx(i,2)=dsdx2
                  dsdx(i,3)=dsdx3
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!        
      DO i=1,ncell_fluid
!
!........Least-square method only available when ngrad(i)>0
!      
         IF(ngrad(i).gt.0)THEN
!
!...........Avoid Least-square method at vertical walls except bottom and top wall when lsindex(i)>=0
!           
            IF(lsindex(i).eq.0) CYCLE !e.g. gori_sbo,rocom9to1,SMALL-POOL,SMALL-POOL-3D
            IF(ndim.eq.2) det_tmp=det_2(i)
            IF(ndim.eq.3) det_tmp=det_3(i)
!            
            IF(det_tmp.ne.1.0)THEN
               DO ix=1,ndim
                  sdx(ix)=ssdx(i,ix)
                  IF(sdx(ix).gt.0.)THEN
                     dsdx(i,ix)=dsdx(i,ix)/sdx(ix)
                  ELSE
                     WRITE(*,*)'### No neighbor cells for least square fitting ###,i=',i
                     WRITE(*,*)'### myrank,global cell=',myrank,jperm(i)
                     WRITE(*,*)'### neigh,neighbor:',i_neigh(i+1)-i_neigh(i),(neigh(j),j=i_neigh(i),i_neigh(i+1)-1)
                     STOP
                  ENDIF
               ENDDO
!               
            ENDIF
!            
         ENDIF
      ENDDO
!
      END SUBROUTINE grad_press_ls

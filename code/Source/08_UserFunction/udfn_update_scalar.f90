!
      SUBROUTINE udfn_update_scalar(ncvgd)
!
!     Define Repeat Criteria for dt_opts
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np
      USE Znum_cell     , ONLY: right_nb_k,istart_nf,                     &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_fluxk1
      USE Zdel_scalar   , ONLY: del_el,del_eg,del_x
      USE Zscalar_coeff , ONLY: sb
      USE Zconst2       , ONLY: dt
      USE Ztimecon      , ONLY: dt_opt,alpha_min
      USE Zbc_index     , ONLY: npb
      USE Zconst1       , ONLY: vv_prob
      USE Zscalar_coeff , ONLY: sfg_nf,sfl_nf,sfd_nf,         &
                                sfg_non_k,sfl_non_k,sfd_non_k
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      INTEGER :: ncvgd(3)
!     REAL(8) eps_vol
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: eng_g,eng_l,             &
                                        alpha_g,alpha_l,alpha_d, &
                                        qu
!.....Local vector arrays
      REAL(8),DIMENSION(nf_fluxk1) :: eng_g_nf,eng_l_nf,alpha_g_nf,alpha_d_nf,qu_nf
!
!     DATA eps_vol/1.d-3/  
!
!.....PAFS-POOL, fluidic_device  
!
!      IF(vv_prob.eq.'PAFS-POOL'.or.vv_prob.eq.'fluidic_device' .or. vv_prob.eq.'DIVA-NEW')THEN 
      IF(vv_prob.eq.'PAFS-POOL'      .or. &
         vv_prob.eq.'fluidic_device' .or. &
         vv_prob.eq.'DIVA-NEW'       .or. &
         vv_prob.eq.'ST2-CT-01'      .or. &
         vv_prob.eq.'ST2-CT-02'      .or. &
         vv_prob.eq.'ST2-CT-03'            )THEN
         ncvgd(:)=0
!
!........Build summation info for non,mcc,inl
!
         nf_number_nb=2
         nf_number_id(-1)=-1
         nf_number_id(0)=0
         nf_number_id(1)=1
         nf_number_id(2)=2
         istart_nfs(-1)=0
         istart_nfs(0)=istart_nfs(-1)+nf_nonk
         istart_nfs(1)=istart_nfs(0) +nf_non
         istart_nfs(2)=istart_nfs(1) +nf_mcc
         lens         =istart_nfs(2) +nf_inl
!
!........Computing Cell
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
            IF(npb(ii).eq.0) THEN
               eng_g_nf(i0)  =sfg_nf(i1,1)*flux_g_nf(i1)+sfl_nf(i1,1)*flux_l_nf(i1)+sfd_nf(i1,1)*flux_d_nf(i1)
               eng_l_nf(i0)  =sfg_nf(i1,2)*flux_g_nf(i1)+sfl_nf(i1,2)*flux_l_nf(i1)+sfd_nf(i1,2)*flux_d_nf(i1)
               alpha_g_nf(i0)=sfg_nf(i1,3)*flux_g_nf(i1)+sfl_nf(i1,3)*flux_l_nf(i1)+sfd_nf(i1,3)*flux_d_nf(i1)
               alpha_d_nf(i0)=sfg_nf(i1,4)*flux_g_nf(i1)+sfl_nf(i1,4)*flux_l_nf(i1)+sfd_nf(i1,4)*flux_d_nf(i1)
               qu_nf(i0)     =sfg_nf(i1,5)*flux_g_nf(i1)+sfl_nf(i1,5)*flux_l_nf(i1)+sfd_nf(i1,5)*flux_d_nf(i1)
            ELSE 
               eng_g_nf(i0)  =0.d0
               eng_l_nf(i0)  =0.d0
               alpha_g_nf(i0)=0.d0
               alpha_d_nf(i0)=0.d0
               qu_nf(i0)     =0.d0
            ENDIF
         ENDDO
!           
         nv=-1 
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            IF(npb(ii).eq.0) THEN
               eng_g_nf(i)  =-(sfg_non_k(i,1)*flux_g_nf(k)+sfl_non_k(i,1)*flux_l_nf(k)+sfd_non_k(i,1)*flux_d_nf(k))
               eng_l_nf(i)  =-(sfg_non_k(i,2)*flux_g_nf(k)+sfl_non_k(i,2)*flux_l_nf(k)+sfd_non_k(i,2)*flux_d_nf(k))
               alpha_g_nf(i)=-(sfg_non_k(i,3)*flux_g_nf(k)+sfl_non_k(i,3)*flux_l_nf(k)+sfd_non_k(i,3)*flux_d_nf(k))
               alpha_d_nf(i)=-(sfg_non_k(i,4)*flux_g_nf(k)+sfl_non_k(i,4)*flux_l_nf(k)+sfd_non_k(i,4)*flux_d_nf(k))
               qu_nf(i)     =-(sfg_non_k(i,5)*flux_g_nf(k)+sfl_non_k(i,5)*flux_l_nf(k)+sfd_non_k(i,5)*flux_d_nf(k))
            ELSE 
               eng_g_nf(i)  =0.d0
               eng_l_nf(i)  =0.d0
               alpha_g_nf(i)=0.d0
               alpha_d_nf(i)=0.d0
               qu_nf(i)     =0.d0
            ENDIF
         ENDDO
!
!........Cells mcc,inl
!
         DO nv=1,2
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               IF(npb(ii).eq.0) THEN
                  eng_g_nf(i0)  =sfg_nf(i1,1)*flux_g_nf(i1)+sfl_nf(i1,1)*flux_l_nf(i1)+sfd_nf(i1,1)*flux_d_nf(i1)
                  eng_l_nf(i0)  =sfg_nf(i1,2)*flux_g_nf(i1)+sfl_nf(i1,2)*flux_l_nf(i1)+sfd_nf(i1,2)*flux_d_nf(i1)
                  alpha_g_nf(i0)=sfg_nf(i1,3)*flux_g_nf(i1)+sfl_nf(i1,3)*flux_l_nf(i1)+sfd_nf(i1,3)*flux_d_nf(i1)
                  alpha_d_nf(i0)=sfg_nf(i1,4)*flux_g_nf(i1)+sfl_nf(i1,4)*flux_l_nf(i1)+sfd_nf(i1,4)*flux_d_nf(i1)
                  qu_nf(i0)     =sfg_nf(i1,5)*flux_g_nf(i1)+sfl_nf(i1,5)*flux_l_nf(i1)+sfd_nf(i1,5)*flux_d_nf(i1)
               ELSE 
                  eng_g_nf(i0)  =0.d0
                  eng_l_nf(i0)  =0.d0
                  alpha_g_nf(i0)=0.d0
                  alpha_d_nf(i0)=0.d0
                  qu_nf(i0)     =0.d0
               ENDIF
            ENDDO
         ENDDO
         CALL sum_nf(0,0,                &
                     eng_g_nf  ,eng_g,   &
                     eng_l_nf  ,eng_l,   &
                     alpha_g_nf,alpha_g, &
                     alpha_d_nf,alpha_d, &
                     qu_nf     ,qu)
!
         DO i=1,ncell_fluid
            IF(npb(i).gt.0) cycle
            eng_g(i)  =eng_g(i)+sb(i,1)
            eng_l(i)  =eng_l(i)+sb(i,2)
            alpha_g(i)=alpha_g(i)+cell%alphag_o(i)+sb(i,3)
            alpha_d(i)=alpha_d(i)+cell%alphad_o(i)+sb(i,4)
            qu(i)     =qu(i)+sb(i,5)
!
            eng_g(i)=cell%eg_o(i)+eng_g(i)
            eng_l(i)=cell%el_o(i)+eng_l(i)        
            alpha_g(i)=MIN(MAX(0.d0,alpha_g(i)),1.d0)   
            alpha_d(i)=MIN(MAX(0.d0,alpha_d(i)),1.d-8)
            alpha_d(i)=0.d0
            alpha_l(i)=MIN(MAX(0.d0,1.d0-alpha_g(i)-alpha_d(i)),1.d0)
            qu(i)=cell%quala_o(i)+qu(i)
         ENDDO
!            
         IF (vv_prob.eq.'fluidic_device') THEN
            DO i=1,ncell_fluid
               IF(npb(i).gt.0) cycle
               IF(qu(i).lt.0.d0)THEN
                  eng_g(i)=cell%eg_o(i)
                  qu(i)=MAX(0.d0,qu(i))
               ENDIF
               IF(qu(i).gt.1.d0)THEN
                  eng_g(i)=cell%eg_o(i)
                  qu(i)=cell%quala_o(i)
               ENDIF   
               IF(alpha_l(i).lt.0.d0)THEN
                  alpha_l(i)=1.d-8
                  alpha_g(i)=1.d0-cell%alphad(i)
                  eng_l(i)=cell%el_o(i)
               ENDIF 
               IF(alpha_l(i).le.0.d0)THEN
                  alpha_l(i)=0.d0
                  alpha_g(i)=1.d0
                  eng_l(i)=cell%el_o(i)
               ENDIF             
               IF(eng_l(i).le.0.d0)THEN                                  
                  eng_l(i)=cell%el_o(i)
               ENDIF
               IF(alpha_g(i).le.1.d-8)THEN
                  eng_g(i)=cell%eg_o(i)
               ENDIF             
               IF(alpha_l(i).le.1.d-8)THEN
                  eng_l(i)=cell%el_o(i)
               ENDIF                      
!
!..............Limit energy change less than eps*old_value
!
               IF(ABS(eng_l(i)-cell%el_o(i)).gt.0.1d0*cell%el_o(i) .or. &
                  ABS(eng_g(i)-cell%eg_o(i)).gt.0.1d0*cell%eg_o(i)) THEN
                  ncvgd(1)=1
               ENDIF
!
!...........limit void fraction error less than eps_vol
!
!           IF (alpha_g(i).lt.-eps_vol  .and. dt.ge.5.0e-4) ncvgd(2)=.TRUE.
!           IF ((alpha_g(i)-1.d0).gt.eps_vol .and. dt.ge.5.0e-4) ncvgd(2)=.TRUE.  
!           IF (eng_l(i)/cell(i)%elsat.gt.1.05d0 .and. cell(i)%regime.eq.13 .and. dt.ge.1.0e-4) ncvgd(3)=.TRUE.    
!
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               IF(npb(i).gt.0) cycle
               IF(qu(i).lt.0.d0)THEN
                  eng_g(i)=cell%egsat(i)
                  qu(i)=MAX(0.d0,qu(i))
               ENDIF
               IF(qu(i).gt.1.d0)THEN
                  eng_g(i)=cell%eg_o(i)
                  qu(i)=cell%quala_o(i)
               ENDIF   
               IF(alpha_l(i).lt.0.d0)THEN
                  alpha_l(i)=1.d-8
                  alpha_g(i)=1.d0-cell%alphad(i)
                  eng_l(i)=cell%el_o(i)
               ENDIF 
               IF(alpha_l(i).le.0.d0)THEN
                  alpha_l(i)=0.d0
                  alpha_g(i)=1.d0
                  eng_l(i)=cell%elsat(i)
               ENDIF             
               IF(eng_l(i).le.0.d0)THEN                                  
                  eng_l(i)=cell%elsat(i)
               ENDIF
!
!..............Limit energy change less than eps*old_value
!
               IF(ABS(eng_l(i)-cell%el_o(i)).gt.0.1d0*cell%el_o(i) .and. dt.ge.1.d-4) ncvgd(1)=1
               IF(ABS(eng_l(i)-cell%el_o(i)).gt.0.1d0*cell%el_o(i) .and. dt.le.1.d-4 .and.  (ncvgd(1).eq.0))THEN
                  IF(cell%alphal_o(i).lt.alpha_min .and. alpha_l(i).ge.alpha_min .and. dt.ge.1.d-6)THEN
                     ncvgd(1)=1
                  ENDIF
               ENDIF
!
!..............limit void fraction error less than eps_vol
!
!           IF (alpha_g(i).lt.-eps_vol  .and. dt.ge.5.0e-4) ncvgd(2)=.TRUE.
!           IF ((alpha_g(i)-1.d0).gt.eps_vol .and. dt.ge.5.0e-4) ncvgd(2)=.TRUE.  
!           IF (eng_l(i)/cell(i)%elsat.gt.1.05d0 .and. cell(i)%regime.eq.13 .and. dt.ge.1.0e-4) ncvgd(3)=.TRUE.    
!
            ENDDO
         ENDIF

         IF(np.gt.1) CALL allreducei_i(ncvgd,3)
         IF(((ncvgd(1).gt.0.or.ncvgd(2).gt.0) .and. dt_opt.eq.2) .or. &
            ( ncvgd(3).gt.0 .and. dt_opt.eq.3)                   .or. & 
            ( ncvgd(1).gt.0 .and. dt_opt.eq.4)) THEN
            RETURN
         ELSE
            IF (vv_prob.eq.'fluidic_device') THEN
               DO i=1,ncell_fluid
                  IF(npb(i).gt.0) cycle
                  cell%eg(i)=eng_g(i)
                  cell%el(i)=eng_l(i)
                  cell%alphag(i)=alpha_g(i)
                  cell%alphal(i)=alpha_l(i)
                  cell%alphad(i)=alpha_d(i)
                  cell%quala(i)=qu(i)
                  IF(cell%quala(i).lt.0.d0)THEN
!                    cell%eg(i)=cell%egsat(i)
                     cell%eg(i)=cell%eg_o(i)
                     cell%quala(i)=MAX(0.d0,cell%quala(i))
                  ENDIF
                  IF(cell%quala(i).gt.1.d0)THEN
                     cell%eg(i)=cell%eg_o(i)
                     cell%quala(i)=cell%quala_o(i)
                  ENDIF   
                  IF(cell%alphal(i).lt.0.d0)THEN
                     cell%alphal(i)=1.d-8
                     cell%alphag(i)=1.d0-cell%alphad(i)
                     cell%el(i)=cell%el_o(i)
                  ENDIF 
                  IF(cell%alphal(i).le.0.d0)THEN
                     cell%alphal(i)=0.d0
                     cell%alphag(i)=1.d0
                     cell%el(i)=cell%elsat(i)
                     cell%eg(i)=cell%eg_o(i)
                  ENDIF             
                  IF(cell%el(i).le.0.d0)THEN 
!                    cell%el(i)=cell%elsat(i)
                     cell%el(i)=cell%el_o(i)
                  ENDIF
               ENDDO
            ELSE
               DO i=1,ncell_fluid
                  IF(npb(i).gt.0) cycle
                     cell%eg(i)=eng_g(i)
                     cell%el(i)=eng_l(i)
                     cell%alphag(i)=alpha_g(i)
                     cell%alphal(i)=alpha_l(i)
                     cell%alphad(i)=alpha_d(i)
                     cell%quala(i)=qu(i)
                     IF(cell%quala(i).lt.0.d0)THEN
                        cell%eg(i)=cell%egsat(i)
                        cell%quala(i)=MAX(0.d0,cell%quala(i))
                     ENDIF
                     IF(cell%quala(i).gt.1.d0)THEN
                        cell%eg(i)=cell%eg_o(i)
                        cell%quala(i)=cell%quala_o(i)
                     ENDIF   
                     IF(cell%alphal(i).lt.0.d0)THEN
                        cell%alphal(i)=1.d-8
                        cell%alphag(i)=1.d0-cell%alphad(i)
                        cell%el(i)=cell%el_o(i)
                     ENDIF 
                     IF(cell%alphal(i).le.0.d0)THEN
                        cell%alphal(i)=0.d0
                        cell%alphag(i)=1.d0
                        cell%el(i)=cell%elsat(i)
                     ENDIF             
                     IF(cell%el(i).le.0.d0)THEN
                        cell%el(i)=cell%elsat(i)
                     ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!.....Update material properties
!
      CALL property_calc(1)
!
!.....Vapor generation rate
!
      CALL set_vapor_generation
!
!.....Save scalar changes
!
      DO i=1,ncell_fluid
         del_eg(i)=cell%eg(i)-cell%eg_o(i)
         del_el(i)=cell%el(i)-cell%el_o(i)
         del_x(i) =cell%quala(i)-cell%quala_o(i)
      ENDDO
!
      END SUBROUTINE udfn_update_scalar

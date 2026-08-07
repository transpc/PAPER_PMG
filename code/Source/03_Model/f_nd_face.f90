!
      SUBROUTINE f_nd_face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non,jneigh_nf
      USE Zconst1      , ONLY: nwlf,ntdf,nlift
      USE Zndforce     , ONLY: clift,f_wl,ctd
      USE Zvector      , ONLY: vg_o,vl_o
      USE Zvector      , ONLY: vl_f_non,vg_f_non
      USE Zvec_geo     , ONLY: fac_non,fac1_non,xn_nf,djir_non
!
      IMPLICIT NONE  
!
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,jj,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: d31,d13,d32,d23,d12,d21
      REAL(8) :: dal,tdf,car,vl1,vl2,vg1,vg2
!.....Local vector arrays
      REAL(8) :: ff_non(nf_non,ndim)
!
!.....Wall Lubrication Force
!
      IF(nwlf.ne.-1)THEN
!
         IF(np.gt.1) CALL communicate_2d(F_wl)
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         IF(ndim.eq.2) THEN
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               ff_non(i,1)=(F_wl(ii,1)*fac1_non(i)+F_wl(kk,1)*fac_non(i))*0.5d0
               ff_non(i,2)=(F_wl(ii,2)*fac1_non(i)+F_wl(kk,2)*fac_non(i))*0.5d0
               vg_f_non(i,1)=vg_f_non(i,1)-ff_non(i,1)
               vg_f_non(i,2)=vg_f_non(i,2)-ff_non(i,2)
               vl_f_non(i,1)=vl_f_non(i,1)+ff_non(i,1)
               vl_f_non(i,2)=vl_f_non(i,2)+ff_non(i,2)
            ENDDO      
         ELSE
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               ff_non(i,1)=(F_wl(ii,1)*fac1_non(i)+F_wl(kk,1)*fac_non(i))*0.5d0
               ff_non(i,2)=(F_wl(ii,2)*fac1_non(i)+F_wl(kk,2)*fac_non(i))*0.5d0
               ff_non(i,3)=(F_wl(ii,3)*fac1_non(i)+F_wl(kk,3)*fac_non(i))*0.5d0
               vg_f_non(i,1)=vg_f_non(i,1)-ff_non(i,1)
               vg_f_non(i,2)=vg_f_non(i,2)-ff_non(i,2)
               vg_f_non(i,3)=vg_f_non(i,3)-ff_non(i,3)
               vl_f_non(i,1)=vl_f_non(i,1)+ff_non(i,1)
               vl_f_non(i,2)=vl_f_non(i,2)+ff_non(i,2)
               vl_f_non(i,3)=vl_f_non(i,3)+ff_non(i,3)
            ENDDO      
         ENDIF
!
      ENDIF
!
!.....Turbulent Dispersion Force
!
      IF(ntdf.ne.-1)THEN
!
         IF(np.gt.1) CALL communicate_1d(ctd)
!
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         IF(ndim.eq.2) THEN
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               jj=jneigh_nf(i1)
               kk=right_non(i)
               dal=cell%alphal_o(kk)-cell%alphal_o(ii)
               tdf=0.5d0*(ctd(ii)+ctd(kk))
               ff_non(i,1)=tdf*(dal*djir_non(i)*xn_nf(i1,1))
               ff_non(i,2)=tdf*(dal*djir_non(i)*xn_nf(i1,2))
               vg_f_non(i,1)=vg_f_non(i,1)-ff_non(i,1)
               vg_f_non(i,2)=vg_f_non(i,2)-ff_non(i,2)
               vl_f_non(i,1)=vl_f_non(i,1)+ff_non(i,1)
               vl_f_non(i,2)=vl_f_non(i,2)+ff_non(i,2)
            ENDDO      
         ELSE
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               jj=jneigh_nf(i1)
               kk=right_non(i)
               dal=cell%alphal_o(kk)-cell%alphal_o(ii)
               tdf=0.5d0*(ctd(ii)+ctd(kk))
               ff_non(i,1)=tdf*(dal*djir_non(i)*xn_nf(i1,1))
               ff_non(i,2)=tdf*(dal*djir_non(i)*xn_nf(i1,2))
               ff_non(i,3)=tdf*(dal*djir_non(i)*xn_nf(i1,3))
               vg_f_non(i,1)=vg_f_non(i,1)-ff_non(i,1)
               vg_f_non(i,2)=vg_f_non(i,2)-ff_non(i,2)
               vg_f_non(i,3)=vg_f_non(i,3)-ff_non(i,3)
               vl_f_non(i,1)=vl_f_non(i,1)+ff_non(i,1)
               vl_f_non(i,2)=vl_f_non(i,2)+ff_non(i,2)
               vl_f_non(i,3)=vl_f_non(i,3)+ff_non(i,3)
            ENDDO      
         ENDIF
!
      ENDIF
!
!.....Lift Force
!
      IF(nlift.ne.-1)THEN
!
         IF(np.gt.1) CALL communicate_1d(Clift,cell%rhol)
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         IF(ndim.eq.2) THEN
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               jj=jneigh_nf(i1)
               kk=right_non(i)
   !
                     Car=0.5d0*(Clift(ii)*cell%alphag(ii)*cell%rhol(ii)+Clift(kk)*cell%alphag(kk)*cell%rhol(kk))
   !
                     d21=(vl_o(kk,2)-vl_o(ii,2))*djir_non(i)*xn_nf(i1,1)
                     d12=(vl_o(kk,1)-vl_o(ii,1))*djir_non(i)*xn_nf(i1,2)
   !
                     vg1=fac1_non(i)*vg_o(ii,2)+fac_non(i)*vg_o(kk,2)
                     vl1=fac1_non(i)*vl_o(ii,2)+fac_non(i)*vl_o(kk,2)
                     ff_non(i,1)=Car*(vg1-vl1)*(d21-d12)
   !
                     vg1=fac1_non(i)*vg_o(ii,1)+fac_non(i)*vg_o(kk,1)
                     vl1=fac1_non(i)*vl_o(ii,1)+fac_non(i)*vl_o(kk,1)
                     ff_non(i,2)=-Car*(vg1-vl1)*(d21-d12)
                     vg_f_non(i,1)=vg_f_non(i,1)-ff_non(i,1)
                     vg_f_non(i,2)=vg_f_non(i,2)-ff_non(i,2)
                     vl_f_non(i,1)=vl_f_non(i,1)+ff_non(i,1)
                     vl_f_non(i,2)=vl_f_non(i,2)+ff_non(i,2)
            ENDDO      
         ELSE
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               jj=jneigh_nf(i1)
               kk=right_non(i)
   !
                     Car=0.5d0*(Clift(ii)*cell%alphag(ii)*cell%rhol(ii)+Clift(kk)*cell%alphag(kk)*cell%rhol(kk))
   !
                     d21=(vl_o(kk,2)-vl_o(ii,2))*djir_non(i)*xn_nf(i1,1)
                     d12=(vl_o(kk,1)-vl_o(ii,1))*djir_non(i)*xn_nf(i1,2)
   !
                     d31=(vl_o(kk,3)-vl_o(ii,3))*djir_non(i)*xn_nf(i1,1)
                     d32=(vl_o(kk,3)-vl_o(ii,3))*djir_non(i)*xn_nf(i1,2)
                     d13=(vl_o(kk,1)-vl_o(ii,1))*djir_non(i)*xn_nf(i1,3)
                     d23=(vl_o(kk,2)-vl_o(ii,2))*djir_non(i)*xn_nf(i1,3)
   !
                     vg1=fac1_non(i)*vg_o(ii,2)+fac_non(i)*vg_o(kk,2)
                     vl1=fac1_non(i)*vl_o(ii,2)+fac_non(i)*vl_o(kk,2)
                     ff_non(i,1)=Car*(vg1-vl1)*(d21-d12)
   !
                     vg1=fac1_non(i)*vg_o(ii,1)+fac_non(i)*vg_o(kk,1)
                     vl1=fac1_non(i)*vl_o(ii,1)+fac_non(i)*vl_o(kk,1)
                     ff_non(i,2)=-Car*(vg1-vl1)*(d21-d12)
   !
                        vg2=fac1_non(i)*vg_o(ii,3)+fac_non(i)*vg_o(kk,3)
                        vl2=fac1_non(i)*vl_o(ii,3)+fac_non(i)*vl_o(kk,3)
                        ff_non(i,1)=ff_non(i,1)+Car*(vg2-vl2)*(d31-d13)
                        vg2=fac1_non(i)*vg_o(ii,3)+fac_non(i)*vg_o(kk,3)
                        vl2=fac1_non(i)*vl_o(ii,3)+fac_non(i)*vl_o(kk,3)
                        ff_non(i,2)=ff_non(i,2)+Car*(vg2-vl2)*(d32-d23)
                        vg1=fac1_non(i)*vg_o(ii,1)+fac_non(i)*vg_o(kk,1)
                        vl1=fac1_non(i)*vl_o(ii,1)+fac_non(i)*vl_o(kk,1)
                        vg2=fac1_non(i)*vg_o(ii,2)+fac_non(i)*vg_o(kk,2)
                        vl2=fac1_non(i)*vl_o(ii,2)+fac_non(i)*vl_o(kk,2)
                        ff_non(i,3)=-Car*(vg1-vl1)*(d31-d13)-Car*(vg2-vl2)*(d32-d23)
                     vg_f_non(i,1)=vg_f_non(i,1)-ff_non(i,1)
                     vg_f_non(i,2)=vg_f_non(i,2)-ff_non(i,2)
                     vg_f_non(i,3)=vg_f_non(i,3)-ff_non(i,3)
                     vl_f_non(i,1)=vl_f_non(i,1)+ff_non(i,1)
                     vl_f_non(i,2)=vl_f_non(i,2)+ff_non(i,2)
                     vl_f_non(i,3)=vl_f_non(i,3)+ff_non(i,3)
            ENDDO      
         ENDIF
!
      ENDIF
!
      ENDSUBROUTINE f_nd_face

!
      SUBROUTINE choke_massflowrate1(ml_nf,mg_nf,md_nf)
      
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_flux
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zvec_geo      , ONLY: f1,f0
      USE Zporous       , ONLY: vfporous
      USE Zare          , ONLY: ar_gas,ar_liq,ar_drp
      USE Znum_cell     , ONLY: istart_nf
      USE Zcoord3       , ONLY: porosity
      USE Zvector       , ONLY: vl_n,vg_n,vd_n
!
      IMPLICIT NONE
!      
      INTEGER :: nf_number,istart,len,i,ii,kk,i1
      REAL(8) :: a,b
      REAL(8),DIMENSION(nf_flux,ndim) :: ml_nf,mg_nf,md_nf

      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(vfporous.eq.0) THEN
         IF(ndim.eq.2)THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)
               b=f0(i)
  
               ml_nf(i1,1)=a*vl_n(ii,1)*ar_liq(ii)+b*vl_n(kk,1)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,1)=a*vg_n(ii,1)*ar_gas(ii)+b*vg_n(kk,1)*ar_gas(kk)
               md_nf(i1,1)=a*vd_n(ii,1)*ar_drp(ii)+b*vd_n(kk,1)*ar_drp(kk)    
               ml_nf(i1,2)=a*vl_n(ii,2)*ar_liq(ii)+b*vl_n(kk,2)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,2)=a*vg_n(ii,2)*ar_gas(ii)+b*vg_n(kk,2)*ar_gas(kk)
               md_nf(i1,2)=a*vd_n(ii,2)*ar_drp(ii)+b*vd_n(kk,2)*ar_drp(kk)
!               
            ENDDO
         ELSE
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)
               b=f0(i)

               ml_nf(i1,1)=a*vl_n(ii,1)*ar_liq(ii)+b*vl_n(kk,1)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,1)=a*vg_n(ii,1)*ar_gas(ii)+b*vg_n(kk,1)*ar_gas(kk)
               md_nf(i1,1)=a*vd_n(ii,1)*ar_drp(ii)+b*vd_n(kk,1)*ar_drp(kk)    
               ml_nf(i1,2)=a*vl_n(ii,2)*ar_liq(ii)+b*vl_n(kk,2)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,2)=a*vg_n(ii,2)*ar_gas(ii)+b*vg_n(kk,2)*ar_gas(kk)
               md_nf(i1,2)=a*vd_n(ii,2)*ar_drp(ii)+b*vd_n(kk,2)*ar_drp(kk)
               ml_nf(i1,3)=a*vl_n(ii,3)*ar_liq(ii)+b*vl_n(kk,3)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,3)=a*vg_n(ii,3)*ar_gas(ii)+b*vg_n(kk,3)*ar_gas(kk)
               md_nf(i1,3)=a*vd_n(ii,3)*ar_drp(ii)+b*vd_n(kk,3)*ar_drp(kk)               
!               
            ENDDO
         ENDIF
      ELSE
         IF(ndim.eq.2)THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)*porosity(ii)
               b=f0(i)*porosity(kk)
 
               ml_nf(i1,1)=a*vl_n(ii,1)*ar_liq(ii)+b*vl_n(kk,1)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,1)=a*vg_n(ii,1)*ar_gas(ii)+b*vg_n(kk,1)*ar_gas(kk)
               md_nf(i1,1)=a*vd_n(ii,1)*ar_drp(ii)+b*vd_n(kk,1)*ar_drp(kk)    
               ml_nf(i1,2)=a*vl_n(ii,2)*ar_liq(ii)+b*vl_n(kk,2)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,2)=a*vg_n(ii,2)*ar_gas(ii)+b*vg_n(kk,2)*ar_gas(kk)
               md_nf(i1,2)=a*vd_n(ii,2)*ar_drp(ii)+b*vd_n(kk,2)*ar_drp(kk)
!                   
            ENDDO
         ELSE
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)*porosity(ii)
               b=f0(i)*porosity(kk)

               ml_nf(i1,1)=a*vl_n(ii,1)*ar_liq(ii)+b*vl_n(kk,1)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,1)=a*vg_n(ii,1)*ar_gas(ii)+b*vg_n(kk,1)*ar_gas(kk)
               md_nf(i1,1)=a*vd_n(ii,1)*ar_drp(ii)+b*vd_n(kk,1)*ar_drp(kk)    
               ml_nf(i1,2)=a*vl_n(ii,2)*ar_liq(ii)+b*vl_n(kk,2)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,2)=a*vg_n(ii,2)*ar_gas(ii)+b*vg_n(kk,2)*ar_gas(kk)
               md_nf(i1,2)=a*vd_n(ii,2)*ar_drp(ii)+b*vd_n(kk,2)*ar_drp(kk)
               ml_nf(i1,3)=a*vl_n(ii,3)*ar_liq(ii)+b*vl_n(kk,3)*ar_liq(kk)  !kg/s/m2
               mg_nf(i1,3)=a*vg_n(ii,3)*ar_gas(ii)+b*vg_n(kk,3)*ar_gas(kk)
               md_nf(i1,3)=a*vd_n(ii,3)*ar_drp(ii)+b*vd_n(kk,3)*ar_drp(kk)    
!                     
            ENDDO
         ENDIF
      ENDIF     
!            
      END SUBROUTINE choke_massflowrate1    
!
      SUBROUTINE choke_massflowrate2(ml_nf,mg_nf,md_nf,mflux_l_nf,mflux_g_nf,mflux_d_nf)
      
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_flux
      USE Zvec_geo      , ONLY: svp_nf,sa_nf
      USE Znum_cell     , ONLY: istart_nf
!
      IMPLICIT NONE
      
      INTEGER :: nf_number,istart,len,i,i1
      REAL(8) :: mli1,mli2,mli3
      REAL(8) :: mgi1,mgi2,mgi3
      REAL(8) :: mdi1,mdi2,mdi3  
      REAL(8),DIMENSION(nf_flux,ndim) :: ml_nf,mg_nf,md_nf
      REAL(8),DIMENSION(nf_flux) :: mflux_l_nf,mflux_g_nf,mflux_d_nf

      IF(ndim.eq.2)THEN
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            mli1=ml_nf(i1,1)
            mli2=ml_nf(i1,2)
            mgi1=mg_nf(i1,1)
            mgi2=mg_nf(i1,2)
            mdi1=md_nf(i1,1)
            mdi2=md_nf(i1,2)               
            mflux_l_nf(i1)=(mli1*svp_nf(i1,1)+mli2*svp_nf(i1,2))/sa_nf(i1)  !kg/s/m2
            mflux_g_nf(i1)=(mgi1*svp_nf(i1,1)+mgi2*svp_nf(i1,2))/sa_nf(i1)
            mflux_d_nf(i1)=(mdi1*svp_nf(i1,1)+mdi2*svp_nf(i1,2))/sa_nf(i1)    
         ENDDO
     ELSE
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            mli1=ml_nf(i1,1)
            mli2=ml_nf(i1,2)
            mli3=ml_nf(i1,3)
            mgi1=mg_nf(i1,1)
            mgi2=mg_nf(i1,2)
            mgi3=mg_nf(i1,3)
            mdi1=md_nf(i1,1)
            mdi2=md_nf(i1,2)
            mdi3=md_nf(i1,3)               
            mflux_l_nf(i1)=(mli1*svp_nf(i1,1)+mli2*svp_nf(i1,2)+mli3*svp_nf(i1,3))/sa_nf(i1) !kg/s/m2
            mflux_g_nf(i1)=(mgi1*svp_nf(i1,1)+mgi2*svp_nf(i1,2)+mgi3*svp_nf(i1,3))/sa_nf(i1)
            mflux_d_nf(i1)=(mdi1*svp_nf(i1,1)+mdi2*svp_nf(i1,2)+mdi3*svp_nf(i1,3))/sa_nf(i1)
         ENDDO
      ENDIF
!            
      END SUBROUTINE choke_massflowrate2    

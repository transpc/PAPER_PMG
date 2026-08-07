!
      SUBROUTINE udfn_model_overwrite
!
!.....Overwrite the coefficient of model by being called int_swap
!
      USE Vol_DATA
      USE Zconst1         , ONLY: vv_prob
      USE Zzone           , ONLY: nzone,ncell_fluid
      USE Zcoord3         , ONLY: porosity
!      
      IMPLICIT NONE
!
      INTEGER i    
      REAL(8) al_min,ag_min,xg1,xl1,y1,xg2,xl2,y2,x,y,slope_g,slope_l 
!
      IF(vv_prob.eq.'opr1000_rv_lbloca'   .or.&
         vv_prob.eq.'opr1000_mc_rv_lbloca' )THEN
         !cell%vfgl(i)=DMAX1(1.0d4,cell%vfgl(i))                       !effective  
         !cell%vfgl(i)=DMAX1(1.0d2,cell%vfgl(i))                       !not effective 
         !IF(nbcon_cell(i).ge.1)cell%vfgl(i)=DMAX1(1.0d4,cell%vfgl(i)) !not effective
         al_min=1.0d-2  !apr1400_lbloca_debug
         ag_min=1.0d-2
         xg1=ag_min
         xl1=al_min
         y1=0.0d0
         xg2=ag_min*0.1d0
         xl2=al_min*0.1d0
         y2=1.0d4
         slope_g=(y2-y1)/(xg2-xg1)
         slope_l=(y2-y1)/(xl2-xl1)
         DO i=1,ncell_fluid
            IF(cell%alphag(i).lt.ag_min)THEN !no effect
              x=cell%alphag(i)
              y=slope_g*(x-xg1) !+y1
              y=DMIN1(1.0d4,DMAX1(0.0d0,y))            
              cell%vfgl(i)=cell%vfgl(i)+y
            ELSEIF(cell%alphal(i).lt.al_min)THEN
              x=cell%alphal(i)
              y=slope_l*(x-xl1) !+y1
              y=DMIN1(1.0d4,DMAX1(0.0d0,y))            
              cell%vfgl(i)=cell%vfgl(i)+y
            ENDIF
         ENDDO
!         
         CALL average_spatially(cell%vfgl)
!         
      ELSEIF(vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN 
!      
        !See SUBROUTINE wall_drag & s_wall_fric
        DO i=1,ncell_fluid
           IF(porosity(i).lt.1.d0)THEN
               cell%vfwg(i)=cell%vfwg(i)*0.3d0 !darcy
               cell%vfwl(i)=cell%vfwl(i)*0.3d0 !darcy
           ELSE
               cell%vfwg(i)=cell%vfwg(i)*0.0d0
               cell%vfwl(i)=cell%vfwl(i)*0.0d0
           ENDIF
        ENDDO
! 
      ELSEIF(vv_prob.eq.'OPR1000_fullvessel_1x1'.or.vv_prob.eq.'opr1000_rv')THEN 
!      
        DO i=1,ncell_fluid
           IF(nzone(i).eq.6)THEN
              cell%vfwg(i)=cell%vfwg(i)*0.1d0 !darcy
              cell%vfwl(i)=cell%vfwl(i)*0.1d0 !darcy
           ELSE
              cell%vfwg(i)=cell%vfwg(i)*0.0d0
              cell%vfwl(i)=cell%vfwl(i)*0.0d0
           ENDIF
        ENDDO   
!
      ENDIF 

      RETURN
      END SUBROUTINE udfn_model_overwrite
!-------------------------------------------------------------------------
      SUBROUTINE average_spatially(ftmp) !siphon,apr1400_lbloca
!      
      USE Zinterface
      USE Zmpi          , ONLY: ncell_fp           
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: myrank,np
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zvec_param    , ONLY: nf_non
      USE Znum_cell     , ONLY: istart_nf, &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_geo      , ONLY: saa_nf,fac_non,fac1_non 
      USE Zio_unit      , ONLY: unit_log
!      
      IMPLICIT NONE
!      
!.....Input
      REAL(8) :: ftmp(ncell_fp)
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,len,i1
      LOGICAL,SAVE :: initial=.true.
!.....Local arrays
      REAL(8) :: ftmp1(ncell_fluid),saa1(ncell_fluid)
!.....Local vector arrays
      REAL(8) :: ftmp_non(nf_non)
!
      IF(initial)THEN
         initial=.false.
         IF(myrank.eq.0) WRITE(*,"(11x,a)") 'Communicating & spatially averaging int_drag in communicate_drag!'
         IF(myrank.eq.0) WRITE(unit_log,"(11x,a)")'Communicating & spatially averaging int_drag in communicate_drag!'
      ENDIF
      IF(np.gt.1) CALL communicate_1d(ftmp)
!
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=0
      nf_number_id(0)=0
      istart_nfs(0)=0
      lens         =istart_nfs(0)+nf_non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         ftmp_non(i)=(ftmp(ii)*fac1_non(i)+ftmp(kk)*fac_non(i))*saa_nf(i1)
      ENDDO
!      
      CALL sum_nf(0,1,            &
                  ftmp_non,ftmp1, &
                  saa_nf  ,saa1)
!      
      DO i=1,ncell_fluid
         IF(saa1(i).gt.0.d0) ftmp(i)=ftmp1(i)/saa1(i)
      ENDDO
!
      END SUBROUTINE average_spatially

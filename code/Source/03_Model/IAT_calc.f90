!
      SUBROUTINE IAT_calc
!
!     This routine calculates iterfacial area for IAT
!     ia(i) = ia(i)+dt*(-ia_conv(i)/Vol(i)+IAT_size(i)+iat_coal(i)+iat_break(i)+iat_nucl(i)+iat_bulk(i))
!     iat== 1; iat dominates,
!     iat== 2; constant bubble diameter--> change to iat=1 after 0.1s
!     iat== 3; constant bubble diameter
!     iat== 0; iac dominates, but call iat_calc for comparing iac and iat
!     iat==-1; iac dominates & never use iat_calc
!
      USE VOL_DATA                 
      USE Zparam       , ONLY: pi
      USE Zb_condition , ONLY: alphab_gas
      USE Zbc_index    , ONLY: nvin
      USE Zcoord3      , ONLY: vol
      USE Zconst1      , ONLY: iat,iturb
      USE Zconst2      , ONLY: dt
      USE Ziat         , ONLY: ia,ia_b,ia_conv,ia_old,iat_break,iat_coal,iat_nucl,iat_size, &
                                dbubble_init,dsm_b
      USE Ztimecon     , ONLY: itim_restart,time
      USE Zturb        , ONLY: turb_dp_o
      USE Zvector      , ONLY: vrel_o
      USE Zzone        , ONLY: ncell_fluid
      USE Zio_unit     , ONLY: unit_log
      
!
      IMPLICIT NONE                         
!                                     
      INTEGER i
!      
      LOGICAL,SAVE:: initial
!      
      REAL(8) IAT_rc_co_i,IAT_TI_BK_i
      REAL(8) drhogdt,psii,dsmi
      REAL(8) epsilon
      REAL(8) ia_source
!
      DATA initial/.true./
!
      IF(initial .and. itim_restart.eq.1)THEN
         initial=.false.
         iat_nucl(1:ncell_fluid)=0.d0
         cell%d1(:)=dbubble_init
         DO i=1,nvin
            ia_b(i)=6.d0*alphab_gas(i)/dbubble_init
            dsm_b(i)=dbubble_init
         ENDDO
         DO i=1,ncell_fluid
            ia(i)=6.d0*cell%alphag(i)/dbubble_init
         ENDDO
      ENDIF
!      
!.....Use constant dbubble_init
!
      IF(iat.eq.2 .or. iat.eq.3)THEN
         DO i=1,ncell_fluid
            ia(i)=6.d0*cell%alphag(i)/dbubble_init
            ia(i)=DMAX1(1.0d-8,ia(i))
         ENDDO
         IF(time.gt.0.1d0 .and. iat.eq.2)iat=1
         RETURN
      ENDIF
!
!.....Convection term
!
      CALL IAT_convection
!
!.....Size variation and condensation for gamma<0
!
      CALL IAT_bubsize_change
!
!.....Calculate ia by using ia_old, convection, source=size, coalescence, breakup, nucleation
!
      DO i=1,ncell_fluid
         psii=1.d0/36.d0*pi
         dsmi=cell%d1(i)
         drhogdt=cell%rhog(i)-cell%rhog_o(i)
         epsilon=turb_dp_o(i) 
!
!........Change epsilon according to input options
!
         IF(iat.eq.3)epsilon=0.d0
         IF(iturb.lt.0)epsilon=0.d0
!
!........Calculate iterfacial area
!
         ia_source=0.d0
         iat_coal(i)=1.0d0*IAT_rc_co_i(vrel_o(i),cell%rhol(i),cell%sigma(i),psii,cell%alphag(i),ia_old(i),epsilon,dsmi)
         iat_break(i)=IAT_TI_BK_i(vrel_o(i),cell%rhol(i),cell%sigma(i),psii,cell%alphag(i),ia_old(i),epsilon,dsmi)
!
         ia_source=ia_source+iat_coal(i)+iat_break(i)
         ia_source=ia_source+iat_nucl(i)
         ia_source=ia_source+iat_size(i)
         ia_source=ia_source-ia_conv(i)/Vol(i)
         ia(i)=ia(i)+dt*ia_source
         ia(i)=DMAX1(1.0d-8,ia(i))
         IF(ISNAN(ia_old(i)).eq.1.or.ISNAN(ia(i)).eq.1)THEN
            WRITE(*,*)'iat1:',iat_coal(i),ia(i)
            WRITE(unit_log,*)'iat1:',iat_coal(i),ia(i)
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE IAT_calc

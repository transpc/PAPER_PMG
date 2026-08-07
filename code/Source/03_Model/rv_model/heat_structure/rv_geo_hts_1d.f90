!
      SUBROUTINE rv_geo_hts_1d
!
      USE Zparam,         ONLY:pi,pi2,pi4
      USE Zrv_hts_1d,     ONLY:ng_hts_1d,nr_1d,ri_1d,ht_geo_1d,vl_1d,vr_1d,sl_1d,sr_1d,slw_1d,srw_1d, &
                                ncell_hts_1d,ig_hts_1d,dz_1d,ht_area_left_1d,ht_area_right_1d,num_rod_1d
      USE Zrv_hts_2d,     ONLY:nz0_2d,height_fuel
!                    
      IMPLICIT NONE
!
      INTEGER :: i,k,g
      REAL(8) :: tmp,dz
      REAL(8) :: dxr(nr_1d-1)
!
      dz_1d=height_fuel/DBLE(nz0_2d)
!
      DO k=1,ng_hts_1d
!
         SELECT CASE(ht_geo_1d(k))
!
         CASE(1) ! Rectangular structure
            DO i=1,nr_1d
               IF(i.eq.1)THEN
                  dxr(i)=ri_1d(k,i+1)-ri_1d(k,i)
                  sr_1d(k,i)=1.0d0/dxr(i)
                  vr_1d(k,i)=0.5d0*dxr(i)
               ELSEIF(i.eq.nr_1d)THEN
                  sl_1d(k,i)=1.0d0/dxr(i-1)
                  vl_1d(k,i)=0.5d0*dxr(i-1)
               ELSE
                  dxr(i)=ri_1d(k,i+1)-ri_1d(k,i)
                  sr_1d(k,i)=1.0d0/dxr(i)
                  vr_1d(k,i)=0.5d0*dxr(i)
                  sl_1d(k,i)=1.0d0/dxr(i-1)
                  vl_1d(k,i)=0.5d0*dxr(i-1)
               ENDIF
            ENDDO
            !slw_1d(k)=1.0d0
            !srw_1d(k)=1.0d0
!
            CASE(2) ! Cylindrical structure
               DO i=1,nr_1d
                  IF(i.eq.1)THEN
                     dxr(i)=ri_1d(k,i+1)-ri_1d(k,i)
                     sr_1d(k,i)=pi2*(ri_1d(k,i)+0.5d0*dxr(i))/dxr(i)
                     vr_1d(k,i)=pi*dxr(i)*(ri_1d(k,i)+0.25d0*dxr(i))
                  ELSEIF(i.eq.nr_1d)THEN
                     sl_1d(k,i)=pi2*(ri_1d(k,i)-0.5d0*dxr(i-1))/dxr(i-1)
                     vl_1d(k,i)=pi*dxr(i-1)*(ri_1d(k,i)-0.25d0*dxr(i-1))
                  ELSE
                     dxr(i)=ri_1d(k,i+1)-ri_1d(k,i)
                     sr_1d(k,i)=pi2*(ri_1d(k,i)+0.5d0*dxr(i))/dxr(i)
                     vr_1d(k,i)=pi*dxr(i)*(ri_1d(k,i)+0.25d0*dxr(i))
                     sl_1d(k,i)=pi2*(ri_1d(k,i)-0.5d0*dxr(i-1))/dxr(i-1)
                     vl_1d(k,i)=pi*dxr(i-1)*(ri_1d(k,i)-0.25d0*dxr(i-1))
                  ENDIF
               ENDDO
               slw_1d(k)=pi2*ri_1d(k,1)
               srw_1d(k)=pi2*ri_1d(k,nr_1d)
!
            CASE(3) ! Spherical structure
               DO i=1,nr_1d
                  IF(i.eq.1)THEN
                     dxr(i)=ri_1d(k,i+1)-ri_1d(k,i)
                     tmp=ri_1d(k,i)+0.5d0*dxr(i)
                     sr_1d(k,i)=pi4*tmp*tmp/dxr(i)
                     vr_1d(k,i)=pi4*(tmp*tmp*tmp-ri_1d(k,i)*ri_1d(k,i)*ri_1d(k,i))/3.0d0
                  ELSEIF(i.eq.nr_1d)THEN
                     tmp=ri_1d(k,i)-0.5d0*dxr(i-1)
                     sl_1d(k,i)=pi4*tmp*tmp/dxr(i-1)
                     vl_1d(k,i)=pi4*(ri_1d(k,i)*ri_1d(k,i)*ri_1d(k,i)-tmp*tmp*tmp)/3.0d0
                  ELSE
                     dxr(i)=ri_1d(k,i+1)-ri_1d(k,i)
                     tmp=ri_1d(k,i)+0.5d0*dxr(i)
                     sr_1d(k,i)=pi4*tmp*tmp/dxr(i-1)
                     vr_1d(k,i)=pi4*(tmp*tmp*tmp-ri_1d(k,i)*ri_1d(k,i)*ri_1d(k,i))/3.0d0
                     tmp=ri_1d(k,i)-0.5d0*dxr(i-1)
                     sl_1d(k,i)=pi4*tmp*tmp/dxr(i-1)
                     vl_1d(k,i)=pi4*(ri_1d(k,i)*ri_1d(k,i)*ri_1d(k,i)-tmp*tmp*tmp)/3.0d0
                  ENDIF
               ENDDO
               slw_1d(k)=pi4*ri_1d(k,1)*ri_1d(k,1)
               srw_1d(k)=pi4*ri_1d(k,nr_1d)*ri_1d(k,nr_1d)
!
         END SELECT
!
      ENDDO
!
!.....Heat transfer area for the left and right surfaces
!
      DO k=1,ncell_hts_1d
         g=ig_hts_1d(k)
!         n=num_rod_1d(g)
!         dz=dz_1d*DBLE(n)
         dz=dz_1d*num_rod_1d(k)
         ht_area_left_1d(k)=slw_1d(g)*dz
         ht_area_right_1d(k)=srw_1d(g)*dz
      ENDDO
!
      END SUBROUTINE rv_geo_hts_1d

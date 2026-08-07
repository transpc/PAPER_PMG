!
      SUBROUTINE rv_geo_hts_2d
!
      USE Zparam      , ONLY: pi,pi2
      USE Zconst1     , ONLY: cplmaster
      USE Zrv_hts_2d  , ONLY: nz0_2d,nr_2d,ri_2d,height_fuel,dz_half_fuel, &
                              azl_2f,arl_2f,vl_2f,                         &
                              azr_2f,arr_2f,vr_2f,arw_2f,                  &
                              dxr_2d,ht_geo_2d,z_fuel,v_2f,nrod_2d
      USE Zrv_hts_2d  , ONLY: dz_fuel00,dz_fuel0
      USE Zrv_ncell   , ONLY: nz_fine
!      
      USE Zrv_hts_2d,      ONLY:rod_azr_2f,rod_arr_2f,rod_vr_2f,rod_azl_2f,           &
                                rod_arl_2f,rod_vl_2f,rod_arw_2f,ri_2d_opt,ri_2d_input
!
      IMPLICIT NONE
!
      INTEGER :: i,k
      INTEGER :: j,nz_ff,n_rod
      REAL*8 :: dz_fuelr(nz0_2d*nz_fine),dz_f
      REAL*8 :: dxr_2dr(nr_2d-1)
!
      nz_ff=nz0_2d*nz_fine
      !Original coding
      !dz_fuel=height_fuel/DBLE(nz0_2d*nz_fine)
      !dz_fuelr=DBLE(nz0_2d*nz_fine)/height_fuel
      !dz_half_fuel=0.5d0*dz_fuel
!
      ! Modified dz, various dz
      dz_f=height_fuel/DBLE(nz_ff)
      dz_fuel0(:)=dz_f
!
!.....MASTER coupling example (OPR1000 MSLB)
!
      IF(cplmaster.ne.0)then
         dz_fuel0(:)=dz_fuel00(:)
      ENDIF   
      
      ALLOCATE(dz_half_fuel(nz_ff))
      dz_half_fuel=0.d0
      DO i=1,nz_ff
         dz_fuelr(i)=1.d0/dz_fuel0(i)
         dz_half_fuel(i)=0.5d0*dz_fuel0(i)
      ENDDO
!      
      IF(ri_2d_opt.eq.1)THEN
         n_rod=nrod_2d
      ELSE
         n_rod=1
      ENDIF  
      ALLOCATE(rod_azr_2f(nz_ff,n_rod,nr_2d))
      ALLOCATE(rod_arr_2f(nz_ff,n_rod,nr_2d))
      ALLOCATE(rod_vr_2f (nz_ff,n_rod,nr_2d))
      ALLOCATE(rod_azl_2f(nz_ff,n_rod,nr_2d))
      ALLOCATE(rod_arl_2f(nz_ff,n_rod,nr_2d))
      ALLOCATE(rod_vl_2f (nz_ff,n_rod,nr_2d))
      ALLOCATE(rod_arw_2f(nz_ff,n_rod))
!         
      DO j=1,n_rod
!     
         ri_2d(:)=ri_2d_input(j,:)
         DO i=1,nr_2d-1
            dxr_2d(i)=ri_2d(i+1)-ri_2d(i)
            dxr_2dr(i)=1.d0/dxr_2d(i)
!           dxl_2d(i+1)=dxr_2d(i)
         ENDDO !i
!
           !Origin
           !z_fuel(k)=dz_half_fuel(k)+dz_fuel0(k)*dble(k-1)
           !Modified
            k=1
            z_fuel(k)=dz_half_fuel(k)
         DO k=2,nz_ff
            z_fuel(k)=z_fuel(k-1)+dz_half_fuel(k-1)+dz_half_fuel(k)
         ENDDO
!
         SELECT CASE(ht_geo_2d)
!     
            CASE(1) ! Rectangular structure
               DO i=1,nr_2d-1
                  DO k=1,nz_ff
                     azr_2f(k,i)=0.5d0*dxr_2d(i) *dz_fuelr(k)
                     arr_2f(k,i)=      dxr_2dr(i)*dz_fuel0(k)
                     vr_2f(k,i) =0.5d0*dxr_2d(i) *dz_fuel0(k)
!     
                     azl_2f(k,i+1)=0.5d0*dxr_2d(i) *dz_fuelr(k)
                     arl_2f(k,i+1)=      dxr_2dr(i)*dz_fuel0(k)
                     vl_2f(k,i+1) =0.5d0*dxr_2d(i) *dz_fuel0(k)
!
                     v_2f(k,i)=vr_2f(k,i)+vl_2f(k,i+1)
                  ENDDO
               ENDDO
!!    
                  i=nr_2d
                  DO k=1,nz_ff
                     arw_2f(k)=dz_fuel0(k)
                     v_2f(k,i)=1.d0            
                  ENDDO
!     
            CASE(2) ! Cylindrical structure
!!    
               DO i=1,nr_2d-1
                  DO k=1,nz_ff
                     azr_2f(k,i)=pi*dxr_2d(i)*(ri_2d(i)+0.25*dxr_2d(i))*dz_fuelr(k)
                     arr_2f(k,i)=pi*(ri_2d(i)+0.5*dxr_2d(i))*dz_fuel0(k)*2.d0*dxr_2dr(i)
                     vr_2f(k,i) =pi*dxr_2d(i)*(ri_2d(i)+0.25*dxr_2d(i))*dz_fuel0(k)
!     
                     azl_2f(k,i+1)=pi*dxr_2d(i)*(ri_2d(i+1)-0.25*dxr_2d(i))*dz_fuelr(k)
                     arl_2f(k,i+1)=pi*(ri_2d(i+1)-0.5*dxr_2d(i))*dz_fuel0(k)*2.d0*dxr_2dr(i)
                     vl_2f(k,i+1) =pi*dxr_2d(i)*(ri_2d(i+1)-0.25*dxr_2d(i))*dz_fuel0(k)
!
                     v_2f(k,i)=vr_2f(k,i)+vl_2f(k,i+1)
                  ENDDO
               ENDDO
!!    
                  i=nr_2d
                  DO k=1,nz_ff
                     arw_2f(k)=pi2*ri_2d(i)*dz_fuel0(k)
                     v_2f(k,i)=1.d0
                  ENDDO
!     
         END SELECT
!
         IF(ri_2d_opt.eq.1)THEN
            DO i=1,nr_2d
               DO k=1,nz_ff
                  rod_azr_2f(k,j,i)=azr_2f(k,i)
                  rod_arr_2f(k,j,i)=arr_2f(k,i)
                  rod_vr_2f (k,j,i)=vr_2f (k,i) 
                  rod_azl_2f(k,j,i)=azl_2f(k,i)
                  rod_arl_2f(k,j,i)=arl_2f(k,i)
                  rod_vl_2f (k,j,i)=vl_2f (k,i)
               ENDDO
            ENDDO
            DO k=1,nz_ff
               rod_arw_2f(k,j)  =arw_2f(k)
            ENDDO
         ENDIF
!      
      ENDDO !j     
!      
      END SUBROUTINE rv_geo_hts_2d
!

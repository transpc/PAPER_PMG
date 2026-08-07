!
      SUBROUTINE udfn_wall_drag(vfwl,vfwg)
!
!     This routine includes user defined wall drag model
!     This code assumes that for each i cell there is only one neighbor
!     with nbcon=-5    
!
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Znum_cell    , ONLY: istart_nf,istart_nb1,istart_nbcon_nf, &
                               ia_nb,icell_nb
      USE Zvec_index   , ONLY: nbcon_nf
      USE Zvector      , ONLY: vl_o
      USE Zb_condition , ONLY: v_wall        
      USE Zvec_geo     , ONLY: saa_nf,djia_nf
!
      IMPLICIT NONE
! 
!.....Output
      REAL(8) :: vfwg(ncell_fp),vfwl(ncell_fp)
!.....Local variables
      INTEGER :: i,ix,nb
      INTEGER :: ii
      INTEGER :: istart3,i3
      INTEGER :: nf_number,istart,len,istart1,i0,i1
      REAL(8) :: vfwl_i
      REAL(8) :: tau_w_film 
      REAL(8) :: vl1
!
      DO i=1,ncell_fluid
         vfwl(i)=0.d0
         vfwg(i)=0.d0
      ENDDO
      nf_number=7
      istart =istart_nf(1,nf_number)
      istart3=istart_nbcon_nf(nf_number)
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len  
         i1=istart1+nb
         ii=icell_nb(i1)
         IF(cell%alphag(ii).gt.0.999d0) cycle
         vfwl_i=0.d0
!DIR$ NOVECTOR
         DO i=ia_nb(i1),ia_nb(i1+1)-1
            i3=istart3+i
            IF(nbcon_nf(i3).ne.-5) cycle
            i0=istart+i
            DO ix=1,1                                ! ndim for general problems
               vl1=vl_o(ii,ix)
               cell%film_thickness(ii)=cell%alphal(ii)*djia_nf(i0)*2.d0
               !tau_w_film=3.d0*cell(i)%lviscosl*(v_wall(ix)-vl1)/cell(i)%film_thickness
               tau_w_film= 3.d0*cell%lviscosl(ii)                                                               &
                          *(v_wall(ix)-vl1*djia_nf(i0)*2.d0/cell%film_thickness(ii))/cell%film_thickness(ii)**2  
               tau_w_film=ABS(tau_w_film)
               IF(tau_w_film.eq.0.d0)THEN
                  vfwl_i=0.d0
               ELSE
                  IF(vl1.ne.0.d0) THEN                        
                     vfwl_i=tau_w_film*saa_nf(i0)                
                     vfwl_i=MAX(100.d0,MIN(vfwl_i,50000.d0))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         vfwl(ii)=vfwl_i
         vfwg(ii)=vfwl_i   ! assume that the wall friction for gas phase is equal to that for liquid phase - need to be modified 
      ENDDO
!
      END SUBROUTINE udfn_wall_drag 

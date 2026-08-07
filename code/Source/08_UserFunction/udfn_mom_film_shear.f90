!
      SUBROUTINE udfn_mom_film_shear
!
!.....This routine calculates wall shear force for laminar liquid film (currently adopted in DIVA problem)
!     This code assumes that for each i cell there is only one neighbor
!     with nbcon=-3
!

      USE VOL_DATA     , ONLY: cell
      USE Zvec_index   , ONLY: nbcon_nf
      USE Znum_cell    , ONLY: istart_nf,istart_nb1,istart_nbcon_nf, &
                               ia_nb,icell_nb
      USE Zb_condition , ONLY: v_wall      
!!!      USE Zsrc         , ONLY: wall_shear
      USE Zvector      , ONLY: vl_o
      USE Zbc_index    , ONLY: npb
      USE Zvec_geo     , ONLY: sap_nf,djia_nf
      USE Zm_src       , ONLY: src_gas,src_liq,src_drp
      USE Zparam       , ONLY: ndim      
      USE Zzone        , ONLY: ncell_fluid
      USE Zmpi         , ONLY: ncell_fp
!
      IMPLICIT NONE
!
!     local variables
      INTEGER :: i,ix,nb
      INTEGER :: ii
      INTEGER :: istart3,i3
      INTEGER :: nf_number,istart,len,istart1,i0,i1
      REAL(8) tau_w_film   
      REAL(8) wall_shear(ncell_fp,ndim)  
!
      nf_number=6
      istart=istart_nf(1,nf_number)
      istart1=istart_nb1(1,nf_number)
      istart3=istart_nbcon_nf(nf_number)
      len   =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
         ii=icell_nb(i1)
         IF(cell%alphag(ii).gt.0.9999d0 .or. npb(ii).gt.0) cycle
!DIR$ NOVECTOR
         DO i=ia_nb(i1),ia_nb(i1+1)-1
            i3=istart3+i
            IF(nbcon_nf(i3).ne.-3) cycle                 ! specify nbcon=3 to define the cells which contain the liquid film 
            i0=istart+i
            ix=1                                    ! specify the directions where the wall shear is applied
            cell%film_thickness(ii)=cell%alphal(ii)*djia_nf(i0)*2.0d0
!           tau_w_film=3.0d0*cell(i)%lviscosl*(v_wall(ix)-vl_o(i,ix))/cell(ii)%film_thickness
            tau_w_film=3.0d0*cell%lviscosl(ii)*(v_wall(ix)-vl_o(i,ix)*djia_nf(i0)*2.0d0/  &
                       cell%film_thickness(ii))/cell%film_thickness(ii)                                
            IF(vl_o(ii,ix).ne.0.0d0) THEN
               wall_shear(ii,ix)=tau_w_film*sap_nf(i0) !*abs(vl_o(i,ix))/vl_o(i,ix)
               wall_shear(ix,i)=MAX(1000.0d0,MIN(wall_shear(ix,i),50000.d0)) ! Set the minimum and maximum values
            ENDIF
!
            ix=3                                    ! specify the directions where the wall shear is applied
            cell%film_thickness(ii)=cell%alphal(ii)*djia_nf(i0)*2.0d0
!           tau_w_film=3.0d0*cell(i)%lviscosl*(v_wall(ix)-vl_o(i,ix))/cell(ii)%film_thickness
            tau_w_film=3.0d0*cell%lviscosl(ii)*(v_wall(ix)-vl_o(i,ix)*djia_nf(i0)*2.0d0/  &
                       cell%film_thickness(ii))/cell%film_thickness(ii)                                
            IF(vl_o(ii,ix).ne.0.0d0) THEN
               wall_shear(ii,ix)=tau_w_film*sap_nf(i0) !*abs(vl_o(i,ix))/vl_o(i,ix)
            ENDIF
         ENDDO
      ENDDO
! 
!.....Add wall shear to momentum source
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            src_gas(i,1)=src_gas(i,1)+wall_shear(i,1)
            src_gas(i,2)=src_gas(i,2)+wall_shear(i,2)
            src_liq(i,1)=src_liq(i,1)+wall_shear(i,1)
            src_liq(i,2)=src_liq(i,2)+wall_shear(i,2)
            src_drp(i,1)=src_drp(i,1)+wall_shear(i,1)
            src_drp(i,2)=src_drp(i,2)+wall_shear(i,2)
         ENDDO
      ELSEIF(ndim.eq.3)THEN
         DO i=1,ncell_fluid
            src_gas(i,1)=src_gas(i,1)+wall_shear(i,1)
            src_gas(i,2)=src_gas(i,2)+wall_shear(i,2)
            src_gas(i,3)=src_gas(i,3)+wall_shear(i,3)
            src_liq(i,1)=src_liq(i,1)+wall_shear(i,1)
            src_liq(i,2)=src_liq(i,2)+wall_shear(i,2)
            src_liq(i,3)=src_liq(i,3)+wall_shear(i,3)
            src_drp(i,1)=src_drp(i,1)+wall_shear(i,1)
            src_drp(i,2)=src_drp(i,2)+wall_shear(i,2)
            src_drp(i,3)=src_drp(i,3)+wall_shear(i,3)
         ENDDO
      ENDIF
!         
      END SUBROUTINE

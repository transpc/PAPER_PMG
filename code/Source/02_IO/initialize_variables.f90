!     
      SUBROUTINE initialize_variables(v0,p0,tl0,tg0,a0_g,quala0,tsol0,cboron0)
!
!     This routine initializes cell variables
!
      USE VOL_DATA     , ONLY: cell            
      USE WALL_DATA    , ONLY: face
      USE SOLID_DATA   , ONLY: solid
      USE Zparam       , ONLY: ndim,nb_max
      USE Znum_cell    , ONLY: i_neigh
      USE Zb_condition , ONLY: vb_lold,vb_gold,vb_gas,vb_liq,v_wall
      USE Zbc_index    , ONLY: nbcon,icell_type,iface_wall
      USE Zpress       , ONLY: p
      USE Zuserdefined , ONLY: udfl_init_variables
      USE Zvector      , ONLY: vg_n,vl_n,vd_n
      USE Zzone        , ONLY: nzone,ncell_fluid,ncell_cond,num_max_zone,nmaterial_c,nzone_c
!
      IMPLICIT NONE
!      
!.....Input
      REAL(8) v0(num_max_zone,ndim),p0(num_max_zone),tl0(num_max_zone)
      REAL(8) tg0(num_max_zone),quala0(num_max_zone),a0_g(num_max_zone)
      REAL(8) tsol0(num_max_zone),cboron0(num_max_zone)
!
!....Local variables
      INTEGER i,ii,j,j0,ix
      INTEGER iokr,iokk
      REAL(8) CpVol,Condu      
!
!.....Initialize cell variables
!
      DO ix=1,ndim
         DO i=1,ncell_fluid
            ii=nzone(i)
            vl_n(i,ix)=v0(ii,ix)
            vd_n(i,ix)=v0(ii,ix)
            vg_n(i,ix)=v0(ii,ix)
         ENDDO
      ENDDO
      DO i=1,ncell_fluid
         ii=nzone(i)
         p(i)=p0(ii)
         cell%tl(i)=tl0(ii)
         cell%tg(i)=tg0(ii)
         cell%alphag(i)=a0_g(ii)
         cell%alphal(i)=1.d0-a0_g(ii)
         cell%quala(i)=quala0(ii)
         cell%cboron(i)=cboron0(ii)
         cell%cboron_o(i)=cboron0(ii)
      ENDDO
!
!!!      DO i=1,ncell_fluid
!!!         IF(npb(i).gt.0) p(i)=pbnd(npb(i))
!!!      ENDDO
!
      DO i=1,ncell_fluid
         cell%alphad(i)=0.d0
         cell%dentr(i)=0.d0
         cell%entr(i)=0.d0
         cell%vfgl(i)=0.d0
         cell%vfgd(i)=0.d0
         cell%eviscosl(i)=0.d0
         cell%eviscosg(i)=0.d0
         cell%vfwg(i)=0.d0
         cell%vfwl(i)=0.d0
         cell%mdiff(i)=0.d0
         cell%film_thickness(i)=0.d0
         cell%film_shear(i)=0.d0     
         cell%wf_vst(i)=0.d0    
!
!........Initialize to solve somaplot.viw problem
!
         cell%regime(i)=0
         cell%ced33(i)=0.d0
      ENDDO
!
!.....Initialize wall variables
!      
      face%wall_fluxl_diff(:)=0.d0
      face%wall_fluxg_diff(:)=0.d0
      face%wall_fluxd_diff(:)=0.d0
      face%ddepartw(:)=0.d0
      face%ratio_evap(:)=0.d0
      face%twall_partition(:)=0.d0
      IF(ncell_cond.gt.0)solid%tsol(:)=0.d0
!
!.....Call & save material properties
!    
      DO i=1,ncell_cond
         solid%tsol(i)=tsol0(nzone_c(i))
         IF(ABS(nmaterial_c(i)).lt.50) THEN
            CALL mat_prop(ABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iokr,iokk)
            solid%rhocps(i)=CpVol
            solid%conds(i)=Condu
         ENDIF
         solid%matnum(i)=ABS(nmaterial_c(i))         
      ENDDO
!
!.....Save boundary velocity into old values
!      
      DO ix=1,ndim
!DIR$ NOVECTOR
         DO i=1,nb_max
            vb_lold(i,ix)=vb_liq(i,ix)
            vb_gold(i,ix)=vb_gas(i,ix)
         ENDDO
      ENDDO
!
!.....Initialize heat partitioning related variables
!       
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1)THEN
            j0=i_neigh(i)-1
            j=iface_wall(i)
            IF(nbcon(j+j0).le.-1)THEN
               face%ratio_evap(i)=0.d0
               face%twall_partition(i)=cell%tg(i)-0.5d0
            ENDIF   
         ENDIF
      ENDDO
!
!.....Initialize wall variables
! 
      cell%twall(:)=0.d0
!DIR$ NOVECTOR
      DO ix=1,ndim
         v_wall(ix)=0.d0
      ENDDO
!
!.....Call user functions
! 
      CALL initialize_topology_criteria
      IF(udfl_init_variables) CALL initialize_specific_variables
!      
      END SUBROUTINE initialize_variables

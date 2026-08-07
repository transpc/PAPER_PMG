!
      SUBROUTINE initialize_topology(itp)
!
!     Define topology map
! 
      USE VOL_DATA    , ONLY: cell
      USE Zzone       , ONLY: ncell_fluid
      USE Zflowregime , ONLY: alphag_bc,gamma_1,gamma_2,alphag_cm
      USE Zvoid       , ONLY: gamma_void 
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,itp
!	  
!.....Calculate gamma (void fraction gradient) and return 
!
      IF(itp.eq.0) THEN
         CALL topology_gamma
         RETURN 
      ENDIF    
!
!.....Determine Topology according to gamma and alphag
!
      DO i = 1, ncell_fluid
         IF(cell%alphag(i).le.alphag_bc)THEN
            IF(gamma_void(i).le.gamma_1)THEN
               cell%regime(i)=11
            ELSEIF(gamma_void(i).lt.gamma_2)THEN
               cell%regime(i)=21
            ELSE
               cell%regime(i)=3
            ENDIF
         ELSEIF(cell%alphag(i).lt.alphag_cm)THEN
            IF(gamma_void(i).le.gamma_1)THEN
               cell%regime(i)=12
            ELSEIF(gamma_void(i).lt.gamma_2)THEN
               cell%regime(i)=22
            ELSE
               cell%regime(i)=3
            ENDIF
         ELSE
            IF(gamma_void(i).le.gamma_1)THEN
               cell%regime(i)=13
            ELSEIF(gamma_void(i).lt.gamma_2)THEN
               cell%regime(i)=23
            ELSE
               cell%regime(i)=3
            ENDIF
         ENDIF
      ENDDO
!
      END SUBROUTINE initialize_topology
!
!----------------------------------------------------------------------
!
      SUBROUTINE topology_gamma
!	  
!     Calculate gamma
!
      USE Zinterface
      USE VOL_DATA    , ONLY: cell
      USE Zzone       , ONLY: ncell_fluid
      USE Zparam      , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_tot
      USE Znum_cell   , ONLY: i_neigh, &
                              nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zcoord3     , ONLY: vol
      USE Zvoid       , ONLY: dagdx,gamma_void
      USE Zconst1     , ONLY: vv_prob
      USE Zvec_geo    , ONLY: sv_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i
      REAL(8) :: char_l
      REAL(8) :: delta1,delta2,delta3
      REAL(8) :: numneigh
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid,ndim) :: delta
!.....Local vector arrays
      REAL(8),DIMENSION(nf_tot,ndim) :: sv_nf0
!        
      CALL grad_scalar(cell%alphag,dagdx,ncell_fluid)
!
      DO i=1,ncell_fluid
         cell%aint3(i)=sqrt(dot_product(dagdx(i,:),dagdx(i,:)))
      ENDDO
!
!.....'check_Hik' defines its own gamma_void to check continuity in HTCs.  
!	  
       IF(vv_prob.eq.'check_Hik') RETURN  
!	  
      IF(ndim.eq.2) THEN
         DO i=1,nf_tot
            sv_nf0(i,1)=abs(sv_nf(i,1))
            sv_nf0(i,2)=abs(sv_nf(i,2))
         ENDDO
      ELSE
         DO i=1,nf_tot
            sv_nf0(i,1)=abs(sv_nf(i,1))
            sv_nf0(i,2)=abs(sv_nf(i,2))
            sv_nf0(i,3)=abs(sv_nf(i,3))
         ENDDO
      ENDIF
!
!.....Build summation info for non
!
      nf_number_nb=8
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      nf_number_id(4)=4
      nf_number_id(5)=5
      nf_number_id(6)=6
      nf_number_id(7)=7
      nf_number_id(8)=8
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      istart_nfs(4)=istart_nfs(3)+nf_out
      istart_nfs(5)=istart_nfs(4)+nf_adw
      istart_nfs(6)=istart_nfs(5)+nf_fsw
      istart_nfs(7)=istart_nfs(6)+nf_ctw
      istart_nfs(8)=istart_nfs(7)+nf_chw
      lens         =istart_nfs(8)+nf_sym
!
      CALL sum_nf_ndim(0,1,ncell_fluid, &
                       sv_nf0,delta)
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            numneigh=real(i_neigh(i+1)-i_neigh(i))
            delta1=numneigh*vol(i)/delta(i,1)
            delta2=numneigh*vol(i)/delta(i,2)
            char_l=delta1*abs(dagdx(i,1))+delta2*abs(dagdx(i,2))
            gamma_void(i)=0.5d0*char_l
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            numneigh=dble(i_neigh(i+1)-i_neigh(i))
            delta1=numneigh*vol(i)/delta(i,1)
            delta2=numneigh*vol(i)/delta(i,2)
            delta3=numneigh*vol(i)/delta(i,3)
            char_l=delta1*abs(dagdx(i,1))+delta2*abs(dagdx(i,2))+delta3*abs(dagdx(i,3))
            gamma_void(i)=0.5d0*char_l
         ENDDO
      ENDIF
!
      END SUBROUTINE topology_gamma

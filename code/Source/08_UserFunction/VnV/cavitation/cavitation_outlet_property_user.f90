!
      SUBROUTINE cavitation_outlet_property_user
!
!     This routine declares flow properties at pressure outlet boundaries
!
      USE VOL_DATA              
      USE Zzone     , ONLY: ncell_fluid
      USE Zbc_index , ONLY: npb
      USE Znum_cell , ONLY: i_neigh,neigh,indexr_sort
      USE Zvector   , ONLY: vl_o,vg_o,vd_o
!
      IMPLICIT NONE
!
      INTEGER i,j,k,j0
!      
!.....Declare flow properties at outlet boundary as follows. 
!
      DO i=1,ncell_fluid
         IF(npb(i).eq.1)THEN
            vl_o(i,2)=0.0d0
            vg_o(i,2)=0.0d0
            vd_o(i,2)=0.0d0
         ENDIF
         IF(npb(i).eq.2)THEN
!..............neighbor have been sorted use indexr_sort to retrieve j=3
!           k=neigh(3,i)
            j0=i_neigh(i)-1
            j=indexr_sort(3+j0)
            k=neigh(j+j0)
            vl_o(i,1)=max(vl_o(k,1),0.0d0)
            vg_o(i,1)=max(vg_o(k,1),0.0d0)
            vd_o(i,1)=max(vd_o(k,1),0.0d0)
            vl_o(i,2)=0.0d0
            vg_o(i,2)=0.0d0
            vd_o(i,2)=0.0d0
            cell%el(i)=cell%el(k)
            cell%eg(i)=cell%eg(k)
            cell%alphal(i)=cell%alphal(k)
            cell%alphag(i)=cell%alphag(k)
            cell%quala(i)=cell%quala(k)
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE cavitation_outlet_property_user

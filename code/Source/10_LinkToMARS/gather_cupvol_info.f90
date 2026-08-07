!
      SUBROUTINE gather_cupvol_info
!      
      USE CMP_DAT !cmp_hd,cmp_da
      USE VOL_DAT !v_hd,v_da
      USE Zcore   !myrank
!      
      IMPLICIT NONE
!      
!DEC$IF defined (MCC_DLL)
      !DEC$ ATTRIBUTES DLLIMPORT:: /c3com/,/c3com_dll/
!DEC$ENDIF      
      INCLUDE 'c3com.h'
      INCLUDE 'contrl.h' !timehy  
!      
      INTEGER i,j 
!
!.....Get MARS time
!     
      IF(nrank.gt.1)CALL broadcast_r(timehy,1)
      c3time_sys=timehy
!
!.....Count the number of CUPVOLs
!      
      IF(myrank.ne.0)GOTO 98701 !mcc-mpi      
      i3nic(2)=0 
      DO i=1,cmp_hd%nComp(2)
         IF(cmp_da(i)%type(2).eq.21)THEN !MCC-148
            i3nic(2)=i3nic(2)+1
         ENDIF 
      ENDDO
      i3nic(3)=i3nic(1)+i3nic(2)      
98701 CONTINUE !mcc-mpi
      IF(nrank.gt.1)CALL broadcast_i(i3nic(3),1) !mcc-mpi
      IF(nrank.gt.1)CALL broadcast_i(i3nic(2),1) !mcc-mpi
!
!.....Search CUPVOL and save its component number 
!
      IF(myrank.ne.0)GOTO 98702 !mcc-mpi      
      j=0
      DO i=1,cmp_hd%nComp(2)
         IF(cmp_da(i)%type(2).eq.21)THEN !MCC-148
		    j=j+1  
            i1Cvoln(1,j)=cmp_da(i)%number(2)*1000000+10000
         ENDIF 
      ENDDO
!
!......Save the CUPID cell index of CUPVOLs
      DO j=1,i3nic(2)
         DO i=1,v_hd%nVols(2) 
            IF(i1Cvoln(1,j).eq.v_da(i)%VolNo(2))THEN
               i3cell(1,j)=-v_da(i)%VolNo(1)-100 !cupid global index
               EXIT
            ENDIF
         ENDDO
      ENDDO
98702 CONTINUE !mcc-mpi
      IF(nrank.gt.1)CALL broadcast_i(i3cell(1,:),72) !mcc-mpi
!
      RETURN 
      ENDSUBROUTINE gather_cupvol_info

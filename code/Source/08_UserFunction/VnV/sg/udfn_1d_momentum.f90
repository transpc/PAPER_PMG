!
      SUBROUTINE udfn_1d_momentum
!
      USE Zsg,     ONLY:n_group,n_1d,vn_1d,rho_1d,pr_flow,ar_tube
!
      IMPLICIT NONE
!
      INTEGER i,j
      REAL(8) fw
! 
      DO i=1,n_group
         fw=pr_flow(i)/ar_tube(i)
         DO j=1,n_1d(i)
            vn_1d(i,j)=fw/rho_1d(i,j)
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE udfn_1d_momentum

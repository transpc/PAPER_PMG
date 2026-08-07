!
      SUBROUTINE udfn_1d_energy
!
      USE Zzone,   ONLY: ncell_fluid,ncell_fluid_all
      USE Zconst2, ONLY: dt
      USE Zsg,     ONLY: n_group,n_1d,en_1d,eo_1d,ein_1d,rho_1d,q_pri,vn_1d,vin_1d,h_tube,  &
                         drde_1d,pr_flow,ar_tube,pin_1d,p_1d,sd_cell
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,m
      REAL(8) :: a,b,fw,dtlj
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: q_pri_tmp
!
!.....Solve 1D Energy equation
!
      ALLOCATE(q_pri_tmp(ncell_fluid_all))
      CALL allgatherv_r(q_pri,q_pri_tmp,ncell_fluid,ncell_fluid_all,0)
!
      DO i=1,n_group
         m=sd_cell(i,1)
         dtlj=2.0d0*dt/h_tube(i,1)
         fw=pr_flow(i)/ar_tube(i)
         a=eo_1d(i,1)*drde_1d(i,1)+rho_1d(i,1)
         b=a+dtlj*fw
         en_1d(i,1)=a*eo_1d(i,1)+dtlj*fw*ein_1d+dt*q_pri_tmp(m)+dtlj*(pin_1d*vin_1d(i)-p_1d(i,1)*vn_1d(i,1))
!         en_1d(i,1)=a*eo_1d(i,1)+dtlj*fw*ein_1d+dt*q_pri_tmp(m)
         en_1d(i,1)=en_1d(i,1)/b
         DO j=2,n_1d(i)
            m=sd_cell(i,j)
            dtlj=2.0d0*dt/(h_tube(i,j)+h_tube(i,j-1))
            a=eo_1d(i,j)*drde_1d(i,j)+rho_1d(i,j)
            b=a+dtlj*fw
            en_1d(i,j)=a*eo_1d(i,j)+dtlj*fw*en_1d(i,j-1)+dt*q_pri_tmp(m)+dtlj*(p_1d(i,j-1)*vn_1d(i,j-1)-p_1d(i,j)*vn_1d(i,j))
!            en_1d(i,j)=a*eo_1d(i,j)+dtlj*fw*en_1d(i,j-1)+dt*q_pri_tmp(m)
            en_1d(i,j)=en_1d(i,j)/b
         ENDDO
         eo_1d(i,:)=en_1d(i,:)
!
      ENDDO
!
      DEALLOCATE(q_pri_tmp)
!
      END SUBROUTINE udfn_1d_energy

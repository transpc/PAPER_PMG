!
      SUBROUTINE udfn_sg_htc
!
      USE Zzone      ,ONLY: ncell_fluid
      USE Zsg        ,ONLY: htc_pr,rho_1d,di_tube,vn_1d,vis_1d,cp_1d,cond_1d,igr,j1d,t_tube,t_1d
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER g,m,i
!
!.....Calculate heat transfer coefficient inside U-tube (primary coolant)
!
      DO i=1,ncell_fluid
         g=igr(i)
         IF(g.eq.0) CYCLE
         m=j1d(i)
         htc_pr(i)=0.23d0*(rho_1d(g,m)*di_tube*vn_1d(g,m)/vis_1d(g,m))**0.8
         htc_pr(i)=htc_pr(i)*(cp_1d(g,m)*vis_1d(g,m)/cond_1d(g,m))**0.333
         htc_pr(i)=htc_pr(i)*cond_1d(g,m)/di_tube
         IF(t_tube(1,g,m).gt.t_1d(g,m)) htc_pr(i)=0.0d0
      ENDDO
!
!.....Calculate heat transfer coefficient inside U-tube (secondary coolant)
!
      CALL udfn_sg_heat_partition
!
      END SUBROUTINE udfn_sg_htc

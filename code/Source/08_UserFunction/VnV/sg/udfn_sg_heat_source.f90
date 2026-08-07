!
      SUBROUTINE udfn_sg_heat_source
!
      USE Zzone        ,ONLY: ncell_fluid
      USE Zsg          ,ONLY: t_1d,mult_cell,mult_3d_cell2, &
                              ih,iavb,izp,ihp,igr,j1d,q_sd,ht_area
      USE Zturb        ,ONLY: yplus,yplusg,utau,tauw,turb_ke,turb_keg
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER g,m,i,k
      REAL(8) qflux
!
      CALL udfn_sg_primary
!
!.....Assign post-processing ouput to unused variables
!
      yplus(:)=0.0d0
      yplusg(:)=0.0d0
      utau(:)=0.0d0
      tauw(:)=0.0d0
      turb_ke(:)=0.0d0
      turb_keg(:)=0.0d0
      DO i=1,ncell_fluid
         g=igr(i)
         IF(g.eq.0) CYCLE
         m=j1d(i)
         qflux=q_sd(g,m)/ht_area(g,m)
         turb_ke(i)=t_1d(g,m)
         turb_keg(i)=qflux
         yplusg(i)=iavb(g,m)
         utau(i)=ih(g,m)
         yplus(i)=izp(g,m)
         tauw(i)=ihp(g,m)
         IF(mult_cell(g,m).gt.0)THEN
            k=mult_3d_cell2(mult_cell(g,m))
            turb_ke(k)=t_1d(g,m)
            turb_keg(k)=qflux
            yplusg(k)=iavb(g,m)
            utau(k)=ih(g,m)
            yplus(k)=izp(g,m)
            tauw(k)=ihp(g,m)
         ENDIF
      ENDDO
!
      END SUBROUTINE udfn_sg_heat_source

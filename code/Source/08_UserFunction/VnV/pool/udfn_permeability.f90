      SUBROUTINE udfn_permeability
!
!.....User defined Permeability Update for PAFS-POOL
!
      USE Znum_cell    , ONLY: ncell,i_neigh_tmp,j_neigh_tmp, &
                               perm_tmp1
      USE Zconst1      , ONLY: vv_prob      
!      
      IMPLICIT NONE
!
!     local variables
      INTEGER i,j,k
!      
      IF(vv_prob.eq.'PAFS-POOL') THEN
          DO i=1,ncell
             DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                k=j_neigh_tmp(j)
                IF(i.eq.1700.and.(k.eq.1811.or.k.eq.1806)) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1702.and.(k.eq.1812.or.k.eq.1807)) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1705.and.(k.eq.1813.or.k.eq.1808)) THEN
                   perm_tmp1(j)=1.d0            
                ELSEIF(i.eq.1709.and.(k.eq.1814.or.k.eq.1809)) THEN
                   perm_tmp1(j)=1.d0            
                ELSEIF(i.eq.1713.and.(k.eq.1815.or.k.eq.1810)) THEN
                   perm_tmp1(j)=1.d0                           
                ELSEIF(i.eq.1683.and.(k.eq.1801.or.k.eq.1796)) THEN
                   perm_tmp1(j)=1.d0                                      
                ELSEIF(i.eq.1685.and.(k.eq.1802.or.k.eq.1797)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1688.and.(k.eq.1803.or.k.eq.1798)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1692.and.(k.eq.1804.or.k.eq.1799)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1696.and.(k.eq.1805.or.k.eq.1800)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1811.and.k.eq.1700) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1812.and.k.eq.1702) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1813.and.k.eq.1705) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1814.and.k.eq.1709) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1815.and.k.eq.1713) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1806.and.k.eq.1700) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1807.and.k.eq.1702) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1808.and.k.eq.1705) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1809.and.k.eq.1709) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1810.and.k.eq.1713) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1801.and.k.eq.1683) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1802.and.k.eq.1685) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1803.and.k.eq.1688) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1804.and.k.eq.1692) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1805.and.k.eq.1696) THEN
                   perm_tmp1(j)=1.d0                        
                ELSEIF(i.eq.1796.and.k.eq.1683) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1797.and.k.eq.1685) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1798.and.k.eq.1688) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1799.and.k.eq.1692) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1800.and.k.eq.1696) THEN
                   perm_tmp1(j)=1.d0                        
                ENDIF
             ENDDO
          ENDDO 
      ENDIF           
!
      RETURN
      END SUBROUTINE udfn_permeability

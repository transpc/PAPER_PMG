!
      SUBROUTINE property_calc(imode)
!
!     This routine calls steamtable 
!
      USE VOL_DATA       , ONLY: cell
      USE STM_TBL_cupid  , ONLY: np
      USE Zncg           , ONLY: tao,cvao_cell,uao_cell,dcva_cell,ra_cell,qn_cell,n_ncg_sp
      USE Zpress         , ONLY: p  
      USE Zzone          , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i,imode
      REAL(8) qn_cell0(n_ncg_sp)
!
      DO i=1,ncell_fluid
         cell%p(i)     =p(i)                                 
         cell%rhog_o(i)=cell%rhog(i)   !pik-2010-03-29-ins                           
         cell%rhol_o(i)=cell%rhol(i)   !pik-boron-2010-04-20-ins
      ENDDO
!
      IF(imode.eq.0)THEN
         DO i=1,ncell_fluid
            qn_cell0(:)=qn_cell(i,:)
            CALL convert_temp2erg(cell%p(i),cell%tl(i),cell%tg(i),cell%quala(i),cell%el(i),cell%eg(i), &
                                  cell%rhol(i),cell%rhog(i),cell%pps_o(i),cell%estm_o(i),              &
                                  tao,cvao_cell(i),uao_cell(i),dcva_cell(i),ra_cell(i),qn_cell0)
         ENDDO   
      ENDIF
      CALL steamtable(imode,np)
!
      DO i=1,ncell_fluid
         IF(cell%alphal(i).le.1.d-8) cell%el(i)=cell%elsat(i)
         IF(cell%alphal(i).le.1.d-8) cell%tl(i)=cell%ts(i)
         cell%ed(i)=cell%el(i)        
      ENDDO  
!
      IF(imode.eq.0)THEN
         DO i=1,ncell_fluid
             cell%eviscosl(i)=cell%lviscosl(i)
             cell%tviscosl(i)=cell%lviscosl(i)
             cell%eviscosg(i)=cell%lviscosg(i)
             cell%tviscosg(i)=cell%lviscosg(i)
         ENDDO
      ENDIF         
!
      END SUBROUTINE property_calc
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
      FUNCTION sat_temp(pi)
!
      USE STM_TBL_cupid  , ONLY: st_tbl,          &
                                 nt,ndxstd,nfluid
!      
      IMPLICIT NONE 
!
      LOGICAL erx
      REAL(8) s(36)
      REAL(8) pi,sat_temp
!      
      s(2)=pi         
      s(9)=0.d0
      IF(nfluid.eq.1)then 
         CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
      ELSEIF(nfluid.eq.2)then 
         CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
      ELSEIF(nfluid.eq.15)then 
         CALL nth2x2_cupid(st_tbl(ndxstd),s,erx) 
      ELSE 
         CALL strpx_cupid(st_tbl(ndxstd),s,erx) 
      ENDIF        
      
      sat_temp=s(1)
!         
      END FUNCTION sat_temp
      

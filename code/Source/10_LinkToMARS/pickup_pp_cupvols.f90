!
      SUBROUTINE pickup_pp_cupvols 
!
!    This routine stores subdomain cupvols' pressure correction, pp,  
!    into global array, ppcup. 
!      
      USE  Zpress , ONLY: pp
      USE  Zmars  , ONLY: marsindex,n_marsbc,ppcup,ppcup_tmp
      USE  Zzone  , ONLY: ncell_fluid
!
      IMPLICIT none
!
!.....Local variables
      INTEGER :: i
!
      DO i=1,ncell_fluid
         IF(marsindex(i).ge.1 .and. marsindex(i).le.n_marsbc) ppcup(marsindex(i))=pp(i) !marsindex(i)=k
      ENDDO
!      
      CALL allreduce_r(ppcup,ppcup_tmp,n_marsbc)
!
      END SUBROUTINE pickup_pp_cupvols       

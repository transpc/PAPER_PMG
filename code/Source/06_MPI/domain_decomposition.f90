!
      SUBROUTINE domain_decomposition(nelem,nelem_sub,ia_sub,ja_sub,celem,celem_max)
!
      USE Zcore     , ONLY: np,myrank
!
!DEC$IF defined (metis_flag) 
!
      USE Zmpi      , ONLY: metis
!
      IMPLICIT NONE
!.....Input
      INTEGER :: nelem,celem_max
      INTEGER :: nelem_sub
      INTEGER :: ia_sub(nelem+1),ja_sub(nelem_sub)
!.....Output
      INTEGER :: celem(nelem)
!.....Local variables
      INTEGER :: moptions(40)
      INTEGER :: edgecut(nelem)
      INTEGER :: ncon(np),vwgt(nelem*1),adjwgt(nelem_sub)
      INTEGER :: celem_offset
!      
      REAL(4) :: tpwgts(np),ubvec(1)
!
!.....METIS pre-requisite procedure
!
!      
!.....Domain partition by METIS
!
      celem=0
      IF (metis.eq.1) THEN
         vwgt(:)=1
         ncon(:)=1
         adjwgt(:)=1
         tpwgts(:)=1.0/FLOAT(np)
         ubvec(:)=1.001
         CALL metis_setdefaultoptions(moptions)
         moptions(17)=1
!        CALL METIS_PartGraphKway(nelem,1,ia_sub,ja_sub,vwgt,0,adjwgt,np,tpwgts,ubvec,moptions,edgecut,celem) ! METIS-5.0,
         CALL METIS_PartGraphrecursive(nelem,1,ia_sub,ja_sub,vwgt,0,adjwgt,np,tpwgts,ubvec,moptions,edgecut,celem) ! METIS-5.0,
         celem_offset=1-MINVAL(celem(1:nelem))
         celem=celem+celem_offset   
         IF(myrank.eq.0) PRINT *,'          METIS domain decomposition finished.'
      ENDIF
!
      celem_max=maxval(celem(:))
!
!DEC$ELSE
      IF(myrank.eq.0)WRITE(*,*)'domain_decomposition is deactivated due to false metis_flag!!!'
      PAUSE
      STOP
!DEC$ENDIF 
!
      END SUBROUTINE 
      

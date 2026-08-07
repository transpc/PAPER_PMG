!
      SUBROUTINE udfn_convert_t2e_liq
!
!     Convert temparature into internal energy
!
      USE STM_TBL_cupid  , ONLY: st_tbl,             &
                                 nt,np,ns,ns2,ndxstd
      USE Zsg            , ONLY: n_group,n_1d,t_1d,p_1d,en_1d
!
      IMPLICIT NONE 
!
      INTEGER :: i,j
      INTEGER :: it    
      LOGICAL :: erx
!
      REAL(8) :: s(36),usubfs
      EQUIVALENCE(s(13),usubfs)
!
      DO i=1,n_group
         DO j=1,n_1d(i)
!     initialize s for sth2x3_cupid
!           s(:)=0.d0
            s(1)=t_1d(i,j)
            s(2)=p_1d(i,j)
            CALL sth2x3_cupid(s,it,erx,                          &
                              st_tbl(ndxstd),                    &
                              st_tbl(ndxstd+nt),                 &
                              st_tbl(ndxstd+nt+np+13*ns+13*ns2))
            IF(erx)then
               print *, '#### ERROR: sth2x3_cupid called from udfn_convert_t2e_liq'
               pause
               stop
            ENDIF
            en_1d(i,j)=usubfs
         ENDDO
      ENDDO 
! 
      END SUBROUTINE udfn_convert_t2e_liq

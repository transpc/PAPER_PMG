!
      SUBROUTINE udfn_prop_in
!
!     Convert temparature into internal energy
!
      USE STM_TBL_cupid , ONLY: st_tbl,             &
                                nt,np,ns,ns2,ndxstd
      USE Zsg      , ONLY: tin_1d,pin_1d,ein_1d,visin_1d,rhoin_1d,cpin_1d
!
      IMPLICIT NONE
!
      INTEGER :: it
      LOGICAL :: erx
!
      REAL(8) :: s(36),usubfs,vsubfs,cpf
      EQUIVALENCE(s(11),vsubfs),  &
                 (s(13),usubfs),  &
                 (s(21),cpf)
!
!     Initialize s  for sth2x3_cupid
!     s(:)=0.d0
      s(1)=tin_1d
      s(2)=pin_1d
      CALL sth2x3_cupid(s,it,erx,                          & 
                        st_tbl(ndxstd),                    &
                        st_tbl(ndxstd+nt),                 &
                        st_tbl(ndxstd+nt+np+13*ns+13*ns2))
      IF(erx)then
         print *, '#### ERROR: sth2x3_cupid called from udfn_prop_in'
         pause
         stop
      ENDIF
      rhoin_1d=1.d0/vsubfs
      ein_1d=usubfs
      cpin_1d=cpf
      CALL udfn_sg_viscos_lw_cupid(tin_1d,rhoin_1d,visin_1d)
!
      RETURN
      END SUBROUTINE udfn_prop_in

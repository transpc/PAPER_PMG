      MODULE c3com_cupid
      SAVE
!      
!.....MCC-jjj-new
!
      INTEGER(4) i3cupid(72),i3invtbl(72),mcdirect(72),mcgdirect(72),nvols_mars
      REAL(8) c3dpv(72,8)
!
!.....Added by C.W.Choi, SCC
!
!      REAL(8) c3vpp(72),c3ngpp(72),c3vt(72),c3lt(72),c3dt(72),c3ngmf(72,10)
!     COMMON/c3com_CUPID/c3dpv,i3cupid,i3invtbl,mcdirect,mcgdirect,nvols_mars
!     c3dpv(i,1): alphag*rhog*eg
!     c3dpv(i,2): alphal*rhol*el
!     c3dpv(i,3): alphag*rhog
!     c3dpv(i,4): alphad*rhol
!     c3dpv(i,5): alphal*rhol
!     c3dpv(i,6): alphag*rhog*x
!     c3dpv(i,7): alphad*rhol*el
!
      ENDMODULE c3com_cupid
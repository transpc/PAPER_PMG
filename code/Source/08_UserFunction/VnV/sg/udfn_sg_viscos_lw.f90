!
      SUBROUTINE udfn_sg_viscos_lw_cupid (tk,rhoc,term)
!
!     SUBROUTINE cond (tk,rhoc,convt) written by WJL on Dec. '97
!     This SUBROUTINE calculates the LIGHT WATER THERMAL CONDUCTIVITY
!     based on ASME '93 Steam Table, AppENDix 7
!     Input quantities are t(K) and rhoc(kg/m3)
!     Output quantities is term(Pa.s)

      IMPLICIT NONE
!
!.....Input
      REAL(8) tk,rhoc,term
!.....Local varaibles
      INTEGER k,ii,jj
      REAL(8) ts1,dt1,r1,dr,sum,temp_,dt1i
!
!.....Data Statements:
!     For Light Water Properties
!
      REAL(8) a(4),bb(6,5) 
      DATA a/0.0181583d0,0.0177624d0,0.0105287d0,-0.0036744d0/ 
      DATA bb/0.501938d0,0.162888d0,-0.130356d0,0.907919d0,-0.551119d0,  &
      0.146543d0,0.235622d0,0.789393d0,0.673665d0,1.207552d0,           &
      0.0670665d0,-0.0843370d0,-0.274637d0,-0.743539d0,-0.959456d0,-    &
      0.687343d0,-0.497089d0,0.195286d0,0.145831d0,0.263129d0,          &
      0.347247d0,0.213486d0,0.100754d0,-0.032932d0,-0.0270448d0,-       &
      0.0253093d0,-0.0267758d0,-0.0822904d0,0.0602253d0,-0.0202595d0/   
!
!.....Test ranges of validity
!
      ts1=647.27d0/tk
      dt1=ts1-1.0d0
      r1=rhoc/317.763d0
      dr=r1-1.0d0
!
      sum=0.d0 
      temp_=1.0d0
      DO  k=1,4 
         sum=sum+a(k)*temp_
         temp_=temp_*ts1
      ENDDO 
      term=1.d-06*sqrt(1.d0/ts1)/sum
!
      IF(tk.eq.647.27d0)  return
!
      sum=0.d0
      dt1i=1.0d0
      DO ii=1,6 
         temp_=1.0d0
         DO jj=1,5 
            sum=sum+bb(ii,jj)*dt1i*temp_
            temp_=temp_*dr
         ENDDO
         dt1i=dt1i*dt1 
      ENDDO
! 
      term=term*exp(r1*sum) 
!
      END SUBROUTINE udfn_sg_viscos_lw_cupid

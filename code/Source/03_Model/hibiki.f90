!
      SUBROUTINE hibiki_bubble_diameter(ag1,ag2,dh,ul,ug,rl,rg,visl,visg,sig, &
      D1,D2,D3,D4,n)
!
!    This routine calculates bubble diameter using Hibiki's bubble diameter model
!
      USE Zparam       , ONLY: ndim
      USE Zvector      , ONLY: vl_o,vg_o
      IMPLICIT NONE
!
!     input
      INTEGER :: n
      REAL(8) :: dh
      REAL(8) :: ag1(*),ag2
      REAL(8) :: rl(*),rg(*),visl(*),visg(*),sig(*)
!     output
      REAL(8) :: ul(*),ug(*)
      REAL(8) :: D1(*),D2(*),D3(*),D4(*)
!     local variables
      INTEGER :: i
      REAL(8) :: ag
      REAL(8) :: al,jg,jl,fil,fig,dpdzl,dpdzg,xtt,phil,dpdz,epsil,lo,Rel,Reg,Reb
      REAL(8) :: t1,t2,t3,t4,t5,t6,t7,t8,xtt2
      REAL(8) :: tag1,tag2,tag3,tag4
!
      ag=0.1d0
      tag2=ag**(0.170d0)
      ag=ag2
      tag3=ag**(0.170d0)
      ag=0.9d0
      tag4=ag**(0.170d0)
!
      IF(ndim.eq.2) THEN
        DO i=1,n
           ul(i)=DSQRT(vl_o(i,1)**2+vl_o(i,2)**2)
           ug(i)=DSQRT(vg_o(i,1)**2+vg_o(i,2)**2)
         ENDDO
      ELSE
        DO i=1,n
           ul(i)=DSQRT(vl_o(i,1)**2+vl_o(i,2)**2+vl_o(i,3)**2)
           ug(i)=DSQRT(vg_o(i,1)**2+vg_o(i,2)**2+vg_o(i,3)**2)
         ENDDO
      ENDIF
!
      DO i=1,n
!
!....Assign parameters of Hibiki's correlations
!     
      lo=DSQRT(sig(i)/9.8d0/DMAX1(1.0d-8,(rl(i)-rg(i))))
      Rel=DMAX1(1.0d0,rl(i)*ul(i)*dh/visl(i))
      Reg=DMAX1(1.0d0,rg(i)*ug(i)*dh/visg(i))
      fil=DMAX1(0.d0,0.046d0*(Rel)**(-0.2))
      fig=DMAX1(0.d0,0.046d0*(Reg)**(-0.2))             
      t1=DEXP(-0.0005839d0*Rel)
      t2=1.d0-t1
      t3=lo*rl(i)/visl(i)
      t4=(lo/dh)**(-0.335d0)
      t6=(rl(i)/rg(i))**(0.138d0) 
      t7=fig/2.0d0*rg(i)*ug(i)**2/dh
      t8=fil/2.0d0*rl(i)*ul(i)**2/dh
!!!!!!!
      ag=ag1(i)
      ag=DMAX1(ag,1.0d-10)
      ag=DMIN1(ag,1.0d0-1.0d-10)
      al=1.0d0-ag
      jg=ag*ug(i)
      jl=al*ul(i)
!     dpdzg=ag**2*(fig/2.0d0*rg(i)*ug(i)**2/dh)
      dpdzg=ag**2*t7
      IF(dpdzg.eq.0.0d0)Then
         dpdz=0.0d0
         epsil=9.81d0*jg*t1
      ELSE
!        dpdzl=al**2*(fil/2.0d0*rl(i)*ul(i)**2/dh)
         dpdzl=al**2*t8
         dpdzl=DMAX1(dpdzl,1.d-10)
         xtt2=dpdzg/dpdzl
         xtt=DSQRT(xtt2)
         phil=(1.0d0+xtt+xtt2)
         dpdz=dpdzl*phil
         epsil=9.81d0*jg*t1+(jg+jl)/(ag*rg(i)+al*rl(i))*dpdz*t2
      ENDIF               
      Reb=(epsil*lo)**(1.0d0/3.0d0)*t3
!
!....Calculate bubble diameter
!     
      IF(Reb.eq.0.0d0)Then
         D1(i)=0.0d0
      ELSE
         t5=Reb**(-0.239d0)
         tag1=ag**(0.170d0)
         D1(i)=lo*1.99d0*t4*t5/1.22d0*tag1*t6
      ENDIF         
!!!!!!!
!!!!!!!
      ag=0.1d0
      al=1.0d0-ag
      jg=ag*ug(i)
      jl=al*ul(i)
!     dpdzg=fig/2.0d0*(ag*rg(i)*ug(i))**2/(rg(i)*dh)
      dpdzg=ag**2*t7
      IF(dpdzg.eq.0.0d0)Then
         dpdz=0.0d0
         epsil=9.81d0*jg*t1
      ELSE
!        dpdzl=fil/2.0d0*(al*rl(i)*ul(i))**2/(rl(i)*dh)
         dpdzl=al**2*t8
         dpdzl=DMAX1(dpdzl,1.d-10)
         xtt2=dpdzg/dpdzl
         xtt=DSQRT(xtt2)
         phil=(1.0d0+xtt+xtt2)
         dpdz=dpdzl*phil
         epsil=9.81d0*jg*t1+(jg+jl)/(ag*rg(i)+al*rl(i))*dpdz*t2
      ENDIF               
      Reb=(epsil*lo)**(1.0d0/3.0d0)*t3
!
!....Calculate bubble diameter
!     
      IF(Reb.eq.0.0d0)Then
         D2(i)=0.0d0
      ELSE
!        t4=(lo/dh)**(-0.335d0)
         t5=Reb**(-0.239d0)
!        tag2=ag**(0.170d0)
!        t6=(rl(i)/rg(i))**(0.138d0) 
         D2(i)=lo*1.99d0*t4*t5/1.22d0*tag2*t6
      ENDIF         
!!!!!!!
!!!!!!!
      ag=ag2
      al=1.0d0-ag
      jg=ag*ug(i)
      jl=al*ul(i)
!     dpdzg=fig/2.0d0*(ag*rg(i)*ug(i))**2/(rg(i)*dh)
      dpdzg=ag**2*t7
      IF(dpdzg.eq.0.0d0)Then
         dpdz=0.0d0
         epsil=9.81d0*jg*t1
      ELSE
!        dpdzl=fil/2.0d0*(al*rl(i)*ul(i))**2/(rl(i)*dh)
         dpdzl=al**2*t8
         dpdzl=DMAX1(dpdzl,1.d-10)
         xtt2=dpdzg/dpdzl
         xtt=DSQRT(xtt2)
         phil=(1.0d0+xtt+xtt2)
         dpdz=dpdzl*phil
         epsil=9.81d0*jg*t1+(jg+jl)/(ag*rg(i)+al*rl(i))*dpdz*t2
      ENDIF               
      Reb=(epsil*lo)**(1.0d0/3.0d0)*t3
!
!....Calculate bubble diameter
!     
      IF(Reb.eq.0.0d0)Then
         D3(i)=0.0d0
      ELSE
!        t4=(lo/dh)**(-0.335d0)
         t5=Reb**(-0.239d0)
!        tag3=ag**(0.170d0)
!        t6=(rl(i)/rg(i))**(0.138d0) 
         D3(i)=lo*1.99d0*t4*t5/1.22d0*tag3*t6
      ENDIF         
!!!!!!!
!!!!!!!
      ag=0.9d0
      al=1.0d0-ag
      jg=ag*ug(i)
      jl=al*ul(i)
!     dpdzg=fig/2.0d0*(ag*rg(i)*ug(i))**2/(rg(i)*dh)
      dpdzg=ag**2*t7
      IF(dpdzg.eq.0.0d0)Then
         dpdz=0.0d0
         epsil=9.81d0*jg*t1
      ELSE
!        dpdzl=fil/2.0d0*(al*rl(i)*ul(i))**2/(rl(i)*dh)
         dpdzl=al**2*t8
         dpdzl=DMAX1(dpdzl,1.d-10)
         xtt2=dpdzg/dpdzl
         xtt=DSQRT(xtt2)
         phil=(1.0d0+xtt+xtt2)
         dpdz=dpdzl*phil
         epsil=9.81d0*jg*t1+(jg+jl)/(ag*rg(i)+al*rl(i))*dpdz*t2
      ENDIF               
      Reb=(epsil*lo)**(1.0d0/3.0d0)*t3
!
!....Calculate bubble diameter
!     
      IF(Reb.eq.0.0d0)Then
         D4(i)=0.0d0
      ELSE
!        t4=(lo/dh)**(-0.335d0)
         t5=Reb**(-0.239d0)
!        tag4=ag**(0.170d0)
!        t6=(rl(i)/rg(i))**(0.138d0) 
         D4(i)=lo*1.99d0*t4*t5/1.22d0*tag4*t6
      ENDIF         
!!!!!!!
      ENDDO 
!
      RETURN
      END SUBROUTINE hibiki_bubble_diameter
 

 ! = = = = = = = = = = = = = = = = = = = = = = = = =
!
      SUBROUTINE P_distance(ndim, jmax, nnodef, nnodec, coordf, coordc, nnz, ia, ja, au)
!
         IMPLICIT NONE
! input
         INTEGER jmax, nnz, ndim
         INTEGER nnodef, nnodec
         INTEGER ia(nnodef+1), ja(nnz)
         REAL(8) coordf(ndim, nnodef), coordc(ndim, nnodec)
! output
         REAL(8) au(nnz)

! free
         INTEGER i, j, k, l, i1, i2, k1, k2, l1, l2, j1, j2, ie, ne, ix
         INTEGER nnd, id, id1, id2, imax
         INTEGER ni(jmax)
         REAL(8) xtmp, xc(ndim), dx(jmax), dd(jmax)
! ---

         ni = 0
         dx = 0.d0
         dd = 0.d0
         
         DO ie = 1, nnodef

            i1 = ia(ie)
            i2 = ia(ie+1)-1
            imax = i2-i1+1

            IF(imax.EQ.1) THEN
                au(i1) = 1.d0
                CYCLE
            ENDIF
!
            xc(1:ndim) = coordf(1:ndim, ie)
            ni(1:imax) = ja(i1:i2)
!            dx(1:imax) = 0.d0

            DO i = 1, imax
               ne = ni(i)
               dx(i) = (xc(1)-coordc(1, ne))**2.d0+(xc(2)-coordc(2, ne))**2.d0
               IF (ndim .EQ. 3) THEN
                  dx(i) = dx(i)+(xc(3)-coordc(3, ne))**2.d0
               END IF
!               dx(i) = (dx(i))
            END DO

            dd(1:imax) = 1.d0
            DO i = 1, imax
               DO j = 1, imax
                  IF (j .eq. i) CYCLE
                  dd(i) = dd(i)*dx(j)
               END DO
            END DO

            xtmp = SUM(dd(1:imax))
            au(i1:i2) = dd(1:imax)/xtmp

         END DO
! -----------------------!
         RETURN

      END

! - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - !
    
    Subroutine P_linear_2D(jmax,nnodef,nnodec,coordf,coordc,nnz,ia,ja,au)
!
      IMPLICIT NONE
! input
      INTEGER(4) jmax,nnz
      INTEGER(4) nnodef,nnodec
      INTEGER(4) ia(nnodef+1),ja(nnz)
      REAL(8) coordf(2,nnodef),coordc(2,nnodec)
! output
      REAL(8) au(nnz)

! free 
      INTEGER(4) i,j,k,l,i1,i2,k1,k2,l1,l2,j1,j2,j3,ie,ne,ix,insd
      INTEGER(4) nnd,id,id1,id2,imax
      INTEGER(4) ni(jmax),mi(jmax)
      REAL(8) xtmp,xc(2),dx(jmax),dd(jmax),x0(2,3) 
      REAL(8) small,shape_f,shape0,shape1,shape2,shape3	  

	  
! 
      small = 1.d-8
      dx = 0.d0
      ni = 0
      mi = 0
! ---	  
	  
      DO ie = 1,nnodef

         xc(1:2) = coordf(1:2,ie)
		  			
          i1 = ia(ie)
		  i2 = ia(ie+1)-1
		  imax = i2-i1+1
          
          ni(1:imax) = ja(i1:i2)
!		  dx = 0.d0
!		  
        IF(imax.eq.1) THEN
          au(i1) = 1.d0
!        
        ELSEIF(imax.eq.2) THEN
          dx(1:2) = (xc(1)-coordc(1,ni(1:2)))**2.d0+(xc(2)-coordc(2,ni(1:2)))**2.d0
          dx(1:2) = DSQRT(dx(1:2))
          au(i1) = dx(2)/(dx(1)+dx(2))
          au(i2) = dx(1)/(dx(1)+dx(2))
!         
        ELSE 
! imax >= 3 ! linear 
          xtmp = 1.d8
		  
          DO i = 1,imax
             ne = ni(i)
             dx(i) = (xc(1)-coordc(1,ne))**2.d0+(xc(2)-coordc(2,ne))**2.d0
             dx(i)=DSQRT(dx(i))
             if(dx(i).lt.xtmp) then
             xtmp = dx(i)
             j1 = i
             endif
			 
          ENDDO

          IF (xtmp.le.small) THEN    ! same position 
            au(i1+j1-1) = 1.d0
		    CYCLE
          ENDIF
		  
! ordering by distance 
          mi(1:imax) = [1:imax]
          CALL bubble_sort_real(imax,mi,dx)
		  
          j1 = ni(mi(1))
          j2 = ni(mi(2))
          x0(1:2,1) = coordc(1:2,j1)
          x0(1:2,2) = coordc(1:2,j2)          
!
          do j = 3,imax
		    j3 = ni(mi(j))
            x0(1:2,3) = coordc(1:2,j3)

!check inside 
            insd = 0
            call check_node(xc(1),xc(2),x0,insd)
            if(insd.eq.1) then
            shape0 = shape_f(x0(1,1),x0(2,1),x0(1,2),x0(2,2),x0(1,3),x0(2,3))
            if(shape0.lt.small) then    
! j1,j2,j3: aline 
             au(i1+mi(1)-1) = dx(mi(2))/(dx(mi(1))+dx(mi(2)))
             au(i1+mi(2)-1) = dx(mi(1))/(dx(mi(1))+dx(mi(2)))
                
            else
            
             shape1 = shape_f(xc(1),xc(2),x0(1,2),x0(2,2),x0(1,3),x0(2,3))
             shape2 = shape_f(xc(1),xc(2),x0(1,1),x0(2,1),x0(1,3),x0(2,3))
             shape1 = shape1/shape0
             shape2 = shape2/shape0
             shape3 = 1.d0 - shape1 - shape2
! 
             au(i1+mi(1)-1) = shape1
             au(i1+mi(2)-1) = shape2
             au(i1+mi(j)-1) = shape3
!
             
            endif
            
            exit			
            endif
			
          end do

! if outside 		  
          if(insd.eq.0) then 
             au(i1+mi(1)-1) = dx(mi(2))/(dx(mi(1))+dx(mi(2)))
             au(i1+mi(2)-1) = dx(mi(1))/(dx(mi(1))+dx(mi(2)))
!   
          endif
		  
        ENDIF
		
      ENDDO
	   
! -----------------------------!
      return
      end
! ------------------------------------------------------------!
! ===========================end sub =========================!
! ------------------------------------------------------------!
! ------------ area shape function in triangle -------
       real*8 function shape_f(x1,y1,x2,y2,x3,y3)
       real*8 x1,y1,x2,y2,x3,y3
       real*8 d1,d2,d3,p
! ---
       d1 = dsqrt((x2-x3)**2.0+(y2-y3)**2.0)
       d2 = dsqrt((x1-x3)**2.0+(y1-y3)**2.0)
       d3 = dsqrt((x1-x2)**2.0+(y1-y2)**2.0)
       p = 0.5d0*(d1+d2+d3)
        shape_f = dsqrt(abs(p*(p-d1)*(p-d2)*(p-d3)))
       end function      
! ---   sub --check_node -----!
      subroutine check_node(x,y,xc,id)
      implicit none
! ---
      integer id
      real*8 x,y,xc(2,3)
! 
      real*8 xm,ym,f_line,xtest,small
! ------------------------------------!
      id = 0 
      small = -1.d-25
	  
! center ---
      xm = (xc(1,1)+xc(1,2)+xc(1,3))/3.d0
      ym = (xc(2,1)+xc(2,2)+xc(2,3))/3.d0
! edge 1-2
      xtest = f_line(xm,ym,xc(1,1),xc(2,1),xc(1,2),xc(2,2))
      xtest = f_line(x,y,xc(1,1),xc(2,1),xc(1,2),xc(2,2))*xtest
      if(xtest.lt.small) goto 10
! edge -2-3
      xtest = f_line(xm,ym,xc(1,2),xc(2,2),xc(1,3),xc(2,3))
      xtest = f_line(x,y,xc(1,2),xc(2,2),xc(1,3),xc(2,3))*xtest
      if(xtest.lt.small) goto 10
! edge -1-3
      xtest = f_line(xm,ym,xc(1,1),xc(2,1),xc(1,3),xc(2,3))
      xtest = f_line(x,y,xc(1,1),xc(2,1),xc(1,3),xc(2,3))*xtest
      if(xtest.lt.small) goto 10
      id = 1
! ----
10    continue
! ----------------------------!
      return
      end
! ------------------------------------------------------------!
! ----------------------function------------------------------!
       real*8 function f_line(x,y,x1,y1,x2,y2)
       real*8 x,y,x1,y1,x2,y2
! ---
       f_line = (y2-y1)*(x-x1) - (x2-x1)*(y-y1)
       end function

!  - - - - - - - - - - - - - - - - - - - - 
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
      Subroutine P_linear_3D(jmax,nnodef,nnodec,coordf,coordc,nnz,ia,ja,au)
!
      IMPLICIT NONE
! input
      INTEGER(4) jmax,nnz
      INTEGER(4) nnodef,nnodec
      INTEGER(4) ia(nnodef+1),ja(nnz)
      REAL(8) coordf(3,nnodef),coordc(3,nnodec)

! output
      REAL(8) au(nnz)

! free 
      INTEGER(4) i,j,k,l,i1,i2,i3,k1,k2,l1,l2,j1,j2,j3,j4,ie,ne,ix,insd
      INTEGER(4) nnd,id,id1,id2,imax
      INTEGER(4) ni(jmax),mi(jmax)
      REAL(8) xtmp,xc(3),dx(jmax),dd(jmax),x0(3,4) 
      REAL(8) small,volm,shape0,shape1,shape2,shape3	  

	  
! 
      small = 1.d-10
      ni = 0
      mi = 0
      dx = 0.d0
! ---	  
	  
      DO ie = 1,nnodef

         xc(1:3) = coordf(1:3,ie)
		  			
          i1 = ia(ie)
		  i2 = ia(ie+1)-1
		  imax = i2-i1+1
          
          ni(1:imax) = ja(i1:i2)
!		  
        IF(imax.eq.1) THEN
          au(i1) = 1.d0
!         
        ELSEIF(imax.eq.2) THEN
          dx(1:2) = (xc(1)-coordc(1,ni(1:2)))**2.d0+(xc(2)-coordc(2,ni(1:2)))**2.d0+(xc(3)-coordc(3,ni(1:2)))**2.d0
          dx(1:2) = DSQRT(dx(1:2))
          au(i1) = dx(2)/(dx(1)+dx(2))
          au(i2) = dx(1)/(dx(1)+dx(2))
!          
        ELSEIF(imax.eq.3) THEN 
          dx(1:3) = (xc(1)-coordc(1,ni(1:3)))**2.d0+(xc(2)-coordc(2,ni(1:3)))**2.d0+(xc(3)-coordc(3,ni(1:3)))**2.d0
          dx(1:3) = DSQRT(dx(1:3))
          xtmp = dx(1)*dx(2)+dx(2)*dx(3)+dx(1)*dx(3)
          
          au(i1) = dx(2)*dx(3)/xtmp
          au(i1+1) = dx(1)*dx(3)/xtmp
          au(i1+2) = dx(1)*dx(2)/xtmp
        ELSE
! imax >= 4 ! linear 
          xtmp = 1.d6
		  
          DO i = 1,imax
             ne = ni(i)
             dx(i) = (xc(1)-coordc(1,ne))**2.d0+(xc(2)-coordc(2,ne))**2.d0+(xc(3)-coordc(3,ne))**2.d0
             dx(i)=DSQRT(dx(i))
             if(dx(i).lt.xtmp) then
             xtmp = dx(i)
             j1 = i
             endif
			 
          ENDDO

          IF (xtmp.le.small) THEN    ! same position 
            au(i1+j1-1) = 1.d0
		    CYCLE
          ENDIF
		  
! ordering by distance 
          mi(1:imax) = [1:imax]
          CALL bubble_sort_real(imax,mi,dx)
		  
          j1 = ni(mi(1))
          j2 = ni(mi(2))
		  j3 = ni(mi(3))
          x0(1:3,1) = coordc(1:3,j1)
          x0(1:3,2) = coordc(1:3,j2)
          x0(1:3,3) = coordc(1:3,j3)          
!
          do j = 4,imax
		    j4 = ni(mi(j))
            x0(1:3,4) = coordc(1:3,j4)

!check inside 
            insd = 0
            call check_node_3D(xc(1),xc(2),xc(3),x0,insd)
            if(insd.eq.1) then
            shape0 = volm(x0(1,1),x0(2,1),x0(3,1),x0(1,2),x0(2,2),x0(3,2),x0(1,3),x0(2,3),x0(3,3),x0(1,4),x0(2,4),x0(3,4))
            if(shape0.lt.small) then    
! j1,j2,j3,j4: on a plan 
             xtmp = dx(mi(1))*dx(mi(2))+dx(mi(2))*dx(mi(3))+dx(mi(1))*dx(mi(3))
             au(i1+mi(1)-1) = dx(mi(2))*dx(mi(3))/xtmp
             au(i1+mi(2)-1) = dx(mi(1))*dx(mi(3))/xtmp
			 au(i1+mi(3)-1) = dx(mi(1))*dx(mi(2))/xtmp
                
            else
            
             shape1 = volm(xc(1),xc(2),xc(3),x0(1,2),x0(2,2),x0(3,2),x0(1,3),x0(2,3),x0(3,3),x0(1,4),x0(2,4),x0(3,4))
             shape2 = volm(xc(1),xc(2),xc(3),x0(1,1),x0(2,1),x0(3,1),x0(1,3),x0(2,3),x0(3,3),x0(1,4),x0(2,4),x0(3,4))
             shape3 = volm(xc(1),xc(2),xc(3),x0(1,1),x0(2,1),x0(3,1),x0(1,2),x0(2,2),x0(3,2),x0(1,4),x0(2,4),x0(3,4))			 
             shape1 = shape1/shape0
             shape2 = shape2/shape0
             shape3 = shape3/shape0			 
   
! 
             au(i1+mi(1)-1) = shape1
             au(i1+mi(2)-1) = shape2
             au(i1+mi(3)-1) = shape3			 
             au(i1+mi(j)-1) = 1.d0 - shape1 - shape2 - shape3
             
            endif
            
            exit			
            endif
			
          end do

! if outside 		  
          if(insd.eq.0) then 
             xtmp = dx(mi(1))*dx(mi(2))+dx(mi(2))*dx(mi(3))+dx(mi(1))*dx(mi(3))
             au(i1+mi(1)-1) = dx(mi(2))*dx(mi(3))/xtmp
             au(i1+mi(2)-1) = dx(mi(1))*dx(mi(3))/xtmp
			 au(i1+mi(3)-1) = dx(mi(1))*dx(mi(2))/xtmp
!            
          endif
		  
        ENDIF
		
      ENDDO
	   
! -----------------------------!
      return
      end
! ------------------------------------------------------------!
! ===========================end sub =========================!
! ------------------------------------------------------------!

! ------------------------------------------------------------ !
			
! ------------------------------------------------------------!
! ===========================end sub =========================!
! ------------------------------------------------------------!
! ------------ area shape function in triangle -------
!       real*8 function volm(x(1,1),x(2,1),x(3,1),x(1,2),x(2,2),x(3,2),x(1,3),x(2,3),x(3,3),x(1,4),x(2,4),x(3,4))
       real*8 function volm(x11,x21,x31,x12,x22,x32,x13,x23,x33,x14,x24,x34)
       real*8 a(3,3),xj,x11,x21,x31,x12,x22,x32,x13,x23,x33,x14,x24,x34
! ---
       
       a(1,1)=x11-x13
       a(1,2)=x12-x13
       a(1,3)=x14-x13
       a(2,1)=x21-x23
       a(2,2)=x22-x23
       a(2,3)=x24-x23
       a(3,1)=x31-x33
       a(3,2)=x32-x33
       a(3,3)=x34-x33
!---
       xj=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(3,1)*a(2,3)+a(1,3)*a(2,1)*a(3,2) &
         -a(1,1)*a(3,2)*a(2,3)-a(1,2)*a(2,1)*a(3,3)-a(1,3)*a(3,1)*a(2,2)

       volm = dabs(xj)/6.d0
       
       end function      
    
    
! ---   sub --check_node -----------------------------------!
      subroutine check_node_3D(x,y,z,xc,id)
      implicit none
! ---
      integer id
      real*8 x,y,z,xc(3,4)
! 
      real*8 xm,ym,zm,f_face,xtest
      real*8 small
! ------------------------------------!
      id = 0 
      small = -1.d-25
! center ---
      xm = (xc(1,1)+xc(1,2)+xc(1,3)+xc(1,4))/4.d0
      ym = (xc(2,1)+xc(2,2)+xc(2,3)+xc(2,4))/4.d0
      zm = (xc(3,1)+xc(3,2)+xc(3,3)+xc(3,4))/4.d0
! face 1-2-3
      xtest = f_face(xm,ym,zm,xc(1,1),xc(2,1),xc(3,1),xc(1,2),xc(2,2),xc(3,2),xc(1,3),xc(2,3),xc(3,3))
      xtest = f_face(x,y,z,xc(1,1),xc(2,1),xc(3,1),xc(1,2),xc(2,2),xc(3,2),xc(1,3),xc(2,3),xc(3,3))*xtest
      if(xtest.lt.small) goto 10
! face -2-3-4
      xtest = f_face(xm,ym,zm,xc(1,2),xc(2,2),xc(3,2),xc(1,3),xc(2,3),xc(3,3),xc(1,4),xc(2,4),xc(3,4))
      xtest = f_face(x,y,z,xc(1,2),xc(2,2),xc(3,2),xc(1,3),xc(2,3),xc(3,3),xc(1,4),xc(2,4),xc(3,4))*xtest
      if(xtest.lt.small) goto 10
! face 1-2-4
      xtest = f_face(xm,ym,zm,xc(1,1),xc(2,1),xc(3,1),xc(1,2),xc(2,2),xc(3,2),xc(1,4),xc(2,4),xc(3,4))
      xtest = f_face(x,y,z,xc(1,1),xc(2,1),xc(3,1),xc(1,2),xc(2,2),xc(3,2),xc(1,4),xc(2,4),xc(3,4))*xtest
      if(xtest.lt.small) goto 10
! face -1-3-4
      xtest = f_face(xm,ym,zm,xc(1,1),xc(2,1),xc(3,1),xc(1,3),xc(2,3),xc(3,3),xc(1,4),xc(2,4),xc(3,4))
      xtest = f_face(x,y,z,xc(1,1),xc(2,1),xc(3,1),xc(1,3),xc(2,3),xc(3,3),xc(1,4),xc(2,4),xc(3,4))*xtest
      if(xtest.lt.small) goto 10
      id = 1
! ----
10    continue
! ----------------------------!
      return
      end
! ------------------------------------------------------------!
! ----------------------function------------------------------!
       real*8 function f_face(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3)
       real*8 x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3
! ---
       f_face = (x-x1)*((y2-y1)*(z3-z1)-(y3-y1)*(z2-z1))-(y-y1)*((x2-x1)*(z3-z1)-(x3-x1)*(z2-z1))+       &
                (z-z1)*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
       
       end function
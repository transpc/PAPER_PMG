      subroutine Relax_GSP(maxit,relax,n,nintf,nnz,ia,ja,ju,au,u,b)
! ---
      implicit none
! ---
      integer n,maxit,nnz,nintf
      integer i,j,k,iter,j1,j2
      integer ia(*)
      integer ja(*)
      integer ju(*)
      real*8 temp,temp1,relax
      real*8 u(n),b(n),au(nnz)    
 !     
! ---
!---Gauss-Seidel method -Parallel
                  
  Do iter=1,maxit
 !  
   Do i=1,nintf
      j1 = ia(i)
	  j2 = ia(i+1)-1
      j = ju(i)
      temp1 = au(j)
! ---
      temp = 0.d0
      do j=j1,j2
      temp = temp + au(j)*u(ja(j))
      enddo

!     temp = dot_product(au(j1:j2),v(ja(j1:j2)))
	 u(i)=(b(i)-temp)/temp1 + u(i)
   End do 
   
  End do 
!  
!  deallocate(v)
      return
    end
    
! - - - - - - - - - - - - - - - - - - - -- 
      subroutine Relax_GS(maxit,relax,n,nnz,ia,ja,ju,au,u,b)
      
      use omp_lib
!      USE MD_OpenMP    
! ---
      implicit none
! ---
      integer n,maxit,nnz
      integer i,j,k,iter,j1,j2
      integer ia(*),ja(*),ju(*)
      real*8 temp,temp1,relax
      real*8 u(*),b(*),au(*)    

! ---
!---Gauss-Seidel method -Parallel
  Do iter=1,maxit

!$omp PARALLEL DO private(j1,j2,temp,temp1,j)
   Do i=1,n
      j1 = ia(i)
	  j2 = ia(i+1)-1
      j = ju(i)
      temp1 = au(j)
! ---
     temp = dot_product(au(j1:j2),u(ja(j1:j2)))
	 u(i)=(b(i)-temp)/temp1 + u(i)
   End do 
       
!$omp end PARALLEL DO
   
  End do 
  
      return
    end
    
! 
! - - - - - - - - - - - - - - - - - - - -- 
! - - - - - - - - - - - - - - - - - - - -- 
      subroutine Relax_GS_MPI(maxit,relax,n,nnz,ia,ja,ju,au,u,b)
      
!      use omp_lib
!      USE MD_OpenMP    
! ---
      implicit none
! ---
      integer n,maxit,nnz
      integer i,j,k,iter,j1,j2
      integer ia(*),ja(*),ju(*)
      real*8 temp,temp1,relax
      real*8 u(*),b(*),au(*)    

! ---
!---Gauss-Seidel method -Parallel
  Do iter=1,maxit

!!$omp PARALLEL DO private(j1,j2,temp,temp1,j)
   Do i=1,n
      j1 = ia(i)
	  j2 = ia(i+1)-1
      j = ju(i)
      temp1 = au(j)
! ---
     temp = dot_product(au(j1:j2),u(ja(j1:j2)))
	 u(i)=(b(i)-temp)/temp1 + u(i)
   End do 
       
!!$omp end PARALLEL DO
   
  End do 
  
      return
    end
    
! 
! - - - - - - - - - - - - - - - - - - - -- 
      subroutine Relax_GS_SYM(id,maxit,relax,n,nnz,ia,ja,ju,au,u,b)
! ---
      implicit none
! ---
      integer id,n,maxit,nnz
      integer i,j,k,iter,j1,j2
      integer ia(*),ja(*),ju(*)
      real*8 temp,temp1,relax
      real*8 u(*),b(*),au(*)    
      INTEGER(4) ista, iend,iste

! ---
!---Gauss-Seidel method
  Do iter=1,maxit
      IF(iter.GT.1) id = 1-id
      IF(id.EQ.0) THEN
          ista = 1
          iend = n
          iste = 1
      ELSE
          ista = n
          iend = 1
          iste = -1
      ENDIF
!
   Do i=ista,iend,iste
      j1 = ia(i)
	  j2 = ia(i+1)-1
      j = ju(i)
      temp1 = au(j)
! ---
     temp = dot_product(au(j1:j2),u(ja(j1:j2)))
	 u(i)=(b(i)-temp)/temp1 + u(i)
   End do 
       
  End do 
  
      return
    end
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
SUBROUTINE Smooth_GS(maxit,ista,iend,b,u,au,ia,ja,ju)
! ---
! note that it is for GS forward 
! maxit = 1
! 
use omp_lib
!USE MD_OpenMP
      
IMPLICIT NONE
! ---
INTEGER(4) maxit,ista,iend
INTEGER(4) ia(*)
INTEGER(4) ja(*)
INTEGER(4) ju(*)
REAL(8) b(*)
REAL(8) u(*)
REAL(8) au(*)
! temp
INTEGER(4) i,j,iter,j1,j2,k
REAL(8) temp,temp1   
! ---
!---Gauss-Seidel method
!  DO iter=1,maxit

!$omp PARALLEL DO private(j1,j2,j,temp1,temp,k) 

   DO i=ista,iend
      j1 = ia(i)
	  j2 = ia(i+1)-1
      j = ju(i)
      temp1 = au(j)
!
      temp = b(i)
      do k=j1,j2
         temp = temp -au(k)*u(ja(k)) 
      enddo
      
! ---
!    temp = DOT_PRODUCT( au(j1:j2),u(ja(j1:j2)) )
    u(i)=temp/temp1 + u(i) !(b(i)-temp)/temp1 + u(i)
   ENDDO 
  
!$omp end PARALLEL DO
!ENDDO 
  
RETURN
    END
! = = = = = = = = = = = = = = = = = = = = 
SUBROUTINE Smooth_GS2(maxit,ista,iend,b,u,au,ia,ja,ju,diagr)
! ---
! note that it is for GS forward 
! maxit = 1
use omp_lib
! 
IMPLICIT NONE
! ---
INTEGER(4) maxit,ista,iend
INTEGER(4) ia(*)
INTEGER(4) ja(*)
INTEGER(4) ju(*)
REAL(8) b(*)
REAL(8) u(*)
REAL(8) au(*)
REAL(8) diagr(*)
! temp
INTEGER(4) i,j,iter,j1,j2,k
REAL(8) temp,temp1   
! ---
!---Gauss-Seidel method
!  DO iter=1,maxit
!$omp PARALLEL DO private(j1,j2,j,temp,k) 
   DO i=ista,iend
      j1 = ia(i)
	  j2 = ia(i+1)-1
!      j = ju(i)
!      temp1 = au(j)
!
      temp = b(i)
      do k=j1,j2
         j = ja(k)
         temp = temp -au(k)*u(j) 
      enddo
      
! ---
!    temp = DOT_PRODUCT( au(j1:j2),u(ja(j1:j2)) )
    u(i)=temp*diagr(i) + u(i) !(b(i)-temp)/temp1 + u(i)
   ENDDO 
     
!$omp end PARALLEL DO
!ENDDO 
  
RETURN
    END
! = = = = = = = = = = 
! = = = = = = = = = = = = = = = = = = = = 
SUBROUTINE Smooth_GS2_MPI(maxit,ista,iend,b,u,au,ia,ja,ju,diagr)
! ---
! note that it is for GS forward 
! maxit = 1
!use omp_lib
! 
IMPLICIT NONE
! ---
INTEGER(4) maxit,ista,iend
INTEGER(4) ia(*)
INTEGER(4) ja(*)
INTEGER(4) ju(*)
REAL(8) b(*)
REAL(8) u(*)
REAL(8) au(*)
REAL(8) diagr(*)
! temp
INTEGER(4) i,j,iter,j1,j2,k
REAL(8) temp,temp1   
! ---
!---Gauss-Seidel method
!  DO iter=1,maxit
!!$omp PARALLEL DO private(j1,j2,j,temp,k) 
   DO i=ista,iend
      j1 = ia(i)
	  j2 = ia(i+1)-1
!      j = ju(i)
!      temp1 = au(j)
!
      temp = b(i)
      do k=j1,j2
         j = ja(k)
         temp = temp -au(k)*u(j) 
      enddo
      
! ---
!    temp = DOT_PRODUCT( au(j1:j2),u(ja(j1:j2)) )
    u(i)=temp*diagr(i) + u(i) !(b(i)-temp)/temp1 + u(i)
   ENDDO 
     
!!$omp end PARALLEL DO
!ENDDO 
  
RETURN
    END
! = = = = = = = = = = 
    
SUBROUTINE Smooth_GS_BW(maxit,ista,iend,b,u,au,ia,ja,ju)
! ---
! note that it is for GS backward
! maxit = 1
! 
IMPLICIT NONE
! ---
INTEGER(4) maxit,ista,iend
INTEGER(4) ia(*)
INTEGER(4) ja(*)
INTEGER(4) ju(*)
REAL(8) b(*)
REAL(8) u(*)
REAL(8) au(*)
! temp
INTEGER(4) i,j,iter,j1,j2,k
REAL(8) temp,temp1   
! ---
!---Gauss-Seidel method
!  DO iter=1,maxit
   DO i=iend,ista,-1
      j1 = ia(i)
	  j2 = ia(i+1)-1
      j = ju(i)
      temp1 = au(j)
!
      temp = b(i)
      do k=j1,j2
         temp = temp -au(k)*u(ja(k)) 
      enddo
      
! ---
!    temp = DOT_PRODUCT( au(j1:j2),u(ja(j1:j2)) )
    u(i)=temp/temp1 + u(i) !(b(i)-temp)/temp1 + u(i)
ENDDO 
       
!ENDDO 
  
RETURN
END
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
      SUBROUTINE smoothing_fine(iter,isth,ndom,relax,nintf,nnode,nnz,ia,ja,ju,           &
                            au,u,b,nnbd,nbdom,spt,rpt,sintf,rintf)
      
      USE MD_matrix, ONLY: diagr,alu
      USE MD_MG_index, ONLY: id_GS_sym,icheb,rcheb
      
      IMPLICIT NONE
	  
      INTEGER(4) isth,ndom,iter
      INTEGER(4) nintf,nnode,nnz,nnbd
      INTEGER(4) ia(nnode+1),ja(nnz),ju(nnode)
      INTEGER(4) spt(nnbd+1),rpt(nnbd+1),sintf(spt(nnbd+1)-1),rintf(rpt(nnbd+1)-1),nbdom(nnbd)
      REAL(8) relax
      REAL(8) au(nnz),u(nnode),b(nnode)

      INTEGER(4) i
	 
! - - - - - - - - - - - 
      IF(isth.EQ.0) THEN
!      CALL Relax_GSP(1,relax,nnode,nintf,nnz,ia,ja,ju,au,u,b)     ! GAuss-seidel smoothing \
!      ELSEIF(isth.EQ.1) THEN

      DO i = 1,iter
          
      IF((MOD(i,2).EQ.1).OR.(id_GS_sym.EQ.0)) THEN
          
      CALL Relax_GS0P(1,relax,nnode,nintf,nnz,ia,ja,au,u,b,diagr)
      ELSE
      CALL Relax_GS0P_BW(1,relax,nnode,nintf,nnz,ia,ja,au,u,b,diagr)
      ENDIF
      
      
      IF(nnbd.NE.0) THEN
      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)
      ENDIF
      
      Enddo

      ELSEIF(isth == 2) THEN
          
      DO i = 1,iter
          
          call Smooth_ILU(ndom,nintf,nnode,nnbd,nnz,ia,ja,ju,       &   ! ILU- smoothing 
                  spt,rpt,sintf,rintf,nbdom,au,alu,b,u) 
          
      IF(nnbd.NE.0) THEN
      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)
      ENDIF
          
     enddo
          
      ELSEIF(isth == 3) then
          
      DO i = 1,iter
          
         call poly_cheb_smooth(ndom,icheb,rcheb,nintf,nnode,nnz,ia,ja,au,b,nnbd,spt,rpt,nbdom,sintf,rintf,u)       ! polynormial chebyshev
             
      IF(nnbd.NE.0) THEN
      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)
      ENDIF
      
     enddo
      ELSE
          
      DO i = 1,iter
    
      CALL Smooth_CG(ndom,nintf,nnode,nnbd,nnz,ia,ja,ju,       &   ! CG- smoothing 
                  spt,rpt,sintf,rintf,nbdom,au,b,u) 
      IF(nnbd.NE.0) THEN
      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)
      ENDIF
      enddo
      
      ENDIF
	 
! 
      RETURN 
      END
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      SUBROUTINE Smooth_CG(ndom,nintf,neq,niut,nnz,ia,ja,ju,     &
                          si,ri,sintf,rintf,iut,au,arhsu,solu)
!
      IMPLICIT NONE
      
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
!
      INTEGER(4) ndom
      INTEGER(4) nintf,neq,niut,nnz
      INTEGER(4) ia(neq+1),ju(neq),ja(nnz)
      INTEGER(4) si(niut+1),ri(niut+1),iut(niut),sintf(si(niut+1)-1),rintf(ri(niut+1)-1)
      REAL(8) au(nnz),arhsu(neq)
! 
      REAL*8 solu(neq)
! 
!  temp.
      INTEGER(4) i,j,n,Allocatestatus,ierr,np
      REAL*8 bknum,bkden,akden,bk,ak,bknum1,akden1
      REAL*8 err0,err1,err,res0,res1
      REAL*8,DIMENSION(:),ALLOCATABLE::r,z
! 

      n=neq
      np = ndom
      Allocate(r(n),z(n),stat=Allocatestatus)
      if(AllocateStatus/=0) stop "**Not enough memory " 

!-----Exchange solu to compute A*solu
      call amux0_PCG(nintf,n,nnz,solu,r,au,ja,ia)
!-----------------------------------------------
        bknum=0.d0
        DO i=1,nintf
        r(i) = arhsu(i)-r(i)	 
	    bknum=bknum+r(i)*r(i)
        ENDDO
!-----
!--	
!DEC$IF defined (mpi_flag)
        call MPI_ALLREDUCE(bknum,bknum1,1,MPI_DOUBLE_PRECISION,   &
                    MPI_SUM,MPI_COMM_WORLD,ierr)
!DEC$ENDIF

      IF(np.GT.1) THEN
      bknum=bknum1
      ELSE
      bknum = bknum
      ENDIF 

!       r = z  ! p0=> r
!-----Exchange p0 to compute A*p0: p0=>r
      CALL send_receive(niut,neq,si,ri,sintf,rintf,iut,r)

      CALL amux0_PCG(nintf,n,nnz,r,z,au,ja,ia) !!A*p0=z
!-----------------------------------------------------

        akden=0.d0
      DO i=1,nintf
         akden=akden+r(i)*z(i)
      ENDDO
		
!DEC$IF defined (mpi_flag)
      CALL MPI_ALLREDUCE(akden,akden1,1,MPI_DOUBLE_PRECISION,  &
                   MPI_SUM,MPI_COMM_WORLD,ierr)
!DEC$ENDIF
        
      IF(np.GT.1) THEN
      akden1=akden1
      ELSE
      akden1 = akden
      ENDIF
      
        ak = bknum/akden1

       DO i=1,nintf
          solu(i) = solu(i) + ak*r(i)
       ENDDO
	  
! 
      DEALLOCATE(r,z)

!
      RETURN
    END
    
! = = = = = = = = = = = 
      SUBROUTINE Relax_GS0P(maxit,relax,n,nintf,nnz,ia,ja,au,u,b,diagr)
      
! ---
      use omp_lib
!      USE MD_OpenMP
      
      IMPLICIT NONE
! ---
      INTEGER(4) maxit,n,nnz,nintf
      INTEGER(4) ia(n+1)
      INTEGER(4) ja(nnz)
      REAL(8) relax
      REAL(8) au(nnz)    
      REAL(8) u(n),b(n),diagr(n)
!
      INTEGER(4) i,j,k,iter
      INTEGER(4) jj,nn
      INTEGER(4) j0,j1,j2,j3,j4,j5,j6,j7,j8
      REAL(8) temp
      
! ---
!---Gauss-Seidel method
  DO iter=1,maxit
      
!$omp PARALLEL DO private(nn,jj,temp,j,j0,j1,j2,j3,j4,j5,j6,j7,j8)
      
   DO i=1,nintf
     nn=ia(i+1)-ia(i)
     jj=ia(i)
     if(nn.eq.1) then
         j0=ja(jj  )
         temp= b(i)           &
              -au(jj  )*u(j0)
     elseif(nn.eq.2) then
         j0=ja(jj  )
         j1=ja(jj+1)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1)
     elseif(nn.eq.3) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2)
     elseif(nn.eq.4) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3)
     elseif(nn.eq.5) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4)
     elseif(nn.eq.6) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5)
     elseif(nn.eq.7) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6)
     elseif(nn.eq.8) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         j7=ja(jj+7)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6) &
              -au(jj+7)*u(j7)
     elseif(nn.eq.9) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         j7=ja(jj+7)
         j8=ja(jj+8)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6) &
              -au(jj+7)*u(j7) &
              -au(jj+8)*u(j8)
     else
      temp = b(i)
      do jj=ia(i),ia(i+1)-1
         j=ja(jj)
         temp=temp-au(jj)*u(j)
      enddo
     endif
     u(i)=temp*diagr(i) + u(i)
   End do
   
  !$omp end PARALLEL DO
   
  End do
!
      return
    end
    
! = = = = = = = = = = = 
      SUBROUTINE Relax_GS0P_BW(maxit,relax,n,nintf,nnz,ia,ja,au,u,b,diagr)
      
! ---
      IMPLICIT NONE
! ---
      INTEGER(4) maxit,n,nnz,nintf
      INTEGER(4) ia(n+1)
      INTEGER(4) ja(nnz)
      REAL(8) relax
      REAL(8) au(nnz)    
      REAL(8) u(n),b(n),diagr(n)
!
      INTEGER(4) i,j,k,iter
      INTEGER(4) jj,nn
      INTEGER(4) j0,j1,j2,j3,j4,j5,j6,j7,j8
      REAL(8) temp
      
! ---
!---Gauss-Seidel method
  DO iter=1,maxit
   DO i=nintf,1,-1
     nn=ia(i+1)-ia(i)
     jj=ia(i)
     if(nn.eq.1) then
         j0=ja(jj  )
         temp= b(i)           &
              -au(jj  )*u(j0)
     elseif(nn.eq.2) then
         j0=ja(jj  )
         j1=ja(jj+1)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1)
     elseif(nn.eq.3) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2)
     elseif(nn.eq.4) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3)
     elseif(nn.eq.5) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4)
     elseif(nn.eq.6) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5)
     elseif(nn.eq.7) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6)
     elseif(nn.eq.8) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         j7=ja(jj+7)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6) &
              -au(jj+7)*u(j7)
     elseif(nn.eq.9) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         j7=ja(jj+7)
         j8=ja(jj+8)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6) &
              -au(jj+7)*u(j7) &
              -au(jj+8)*u(j8)
     else
      temp = b(i)
      do jj=ia(i),ia(i+1)-1
         j=ja(jj)
         temp=temp-au(jj)*u(j)
      enddo
     endif
     u(i)=temp*diagr(i) + u(i)
   End do
  End do
!
      return
      end
    

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      SUBROUTINE Smooth_ILU(ndom,nintf,neq,niut,nnz,ia,ja,ju,     &
                          si,ri,sintf,rintf,iut,au,alu,arhsu,solu)
!
      IMPLICIT NONE
      
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
!
      INTEGER(4) ndom
      INTEGER(4) nintf,neq,niut,nnz
      INTEGER(4) ia(neq+1),ju(neq),ja(nnz)
      INTEGER(4) si(niut+1),ri(niut+1),iut(niut),sintf(si(niut+1)-1),rintf(ri(niut+1)-1)
      REAL(8) au(nnz),alu(nnz),arhsu(neq)
! 
      REAL(8) solu(neq)
! 
!  temp.
      INTEGER(4) i,j,n,np,Allocatestatus,ierr
      REAL*8,DIMENSION(:),ALLOCATABLE::r,z
      real(8) bknum,bknum1
! 

      n=neq
      np = ndom
      Allocate(r(n),z(n),stat=Allocatestatus)
      if(AllocateStatus/=0) stop "**Not enough memory " 

! step 1: r = b-A*x
!-----Exchange solu to compute A*solu
      
      CALL send_receive(niut,neq,si,ri,sintf,rintf,iut,solu)
      
      call amux0_PCG(nintf,n,nnz,solu,r,au,ja,ia)
!-----------------------------------------------
      
        DO i=1,nintf
        r(i) = arhsu(i)-r(i)	 
!	    bknum=bknum+r(i)*r(i)
        ENDDO
        
! step 2: z = ILU(A)*r
        
!.......ILU PC
!
      
          CALL lusol0P(nintf,n,nnz,r,z,alu,ja,ia,ju)
          
! step 3: x = x+ z
        
        DO i=1,nintf
        solu(i) = solu(i) + z(i)
        ENDDO
      
! 
      DEALLOCATE(r,z)

!
      RETURN
    END
    
! = = = = = = = = = = = 
        




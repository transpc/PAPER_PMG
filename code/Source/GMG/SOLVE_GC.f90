! - - - - - - - - - 
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
  SUBROUTINE SOLVE_GC_all(nintf,nnode,maxit,err_crt,r,e)
          
      USE MD_MG_Global_C, ONLY: i_dir,nlv_glo,nnodeG,nnzG,eG,rG,rG0,imapG,iaG,jaG,juG,auG,      &
                                      igather, nsengatR, irevgatR, idispR , imapgatR
      USE MD_MPI, ONLY: myrank
      USE MD_parameter, ONLY: ndom
       
      use omp_lib
      
      IMPLICIT NONE
          
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
!
! input: 
  INTEGER nintf,nnode,maxit
  REAL*8 err_crt
  REAL*8 r(nnode)
! out:
  REAL*8 e(nnode)
! temp:
  INTEGER i,j,k,ierr
  REAL(8) tmp
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!$omp PARALLEL DO
   DO i=1,nnodeG
   rG0(i) = 0.d0
   rG(i) = 0.d0
   eG(i) = 0.d0
   ENDDO
!$omp end PARALLEL DO 


!step-1: transfer residual from local to global
       IF(igather.EQ.0) THEN
         DO i=1,nintf
            j = imapG(i)
                rG0(j) = r(i)
         ENDDO
       ENDIF
!
    IF(ndom.EQ.1) THEN
      IF(igather.EQ.0) THEN        
      rG = rG0
      ELSE 
      !$omp PARALLEL DO	  
	  DO i=1,nintf
      rG(imapgatR(i)) = r(i)
	  ENDDO
      !$omp end PARALLEL DO 
      ENDIF
    
    ELSE
!DEC$IF defined (mpi_flag)

    IF(igather.EQ.0) THEN
! step-2: S&R to all processors
        CALL  MPI_ALLREDUCE(rG0,rG,nnodeG,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
 
    ELSE
             
    CALL MPI_ALLGATHERV(r,nsengatR,mpi_double_precision,rG0,irevgatR,idispR,mpi_double_precision,mpi_comm_world, ierr)

      !$omp PARALLEL DO
    
        DO j = 1,nnodeG
            rG(imapgatR(j)) = rG0(j)
        ENDDO
            
      !$omp end PARALLEL DO
        
    ENDIF
    
!DEC$ENDIF    
    ENDIF
         
! step-3: SOLVE EXAC on coarsest level
!
      
! 
      
      IF(nlv_glo.EQ.0) THEN
          CALL SOLVE_EXACT(i_dir,nnodeG,eG,rG,nnzG,iaG,jaG,juG,auG)
      ELSE
          CALL SOLVE_COARSE(i_dir,nnodeG,eG,rG,nnzG,iaG,jaG,juG,auG)
          
      ENDIF
!
      
       
!!DEC$IF defined (mpi_flag)

!!DEC$ENDIF
! step -4: transfer to local the error
      
      !$omp PARALLEL DO private(j)
         DO i=1,nnode
            j = imapG(i)
            e(i) = eG(j)
         ENDDO
      !$omp end PARALLEL DO

!
      RETURN
    END
! - - - - - - - - - 
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! - - - - - - - - - 
    
    
! = = = = = = = = = = = = = = = = = = = = = = = =     
    
! = = = = = = = = = = = = = = = = = = = = = = = = 
      SUBROUTINE SOLVE_COARSE(i_dir,nnodeG,eG,rG,nnzG,iaG,jaG,juG,auG)
! 
      USE MD_MG_index, ONLY: iter_mg
	  USE MD_MG_Global_C,ONLY: nlv_glo,                                                   &
	                      nnziG,XintpG,XrestG,iaiG,jaiG,jarG,iarG,                        &
						  nnodeGC,nnzGC,iaGC,jaGC,juGC,auGC,                              &
						  eGC,rGC,eGt,rGt,                                                &
						  nnziGC,iaiGC,jaiGC,iarGC,jarGC,XrestGC,XintpGC,   &
                                                  nnodecm,nnzicm,nnzcm

	  
	  
      USE omp_lib
!---------------------------------------------------------------------!
      implicit none

! ---input
      INTEGER (4) i_dir
	  INTEGER (4) nnodeG,nnzG
	  INTEGER (4) iaG(nnodeG+1),jaG(nnzG),juG(nnodeG)
	  REAL(8) eG(nnodeG),rG(nnodeG),auG(nnzG)
! --- temp
      INTEGER (4) ilv,i,j,id_GS,iter_mgc
	  INTEGER (4) nnode1,nnode2,nnz,nnzi
      
      REAL(8) r(nnodecm),e(nnodecm),rt(nnodecm),et(nnodecm)
      INTEGER (4)ia(nnodecm+1),ja(nnzcm),ju(nnodecm)
      INTEGER (4) iai(nnodecm+1),jai(nnzicm),iar(nnodecm+1),jar(nnzicm)
      REAL (8) Xrest(nnzicm),Xintp(nnzicm),au(nnzcm)

! --- 
      nnode1 = nnodecm !nnodeGC(1)
	  nnz = nnzcm!nnzGC(1)
	  nnzi = nnzicm !MAXVAL(nnziGC)
!/
      iter_mgc = max(iter_mg,2)
      
! set it is allway to 2 because no communication overhead!
!/
	  
! ====================================================================!
! ---   Starting M-G Iteration method --- 
! ====================================================================!    
! --------
    !$omp PARALLEL DO
 	
	  DO i=1,nnodecm
      r(i) = 0.d0
	  e(i) = 0.d0
	  rt(i) = 0.d0
	  ENDDO
	  
    !$omp end PARALLEL DO
        
! step 1: smoothing on level AuG
! step 2: restriction on 
! IF: nlv_glo = 1 -> goto diret solver directly
! IF: nlv_glo >1 ->  some levels more and goto direct solver. 
! 
	  CALL Relax_GS(ITER_MGC,1.d0,nnodeG,nnzG,iaG,jaG,juG,auG,eG,rG)
      

      CALL resi_GC(nnodeG,eG,rG,rGt,auG,jaG,iaG)
! ----restriction 
      ilv = 1
      nnode2 = nnodeGC(ilv)
      call mt_amux2(nnode2,nnodeG,nnziG,rGt,r,XrestG,jarG,iarG)
      
! here: nnode2: fine, nnode1: coarse
	  DO ilv = 1,nlv_glo-1
	  
! smoothing
! 
      nnz = nnzGC(ilv)
    !$omp PARALLEL
	
    !$omp DO	
	  DO i=1,nnode2+1
	  ia(i) = iaGC(i,ilv)
	  ENDDO
    !$omp END DO

    !$omp DO	
	  DO i=1,nnz
      ja(i) = jaGC(i,ilv)
	  au(i) = auGC(i,ilv)
	  ENDDO
    !$omp END DO
! 
    !$omp DO
	  DO i=1,nnode2
	  ju(i) = juGC(i,ilv)
 	  e(i) = 0.d0
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL	
!
	  CALL Relax_GS(ITER_MGC,1.d0,nnode2,nnz,ia,ja,ju,au,e,r)
    
	  
! store
    !$omp PARALLEL DO
	  DO i=1,nnode2
      eGC(i,ilv) = e(i)	  
      rGC(i,ilv) = r(i)
	  ENDDO
    !$omp END PARALLEL DO		  
!
      CALL resi_GC(nnode2,e,r,rt,au,ja,ia)
! restiction
      nnode1 = nnodeGC(ilv+1)
!
	  nnzi = nnziGC(ilv+1)
!
    !$omp PARALLEL
	
    !$omp DO
	  DO i=1,nnzi
      Xrest(i) = XrestGC(i,ilv+1)
      jar(i) = jarGC(i,ilv+1)
	  ENDDO
    !$omp END DO	  
	  
    !$omp DO
	  DO i=1,nnode1+1
	  iar(i) = iarGC(i,ilv+1)
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL
!
      call mt_amux2(nnode1,nnode2,nnzi,rt,r,Xrest,jar,iar)	
! for next level: (coarser)

      nnode2 = nnode1
 
      ENDDO
	  
! = = = = = = = = = = = = = = = = = = 
!    Direct solver at nle_glo

!	  => e2 = A^-1*r2
!
     IF(i_dir.EQ.0) THEN
      ilv = nlv_glo
      nnz = nnzGC(ilv)
	  ia(1:nnode2+1) = iaGC(1:nnode2+1,ilv)
      ja(1:nnz) = jaGC(1:nnz,ilv)
	  ju(1:nnode2) = juGC(1:nnode2,ilv)
	  au(1:nnz) = auGC(1:nnz,ilv)
	 ENDIF
	 
	  CALL SOLVE_EXACT(i_dir,nnode2,e,r,nnz,ia,ja,ju,au)
	  
! = = = = = = = = = = = = = = = = == = 
	  
! notes: nnode1: fine; nnode2: coarse

	  DO ilv = nlv_glo-1,1,-1

      nnode1 = nnodeGC(ilv)
!----interpolation
	  nnzi = nnziGC(ilv+1)
      
    !$omp PARALLEL
	
    !$omp DO
	  DO i=1,nnzi
      Xintp(i) = XintpGC(i,ilv+1)
      jai(i) = jaiGC(i,ilv+1)
	  ENDDO
    !$omp END DO	  
	  
    !$omp DO
	  DO i=1,nnode1+1	  
	  iai(i) = iaiGC(i,ilv+1)
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL
!
      call mt_amux2(nnode1,nnode2,nnzi,e,et,Xintp,jai,iai)
	  
      nnz = nnzGC(ilv)
!
    !$omp PARALLEL
	
    !$omp DO
	  DO i=1,nnode1	  
      e(i) = eGC(i,ilv) + et(i)
	  r(i) = rGC(i,ilv)
	  ju(i) = juGC(i,ilv)
	  ENDDO
    !$omp END DO

    !$omp DO
	  DO i=1,nnode1+1	
	  ia(i) = iaGC(i,ilv)
	  ENDDO
    !$omp END DO
	
    !$omp DO
	  DO i=1,nnz	  
      ja(i) = jaGC(i,ilv)
	  au(i) = auGC(i,ilv)
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL
!
	  CALL Relax_GS(ITER_MGC,1.d0,nnode1,nnz,ia,ja,ju,au,e,r)
     
	  
! for next level: (finer)
      nnode2 = nnode1
	  
	  ENDDO
	  
! level AuG
      CALL mt_amux2(nnodeG,nnode2,nnziG,e,eGt,XintpG,jaiG,iaiG)	  
	  
    !$omp PARALLEL DO
	  DO i=1,nnodeG
	  eG(i) = eG(i) +  eGt(i)
	  ENDDO
    !$omp END PARALLEL DO
!	  
	  CALL Relax_GS(ITER_MGC,1.d0,nnodeG,nnzG,iaG,jaG,juG,auG,eG,rG)  
     
! 
! - - - - - 
	  
      RETURN
      END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = 
subroutine resi_GC(n,x,b,r,a,ja,ia)

      use omp_lib
      
implicit none
!
integer :: n,i,k1,k2,j
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8  :: x(*)
real*8  :: r(*)
real*8  :: b(*)
real*8 temp
!  ...

      !$omp PARALLEL DO private(k1,k2,temp,j)
do i= 1,n
   k1 = ia(i)
   k2 = ia(i+1)-1
   temp = b(i)
   do j=k1,k2
    temp = temp -a(j)*x(ja(j))   
   enddo
   
   r(i) = temp  !b(i)-temp
enddo

      !$omp end PARALLEL DO
!=====
return
      
    end subroutine

! = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = 

! = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    SUBROUTINE SOLVE_EXACT(i_dir,nnode,e,r,nnz,ia,ja,ju,au)
!
    USE MD_MG_Global_C, ONLY: Ainv, aluG
    USE MD_MG_index, ONLY: crit_1,maxit_1
    USE MD_parameter, ONLY: ndom
    USE MD_MPI, ONLY: myrank
!
    use omp_lib
          
    IMPLICIT NONE 
          
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
! 
       integer::tag,ierr
!DEC$IF defined (mpi_flag)
       integer::status(mpi_status_size)
!DEC$ENDIF
       
    INTEGER (4) i_dir
	INTEGER (4) nnode,nnz
	INTEGER (4) ia(*),ja(*),ju(*)
    REAL(8) e(*),r(*), au(*)
    REAL(8) e0(nnode)
! 
    INTEGER(4) i,j,i1,i2
	REAL(8) tmp
!
	
    IF(i_dir.EQ.1) THEN
        
      !$omp PARALLEL DO private(tmp,j)
     DO i=1,nnode
         e(i) =  DOT_PRODUCT(Ainv(i,1:nnode),r(1:nnode))
     ENDDO
     
      !$omp end PARALLEL DO
     
    ELSEIF(i_dir.EQ.2) THEN
        
        j = INT(nnode/ndom-0.5)+1
        IF(j.LT.1) j = 1
        IF(j.GT.nnode) j = nnode
        e0 = 0.d0
        
        i1 = myrank*j+1
        i2 = (myrank+1)*j
        i2 = MIN(i2,nnode)
        
        
     DO i=i1,i2
         e0(i) = DOT_PRODUCT(Ainv(i,1:nnode),r(1:nnode))
     ENDDO    
!
!DEC$IF defined (mpi_flag)
        CALL  MPI_ALLREDUCE(e0,e,nnode,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
!DEC$ENDIF     
!
	ENDIF
	
	RETURN
	
    END

! = = = = = = = = = = = = = = = = = = = = = = = = 
! = = = = = = = = = = = = = = = = = = = = = = = = = = = !

! = = = = = = = = = = = = = = = = = = = = = = = = 

! = = = = = = = = = = = = = = = = = = = = = = = = = = = 
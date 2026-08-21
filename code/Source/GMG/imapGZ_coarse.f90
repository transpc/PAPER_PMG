
      
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Coarsest global:

      SUBROUTINE imapGZ_coarse(nintf,nnode,nnz,ia,ja)
      
       USE MD_MG_Global_C, ONLY: iaG, jaG, imapG, imapGZ
       IMPLICIT NONE
           
! inlet: 
       INTEGER(4) nintf,nnode,nnz
       INTEGER(4) ia(nnode+1),ja(nnz)
           
! 
       INTEGER(4) i,j,k,i1,i2,j1,j2,k1,m 
           
       ALLOCATE(imapGZ(nnz))
       imapGZ = 0
          
      DO i=1,nintf
            i1 = ia(i)
                i2 = ia(i+1)
                
                j = imapG(i)
                j1 = iaG(j)
                j2 = iaG(j+1)
                
                DO k=i1,i2-1
                   k1 = imapG(ja(k))
                   Do m = j1,j2-1
                      IF(k1.EQ.jaG(m)) GOTO 10
                   ENDDO
                   
                   WRITE(*,*)'PMG error in imapGZ'
                   
10         CONTINUE
           imapGZ(k) = m
          ENDDO
          ENDDO
!          
          RETURN
    END

    ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! for gather:

      SUBROUTINE imap_GATHER(np,myrank,nintf,nnz,ia)
      
       USE MD_MG_Global_C, ONLY:  nnodeG, nnzG, imapGZ, imapG,                      &
                          nsengatA, irevgatA, idispA , imapgatA,       &
                          nsengatR, irevgatR, idispR , imapgatR      
       IMPLICIT NONE
           
!DEC$IF defined (mpi_flag)
   INCLUDE 'mpif.h'
!DEC$ENDIF
! inlet: 
       INTEGER(4) np,myrank,nintf,nnz
       INTEGER(4) ia(*)
           
! 
       INTEGER(4) i,j,k,i1,i2,j1,j2,k1,m,nd ,ip
       integer::tag,ierr
!DEC$IF defined (mpi_flag)          
       integer::status(mpi_status_size)
!DEC$ENDIF
       INTEGER(4),DIMENSION(:),ALLOCATABLE:: imark
           
!
       
      ALLOCATE(irevgatA(np), idispA(np),imapgatA(nnzG))
      ALLOCATE(irevgatR(np), idispR(np),imapgatR(nnodeG))
      
      idispA = 0
      imapgatA = 0


      idispR = 0
      imapgatR = 0
      
      nsengatA = ia(nintf+1)-1
      nsengatR = nintf   
      
! initial and used for np = 1
      
      irevgatA = nsengatA
      irevgatR = nsengatR  
      
      IF(np.NE.1) THEN
!DEC$IF defined (mpi_flag)      
      CALL mpi_barrier(mpi_comm_world,ierr) 

      CALL MPI_ALLGATHER(nsengatA,1,MPI_INTEGER,irevgatA,1,MPI_INTEGER,mpi_comm_world,ierr)

      CALL MPI_ALLGATHER(nsengatR,1,MPI_INTEGER,irevgatR,1,MPI_INTEGER,mpi_comm_world,ierr)      

!DEC$ENDIF
      ENDIF
      
!

      DO i=2,np
         idispA(i) = idispA(i-1) +  irevgatA(i-1)
         idispR(i) = idispR(i-1) +  irevgatR(i-1)    
      ENDDO
      
! set imap:

        
        
!
        
! 
      IF(np.EQ.1) THEN
      imapgatA(1:nnzG) =  imapGZ(1:nnzG)
      imapgatR(1:nnodeG) = imapG(1:nnodeG)
      ELSE

!DEC$IF defined (mpi_flag)
      CALL mpi_barrier(mpi_comm_world,ierr) 
      CALL MPI_ALLGATHERV(imapGZ,nsengatA,MPI_INTEGER,imapgatA,irevgatA,idispA,MPI_INTEGER,mpi_comm_world, ierr)
        
      CALL MPI_ALLGATHERV(imapG,nsengatR,MPI_INTEGER,imapgatR,irevgatR,idispR,MPI_INTEGER,mpi_comm_world, ierr)
  
!DEC$ENDIF
      ENDIF
      
! test validity 
        ALLOCATE(imark(nnzG))
        imark = 0
!                
          
            DO j = 1,nnzG
        
            k = imapgatA(j)
            IF(imark(k).EQ.1) THEN
                write(*,*)'error in imap_GAT'
            ENDIF
            imark(k)=1
            ENDDO
        
        DEALLOCATE(imark)
        
        !          
          RETURN
    END

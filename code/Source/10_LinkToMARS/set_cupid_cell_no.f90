!   
      SUBROUTINE set_cupid_cell_no
!                                                                      
!     Make an inverse table of CUPID local cell index
!     by saving CUPID local cell index of CUPVOLs
!
      USE Zmpi       , ONLY: jperm
      USE Zzone      , only: ncell_fluid
      USE Zcore      , only: myrank, np
      USE Znum_cell  , only: i_neigh
      USE Zbc_index  , only: nbcon
      USE Zparam     , ONLY: nb_mars
      USE Zmars      , only: ppcup,pcup,egcup,elcup,alphagcup,cboroncup,qualacup,rholcup,rhogcup,&
                             ppcup_tmp,pcup_tmp,egcup_tmp,elcup_tmp,alphagcup_tmp,cboroncup_tmp,&
                             qualacup_tmp,rholcup_tmp,marsindex,i3invtbl_tmp,n_marsbc,rhogcup_tmp,&
                             i3cupid_loc,&
                             tlcup,tgcup,tlcup_tmp,tgcup_tmp
      USE c3com_cupid, only: i3cupid,i3invtbl,nvols_mars
!
      IMPLICIT NONE
!      
!DEC$IF defined (MCC)      
      INCLUDE 'c3com.h'
!DEC$ELSEIF defined (MCC_DLL)          
      INCLUDE 'c3com.h'
      !dec$ attributes dllexport :: set_cupid_cell_no
!DEC$ELSEIF defined (SPACE)          
      INCLUDE 'c3com.h'
      !dec$ attributes dllexport :: set_cupid_cell_no
!DEC$ENDIF 
!
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
      INTEGER ierr
      INTEGER, ALLOCATABLE::rcounts_tmp(:),displs_tmp(:),rcounts(:),displs(:)       
!DEC$ENDIF      
      INTEGER, ALLOCATABLE::i3invtbl2(:),idxcup(:),idxcup_tmp(:)
!
      INTEGER i,j,k,idx_tmp,idx_,loop,idx_idxcup
!      
      logical,save:: initial_p
!      
      DATA initial_p /.true./
!
      ALLOCATE(marsindex(ncell_fluid))
      ALLOCATE(i3invtbl_tmp(n_marsbc))
      ALLOCATE(ppcup(n_marsbc),pcup(n_marsbc),egcup(n_marsbc),elcup(n_marsbc),&
                alphagcup(n_marsbc),cboroncup(n_marsbc),qualacup(n_marsbc),&
                rholcup(n_marsbc),rhogcup(n_marsbc))
      ALLOCATE(ppcup_tmp(n_marsbc),pcup_tmp(n_marsbc),egcup_tmp(n_marsbc),elcup_tmp(n_marsbc),&
                alphagcup_tmp(n_marsbc),cboroncup_tmp(n_marsbc),qualacup_tmp(n_marsbc),&
                rholcup_tmp(n_marsbc),rhogcup_tmp(n_marsbc))
      ALLOCATE(tlcup(n_marsbc),tgcup(n_marsbc),tlcup_tmp(n_marsbc),tgcup_tmp(n_marsbc))                
      ALLOCATE(idxcup(n_marsbc),idxcup_tmp(n_marsbc))                
!      
      marsindex(:)=0

      i3invtbl(:)=0
      ppcup(:)=0.0d0
      pcup(:)=0.0d0
      egcup(:)=0.0d0
      elcup(:)=0.0d0
      alphagcup(:)=0.0d0
      cboroncup(:)=0.0d0
      qualacup(:)=0.0d0
      rholcup(:)=0.0d0
      rhogcup(:)=0.0d0
     
      i3invtbl_tmp(:)=0
      ppcup_tmp(:)=0.0d0
      pcup_tmp(:)=0.0d0
      egcup_tmp(:)=0.0d0
      elcup_tmp(:)=0.0d0
      alphagcup_tmp(:)=0.0d0
      cboroncup_tmp(:)=0.0d0
      qualacup_tmp(:)=0.0d0
      rholcup_tmp(:)=0.0d0
      rhogcup_tmp(:)=0.0d0      
!DEC$IF defined (MCC)      
      i3cupid(:)=i3cell(1,:)
!DEC$ELSEIF defined (MCC_DLL)          
      i3cupid(:)=i3cell(1,:)
!DEC$ELSEIF defined (SPACE)          
      i3cupid(:)=i3cell(1,:)
!DEC$ENDIF 
!
!...Find the inverse table of the CUPID cell index, next
!
      idx_=0
      i3invtbl(:)=0
      ALLOCATE(i3invtbl2(n_marsbc))
      ALLOCATE(i3cupid_loc(n_marsbc))
      i3cupid_loc(:)=0
      i3invtbl2(:)=0
      DO i=1,ncell_fluid
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ge.nb_mars)THEN
               idx_=idx_+1
!!!               idx_idxcup=(nbcon(j)-200)/10)
!!!               WRITE(*,*)'idx_idxcup=',idx_idxcup,nbcon(j)
!!!               idxcup(idx_idxcup)=jperm(i)
               loop=0
               DO k=1,i3nic(2)
                  IF(jperm(i).eq.i3cupid(k))THEN 
                     i3cupid_loc(k)=i
                     i3invtbl(idx_)=k !connect cupid index of cupvols(idx_) to mars index of cupvols (k)
                     i3invtbl2(k)=k   !to make i3invtbl_tmp
                     marsindex(i)=k   !connect local cell number of cupvols(i)  to mars index of cupvol (k)
!                     WRITE(*,"(11x,a,1i3,1i3,1i10,1i3)")'## myrank,idx_,jperm(i),k=',myrank,idx_,jperm(i),k 
                     loop=loop+1
                     EXIT
                  ENDIF
               ENDDO
               IF(loop.ne.1)THEN
                  WRITE(*,"(11x,a,i)")'## error! in set_cupid_cell_no, loop=',loop
                  PAUSE
                  STOP
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
      CALL allreduce_i1(idx_,idx_tmp)
      IF(1)then
         CALL allreduce_i(i3invtbl2,i3invtbl_tmp,n_marsbc)
      ENDIF
      IF(np.gt.1)then
          IF(1)then
!            WRITE(*,"(11x,a,1i3,1x,8i3)")'## myrank,i3invtbl2(1:8)=',myrank,i3invtbl2(1:n_marsbc) 
!            IF(myrank.eq.0)WRITE(*,"(11x,a,1i3,1x,8i3)")'## myrank,i3invtbl2(1:8)=',myrank,i3invtbl_tmp(1:n_marsbc)
            DEALLOCATE(i3invtbl2)
          ELSE
!DEC$IF defined (mpi_flag)
            CALL barrier_mpi
!            WRITE(*,"(11x,a,1i3,1x,8i3)")'## myrank,i3invtbl(1:8)=',myrank,i3invtbl(1:n_marsbc)
            allocate(rcounts(0:np-1),rcounts_tmp(0:np-1),displs(0:np-1),displs_tmp(0:np-1))
!!!!!!!WHY MAKE IT COMPLICATED ???????
!           rcounts=0
!           rcounts_tmp=0
!           displs=0
!           displs_tmp=0
!           rcounts(myrank)=idx_
!           CALL barrier_mpi
!           CALL allreduce_i(rcounts,rcounts_tmp,np)
!           DO i=0,np-1
!              IF(i.lt.myrank)displs(myrank)=displs(myrank)+rcounts_tmp(i)
!           ENDDO
!           CALL allreduce_i(displs,displs_tmp,np)        
!           rcounts(:)=rcounts_tmp(:)
!           displs(:)=displs_tmp(:)
!           CALL MPI_GatherV(i3invtbl,idx_,MPI_INTEGER,i3invtbl_tmp,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!           !CALL MPI_Gather(i3invtbl,2,MPI_INTEGER,i3invtbl_tmp,2,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!           CALL MPI_BCAST(i3invtbl_tmp,n_marsbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
!           CALL barrier_mpi
!!!
            CALL allgather_i(idx_,rcounts)
            i=0
            displs(i)=0
            DO i=1,np-1
               displs(i)=displs(i-1)+rcounts(i-1)
            ENDDO
            CALL allgather_vec_i(i3invtbl,idx_,i3invtbl_tmp,n_marsbc,rcounts,displs)
!            WRITE(*,"(11x,a,1i3,1x,8i3)")'## myrank,i3invtbl_tmp(1:8)=',myrank,i3invtbl_tmp(1:n_marsbc)
!DEC$ENDIF            
!!!            CALL allreduce_i(idxcup,idxcup_tmp,n_marsbc)           
!!!            DO i=1,n_marsbc
!!!               write(97,"(11x,a,1i3,1x,1i10)")'##nbcon,jperm(i),k=',i*10+200,idxcup_tmp(i)
!!!            ENDDO   
         ENDIF
      ELSE
         i3invtbl_tmp(1:n_marsbc)=i3invtbl(1:n_marsbc)
      ENDIF
!          write(*,"(11x,a,3i3)")'## myrank,n_marsbc,idx_tmp=',myrank,n_marsbc,idx_tmp
!          write(*,"(11x,a,2i3)")'## myrank,n_marsbc',myrank,n_marsbc  
!          write(*,"(11x,a,8i7)")'## i3cupid(j)=',i3cupid(1:n_marsbc)      
!
!...Check and print the mars boundary into 'event.dat'
!
      IF(myrank.eq.0) then
         IF(n_marsbc.ne.idx_tmp)then
            write(97,"(11x,a)")'## MARS-CUPID boundary number is not coherent!'
            write(97,"(11x,a,2i3)")'## n_marsbc,idx_tmp in CUPID=',n_marsbc,idx_tmp
            write(*,"(11x,a)") '## MARS-CUPID boundary number is not coherent!'
            write(*,"(11x,a,2i3)") '## n_marsbc,idx_tmp in CUPID=',n_marsbc,idx_tmp
            pause
            stop 
         ELSE
            write(97,"(11x,a,1i3)")'## MARS-CUPID boundary number is coherent as',n_marsbc        
            write(*,"(11x,a,1i3)")'## MARS-CUPID boundary number is coherent as',n_marsbc        
            DO i=1, idx_tmp
               IF(myrank.eq.0)write(97,1000)i3invtbl_tmp(i),i3cupid(i3invtbl_tmp(i))
               write(*,1000)i3invtbl_tmp(i),i3cupid(i3invtbl_tmp(i))
            ENDDO            
         ENDIF
      ENDIF
!      
1000  FORMAT(1x,2I10)    
!
      END SUBROUTINE set_cupid_cell_no    

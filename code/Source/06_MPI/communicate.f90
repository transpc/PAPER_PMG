!
!     Friendly user inteface for MPI communication
!     It uses variable argument list
!
!     ==>Fluid processing
!     communicate_1d          communicate arrays 1d of the form a(ncell_fp) max 8 arguments
!     communicate_2d          communicate arrays 2d of the form a(ncell_fp,ndim)  max 3 arguments
!     communicate_3d          communicate arrays 3d of the form a(ncell_fp,ndim,ndim) max 3 arguments
!
!     communicate_1d_int      communicate arrays 1d of the form a(ncell_fp) integer form max 3 arguments
!
!     communicate_1d_csr      communicate 1 array  1d of the form CSR
!     communicate_2d_csr      communicate 1 array  2d of the form CSR in first dimension
!
!     communicate_1d_csr_int  communicate 1 array  1d of the form CSR integer form
!
!     ==>Solid processing
!     communicate_1d_c        communicate 1 array 1d of the form a(ncell_ps)
!     communicate_2d_c        communicate 1 array 2d of the form a(ncell_ps,ndim)
!
!     communicate_1d_c_int    communicate 1 array 1d of the form a(ncell_ps) integer form
!
!     communicate_1d_c_csr    communicate 1 array  1d of the form CSR
!
!     ==>RV processing
!     communicate_rv_2d       communicate  arrays 2d of the form  a(ncell_fuel_rod_p,nr_2d) max 2 arguments
!
!     ==>Initialization processing
!     communicate_allb        communicate  specific arrays (→ communicate_allb.f90 로 분리, C010-1)
!
!     ==>Solver processing
!     communicate             switching driver to enable solid or fluid communication via izone flag
!
!     ==>Internal routines handling send-receive-wait: THOSE ROUTINES SHOULD NEVER BE CALLED by THE USER
!     communicate_nb          communicate nb>1 blocks
!     communicate_nb_int      communicate nb>1 blocks integer form
!     mwait                   end-receive wait
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d(a1,a2,a3,a4,a5,a6,a7,a8)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi     , ONLY: ncell_fp,niut,iut,si,ri,sintf,rintf
      USE Zzone    , ONLY: ncell_fluid
!
      IMPLICIT NONE
!.....Input
      REAL(8),INTENT(INOUT) :: a1(ncell_fp)
      REAL(8),DIMENSION(ncell_fp),OPTIONAL :: a2,a3,a4,a5,a6,a7,a8
!.....Local variables
      INTEGER :: i,j,i1,ib
      INTEGER :: ip,ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8
      INTEGER :: nb,lda
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: svar,rvar
!
      nb=1
      IF(.not.PRESENT(a2)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a3)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a4)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a5)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a6)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a7)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a8)) GOTO 100
      nb=nb+1
100   CONTINUE
!
      ALLOCATE(svar((si(niut+1)-1)*nb),rvar((ri(niut+1)-1)*nb))
      IF(nb.eq.1) THEN
         DO i=1,si(niut+1)-1
            j=sintf(i)
            svar(i)=a1(j)
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+  lda
            ip3=ip+2*lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
               svar(i1+ip3)=a3(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.4) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+  lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
               svar(i1+ip3)=a3(j)
               svar(i1+ip4)=a4(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.5) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+  lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
               svar(i1+ip3)=a3(j)
               svar(i1+ip4)=a4(j)
               svar(i1+ip5)=a5(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.6) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+  lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            ip6=ip+5*lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
               svar(i1+ip3)=a3(j)
               svar(i1+ip4)=a4(j)
               svar(i1+ip5)=a5(j)
               svar(i1+ip6)=a6(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.7) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+  lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            ip6=ip+5*lda
            ip7=ip+6*lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
               svar(i1+ip3)=a3(j)
               svar(i1+ip4)=a4(j)
               svar(i1+ip5)=a5(j)
               svar(i1+ip6)=a6(j)
               svar(i1+ip7)=a7(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.8) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+  lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            ip6=ip+5*lda
            ip7=ip+6*lda
            ip8=ip+7*lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
               svar(i1+ip3)=a3(j)
               svar(i1+ip4)=a4(j)
               svar(i1+ip5)=a5(j)
               svar(i1+ip6)=a6(j)
               svar(i1+ip7)=a7(j)
               svar(i1+ip8)=a8(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ENDIF
!
      CALL communicate_nb(niut,iut,si,ri,svar,rvar,nb)
!
      IF(nb.eq.1) THEN
         DO i=1,ncell_fp-ncell_fluid
            j=rintf(i)
            a1(j)=rvar(i)
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            ip3=ip+2*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
               a3(j)=rvar(i1+ip3)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.4) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
               a3(j)=rvar(i1+ip3)
               a4(j)=rvar(i1+ip4)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.5) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
               a3(j)=rvar(i1+ip3)
               a4(j)=rvar(i1+ip4)
               a5(j)=rvar(i1+ip5)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.6) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            ip6=ip+5*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
               a3(j)=rvar(i1+ip3)
               a4(j)=rvar(i1+ip4)
               a5(j)=rvar(i1+ip5)
               a6(j)=rvar(i1+ip6)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.7) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            ip6=ip+5*lda
            ip7=ip+6*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
               a3(j)=rvar(i1+ip3)
               a4(j)=rvar(i1+ip4)
               a5(j)=rvar(i1+ip5)
               a6(j)=rvar(i1+ip6)
               a7(j)=rvar(i1+ip7)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.8) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            ip3=ip+2*lda
            ip4=ip+3*lda
            ip5=ip+4*lda
            ip6=ip+5*lda
            ip7=ip+6*lda
            ip8=ip+7*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
               a3(j)=rvar(i1+ip3)
               a4(j)=rvar(i1+ip4)
               a5(j)=rvar(i1+ip5)
               a6(j)=rvar(i1+ip6)
               a7(j)=rvar(i1+ip7)
               a8(j)=rvar(i1+ip8)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ENDIF
      DEALLOCATE(svar,rvar)
!
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_2d(a1,a2,a3)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_fp,niut,iut,si,ri,sintf,rintf
      USE Zparam       , ONLY: ndim
!
      IMPLICIT NONE
!.....Input
      REAL(8),INTENT(INOUT) :: a1(ncell_fp,ndim)
      REAL(8),DIMENSION(ncell_fp,ndim),OPTIONAL :: a2,a3
!.....Local variables
      INTEGER :: i,j,i1,ib,ix
      INTEGER :: ip,ip1,ip2,ip3
      INTEGER :: nb,nb1,lda
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: svar,rvar
!       
      nb=1
      IF(.not.PRESENT(a2)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a3)) GOTO 100
      nb=nb+1
100   CONTINUE
!
      nb1=nb*ndim
      ALLOCATE(svar((si(niut+1)-1)*nb1),rvar((ri(niut+1)-1)*nb1))
!            
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            DO ix=1,ndim
               ip1=ip+(ix-1)*lda
               DO i=si(ib),si(ib+1)-1
                  i1=i-si(ib)
                  j=sintf(i)
                  svar(i1+ip1)=a1(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            DO ix=1,ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim
               DO i=si(ib),si(ib+1)-1
                  i1=i-si(ib)
                  j=sintf(i)
                  svar(i1+ip1)=a1(j,ix)
                  svar(i1+ip2)=a2(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            DO ix=1,ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim
               ip3=ip2+lda*ndim
               DO i=si(ib),si(ib+1)-1
                  i1=i-si(ib)
                  j=sintf(i)
                  svar(i1+ip1)=a1(j,ix)
                  svar(i1+ip2)=a2(j,ix)
                  svar(i1+ip3)=a3(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ENDIF
      CALL communicate_nb(niut,iut,si,ri,svar,rvar,nb1)
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            DO ix=1,ndim
               ip1=ip+(ix-1)*lda
               DO i=ri(ib),ri(ib+1)-1
                  i1=i-ri(ib)
                  j=rintf(i)
                  a1(j,ix)=rvar(i1+ip1)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            DO ix=1,ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim
               DO i=ri(ib),ri(ib+1)-1
                  i1=i-ri(ib)
                  j=rintf(i)
                  a1(j,ix)=rvar(i1+ip1)
                  a2(j,ix)=rvar(i1+ip2)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            DO ix=1,ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim
               ip3=ip2+lda*ndim
               DO i=ri(ib),ri(ib+1)-1
                  i1=i-ri(ib)
                  j=rintf(i)
                  a1(j,ix)=rvar(i1+ip1)
                  a2(j,ix)=rvar(i1+ip2)
                  a3(j,ix)=rvar(i1+ip3)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ENDIF
      DEALLOCATE(svar,rvar)
!DEC$ENDIF
!
      END SUBROUTINE communicate_2d
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_3d(a1,a2,a3)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_fp,niut,iut,si,ri,sintf,rintf
      USE Zparam       , ONLY: ndim
!
      IMPLICIT NONE
!.....Input
      REAL(8),INTENT(INOUT) :: a1(ncell_fp,ndim*ndim)
      REAL(8),DIMENSION(ncell_fp,ndim*ndim),OPTIONAL :: a2,a3
!.....Local variables
      INTEGER :: i,j,i1,ib,ix
      INTEGER :: ip,ip1,ip2,ip3
      INTEGER :: nb,nb1,lda
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: svar,rvar
!       
      nb=1
      IF(.not.PRESENT(a2)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a3)) GOTO 100
      nb=nb+1
100   CONTINUE
!
      nb1=nb*ndim*ndim
      ALLOCATE(svar((si(niut+1)-1)*nb1),rvar((ri(niut+1)-1)*nb1))
!
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            DO ix=1,ndim*ndim
               ip1=ip+(ix-1)*lda
               DO i=si(ib),si(ib+1)-1
                  i1=i-si(ib)
                  j=sintf(i)
                  svar(i1+ip1)=a1(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            DO ix=1,ndim*ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim*ndim
               DO i=si(ib),si(ib+1)-1
                  i1=i-si(ib)
                  j=sintf(i)
                  svar(i1+ip1)=a1(j,ix)
                  svar(i1+ip2)=a2(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            DO ix=1,ndim*ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim*ndim
               ip3=ip2+lda*ndim*ndim
               DO i=si(ib),si(ib+1)-1
                  i1=i-si(ib)
                  j=sintf(i)
                  svar(i1+ip1)=a1(j,ix)
                  svar(i1+ip2)=a2(j,ix)
                  svar(i1+ip3)=a3(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ENDIF
      CALL communicate_nb(niut,iut,si,ri,svar,rvar,nb1)
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            DO ix=1,ndim*ndim
               ip1=ip+(ix-1)*lda
               DO i=ri(ib),ri(ib+1)-1
                  i1=i-ri(ib)
                  j=rintf(i)
                  a1(j,ix)=rvar(i1+ip1)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            DO ix=1,ndim*ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim*ndim
               DO i=ri(ib),ri(ib+1)-1
                  i1=i-ri(ib)
                  j=rintf(i)
                  a1(j,ix)=rvar(i1+ip1)
                  a2(j,ix)=rvar(i1+ip2)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            DO ix=1,ndim*ndim
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*ndim*ndim
               ip3=ip2+lda*ndim*ndim
               DO i=ri(ib),ri(ib+1)-1
                  i1=i-ri(ib)
                  j=rintf(i)
                  a1(j,ix)=rvar(i1+ip1)
                  a2(j,ix)=rvar(i1+ip2)
                  a3(j,ix)=rvar(i1+ip3)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ENDIF
      DEALLOCATE(svar,rvar)
!DEC$ENDIF
!      
    END SUBROUTINE communicate_3d
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d_int0(b1,a1)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_fp,niut,iut,si,ri,rintf
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: b1(*)
      INTEGER :: a1(ncell_fp)
!.....Local variables
      INTEGER :: i,j,i1,ib
      INTEGER :: ip,ip1
      INTEGER :: nb,lda
!.....Local arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: svar,rvar
      ALLOCATE(svar((si(niut+1)-1)),rvar((ri(niut+1)-1)))
!
      nb=1
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               svar(i1+ip1)=b1(i1+ip1)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      CALL communicate_nb_int(niut,iut,si,ri,svar,rvar,nb)
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      DEALLOCATE(svar,rvar)
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d_int0      
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d_int(a1,a2,a3)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_fp,niut,iut,si,ri,sintf,rintf
!
      IMPLICIT NONE
!
!.....Input
      INTEGER,INTENT(INOUT) :: a1(ncell_fp)
      INTEGER,DIMENSION(ncell_fp),OPTIONAL :: a2,a3
!.....Local variables
      INTEGER :: i,j,i1,ib
      INTEGER :: ip,ip1,ip2,ip3
      INTEGER :: nb,lda
!.....Local arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: svar,rvar
      nb=1
      IF(.not.PRESENT(a2)) GOTO 100
      nb=nb+1
      IF(.not.PRESENT(a3)) GOTO 100
      nb=nb+1
100   CONTINUE
!
      ALLOCATE(svar((si(niut+1)-1)*nb),rvar((ri(niut+1)-1)*nb))
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=si(ib+1)-si(ib)
            ip1=ip
            ip2=ip+  lda
            ip3=ip+2*lda
            DO i=si(ib),si(ib+1)-1
               i1=i-si(ib)
               j=sintf(i)
               svar(i1+ip1)=a1(j)
               svar(i1+ip2)=a2(j)
               svar(i1+ip3)=a3(j)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ENDIF
      CALL communicate_nb_int(niut,iut,si,ri,svar,rvar,nb)
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ELSEIF(nb.eq.3) THEN
         ip=1
         DO ib=1,niut
            lda=ri(ib+1)-ri(ib)
            ip1=ip
            ip2=ip+lda
            ip3=ip+2*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=i-ri(ib)
               j=rintf(i)
               a1(j)=rvar(i1+ip1)
               a2(j)=rvar(i1+ip2)
               a3(j)=rvar(i1+ip3)
            ENDDO
            ip=ip+lda*nb
         ENDDO
      ENDIF
      DEALLOCATE(svar,rvar)
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d_int
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d_csr(a,i_neigh)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_fp,maxmt_fp,niut,iut,si,ri,sintf,rintf
      USE Zparam       , ONLY: ns
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: i_neigh(ncell_fp+1)
      REAL(8) :: a(maxmt_fp)
!.....Local variables
      INTEGER :: i,j,j0,nb
      INTEGER :: ib,lda,ip,i1,j1,ip1
!.....Local arrays
      REAL(8) :: svar(ns*(si(niut+1)-1)),rvar(ns*(ri(niut+1)-1))
!
      nb=ns
!
      ip=0
      DO ib=1,niut
         lda=si(ib+1)-si(ib)
         ip1=ip
         DO i=si(ib),si(ib+1)-1
            i1=ns*(i-si(ib))
            j1=sintf(i)
            j0=i_neigh(j1)-1
            DO j=i_neigh(j1),i_neigh(j1+1)-1
               svar(j-j0+i1+ip1)=a(j)
            ENDDO
            DO j=i_neigh(j1+1),i_neigh(j1)+ns-1
               svar(j-j0+i1+ip1)=0.d0
            ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
      CALL communicate_nb(niut,iut,si,ri,svar,rvar,nb)
      ip=0
      DO ib=1,niut
         lda=ri(ib+1)-ri(ib)
         ip1=ip
         DO i=ri(ib),ri(ib+1)-1
            i1=ns*(i-ri(ib))
            j1=rintf(i)
            j0=i_neigh(j1)-1
            DO j=i_neigh(j1),i_neigh(j1+1)-1
               a(j)=rvar(j-j0+i1+ip1)
            ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d_csr
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_2d_csr(a,i_neigh)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_fp,maxmt_fp,niut,iut,si,ri,sintf,rintf
      USE Zparam       , ONLY: ndim,ns
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: i_neigh(ncell_fp+1)
      REAL(8) :: a(maxmt_fp,ndim)
!.....Local variables
      INTEGER :: i,j,j0,ix,nb,nb1
      INTEGER :: ib,lda,ip,i1,j1,ip1
!.....Local arrays
      REAL(8) :: svar(ns*(si(niut+1)-1)*ndim),rvar(ns*(ri(niut+1)-1)*ndim)
!
      nb=ns*ndim
!
      nb1=ndim
      ip=0
      DO ib=1,niut
         lda=si(ib+1)-si(ib)
         DO ix=1,nb1
         ip1=ip+ns*(ix-1)*lda
         DO i=si(ib),si(ib+1)-1
            i1=ns*(i-si(ib))
            j1=sintf(i)
            j0=i_neigh(j1)-1
            DO j=i_neigh(j1),i_neigh(j1+1)-1
               svar(j-j0+i1+ip1)=a(j,ix)
            ENDDO
            DO j=i_neigh(j1+1),i_neigh(j1)+ns-1
               svar(j-j0+i1+ip1)=0.d0
            ENDDO
         ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
      CALL communicate_nb(niut,iut,si,ri,svar,rvar,nb)
      ip=0
      DO ib=1,niut
         lda=ri(ib+1)-ri(ib)
         DO ix=1,nb1
            ip1=ip+ns*(ix-1)*lda
            DO i=ri(ib),ri(ib+1)-1
               i1=ns*(i-ri(ib))
               j1=rintf(i)
               j0=i_neigh(j1)-1
               DO j=i_neigh(j1),i_neigh(j1+1)-1
                  a(j,ix)=rvar(j-j0+i1+ip1)
               ENDDO
            ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_2d_csr
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d_csr_int(a,i_neigh)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_fp,maxmt_fp,niut,iut,si,ri,sintf,rintf
      USE Zparam       , ONLY: ns
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: i_neigh(ncell_fp+1)
      INTEGER :: a(maxmt_fp)
!.....Local variables
      INTEGER :: i,j,j0,nb
      INTEGER :: ib,lda,ip,i1,j1,ip1
!.....Local arrays
      INTEGER :: svar(ns*(si(niut+1)-1)),rvar(ns*(ri(niut+1)-1))
!
      nb=ns
!
      ip=0
      DO ib=1,niut
         lda=si(ib+1)-si(ib)
         ip1=ip
         DO i=si(ib),si(ib+1)-1
            i1=ns*(i-si(ib))
            j1=sintf(i)
            j0=i_neigh(j1)-1
            DO j=i_neigh(j1),i_neigh(j1+1)-1
               svar(j-j0+i1+ip1)=a(j)
            ENDDO
            DO j=i_neigh(j1+1),i_neigh(j1)+ns-1
               svar(j-j0+i1+ip1)=0.d0
            ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
      CALL communicate_nb_int(niut,iut,si,ri,svar,rvar,nb)
      ip=0
      DO ib=1,niut
         lda=ri(ib+1)-ri(ib)
         ip1=ip
         DO i=ri(ib),ri(ib+1)-1
            i1=ns*(i-ri(ib))
            j1=rintf(i)
            j0=i_neigh(j1)-1
            DO j=i_neigh(j1),i_neigh(j1+1)-1
               a(j)=rvar(j-j0+i1+ip1)
            ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d_csr_int
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d_c(a)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_ps,niut_c,iut_c,si_c,ri_c,sintf_c,rintf_c
      USE Zzone        , ONLY: ncell_cond
!
      IMPLICIT NONE
!
!     input
      REAL(8) a(ncell_ps)
!     local variables
      INTEGER :: i,j,nb
!     local arrays
      REAL(8) :: svar(si_c(niut_c+1)-1),rvar(ri_c(niut_c+1)-1)
!       
      nb=1
!            
      DO i=1,si_c(niut_c+1)-1
         j=sintf_c(i)
         svar(i)=a(j)
      ENDDO
      CALL communicate_nb(niut_c,iut_c,si_c,ri_c,svar,rvar,nb)
      DO i=1,ncell_ps-ncell_cond
         j=rintf_c(i)
         a(j)=rvar(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d_c
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_2d_c(a)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_ps,niut_c,iut_c,si_c,ri_c,sintf_c,rintf_c
      USE Zparam       , ONLY: ndim
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a(ncell_ps,ndim)
!.....Local variables
      INTEGER :: i,j,ix,nb
      INTEGER :: ib,lda,i1,ip,ip1
!.....Local arrays
      REAL(8) :: svar((si_c(niut_c+1)-1)*ndim),rvar((ri_c(niut_c+1)-1)*ndim)
!       
      nb=ndim
!            
      ip=1
      DO ib=1,niut_c
         lda=si_c(ib+1)-si_c(ib)
         DO ix=1,nb
            ip1=ip+(ix-1)*lda
            DO i=si_c(ib),si_c(ib+1)-1
               i1=i-si_c(ib)
               j=sintf_c(i)
               svar(i1+ip1)=a(j,ix)
            ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
      CALL communicate_nb(niut_c,iut_c,si_c,ri_c,svar,rvar,nb)
      ip=1
      DO ib=1,niut_c
         lda=ri_c(ib+1)-ri_c(ib)
         DO ix=1,nb
            ip1=ip+(ix-1)*lda
            DO i=ri_c(ib),ri_c(ib+1)-1
               i1=i-ri_c(ib)
               j=rintf_c(i)
               a(j,ix)=rvar(i1+ip1)
            ENDDO
         ENDDO
         ip=ip+lda*nb
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_2d_c
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d_c_int(a)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: ncell_ps,niut_c,iut_c,si_c,ri_c,sintf_c,rintf_c
      USE Zzone        , ONLY: ncell_cond
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: a(ncell_ps)
!.....Local variables
      INTEGER :: i,j,nb
!.....Local arrays
      INTEGER :: svar(si_c(niut_c+1)-1),rvar(ri_c(niut_c+1)-1)
!       
      nb=1
!            
      DO i=1,si_c(niut_c+1)-1
         j=sintf_c(i)
         svar(i)=a(j)
      ENDDO
      CALL communicate_nb_int(niut_c,iut_c,si_c,ri_c,svar,rvar,nb)
      DO i=1,ncell_ps-ncell_cond
         j=rintf_c(i)
         a(j)=rvar(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d_c_int
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_1d_c_csr(a,i_neigh_c)
!DEC$IF defined (mpi_flag)
!
      USE Zmpi         , ONLY: maxmt_ps,ncell_ps,niut_c,iut_c,si_c,ri_c,sintf_c,rintf_c
      USE Zparam       , ONLY: ns
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: i_neigh_c(ncell_ps+1)
      REAL(8) :: a(maxmt_ps)
!.....Local variables
      INTEGER :: i,j,j0,nb
      INTEGER :: ib,lda,ip,i1,j1,ip1
!.....Local arrays
      REAL(8) :: svar(ns*(si_c(niut_c+1)-1)),rvar(ns*(ri_c(niut_c+1)-1))
!
      nb=ns
!
      ip=0
      DO ib=1,niut_c
         lda=si_c(ib+1)-si_c(ib)
            ip1=ip
            DO i=si_c(ib),si_c(ib+1)-1
               i1=ns*(i-si_c(ib))
               j1=sintf_c(i)
               j0=i_neigh_c(j1)-1
               DO j=i_neigh_c(j1),i_neigh_c(j1+1)-1
                  svar(j-j0+i1+ip1)=a(j)
               ENDDO
               DO j=i_neigh_c(j1+1),i_neigh_c(j1)+ns-1
                  svar(j-j0+i1+ip1)=0.d0
               ENDDO
            ENDDO
         ip=ip+lda*nb
      ENDDO
      CALL communicate_nb(niut_c,iut_c,si_c,ri_c,svar,rvar,nb)
      ip=0
      DO ib=1,niut_c
         lda=ri_c(ib+1)-ri_c(ib)
            ip1=ip
            DO i=ri_c(ib),ri_c(ib+1)-1
               i1=ns*(i-ri_c(ib))
               j1=rintf_c(i)
               j0=i_neigh_c(j1)-1
               DO j=i_neigh_c(j1),i_neigh_c(j1+1)-1
                  a(j)=rvar(j-j0+i1+ip1)
               ENDDO
            ENDDO
         ip=ip+lda*nb
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_1d_c_csr
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_rv_2d(a1,a2)
!DEC$IF defined (mpi_flag)
!
      USE Zrv_mpi      , ONLY: niut_fuel_rod,iut_fuel_rod,si_fuel_rod,ri_fuel_rod,sintf_fuel_rod, &
                               rintf_fuel_rod,ncell_fuel_rod_p
      USE Zrv_hts_2d   , ONLY: nr_2d
!
      IMPLICIT NONE
!.....Input
      REAL(8),INTENT(INOUT) :: a1(ncell_fuel_rod_p,nr_2d)
      REAL(8),DIMENSION(ncell_fuel_rod_p,nr_2d),OPTIONAL :: a2
!.....Local variables
      INTEGER :: i,j,i1,ib,ix
      INTEGER :: ip,ip1,ip2
      INTEGER :: nb,nb1,lda
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: svar,rvar
!
      nb=1
      IF(.not.PRESENT(a2)) GOTO 100
      nb=nb+1
100   CONTINUE
!
      nb1=nb*nr_2d
      ALLOCATE(svar((si_fuel_rod(niut_fuel_rod+1)-1)*nr_2d*nb), &
               rvar((ri_fuel_rod(niut_fuel_rod+1)-1)*nr_2d*nb))
!       
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut_fuel_rod
            lda=si_fuel_rod(ib+1)-si_fuel_rod(ib)
            DO ix=1,nr_2d
               ip1=ip+(ix-1)*lda
               DO i=si_fuel_rod(ib),si_fuel_rod(ib+1)-1
                  i1=i-si_fuel_rod(ib)
                  j=sintf_fuel_rod(i)
                  svar(i1+ip1)=a1(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut_fuel_rod
            lda=si_fuel_rod(ib+1)-si_fuel_rod(ib)
            DO ix=1,nr_2d
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*nr_2d
               DO i=si_fuel_rod(ib),si_fuel_rod(ib+1)-1
                  i1=i-si_fuel_rod(ib)
                  j=sintf_fuel_rod(i)
                  svar(i1+ip1)=a1(j,ix)
                  svar(i1+ip2)=a2(j,ix)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ENDIF
      CALL communicate_nb(niut_fuel_rod,iut_fuel_rod,si_fuel_rod,ri_fuel_rod,svar,rvar,nb1)
      IF(nb.eq.1) THEN
         ip=1
         DO ib=1,niut_fuel_rod
            lda=ri_fuel_rod(ib+1)-ri_fuel_rod(ib)
            DO ix=1,nr_2d
               ip1=ip+(ix-1)*lda
               DO i=ri_fuel_rod(ib),ri_fuel_rod(ib+1)-1
                  i1=i-ri_fuel_rod(ib)
                  j=rintf_fuel_rod(i)
                  a1(j,ix)=rvar(i1+ip1)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ELSEIF(nb.eq.2) THEN
         ip=1
         DO ib=1,niut_fuel_rod
            lda=ri_fuel_rod(ib+1)-ri_fuel_rod(ib)
            DO ix=1,nr_2d
               ip1=ip+(ix-1)*lda
               ip2=ip1+lda*nr_2d
               DO i=ri_fuel_rod(ib),ri_fuel_rod(ib+1)-1
                  i1=i-ri_fuel_rod(ib)
                  j=rintf_fuel_rod(i)
                  a1(j,ix)=rvar(i1+ip1)
                  a2(j,ix)=rvar(i1+ip2)
               ENDDO
            ENDDO
            ip=ip+lda*nb1
         ENDDO
      ENDIF
      DEALLOCATE(svar,rvar)
!DEC$ENDIF
!
      END SUBROUTINE communicate_rv_2d
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate(a,izone)
!DEC$IF defined (mpi_flag)
      USE Zinterface
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: izone
      REAL(8) :: a(*)
!
      IF(izone.eq.0)THEN
         CALL communicate_1d(a)
      ELSE
         CALL communicate_1d_c(a)
      ENDIF
!DEC$ENDIF
!
      END SUBROUTINE communicate
!
!------------------------------------------------------------------------------
!     
      SUBROUTINE communicate_nb(niut,iut,si,ri,svar,rvar,nb)
!DEC$IF defined (mpi_flag)
!        
!     This routine is a general "communicate" function for MPI
!     
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: nb
      INTEGER :: niut
      INTEGER :: iut(niut)
      INTEGER :: si(niut+1),ri(niut+1)
      REAL(8) :: svar(*)
!.....Output
      REAL(8) :: rvar(*)
!.....Local variables
      INTEGER :: ib,ldas,ips,ldar,ipr
      INTEGER :: ierr,tag
!.....Local arrays
      INTEGER :: request(2*niut)
      INTEGER :: status(MPI_STATUS_SIZE)
!     
      tag=1
      ips=1
      ipr=1
      DO ib=1,niut
         ldas=(si(ib+1)-si(ib))*nb
         ldar=(ri(ib+1)-ri(ib))*nb
         CALL MPI_ISEND(svar(ips),ldas,             &
              MPI_DOUBLE_PRECISION,iut(ib),tag,     &
              MPI_COMM_WORLD,request(ib),ierr)
         
         CALL MPI_IRECV(rvar(ipr),ldar,             &
              MPI_DOUBLE_PRECISION,iut(ib),tag,     &
              MPI_COMM_WORLD,request(ib+niut),ierr)
         ips=ips+ldas
         ipr=ipr+ldar
      ENDDO
!
      DO ib=1,niut
         CALL MPI_WAIT(request(ib     ),status,ierr)
         CALL MPI_WAIT(request(ib+niut),status,ierr)
      ENDDO
!
!DEC$ENDIF
!
      END SUBROUTINE communicate_nb
!
!------------------------------------------------------------------------------
!     
      SUBROUTINE communicate_nb_int(niut,iut,si,ri,svar,rvar,nb)
!DEC$IF defined (mpi_flag)
!        
!     This routine is a general "communicate" function for MPI
!     
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: nb
      INTEGER :: niut
      INTEGER :: iut(niut)
      INTEGER :: si(niut+1),ri(niut+1)
      INTEGER :: svar(*)
!.....Output
      INTEGER :: rvar(*)
!.....Local variables
      INTEGER :: ib,ldas,ips,ldar,ipr
      INTEGER :: ierr,tag
!.....Local arrays
      INTEGER :: request(2*niut)
      INTEGER :: status(MPI_STATUS_SIZE)
!     
      tag = 1
      ips=1
      ipr=1
      DO ib=1,niut
         ldas=si(ib+1)-si(ib)
         ldar=ri(ib+1)-ri(ib)
         CALL MPI_ISEND(svar(ips),ldas*nb,         &
              MPI_INTEGER,iut(ib),tag,              &
              MPI_COMM_WORLD,request(ib),ierr)
         
         CALL MPI_IRECV(rvar(ipr),ldar*nb,         &
              MPI_INTEGER,iut(ib),tag,              &
              MPI_COMM_WORLD,request(ib+niut),ierr)
         ips=ips+ldas*nb
         ipr=ipr+ldar*nb
      ENDDO
!
      DO ib=1,niut
         CALL MPI_WAIT(request(ib     ),status,ierr)
         CALL MPI_WAIT(request(ib+niut),status,ierr)
      ENDDO
!
!DEC$ENDIF
!
      END SUBROUTINE communicate_nb_int
!
!------------------------------------------------------------------------------
!
      SUBROUTINE mwait(niut,request)
!DEC$IF defined (mpi_flag)
!
!     This routine is a general "wait" function for MPI 
!
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
      Integer niut,i
      INTEGER status(MPI_STATUS_SIZE),ierr
      INTEGER request(2*niut)
!
      DO i=1,niut
         CALL MPI_WAIT(request(i),status,ierr)
      ENDDO  
!
!DEC$ENDIF
!
      END SUBROUTINE mwait

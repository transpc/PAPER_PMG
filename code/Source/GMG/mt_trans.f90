! ====================================================================!
! ---CSR format of transpose matrix -----------------------
! ====================================================================!
      subroutine mt_trans(nr,nc,nnz,ia,ja,ia1,ja1,a,a1)
!      use md_function
! ---
      implicit none
! ---inlet
      integer nr,nc,nnz
      integer ia(nr+1),ja(nnz)
      real*8 a(nnz)
! ---outlet
      integer ia1(nc+1),ja1(nnz)
      real*8 a1(nnz)
! ---temp
      integer i,j,k,jj,i1,i2,init, j1, j2, nn
!      integer irow,icsr
      integer,dimension(:),allocatable::iwk, ia2
! --- 
      allocate(iwk(nc))
! ------------------------
    
!   iwk for each colume
        iwk = 0
        DO k=1, nnz
            j = ja(k)
            iwk(j) = iwk(j) + 1
        ENDDO
        
! ia1:
        ia1(1) = 1
        DO i = 1,nc
            ia1(i+1) = ia1(i) + iwk(i)
        ENDDO
! 
        IF(ia1(nc+1).NE.(nnz+1)) THEN
            WRITE(*,*)'PMG error in matrix-transpose, ia1(nc+1)=/nnz+1'
            STOP
        ENDIF
!
        DEALLOCATE(iwk)
        ALLOCATE(ia2(nc+1))

! ja1: using temp ia2
        ia2 = ia1
! 
        DO i=1,nr
            j1 = ia(i)
            j2 = ia(i+1)-1
            DO k = j1,j2
                j = ja(k)
                nn = ia2(j)
                ja1(nn) = i
                a1(nn) = a(k)
                ia2(j) = nn+1
            ENDDO
        ENDDO
 
! test
        DO i = 1,nc
            IF(ia2(i).NE.ia1(i+1)) THEN
                WRITE(*,*)'error in ia2, mt_trans'
                STOP
            ENDIF
        ENDDO
!
! set values: a1
 !       DO i = 1,nnz
 !           a1(i) = a(i)
 !       ENDDO
        
! 
        DEALLOCATE(ia2)
! --- 
      return
    End

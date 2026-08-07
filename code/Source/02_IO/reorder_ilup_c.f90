      SUBROUTINE reorder_ilup_c
!
      USE Zmpi         , ONLY: maxmt_c,maxmt_lu0_c,maxmt_lu1_c,                 &
                               ia_a_c,ja_a_c,ju_a_c,iend_c,                     &
                               diag_lu_c,alu0_c,alu1_c,ia0_c,ia1_c,ja0_c,ja1_c, &
                               levt_c,lev_typet_c,                              &
                               maxmt2_c,maxmt_r_c,au_r_c,ia_r_c,ja_r_c,ju_r_c,  &
                               perm_r_c,permi_r_c,index_r_c
      USE Zzone        , ONLY: ncell_cond
      USE Zbicg        , ONLY: lev_type_c,levmpi_type_c,lev_c,levmpi_c
      USE Zcore        , ONLY: np,myrank
!
      IMPLICIT NONE
!
!.....Input
!.....Output
!.....Local variables
      INTEGER :: i,ip
      INTEGER :: jj,l,len,isize
!.....Local arrays
      INTEGER :: iend_r_c(ncell_cond)
      INTEGER :: iar_c(ncell_cond+1)
      INTEGER,DIMENSION(:),ALLOCATABLE :: jar_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: ja_w_c,ia_w_c,idiag_w_c,lev_w_c
!
      If(np.eq.1) THEN
         levt_c=lev_c
         lev_typet_c=lev_type_c
      ELSE
         levt_c=levmpi_c
         lev_typet_c=levmpi_type_c
      ENDIF
!
      If(np.eq.1) THEN
!........Check for tridiagonal
            i=1
            IF(ia_a_c(i+1)-ia_a_c(i).ne.2) GOTO 101
            IF(ja_a_c(1).ne.i  ) GOTO 101
            IF(ja_a_c(2).ne.i+1) GOTO 101
         DO i=2,ncell_cond-1
            IF(ia_a_c(i+1)-ia_a_c(i).ne.3) GOTO 101
            jj=ia_a_c(i)
            IF(ja_a_c(jj  ).ne.i-1) GOTO 101
            IF(ja_a_c(jj+1).ne.i  ) GOTO 101
            IF(ja_a_c(jj+2).ne.i+1) GOTO 101
         ENDDO
            i=ncell_cond
            IF(ia_a_c(i+1)-ia_a_c(i).ne.2) GOTO 101
            jj=ia_a_c(i)
            IF(ja_a_c(jj  ).ne.i-1) GOTO 101
            IF(ja_a_c(jj+1).ne.i  ) GOTO 101
         levt_c=-2
!        write(*,*) 'tridiagonal'
         GOTO 102
101      CONTINUE
!        write(*,*) 'not tridiagonal',i
102      CONTINUE
      ENDIF
!     stop 10
!
!     write(*,*) 'ncell_cond==>',ncell_cond
!     stop 99
      IF(myrank.eq.0) THEN
         IF(levt_c     .eq.-1) WRITE(*,*) '          Requested for solid: direct solver.'
         IF(levt_c     .eq. 0) WRITE(*,*) '          Requested for solid: ilu0 preconditionner.'
         IF(levt_c     .eq. 1) WRITE(*,*) '          Requested for solid: ilu1 preconditionner.'
         IF(levt_c     .eq. 2) WRITE(*,*) '          Requested for solid: ilu2 preconditionner.'
         IF(lev_typet_c.eq. 0) WRITE(*,*) '          Requested for solid: no reordering.'
         IF(lev_typet_c.eq. 1) WRITE(*,*) '          Requested for solid: METIS reordering.'
         IF(lev_typet_c.eq. 2) WRITE(*,*) '          Requested for solid: CUTHILL reordering.'
      ENDIF
!
!     write(*,*) '====>lev_typet_c',levt_c,lev_typet_c
!
!.....Must always allocate used to call csr_cg_solver
!
      ALLOCATE(perm_r_c(ncell_cond),permi_r_c(ncell_cond))
      IF(levt_c.eq.0) THEN
         IF(lev_typet_c.eq.0) THEN
!
!...........No reordering, lu0 default
!
            maxmt_lu0_c=0
            maxmt_lu1_c=0
            DO i=1,ncell_cond
               l=ju_a_c(i)-ia_a_c(i)
               maxmt_lu0_c=maxmt_lu0_c+l
               l=iend_c(i)-ju_a_c(i)
               maxmt_lu1_c=maxmt_lu1_c+l
            ENDDO
         ELSE
!
!...........Build a till iend without diagonal
!
            iar_c(1)=1
            DO i=1,ncell_cond
               l=iend_c(i)-ia_a_c(i)
               iar_c(i+1)=iar_c(i)+l
            ENDDO
            len=iar_c(ncell_cond+1)-1
            ALLOCATE(jar_c(len))
            maxmt2_c=len+ncell_cond
            ip=1
            DO i=1,ncell_cond
               DO jj=ia_a_c(i),ju_a_c(i)-1
                  jar_c(ip)=ja_a_c(jj)
                  ip=ip+1
               ENDDO
               DO jj=ju_a_c(i)+1,iend_c(i)
                  jar_c(ip)=ja_a_c(jj)
                  ip=ip+1
               ENDDO
            ENDDO
!
!...........Choose reordering scheme
!
            write(*,*) 'lev_typet_c',lev_typet_c
!           IF(nnp.eq.0) THEN
               IF(lev_typet_c.eq.1) THEN
                  CALL reorder_metis(ncell_cond,iar_c,jar_c, &
                                     perm_r_c,permi_r_c,len)
               ELSEIF(lev_typet_c.eq.2) THEN
                  CALL reorder_cuthill(ncell_cond,iar_c,jar_c, &
                                       perm_r_c,permi_r_c,len)
               ENDIF
!              write(*,*) '9==>',maxmt,maxmt2_c
!           ELSE
!           ENDIF
            ALLOCATE(ia_r_c(ncell_cond+1))
            ALLOCATE(ju_r_c(ncell_cond))
            ALLOCATE(ja_r_c(maxmt2_c))
            ALLOCATE(index_r_c(maxmt2_c))
            CALL final_r(ncell_cond,maxmt_c,maxmt2_c, &
                         ia_a_c,iend_c,ja_a_c,        &
                         ia_r_c,ju_r_c,ja_r_c,        &
                         perm_r_c,permi_r_c,index_r_c)
            maxmt_lu0_c=0
            maxmt_lu1_c=0
            DO i=1,ncell_cond
               l=ju_r_c(i)-ia_r_c(i)
               maxmt_lu0_c=maxmt_lu0_c+l
               l=ia_r_c(i+1)-1-ju_r_c(i)
               maxmt_lu1_c=maxmt_lu1_c+l
            ENDDO
         ENDIF
      ELSEIF(levt_c.ge.-1) THEN
         IF(lev_typet_c.eq.0) THEN
!
!...........No reordering, lup ia_r,ju_r,ja_r array till iend
!
            ALLOCATE(ia_r_c(ncell_cond+1))
            ia_r_c(1)=1
            DO i=1,ncell_cond
               l=iend_c(i)-ia_a_c(i)+1
               ia_r_c(i+1)=ia_r_c(i)+l
            ENDDO
            maxmt2_c=ia_r_c(ncell_cond+1)-1
            ALLOCATE(ju_r_c(ncell_cond))
            ALLOCATE(ja_r_c(maxmt2_c))
            ip=1
            DO i=1,ncell_cond
               DO jj=ia_a_c(i),iend_c(i)
                  ja_r_c(ip)=ja_a_c(jj)
                  IF(ja_a_c(jj).eq.i) THEN
                     ju_r_c(i)=ip
                     ip=ip+1
                     EXIT
                  ENDIF
                  ip=ip+1
               ENDDO
               DO jj=ju_a_c(i)+1,iend_c(i)
                  ja_r_c(ip)=ja_a_c(jj)
                  ip=ip+1
               ENDDO
            ENDDO
         ELSE
!
!...........Build a till iend without diagonal
!
            iar_c(1)=1
            DO i=1,ncell_cond
               l=iend_c(i)-ia_a_c(i)
               iar_c(i+1)=iar_c(i)+l
            ENDDO
            len=iar_c(ncell_cond+1)-1
            ALLOCATE(jar_c(len))
            maxmt2_c=len+ncell_cond
            ip=1
            DO i=1,ncell_cond
               DO jj=ia_a_c(i),ju_a_c(i)-1
                  jar_c(ip)=ja_a_c(jj)
                  ip=ip+1
               ENDDO
               DO jj=ju_a_c(i)+1,iend_c(i)
                  jar_c(ip)=ja_a_c(jj)
                  ip=ip+1
               ENDDO
            ENDDO
!
!...........Choose reordering scheme
!
!           write(*,*) 'lev_typet_c',lev_typet_c
!           IF(nnp.eq.0) THEN
               IF(lev_typet_c.eq.1) THEN
                  CALL reorder_metis(ncell_cond,iar_c,jar_c, &
                                     perm_r_c,permi_r_c,len)
               ELSEIF(lev_typet_c.eq.2) THEN
                  CALL reorder_cuthill(ncell_cond,iar_c,jar_c, &
                                       perm_r_c,permi_r_c,len)
               ENDIF
!              write(*,*) '9==>',maxmt,maxmt2_c
!           ELSE
!           ENDIF
            ALLOCATE(ia_r_c(ncell_cond+1))
            ALLOCATE(ju_r_c(ncell_cond))
            ALLOCATE(ja_r_c(maxmt2_c))
            ALLOCATE(index_r_c(maxmt2_c))
            CALL final_r(ncell_cond,maxmt_c,maxmt2_c,  &
                         ia_a_c,iend_c,ja_a_c,         &
                         ia_r_c,ju_r_c,ja_r_c,         &
                         perm_r_c,permi_r_c,index_r_c)
!
         ENDIF
!
!........Do symbolic factorisation in ia_w,idiag_w,ja_w
!
401      format(11x,a,f5.2)
            isize=150*maxmt2_c
            write(*,*) '1-==>',isize,maxmt2_c
            ALLOCATE(lev_w_c(isize))
            ALLOCATE(ia_w_c(ncell_cond+1))
            ALLOCATE(ja_w_c(isize))
            ALLOCATE(idiag_w_c(ncell_cond))
            IF(levt_c.eq.-1) THEN
               CALL dgefs2s(ncell_cond,maxmt2_c,isize, &
                            ia_r_c,ja_r_c,             &
                            ia_w_c,ja_w_c,idiag_w_c)
            ELSE
               CALL dgefs2s_lev(ncell_cond,maxmt2_c,isize,              &
                                ia_r_c,ja_r_c,                          &
                                ia_w_c,ja_w_c,idiag_w_c,lev_w_c,levt_c)
            ENDIF
            IF(myrank.eq.0) THEN
               IF(ncell_cond.gt.0) THEN
                  WRITE(*,401) 'fill in ratioc: ',dble(ia_w_c(ncell_cond+1)-1)/(ia_r_c(ncell_cond+1)-1)
               ENDIF
            ENDIF
!
!...........Size with fill-in computed put back in ia_r,ju_r,ha_r
!
            maxmt_r_c=ia_w_c(ncell_cond+1)-1
            DEALLOCATE(ja_r_c)
            ALLOCATE(ja_r_c(maxmt_r_c))
            ALLOCATE(au_r_c(maxmt_r_c))
            do i=1,ncell_cond
               ju_r_c(i)=idiag_w_c(i)
               iend_r_c(i)=ia_w_c(i+1)-1
            enddo
            do i=1,ncell_cond+1
               ia_r_c(i)=ia_w_c(i)
            enddo
            do i=1,maxmt_r_c
               ja_r_c(i)=ja_w_c(i)
            enddo
            maxmt_lu0_c=0
            maxmt_lu1_c=0
            DO i=1,ncell_cond
               l=ju_r_c(i)-ia_r_c(i)
               maxmt_lu0_c=maxmt_lu0_c+l
               l=iend_r_c(i)-ju_r_c(i)
               maxmt_lu1_c=maxmt_lu1_c+l
            ENDDO
      ELSEIF(levt_c.eq.-2) THEN
         RETURN
      ENDIF
!
      ALLOCATE(diag_lu_c(ncell_cond))
      ALLOCATE(alu0_c(maxmt_lu0_c),alu1_c(maxmt_lu1_c))
      ALLOCATE(ja0_c(maxmt_lu0_c),ja1_c(maxmt_lu1_c))
      ALLOCATE(ia0_c(ncell_cond+1),ia1_c(ncell_cond+1))
!
      END SUBROUTINE reorder_ilup_c

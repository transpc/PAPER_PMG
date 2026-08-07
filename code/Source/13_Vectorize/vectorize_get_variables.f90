!
!------------------------------------------------------------------------------
!
      SUBROUTINE get_vector_disp(j,i,i1)
!
!     This routine maps the j,i to offset i1 in vector space
!     It returns negative offset for k < i
!
      USE Znum_cell    , ONLY: i_neigh,neigh, &
                               istart_nf,istart_nb1, &
                               ia_nb,icell_nb,right_nb_k,iptr_nb_k
      USE Zvec_index   , ONLY: jneigh_nf
      USE Zbc_index    , ONLY: nbcon
!
      IMPLICIT NONE
!
!.....External function
      INTEGER :: get_nf_number
!.....Input
      INTEGER :: i,j,j0
!.....Output
      INTEGER :: i1
!.....Local variables
      INTEGER :: m,ii,k,nb
      INTEGER :: nf_number,istart,len,istart1,i0
!
         j0=i_neigh(i)-1
            k=neigh(j+j0)
            nf_number=get_nf_number(nbcon(j+j0))
            IF(nf_number.eq.0 .and. k.lt.i) THEN
!
!..............Cells non_k
!
               nb=iptr_nb_k(i)
               m=ia_nb(nb)+j-1
               i1=-right_nb_k(m)
            ELSE
!
!..............The rest
!
               istart=istart_nf(1,nf_number)
               istart1=istart_nb1(1,nf_number)
               len   =istart_nb1(2,nf_number)
               DO nb=1,len
                  i0=istart1+nb
                  ii=icell_nb(i0)
                  IF(ii.eq.i) THEN
                     DO m=ia_nb(i0),ia_nb(i0+1)-1
                        i1=istart+m
                        IF(jneigh_nf(i1).eq.j) goto 100
                     ENDDO
                  ENDIF
               ENDDO
100            CONTINUE
            ENDIF
!
      END SUBROUTINE get_vector_disp
!
!------------------------------------------------------------------------------
!
      SUBROUTINE get_scalar_variable_p_ndim(x_nf,x)
!
      USE Zmpi         , ONLY: maxmt_fluid
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_tot
      USE Znum_cell    , ONLY: nb_nf,i_neigh,istart_nf
      USE Znum_cell    , ONLY: istart_nb1,                          &
                               ia_nb,icell_nb,iptr_nb_k,right_nb_k
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: x_nf(nf_tot,ndim)
!.....Output
      REAL(8) :: x(maxmt_fluid,ndim)
!.....Local variables
      INTEGER :: i,j,k,ix
      INTEGER :: ii,ip
      INTEGER :: nf_number,istart,len,istart1,i0,i1
!.....Local arrays
      INTEGER :: ip1(0:nb_nf)
!
      DO nf_number=0,nb_nf
         ip1(nf_number)=1
      ENDDO
      DO i=1,ncell_fluid
         ii=iptr_nb_k(i)
         j=i_neigh(i)-1
         IF(ii.gt.0) then
            DO ip=ia_nb(ii),ia_nb(ii+1)-1
               k=right_nb_k(ip)
               j=j+1
               DO ix=1,ndim
                  x(j,ix)=x_nf(k,ix)
               ENDDO
            ENDDO
         ENDIF
         DO nf_number=0,nb_nf
            ii=ip1(nf_number)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            IF(ii.gt.len) cycle
            i1=istart1+ii
            IF(icell_nb(i1).eq.i) THEN
               istart =istart_nf(1,nf_number)
               DO ip=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+ip
                  j=j+1
                  DO ix=1,ndim
                     x(j,ix)=x_nf(i0,ix)
                  ENDDO
               ENDDO
               ip1(nf_number)=ip1(nf_number)+1
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE get_scalar_variable_p_ndim
!
!------------------------------------------------------------------------------
!
      SUBROUTINE get_scalar_variable_n_ndim(x_nf,x)
!
      USE Zmpi         , ONLY: maxmt_fluid
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_tot
      USE Znum_cell    , ONLY: nb_nf,i_neigh,istart_nf,             &
                               istart_nb1,                          &
                               ia_nb,icell_nb,iptr_nb_k,right_nb_k
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: x_nf(nf_tot,ndim)
      REAL(8) :: x(maxmt_fluid,ndim)
!.....Local variables
      INTEGER :: i,j,k,ix
      INTEGER :: ii,ip
      INTEGER :: nf_number,istart,len,istart1,i0,i1
!.....Local arrays
      INTEGER :: ip1(0:nb_nf)
!
      DO nf_number=0,nb_nf
         ip1(nf_number)=1
      ENDDO
      DO i=1,ncell_fluid
         ii=iptr_nb_k(i)
         j=i_neigh(i)-1
         IF(ii.gt.0) then
            DO ip=ia_nb(ii),ia_nb(ii+1)-1
               k=right_nb_k(ip)
               j=j+1
               DO ix=1,ndim
                  x(j,ix)=-x_nf(k,ix)
               ENDDO
            ENDDO
         ENDIF
         DO nf_number=0,nb_nf
            ii=ip1(nf_number)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            IF(ii.gt.len) cycle
            i1=istart1+ii
            IF(icell_nb(i1).eq.i) THEN
               istart =istart_nf(1,nf_number)
               DO ip=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+ip
                  j=j+1
                  DO ix=1,ndim
                     x(j,ix)=x_nf(i0,ix)
                  ENDDO
               ENDDO
               ip1(nf_number)=ip1(nf_number)+1
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE get_scalar_variable_n_ndim
!
!------------------------------------------------------------------------------
!
      SUBROUTINE get_scalar_variable_p(x_nf,x)
!
      USE Zmpi         , ONLY: maxmt_fluid
      USE Zzone        , ONLY: ncell_fluid
      USE Zvec_param   , ONLY: nf_tot
      USE Znum_cell    , ONLY: nb_nf,i_neigh,istart_nf,             &
                               istart_nb1,                          &
                               ia_nb,icell_nb,iptr_nb_k,right_nb_k
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: x_nf(nf_tot)
      REAL(8) :: x(maxmt_fluid)
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: ii,ip
      INTEGER :: nf_number,istart,len,istart1,i0,i1
!.....Local arrays
      INTEGER :: ip1(0:nb_nf)
!
      DO nf_number=0,nb_nf
         ip1(nf_number)=1
      ENDDO
      DO i=1,ncell_fluid
         ii=iptr_nb_k(i)
         j=i_neigh(i)-1
         IF(ii.gt.0) then
            DO ip=ia_nb(ii),ia_nb(ii+1)-1
               k=right_nb_k(ip)
               j=j+1
               x(j)=x_nf(k)
            ENDDO
         ENDIF
         DO nf_number=0,nb_nf
            ii=ip1(nf_number)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            IF(ii.gt.len) cycle
            i1=istart1+ii
            IF(icell_nb(i1).eq.i) THEN
               istart =istart_nf(1,nf_number)
               DO ip=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+ip
                  j=j+1
                  x(j)=x_nf(i0)
               ENDDO
               ip1(nf_number)=ip1(nf_number)+1
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE get_scalar_variable_p
!
!------------------------------------------------------------------------------
!
      SUBROUTINE get_scalar_variable_n(x_nf,x)
!
      USE Zmpi         , ONLY: maxmt_fluid
      USE Zzone        , ONLY: ncell_fluid
      USE Zvec_param   , ONLY: nf_tot
      USE Znum_cell    , ONLY: nb_nf,i_neigh,istart_nf,             &
                               istart_nb1,                          &
                               ia_nb,icell_nb,iptr_nb_k,right_nb_k
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: x_nf(nf_tot)
      REAL(8) :: x(maxmt_fluid)
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: ii,ip
      INTEGER :: nf_number,istart,len,istart1,i0,i1
!.....Local arrays
      INTEGER :: ip1(0:nb_nf)
!
      DO nf_number=0,nb_nf
         ip1(nf_number)=1
      ENDDO
      DO i=1,ncell_fluid
         ii=iptr_nb_k(i)
         j=i_neigh(i)-1
         IF(ii.gt.0) then
            DO ip=ia_nb(ii),ia_nb(ii+1)-1
               k=right_nb_k(ip)
               j=j+1
               x(j)=-x_nf(k)
            ENDDO
         ENDIF
         DO nf_number=0,nb_nf
            ii=ip1(nf_number)
            istart1=istart_nb1(1,nf_number)
            len   =istart_nb1(2,nf_number)
            IF(ii.gt.len) cycle
            i1=istart1+ii
            IF(icell_nb(i1).eq.i) THEN
               istart =istart_nf(1,nf_number)
               DO ip=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+ip
                  j=j+1
                  x(j)=x_nf(i0)
               ENDDO
               ip1(nf_number)=ip1(nf_number)+1
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE get_scalar_variable_n
!
!------------------------------------------------------------------------------
!
      SUBROUTINE get_scalar_variable_n_i_ndim(x_nf,x,i,ix)
!
      USE Zparam       , ONLY: ndim,ns
      USE Zvec_param   , ONLY: nf_tot
      USE Znum_cell    , ONLY: nb_nf,istart_nf,                     &
                               istart_nb1,                          &
                               ia_nb,icell_nb,iptr_nb_k,right_nb_k
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: ix
      REAL(8) :: x_nf(nf_tot,ndim)
      REAL(8) :: x(ns)
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: ii,ip
      INTEGER :: nf_number,istart,len,istart1,i0,i1
!.....Local arrays
      INTEGER :: ip1(0:nb_nf)
!
      DO nf_number=0,nb_nf
         ip1(nf_number)=1
      ENDDO
         ii=iptr_nb_k(i)
         j=0
         IF(ii.gt.0) then
            DO ip=ia_nb(ii),ia_nb(ii+1)-1
               k=right_nb_k(ip)
               j=j+1
               x(j)=-x_nf(k,ix)
            ENDDO
         ENDIF
         DO nf_number=0,nb_nf
            ii=ip1(nf_number)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            IF(ii.gt.len) cycle
            i1=istart1+ii
100         CONTINUE
            IF(icell_nb(i1).lt.i) THEN
               ii=ii+1
               IF(ii.gt.len) goto 110
               i1=istart1+ii
               ip1(nf_number)=ii
               GOTO 100
            ENDIF
            IF(icell_nb(i1).eq.i) THEN
               istart =istart_nf(1,nf_number)
               DO ip=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+ip
                  j=j+1
                  x(j)=x_nf(i0,ix)
               ENDDO
               ip1(nf_number)=ip1(nf_number)+1
            ENDIF
110         CONTINUE
         ENDDO
!
      END SUBROUTINE get_scalar_variable_n_i_ndim
!
!------------------------------------------------------------------------------
!
      FUNCTION get_nf_number(nbcon)
!
!     This function return the nf_number associated with neighbor (j,i)
!
      USE Zparam       , ONLY: nin_max,nb_max
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: nbcon
!.....Output
      INTEGER :: get_nf_number
!
      IF(nbcon.eq.0)THEN
         get_nf_number=0
      ELSEIF(nbcon.ge.201)THEN
         get_nf_number=1
      ELSEIF(nbcon.gt.0.and.nbcon.le.nin_max)THEN
         get_nf_number=2
      ELSEIF(nbcon.gt.nin_max.and.nbcon.le.nb_max)THEN
         get_nf_number=3
      ELSEIF(nbcon.eq.-1)THEN
         get_nf_number=4
      ELSEIF(nbcon.eq.-2)THEN
         get_nf_number=5
      ELSEIF(nbcon.eq.-3.or.nbcon.eq.-4)THEN
         get_nf_number=6
      ELSEIF(nbcon.eq.-5.or.nbcon.eq.-6)THEN
         get_nf_number=7
      ELSEIF(nbcon.eq.101)THEN
         get_nf_number=8
      ELSEIF(nbcon.ge.-39 .and. nbcon.le.-31)THEN
         get_nf_number=9
      ELSEIF(nbcon.ge.-49 .and. nbcon.le.-41)THEN
         get_nf_number=10
      ELSEIF(nbcon.ge.-59 .and. nbcon.le.-51)THEN
         get_nf_number=11
      ENDIF
!
      END FUNCTION get_nf_number
!
      FUNCTION get_nf_number_j(j)
!
!     This function return the nf_number associated with neighbor (j,i)
!
      USE Zparam       , ONLY: nin_max,nb_max
      USE Zbc_index    , ONLY: nbcon
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: j
!.....Output
      INTEGER :: get_nf_number_j
!.....Local
      INTEGER nbcon0      
!
      nbcon0=nbcon(j)
      IF(nbcon0.eq.0)THEN
         get_nf_number_j=0
      ELSEIF(nbcon0.ge.201)THEN
         get_nf_number_j=1
      ELSEIF(nbcon0.gt.0.and.nbcon0.le.nin_max)THEN
         get_nf_number_j=2
      ELSEIF(nbcon0.gt.nin_max.and.nbcon0.le.nb_max)THEN
         get_nf_number_j=3
      ELSEIF(nbcon0.eq.-1)THEN
         get_nf_number_j=4
      ELSEIF(nbcon0.eq.-2)THEN
         get_nf_number_j=5
      ELSEIF(nbcon0.eq.-3.or.nbcon0.eq.-4)THEN
         get_nf_number_j=6
      ELSEIF(nbcon0.eq.-5.or.nbcon0.eq.-6)THEN
         get_nf_number_j=7
      ELSEIF(nbcon0.eq.101)THEN
         get_nf_number_j=8
      ELSEIF(nbcon(j).ge.-39 .and. nbcon(j).le.-31)THEN
         get_nf_number_j=9
      ELSEIF(nbcon(j).ge.-49 .and. nbcon(j).le.-41)THEN
         get_nf_number_j=10
      ELSEIF(nbcon(j).ge.-59 .and. nbcon(j).le.-51)THEN
         get_nf_number_j=11
      ENDIF
!
      END FUNCTION get_nf_number_j
!
      FUNCTION get_global_cell(i)
!
!     This function return the nf_number associated with neighbor (j,i)
!
      USE Zmpi   , ONLY: jperm,mapping_ext
      USE Zzone  , ONLY: ncell_fluid
      IMPLICIT NONE
!
!.....Input
      INTEGER :: i
!.....Output
      INTEGER :: get_global_cell
!.....Local
!      INTEGER nbcon0      
!
      If(i.le.ncell_fluid) THEN
         get_global_cell=jperm(i)
      ELSE
         get_global_cell=mapping_ext(i)
      ENDIF
!
      END FUNCTION get_global_cell       

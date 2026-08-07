!
      SUBROUTINE find_icelltype
!
!.....This routine find icell_type, wall_cell, xfc/xn of wallcell
!
      USE Zinterface
      USE Zmpi         , ONLY: jperm
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: nn,ndim,nin_max,nb_max,nb_sym,nb_mars
      USE Znum_cell    , ONLY: i_neigh
      USE Zbc_index    , ONLY: nbcon,npb,icell_type,icell_type_tmp,iface_wall_tmp, &
                               iface_wall,num_wallcells,wall_cell
      USE Zcoord2      , ONLY: xfc,xfc_wallcell
      USE Zcoord4      , ONLY: dji,sa
      USE Znormal      , ONLY: xn,wall_cell_l,num_wallcells_l,xn_wallcell_l,sa_wallcell_l,xn_wallcell
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,ii,i0,j0,j1,k,ix
      INTEGER :: ip,itype
      REAL(8) :: dr_min
!.....Local arrays
      INTEGER :: num_wallcells_p(np),num_wallcells_p_dsp(np)
      REAL(8), ALLOCATABLE::xfc_wallcell_l(:,:)
      INTEGER, ALLOCATABLE::list(:),indexo_l(:)
      INTEGER, ALLOCATABLE::windexo_l(:)
      REAL(8), ALLOCATABLE::wxfc_wallcell(:,:)
      REAL(8), ALLOCATABLE::wxn_wallcell(:,:)
!
      ALLOCATE(icell_type_tmp(nn),iface_wall_tmp(nn))
!
!........Define cell type depending on face condition
!        icell_type(i) = :
!        0 , If the cell has only computing faces
!        1 , If the cell has has wall boundary faces
!        2 , If the cell has has inlet boundary faces
!        3 , If the cell has has outlet boundary faces
!        4 , If the cell has has symmetry boundary faces
!        5 , If the cell has has mcc boundary faces
!
      DO i=1,ncell_fluid
         j1=i_neigh(i)-1
         itype=0
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).lt.0) THEN
               itype=1
            ELSE
               IF(nbcon(j).gt.0.and.nbcon(j).le.nin_max)THEN
                  IF(itype.eq.0.or.itype.ge.2)itype=2
               ELSEIF(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)THEN
                  IF(itype.eq.0.or.itype.ge.3)itype=3
               ELSEIF(nbcon(j).eq.nb_sym)THEN
                  IF(itype.eq.0.or.itype.ge.4)itype=4
               ELSEIF(nbcon(j).eq.nb_mars)THEN
                  IF(itype.eq.0.or.itype.ge.5)itype=5
               ENDIF
            ENDIF
         ENDDO
!         
!........Save the face index and calculate distance from the wall when a cell has wall boundary face
!
         dr_min=huge(0.d0)
         j0=0
         DO j=i_neigh(i),i_neigh(i+1)-1
             IF((itype.eq.1 .and.  nbcon(j).lt.0) .or. &
                (itype.eq.2 .and. (nbcon(j).gt.0.and.nbcon(j).le.nin_max))     .or. &
                (itype.eq.3 .and. (nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)).or. &
                (itype.eq.4 .and.  nbcon(j).eq.nb_sym))THEN
               IF(dji(j).lt.dr_min) THEN
                  dr_min=dji(j)
                  j0=j-j1
               ENDIF
             ENDIF
         ENDDO
         icell_type(i)=itype
         iface_wall(i)=j0
      ENDDO  
!
      num_wallcells_l=0
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1)THEN
            num_wallcells_l=num_wallcells_l+1
         ENDIF
      ENDDO
!
      ALLOCATE(wall_cell_l(num_wallcells_l))
      ALLOCATE(list(num_wallcells_l))
      ALLOCATE(indexo_l(num_wallcells_l))
      k=0
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1)THEN
            k=k+1
            wall_cell_l(k)=i
         ENDIF
      ENDDO
!
!     The num_wallcells list was originally built on the global
!     If we build on the local the sequence would be different
!     Then we need a mapping scheme obtained via a search of jperm in
!     num_wallcells list
!
      CALL allgatherv_i(icell_type,icell_type_tmp,ncell_fluid,ncell_fluid_all,0)
      CALL allgatherv_i(iface_wall,iface_wall_tmp,ncell_fluid,ncell_fluid_all,0)
!
      num_wallcells=0
      DO i=1,ncell_fluid_all
         IF(icell_type_tmp(i).eq.1)THEN
            num_wallcells=num_wallcells+1
         ENDIF         
      ENDDO
!
!.....Build the wall_cell in the global
!
      ALLOCATE(wall_cell(num_wallcells))
      k=0
      DO i=1,ncell_fluid_all
         IF(icell_type_tmp(i).eq.1)THEN
            k=k+1
            wall_cell(k)=i
         ENDIF         
      ENDDO
!
!.....Build the list to be searched
!
      DO i=1,num_wallcells_l
         i0=wall_cell_l(i)
         list(i)=jperm(i0)
      ENDDO
!
      CALL searchx_list(num_wallcells,wall_cell,num_wallcells_l,list,indexo_l)
!
      IF(num_wallcells_l.gt.0) THEN
         ALLOCATE(sa_wallcell_l(num_wallcells_l))
         ALLOCATE(xfc_wallcell_l(num_wallcells_l,ndim))
         ALLOCATE(xn_wallcell_l(num_wallcells_l,ndim))
      ELSE
         ALLOCATE(sa_wallcell_l(1))
         ALLOCATE(xfc_wallcell_l(1,ndim))
         ALLOCATE(xn_wallcell_l(1,ndim))
      ENDIF
      IF(ndim.eq.2) THEN
         DO ii=1, num_wallcells_l
            i=wall_cell_l(ii)
            j=iface_wall(i)
            k=j
            j1=i_neigh(i)-1
            sa_wallcell_l(ii)   =sa(k+j1)
            xfc_wallcell_l(ii,1)=xfc(k+j1,1)
            xfc_wallcell_l(ii,2)=xfc(k+j1,2)
            xn_wallcell_l(ii,1) =xn(k+j1,1)
            xn_wallcell_l(ii,2) =xn(k+j1,2)
         ENDDO
      ELSE
         DO ii=1,num_wallcells_l
            i=wall_cell_l(ii)
            j=iface_wall(i)
            k=j
            j1=i_neigh(i)-1
            sa_wallcell_l(ii)   =sa(k+j1)
            xfc_wallcell_l(ii,1)=xfc(k+j1,1)
            xfc_wallcell_l(ii,2)=xfc(k+j1,2)
            xfc_wallcell_l(ii,3)=xfc(k+j1,3)
            xn_wallcell_l(ii,1) =xn(k+j1,1)
            xn_wallcell_l(ii,2) =xn(k+j1,2)
            xn_wallcell_l(ii,3) =xn(k+j1,3)
         ENDDO
      ENDIF
!
!.....Get size and displacement to call ALLGATHERV
!
      CALL allgather_i(num_wallcells_l,num_wallcells_p)
      ip=1
      num_wallcells_p_dsp(ip)=0
      DO ip=2,np
         num_wallcells_p_dsp(ip)=num_wallcells_p_dsp(ip-1)+num_wallcells_p(ip-1)
      ENDDO
      ALLOCATE(windexo_l(num_wallcells))
      ALLOCATE(wxfc_wallcell(num_wallcells,ndim))
      ALLOCATE(xfc_wallcell(num_wallcells,ndim))
      ALLOCATE(wxn_wallcell(num_wallcells,ndim))
      ALLOCATE(xn_wallcell(num_wallcells,ndim))
      IF(np.gt.1)THEN
         CALL allgather_vec_i(indexo_l,num_wallcells_l,windexo_l,num_wallcells, &
                              num_wallcells_p,num_wallcells_p_dsp)
      ELSE
         DO i=1,num_wallcells
            windexo_l(i)=indexo_l(i)
         ENDDO
      ENDIF
      IF(np.gt.1)THEN
         DO ix=1,ndim
            CALL allgather_vec_r(xfc_wallcell_l(1,ix),num_wallcells_l,wxfc_wallcell(1,ix),num_wallcells, &
                                 num_wallcells_p,num_wallcells_p_dsp)
            DO i=1,num_wallcells
               k=windexo_l(i)
               xfc_wallcell(k,ix)=wxfc_wallcell(i,ix) 
            ENDDO
            CALL allgather_vec_r(xn_wallcell_l(1,ix),num_wallcells_l,wxn_wallcell(1,ix),num_wallcells, &
                                 num_wallcells_p,num_wallcells_p_dsp)
            DO i=1,num_wallcells
               k=windexo_l(i)
               xn_wallcell(k,ix)=wxn_wallcell(i,ix) 
            ENDDO
         ENDDO
      ELSE
         DO ix=1,ndim
            DO i=1,num_wallcells
               wxfc_wallcell(i,ix)=xfc_wallcell_l(i,ix)
            ENDDO
            DO i=1,num_wallcells
               k=windexo_l(i)
               xfc_wallcell(k,ix)=wxfc_wallcell(i,ix) 
            ENDDO
            DO i=1,num_wallcells
               wxn_wallcell(i,ix)=xn_wallcell_l(i,ix)
            ENDDO
            DO i=1,num_wallcells
               k=windexo_l(i)
               xn_wallcell(k,ix)=wxn_wallcell(i,ix) 
            ENDDO
         ENDDO
      ENDIF
!
!.....Deallocate temporary variables
!
      IF(ALLOCATED(icell_type_tmp)) DEALLOCATE(icell_type_tmp)
      IF(ALLOCATED(iface_wall_tmp)) DEALLOCATE(iface_wall_tmp)
!
!.....Communicate icell_type,npb
!
      IF(np.gt.1) CALL communicate_1d_int(icell_type, &
                                             npb)
!
      END SUBROUTINE find_icelltype 

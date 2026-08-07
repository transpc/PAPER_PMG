!
      SUBROUTINE frink_weight(neigh_face,nd,node_face,n_face,xloc_tmp,  &
                              num_cell_node_tmp,cell_node_tmp,rwcn_tmp)
!
      USE Zmpi         , ONLY: maxmt_cell
      USE Znum_cell    , ONLY: i_neigh_tmp
      USE Zparam       , ONLY: nn,ndim
      USE Znode        , ONLY: nmax_vertex,nd_max,xnode,swcn
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n_face
      INTEGER :: neigh_face(maxmt_cell),nd(n_face),node_face(nmax_vertex,n_face)
      REAL(8) :: xloc_tmp(nn,ndim)
!.....Output
      REAL(8) :: rwcn_tmp(nd_max,nn)
!.....Local variables
      INTEGER :: i,j,ix,k,m,n,i1,ii,n1
      INTEGER :: num_cell_node_tmp(nn),cell_node_tmp(nd_max,nn)
      REAL(8) :: ds1
!
      DO i=1,nn
         n1=0
         DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
            ii=neigh_face(j)
            DO m=1,nd(ii)
               n=node_face(m,ii)
               DO i1=1,nd_max
                  IF(n.eq.cell_node_tmp(i1,i)) goto 110
               ENDDO
               n1=n1+1
               cell_node_tmp(n1,i)=n
110            CONTINUE
            ENDDO
         ENDDO
         num_cell_node_tmp(i)=n1
      ENDDO
!
!.....Inverse distance weightting method
!
      DO i=1,nn
         DO j=1,num_cell_node_tmp(i)
            k=cell_node_tmp(j,i)
            rwcn_tmp(j,i)=0.d0
            DO ix=1,ndim
               ds1=xloc_tmp(i,ix)-xnode(k,ix)
               rwcn_tmp(j,i)=rwcn_tmp(j,i)+ds1**2
            ENDDO
            rwcn_tmp(j,i)=1.0d0/DSQRT(rwcn_tmp(j,i))
            swcn(k)=swcn(k)+rwcn_tmp(j,i)
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE frink_weight
!
!
      SUBROUTINE frink_weight_lap(num_neigh_tmp,neigh_face,nd,node_face,n_face,xloc_tmp,  &
                                  num_cell_node_tmp,cell_node_tmp,rwcn_tmp)
!
      USE Zmpi         , ONLY: maxmt_cell
      USE Znum_cell    , ONLY: i_neigh_tmp
      USE Zparam       , ONLY: nn,ndim
      USE Znode        , ONLY: nmax_vertex,nd_max,xnode,n_node,swcn
!
      IMPLICIT NONE
!.....Input
      INTEGER :: n_face
      INTEGER :: num_neigh_tmp(nn),neigh_face(maxmt_cell),nd(n_face),node_face(nmax_vertex,n_face)
      REAL(8) :: xloc_tmp(nn,ndim)
!.....Output
      REAL(8) :: rwcn_tmp(nd_max,nn)
!.....Local variables
      INTEGER :: i,j,ix,k,m,n,i1,i3,j0
      INTEGER :: num_cell_node_tmp(nn),cell_node_tmp(nd_max,nn)
      REAL(8) :: ds1,ds2,c11,c22,c33,c12,c13,c14,den
!
      INTEGER, ALLOCATABLE::lxy(:,:)
      REAL(8), ALLOCATABLE::rxx(:,:),ixx(:,:),ixy(:,:),lam(:,:)
!
      DO i=1,nn
         j0=i_neigh_tmp(i)-1
         DO j=1,num_neigh_tmp(i)
            k=neigh_face(j+j0)
            DO m=1,nd(k)
               n=node_face(m,k)
               i3=0
               DO i1=1,nd_max
                  IF(n.eq.cell_node_tmp(i1,i)) i3=1
               ENDDO
               IF(i3.eq.0) THEN
                  num_cell_node_tmp(i)=num_cell_node_tmp(i)+1
                  cell_node_tmp(num_cell_node_tmp(i),i)=n
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      ALLOCATE(rxx(n_node,ndim),ixx(n_node,ndim),ixy(n_node,ndim),lam(n_node,ndim),lxy(2,ndim))
      rxx(:,:)=0.d0
      ixx(:,:)=0.d0
      ixy(:,:)=0.d0
      lxy(:,:)=0
!
      DO ix=1,ndim
         lxy(1,ix)=ix
         lxy(2,ix)=ix+1
         IF(ndim.eq.3.and.ix.eq.ndim)lxy(2,ix)=1
      ENDDO
!
!.....Obtain geometric centroid and moments
!
      IF(ndim.eq.2)then
         DO i=1,nn
            DO j=1,num_cell_node_tmp(i)
               k=cell_node_tmp(j,i)
               DO ix=1,ndim
                  ds1=xloc_tmp(i,ix)-xnode(k,ix)
                  rxx(k,ix)=rxx(k,ix)+ds1
                  ixx(k,ix)=ixx(k,ix)+ds1**2
               ENDDO
               ds1=xloc_tmp(i,1)-xnode(k,1)
               ds2=xloc_tmp(i,2)-xnode(k,2)
               ixy(k,1)=ixy(k,1)+ds1*ds2
            ENDDO
         ENDDO
      ELSEIF(ndim.eq.3)then
         DO i=1,nn
            DO j=1,num_cell_node_tmp(i)
               k=cell_node_tmp(j,i)
               DO ix=1,ndim
                  ds1=xloc_tmp(i,ix)-xnode(k,ix)
                  rxx(k,ix)=rxx(k,ix)+ds1
                  ixx(k,ix)=ixx(k,ix)+ds1**2
               ENDDO
               DO ix=1,ndim
                  ds1=xloc_tmp(i,lxy(1,ix))-xnode(k,lxy(1,ix))
                  ds2=xloc_tmp(i,lxy(2,ix))-xnode(k,lxy(2,ix))
                  ixy(k,ix)=ixy(k,ix)+ds1*ds2
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!.....Obtain Lagrangian multiplier
!
      IF(ndim.eq.2)then
         DO i=1,n_node
            den=ixx(i,1)*ixx(i,2)-ixy(i,1)*ixy(i,1)
            IF(dabs(den).lt.1.d-15)then
               WRITE(*,*)'++ Warning ++small denominator in obtaining',' Lagrangian multiplier!!',i,den
               WRITE(*,*)ixx(i,1:2),ixy(i,1)
            ELSE
               den=1.d0/den
            ENDIF
            lam(i,1)=den*(ixy(i,1)*rxx(i,2)-ixx(i,2)*rxx(i,1))
            lam(i,2)=den*(ixy(i,1)*rxx(i,1)-ixx(i,1)*rxx(i,2))
         ENDDO
      ELSEIF(ndim.eq.3)then
         DO i=1,n_node
            c11=ixx(i,1)*ixx(i,2)-ixy(i,1)*ixy(i,1)
            c22=ixx(i,2)*ixx(i,3)-ixy(i,2)*ixy(i,2)
            c33=ixx(i,3)*ixx(i,1)-ixy(i,3)*ixy(i,3)
!
            c12=ixy(i,1)*ixy(i,2)-ixx(i,2)*ixy(i,3)
            c13=ixy(i,1)*ixx(i,3)-ixy(i,3)*ixy(i,2)
            c14=ixx(i,1)*ixy(2,i)-ixy(i,1)*ixy(i,3)
!
            den=(ixx(i,1)*c22-ixy(i,1)*c13-ixy(i,3)*c12)
            IF(dabs(den).lt.1.d-15)then
               WRITE(*,*)'++ Warning ++small denominator in obtaining',' Lagrangian multiplier!!',i,den
            ELSE
               den=1.d0/den
            ENDIF
            lam(i,1)=den*(-rxx(i,1)*c22+rxx(i,2)*c13-rxx(i,3)*c12)
            lam(i,2)=den*(+rxx(i,1)*c13-rxx(i,2)*c33+rxx(i,3)*c14)
            lam(i,3)=den*(-rxx(i,1)*c12+rxx(i,2)*c14-rxx(i,3)*c11)
         ENDDO
      ENDIF
!
!.....Obtain geometric weights based on Lagrangian multiplier
!
      DO i=1,nn
         DO j=1,num_cell_node_tmp(i)
            k=cell_node_tmp(j,i)
            rwcn_tmp(j,i)=1.d0
            DO ix=1,ndim
               rwcn_tmp(j,i)=rwcn_tmp(j,i)+lam(k,ix)*(xloc_tmp(i,ix)-xnode(k,ix))
            ENDDO
!
!.....Test for the inverse distance weightting method
!
            IF(rwcn_tmp(j,i).lt.1.d-10.or.rwcn_tmp(j,i).gt.2.d0)then
               !--> return to inverse-distance weighting
               rwcn_tmp(j,i)=0.d0
               DO ix=1,ndim
                  ds1=xloc_tmp(i,ix)-xnode(k,ix)
                  rwcn_tmp(j,i)=rwcn_tmp(j,i)+ds1**2
               ENDDO
               rwcn_tmp(j,i)=1.0/DSQRT(rwcn_tmp(j,i))
            ENDIF
            swcn(k)=swcn(k)+rwcn_tmp(j,i)
         ENDDO
      ENDDO
!
      DEALLOCATE(rxx,ixx,ixy,lam,lxy)
!
      RETURN
      END SUBROUTINE frink_weight_lap

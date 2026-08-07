!
      SUBROUTINE outlet_horizontal_user
!
!     Treatment of gravity in horizontal outlet boundaries having several outlet cells 
!     (please, same region as same npin)
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Znum_cell    , ONLY: i_neigh
      USE Zb_condition , ONLY: pbnd
      USE Zbc_index    , ONLY: npb,npin,nbcon
      USE Zconst2      , ONLY: grav
      USE Zcoord1      , ONLY: xloc
      USE Zpress       , ONLY: p
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER npb_loop,idex,i,j
      REAL(8) outlet_sum,outlet_mid
!      
!.....Local arrays
      REAL(8) elev(1000)
!
!.....Set the pressures on P.B. to INITIAL values
!
      DO i=1,ncell_fluid
        IF(npb(i).gt.0) p(i)=pbnd(npb(i))
      ENDDO
!
!.....Set the pressures on P.B. considering gravity 
!
      DO npb_loop=5,5+npin-1
         idex=0
         elev(:)=-9.8d0*(-1000.d0) !-1000m
         outlet_sum=0.0d0
         outlet_mid=0.0d0            
!
         DO i=1,ncell_fluid
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).ne.npb_loop)CYCLE
               idex=idex+1
               elev(idex)=dot_product(xloc(i,:),grav(:))
               outlet_sum=outlet_sum+elev(idex)
            ENDDO
         ENDDO
!
         IF(idex.gt.0) outlet_mid=outlet_sum/idex  
         IF(np.gt.1) CALL allreducei_r1(outlet_mid)
         idex=0
!
         DO i=1,ncell_fluid
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).ne.npb_loop)CYCLE
               idex=idex+1
               p(i)=p(i)+(elev(idex)-outlet_mid)*cell%rhom(i)
            ENDDO
         ENDDO
!         
      ENDDO
!      
      END SUBROUTINE outlet_horizontal_user 

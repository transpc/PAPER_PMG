!
      SUBROUTINE natural_convection_user
!
!     Temperature boundary condition
!     F. Ampofo, T.G. Karayiannis, 2003, Experimental benchmark data for turbulent natural convection in an air filled square cavity,
!     International Journal of Heat and Mass Transfer, 46, 3551-3572.
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ns
      USE Znum_cell    , ONLY: i_neigh
      USE Zbc_index    , ONLY: nbcon
      USE Znormal      , ONLY: wall_cell_l,num_wallcells_l
      USE Zcoord1      , ONLY: xloc
      USE Zvec_geo     , ONLY: xn_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,j,k,j0
      INTEGER ii
      LOGICAL, SAVE :: INITIAL=.true.
      REAL(8) dx,dy_twall,dy_bwall
!.....Local arrays
      REAL(8):: yn(ns)
      REAL(8) xpos(21),Ttopwall(21),Tbotwall(21)
!
      DATA xpos /0, 0.002, 0.00667, 0.0133, 0.0267, 0.0533, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, &
                 0.7, 0.8, 0.9, 0.9467, 0.9733, 0.9867, 0.9933, 0.998, 1/ 
      DATA Ttopwall /1, 0.949, 0.9333, 0.9342, 0.9183, 0.8782, 0.821, 0.7597, 0.7107, 0.6779, 0.6393, &
                     0.6135, 0.5578, 0.488, 0.3372, 0.2338, 0.1409, 0.1334, 0.1234, 0.0967, 0/     
      DATA Tbotwall/1, 0.9184, 0.9038, 0.8862, 0.8639, 0.7733, 0.6608, 0.5263, 0.452, 0.396, 0.3503, &
                    0.3116, 0.2722, 0.2214, 0.149, 0.0938, 0.0445, 0.0352, 0.0272, 0.0185, 0/
!
      IF(initial) THEN
!      
!........Dimensionalize      
!
         xpos=xpos*0.75d0
         Ttopwall=Ttopwall*40.d0+283.15d0
         Tbotwall=Tbotwall*40.d0+283.15d0
!      
!........Initialize
!
         DO i=1,ncell_fluid
            cell%T_top(i)=0.d0
            cell%T_bot(i)=0.d0
         ENDDO    
!         
         DO ii=1, num_wallcells_l
            i=wall_cell_l(ii)
!...........Get all the xn(j,2) for cell m1
            CALL get_scalar_variable_n_i_ndim(xn_nf,yn,i,2)
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).ne.-1) cycle
!            
!..............Interpolate Temperature at xloc position from the paper position(xpos)
!
               DO k=1,21-1
                  IF(xloc(i,1).ge.xpos(k).and.xloc(i,1).lt.xpos(k+1)) THEN
                     dx=xpos(k+1)-xpos(k)
                     IF(yn(j-j0).gt.0.5d0) THEN
!                    
!......................Top face: only defined on cell variable but it's meaning is face value.
!
                        dy_twall=Ttopwall(k+1)-Ttopwall(k)
                        cell%T_top(i)=Ttopwall(k)+dy_twall/dx*(xloc(i,1)-xpos(k))
                     ELSE
!
!......................Bottom face: only defined on cell variable but it's meaning is face value.
!
                       dy_bwall=Tbotwall(k+1)-Tbotwall(k)
                       cell%T_bot(i)=Tbotwall(k)+dy_bwall/dx*(xloc(i,1)-xpos(k)) 
                     ENDIF
                     EXIT 
                  ENDIF
               ENDDO 
               EXIT 
            ENDDO   
         ENDDO    
         initial=.false.
      ENDIF  
!     
      END SUBROUTINE natural_convection_user
!

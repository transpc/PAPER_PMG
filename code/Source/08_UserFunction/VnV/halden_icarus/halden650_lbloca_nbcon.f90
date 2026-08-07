!
      SUBROUTINE halden650_lbloca_nbcon(iflag) !pik-halden
!
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: myrank
      USE Zparam          , ONLY: nin_max,nb_max
      USE Zbc_index       , ONLY: nbcon,npb
      USE Znum_cell       , ONLY: i_neigh
      USE Zvec_major      , ONLY: flux_g_nf,flux_l_nf,flux_d_nf
!            
      IMPLICIT NONE
!      
!.....Input
      INTEGER iflag
!.....Local variables
      INTEGER i1
      INTEGER i,j,j0
!
      IF(iflag.eq.-1)THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'just close outlet'
!      
!........close 1 outlet    
         CALL nbcon_change_start    
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)THEN
                  CALL get_vector_disp(j-j0,i,i1)
                  i1=abs(i1)
                  nbcon(j)=-1
                  npb(i)=0
                  flux_g_nf(i1)=0.d0
                  flux_l_nf(i1)=0.d0
                  flux_d_nf(i1)=0.d0
               ENDIF   
            ENDDO
         ENDDO
         CALL halden650_lbloca_clean_bc_user
         CALL nbcon_change_end    
!
      ELSE
!
         IF(myrank.eq.0)WRITE(*,*)'iflag of halden650_lbloca_nbcon should be -1 !!!'
         PAUSE
         STOP
!         
      ENDIF   
!         
      END SUBROUTINE halden650_lbloca_nbcon     

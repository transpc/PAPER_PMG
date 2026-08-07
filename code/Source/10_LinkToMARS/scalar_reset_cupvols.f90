!
      SUBROUTINE scalar_reset_cupvols 
!
!     This routine resets scalar values to that of previous time step
!
      USE Zmars      , ONLY: n_marsbc 
!      
      IMPLICIT NONE
!      
      INCLUDE 'c3com.h' 
!
      INTEGER i,j
!      
!.....MARS interface 
!
      DO i=1,n_marsbc
            c3vg(1,i)=    c3vg_o(1,i)   
            c3vl(1,i)=    c3vl_o(1,i)   
         c3alphg(1,i)= c3alphg_o(1,i)   
         c3alphf(1,i)= c3alphf_o(1,i)   
         c3betag(1,i)= c3betag_o(1,i)   
         c3betaf(1,i)= c3betaf_o(1,i)   
            c3xi(1,i)=    c3xi_o(1,i)   
          c3area(1,i)=  c3area_o(1,i)   
          c3delp(1,i)=  c3delp_o(1,i)   
            c3pa(1,i)=    c3pa_o(1,i)   
            c3ug(1,i)=    c3ug_o(1,i)   
            c3uf(1,i)=    c3uf_o(1,i)   
            c3al(1,i)=    c3al_o(1,i)   
          c3arxq(1,i)=  c3arxq_o(1,i)   
          c3rhof(1,i)=  c3rhof_o(1,i)   
          c3rhog(1,i)=  c3rhog_o(1,i)   
           c3brn(1,i)=   c3brn_o(1,i)   
         c3vpgno(1,i)= c3vpgno_o(1,i)   
      ENDDO           
      DO i=1,n_marsbc
        DO j=1,n_marsbc
           c3yeta(1,i,j)=c3yeta_o(1,i,j)
        ENDDO 
      ENDDO         
      DO i=1,n_marsbc
        DO j=1,15
           c3rtp(1,i,j)=c3rtp_o(1,i,j)
        ENDDO 
      ENDDO           
!      
!DEC$IF defined (SPACE)
      CALL c3com_copy_C2S
!DEC$ENDIF      
!
      RETURN
      END SUBROUTINE scalar_reset_cupvols      
!--------------------------------------------------------------
      SUBROUTINE shift_solutions_cupvols
!
!     This routine shift solutions fo time-marching.
!
      USE Zmars      , ONLY: n_marsbc 
!      
      IMPLICIT NONE
!      
      INCLUDE 'c3com.h' 
!
      INTEGER i,j
!      
!.....MARS interface 
!
      DO i=1,n_marsbc
         c3vg_o(1,i)=    c3vg(1,i) 
         c3vl_o(1,i)=    c3vl(1,i) 
      c3alphg_o(1,i)= c3alphg(1,i) 
      c3alphf_o(1,i)= c3alphf(1,i) 
      c3betag_o(1,i)= c3betag(1,i) 
      c3betaf_o(1,i)= c3betaf(1,i) 
         c3xi_o(1,i)=    c3xi(1,i) 
       c3area_o(1,i)=  c3area(1,i) 
       c3delp_o(1,i)=  c3delp(1,i) 
         c3pa_o(1,i)=    c3pa(1,i) 
         c3ug_o(1,i)=    c3ug(1,i) 
         c3uf_o(1,i)=    c3uf(1,i) 
         c3al_o(1,i)=    c3al(1,i) 
       c3arxq_o(1,i)=  c3arxq(1,i) 
       c3rhof_o(1,i)=  c3rhof(1,i) 
       c3rhog_o(1,i)=  c3rhog(1,i) 
        c3brn_o(1,i)=   c3brn(1,i) 
      c3vpgno_o(1,i)= c3vpgno(1,i) 
      ENDDO           
      DO i=1,n_marsbc
        DO j=1,n_marsbc
           c3yeta_o(1,i,j)=c3yeta(1,i,j)
        ENDDO 
      ENDDO         
      DO i=1,n_marsbc
        DO j=1,15
           c3rtp_o(1,i,j)=c3rtp(1,i,j)
        ENDDO 
      ENDDO           
!
      RETURN
      END SUBROUTINE shift_solutions_cupvols

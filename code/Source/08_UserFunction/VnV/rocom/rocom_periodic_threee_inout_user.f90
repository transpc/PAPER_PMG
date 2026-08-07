!
      SUBROUTINE rocom_periodic_threee_inout_user
!
!     Define inlet boundaries to the ROCOM problem
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: myrank      
      USE Zparam       , ONLY: nb_max
      USE Znum_cell    , ONLY: i_neigh
      USE Zbc_index    , ONLY: nbcon
      USE Zb_condition , ONLY: eb_liq,rhob_liq,rhob_gas,alphab_liq,alphab_gas,cb_pl,cb_pg,cb_pd
      USE Zboron       , ONLY: cboronb_liq
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i,j,k
      LOGICAL,SAVE :: INITIAL=.TRUE.
      REAL(8) :: cvm,c_vm,arhob_gas,arhob_liq,arhob_drp
      REAL(8) :: rhob_m,rhob_vm,denom  
!.....Local arrays
      INTEGER :: number(nb_max)
      REAL(8) :: rholeg_liq(nb_max),eleg_liq(nb_max),cbeg_liq(nb_max)
!      
      IF(INITIAL)THEN
         IF(myrank.eq.0)THEN
            WRITE(*,*) '          Periodic boundary are USEd now!'
            WRITE(unit_log,*)'          Periodic boundary are USEd now!'
         ENDIF   
         INITIAL=.FALSE.
      ENDIF  
!
!.....Define inlet conditions
!
      DO i=1,nb_max
         number(i)    =0
         rholeg_liq(i)=0.d0
         eleg_liq(i)  =0.d0
         cbeg_liq(i)  =0.d0
      ENDDO
!
      DO i=1,ncell_fluid
         DO j=i_neigh(i),i_neigh(i+1)-1
            k=nbcon(j)
            IF(k.ge.1 .and. k.le.nb_max)THEN
               rholeg_liq(k)=rholeg_liq(k)+cell%rhol(i)
               eleg_liq(k)  =eleg_liq(k)  +cell%el(i)
               cbeg_liq(k)  =cbeg_liq(k)  +cell%cboron(i)
               number(k)=number(k)+1
            ENDIF
         ENDDO
      ENDDO
!                             
      DO i=1,nb_max
         IF(number(i).ge.1) rholeg_liq(i)=rholeg_liq(i)/number(i)
         IF(number(i).ge.1) eleg_liq(i)  =eleg_liq(i)  /number(i)
         IF(number(i).ge.1) cbeg_liq(i)  =cbeg_liq(i)  /number(i)         
      ENDDO         
!
!.....Inlet #2,3,4
!
      DO i=2,4
         eb_liq(i)     =eleg_liq(i+4)
         rhob_liq(i)   =rholeg_liq(i+4)
         cboronb_liq(i)=cbeg_liq(i+4)
      ENDDO
!
!.....Miscellaneous coefficients
!
      DO i=2,4
         cvm=c_vm(alphab_gas(i))
         arhob_liq=alphab_liq(i)*rhob_liq(i)
         arhob_liq=(1.d0-alphab_gas(i))*rhob_liq(i)
         arhob_drp=0.0d0
         arhob_gas=alphab_gas(i)*rhob_gas(i)
         rhob_m   =arhob_liq+arhob_gas
         rhob_vm  =cvm*(arhob_liq+arhob_gas)
!         
         denom=rhob_gas(i)*rhob_liq(i)+rhob_m*rhob_vm
         IF(denom.eq.0.d0) denom=1.d0
!         
         cb_pl(i)=(rhob_gas(i)+rhob_vm)/denom
         cb_pg(i)=(rhob_liq(i)+rhob_vm)/denom                !JJJ, vg=vd at boundary
         cb_pd(i)=cb_pg(i)
!            
      ENDDO
!                   
      END SUBROUTINE rocom_periodic_threee_inout_user

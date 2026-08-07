!
    SUBROUTINE v1d3d_cupid 
!     
!     This routine transfers the dP of Cell C(i), and 
!     calculates the phasic velocities at the 1D/3D interfaces      
!     by backsubstitution 
! 
      USE Zmars, only:ppcup_tmp
      USE c3com_cupid, only:i3cupid
!
      USE Zmars   , only:marsindex,n_marsbc,pcup,egcup,elcup,&
                         alphagcup,rholcup,rhogcup,cboroncup,qualacup,&
                         pcup_tmp,egcup_tmp,elcup_tmp,alphagcup_tmp,&
                         rhogcup_tmp,rholcup_tmp,cboroncup_tmp,qualacup_tmp,&
                         tlcup,tgcup,tlcup_tmp,tgcup_tmp      
      USE Zcore   , only:np,myrank
      USE Ztimecon, ONLY:time                
!
      IMPLICIT none
!
      INCLUDE 'c3com.h' 
!
!DEC$IF defined (SPACE)
      INCLUDE 'c3com_space.h'
!DEC$ENDIF           
!                                                                 
 
      INTEGER i,j,k,ic
!      
      REAL(8) c3dv
!
      INTEGER,SAVE:: loop
      DATA loop/0/      
!
!...Transfer the dP of Cell C(i) 
!
      IF(i3nic(2).eq.0)RETURN
      DO j=1,i3nic(2)
         ic=i3cupid(j)
         c3delp(1,j)=ppcup_tmp(j) !!mcc-pik-mpi-2013-06-27 pp(ic) ! from cupid into mars
      ENDDO
!
!...Calculate the phasic velocity at the interface
!
      i=1
      DO j=1,i3nic(2) 
	     c3dv=0.0d0 
         DO k=1,i3nic(2) 
            ic=i3cupid(k)
            c3dv=c3dv+c3yeta(i,j,k)*c3delp(1,k)
         ENDDO 
!                                                                    
!        c3vg(i,j)=c3alphg(i,j)+c3betag(i,j)*(pp(ic)-c3xi(i,j)-c3dv)                                                           
!        c3vl(i,j)=c3alphf(i,j)+c3betaf(i,j)*(pp(ic)-c3xi(i,j)-c3dv)                                                           
         c3vg(i,j)=c3alphg(i,j)+c3betag(i,j)*(c3delp(1,j)-c3xi(i,j)-c3dv)                                                           
         c3vl(i,j)=c3alphf(i,j)+c3betaf(i,j)*(c3delp(1,j)-c3xi(i,j)-c3dv)                                                           
      ENDDO 
! 
!DEC$IF defined (SPACE)
      s3delp(:) =c3delp(1,:)           
      s3vg(:)   =c3vg(1,:)                
      s3vl(:)   =c3vl(1,:)   
!DEC$ENDIF           
!                                                                 
      RETURN 
      END SUBROUTINE v1d3d_cupid                          

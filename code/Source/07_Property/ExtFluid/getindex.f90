      SUBROUTINE getindex_cupid(tables,z,pointer,lower,upper,index,err) 
!define win32dvf                                                        
!define erf                                                             
!define fourbyt                                                         
!define hconden                                                         
!define impnon                                                          
!define in32                                                            
!define newnrc                                                          
!define ploc                                                            
!define sphaccm                                                         
!define unix                                                            
!define noselap                                                         
!define noextvol                                                        
!define noextv20                                                        
!define noextsys                                                        
!define noextjun                                                        
!define noextj20                                                        
!define noparcs                                                         
!define nonpa                                                           
!define nomap                                                           
!define logp                                                            
!deck getindex                                                          
!                                                                       
!  $Id: getindex.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $        
!                                                                       
!  Find index in table 1 or 2 of the nearest (lesser) temperature       
!  or pressure than the input temperature or pressure, z                
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
!                                                                       
      REAL(8) z,tables(1) 
      INTEGER i,lower,upper,index,pointer 
      LOGICAL err 
!                                                                       
!      call timstart ('getindex')                                       
!                                                                       
      err=.false. 
!                                                                       
!  Test to make sure that the lower bound is less than or               
!  equal to the upper bound                                             
!                                                                       
      IF(upper.lt.lower)then 
         err=.true. 
         GOTO 999 
      ENDIF 
!                                                                       
!  Test to make sure that the input temperature or pressure is          
!  within the requested range                                           
!                                                                       
      IF((z.lt.tables(pointer+lower)).or.(z.gt.tables(pointer+upper)))  &
      then                                                              
         err=.true. 
         GOTO 999 
      ENDIF 
!                                                                       
      DO 10 i=lower+1,upper,1 
         IF(z.lt.tables(pointer+i))then 
            index=i-1 
            GOTO 20 
         ENDIF 
   10 END DO 
!                                                                       
!  z is equal to maximum value in table                                 
!                                                                       
      index=upper-1 
!                                                                       
   20 CONTINUE 
!                                                                       
  999 CONTINUE 
!                                                                       
!      call timstop ('getindex')                                        
!                                                                       
      RETURN 
      END SUBROUTINE getindex_cupid                       

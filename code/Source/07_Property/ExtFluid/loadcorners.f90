      SUBROUTINE loadcorners_cupid(tables,a,itleft,itright,iplo,ptableprop,   &
      stableprop,ntableprop,ofirstpres,ptr)                             
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
!deck loadcorners                                                       
!                                                                       
!  $Id: loadcorners.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $     
!                                                                       
!  Loads the information needed at the four corner points for           
!  the four-point (2D) Hermite interpolation formula for                
!  liquid                                                               
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
      INCLUDE 'magic.h' 
!                                                                       
      REAL(8) a(11,4),tables(1) 
!      real f(36)                                                       
      INTEGER p1,p2,t1,t2 
      INTEGER c(4),iplo,itleft,ntableprop,ofirstpres,ptableprop,        &
      stableprop,ptr,itright                                            
!if def,swesty                                                          
!      integer i,j                                                      
!endif                                                                  
!                                                                       
!      call timstart ('loadcorners')                                    
!                                                                       
!  Load pressure and temperature indices                                
!                                                                       
      t1=ptable1+itleft 
      t2=ptable1+itright 
      p1=ptable2+iplo 
      p2=ptable2+iplo+1 
!                                                                       
!   Load corner point info for the 4-point formula                      
!   c(1) is lower left hand corner                                      
!   c(2) is lower right hand corner                                     
!   c(3) is upper left hand corner                                      
!   c(4) is upper right hand corner                                     
!                                                                       
      c(1)=ptableprop+(iplo-ofirstpres)*stableprop+(itleft-1)*          &
      ntableprop                                                        
!      ptr = c(1)+11                                                    
      ptr=c(1)+1 
!if def,old                                                             
!      c(2) = c(1) + ntableprop                                         
!      c(3) = c(1) + stableprop                                         
!      c(4) = c(3) + ntableprop                                         
!endif                                                                  
!                                                                       
!  Load up the a array for the 4-point formula                          
!                                                                       
      a(1,1)=tables(t1) 
      a(1,2)=tables(t2) 
      a(1,3)=tables(t1) 
      a(1,4)=tables(t2) 
      a(2,1)=tables(p1) 
      a(2,2)=tables(p1) 
      a(2,3)=tables(p2) 
      a(2,4)=tables(p2) 
!                                                                       
!if def,swesty                                                          
!      do 100 i = 2,10                                                  
!        do 200 j = 1,4                                                 
!          a(i+1,j) = tables(c(j) + i)                                  
!200     continue                                                       
!100   continue                                                         
!endif                                                                  
!                                                                       
!        do 10 i = 1,36                                                 
!          f(i) = tables(c(1)+11+i)                                     
!          if (f(i) .eq. magic) stop                                    
!10      continue                                                       
!                                                                       
!      call timstop ('loadcorners')                                     
!                                                                       
      RETURN 
      END SUBROUTINE loadcorners_cupid                    

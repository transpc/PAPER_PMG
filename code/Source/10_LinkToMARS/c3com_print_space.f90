!
      SUBROUTINE c3com_print_space
!
!     This routine prints c3com data
!
      USE Zcore,      ONLY: myrank
      USE Ztimecon,   ONLY: itim
      USE c3com_cupid
!
      IMPLICIT none
!      
      INCLUDE 'c3com.h' 
      INCLUDE 'c3com_space.h' 
!
      LOGICAL,SAVE::initial
      DATA initial/.true./
!      
      IF(itim.gt.10)RETURN
      IF(myrank.ne.0)RETURN  
      IF(initial)THEN
         initial=.FALSE.
         OPEN(201,file='c3com_check_space.dat')        
      ENDIF   
   
      WRITE(201,15) ' ' 
      WRITE(201,"(5x,a,1i5)") ' c3com variable check : itim   = ',itim   
      WRITE(201,"(5x,a,1i5)") ' c3com variable check : iwhere = ',i_where   
      WRITE(201,15) 's3dt_super          =',s3dt_super       ,'dt_super            =',dt_super       
      WRITE(201,15) 's3dt_relap          =',s3dt_relap       ,'dt_relap            =',dt_relap       
WRITE(201,15) 's3dt_cobra          =',s3dt_cobra       ,'dt_cobra            =',dt_cobra       
      WRITE(201,15) 's3vg(1)             =',s3vg(1)          ,'c3vg(1)             =',c3vg(1,1)          
      WRITE(201,15) 's3vl(1)             =',s3vl(1)          ,'c3vl(1)             =',c3vl(1,1)          
      WRITE(201,15) 's3alphg(1)          =',s3alphg(1)       ,'c3alphg(1)          =',c3alphg(1,1)       
      WRITE(201,15) 's3alphf(1)          =',s3alphf(1)       ,'c3alphf(1)          =',c3alphf(1,1)       
      WRITE(201,15) 's3betag(1)          =',s3betag(1)       ,'c3betag(1)          =',c3betag(1,1)      
      WRITE(201,15) 's3betaf(1)          =',s3betaf(1)       ,'c3betaf(1)          =',c3betaf(1,1)      
      WRITE(201,15) 's3xi(1)             =',s3xi(1)          ,'c3xi(1)             =',c3xi(1,1)         
      WRITE(201,15) 's3yeta(1,1)         =',s3yeta(1,1)      ,'c3yeta(1,1)         =',c3yeta(1,1,1)    
      WRITE(201,15) 's3area(1)           =',s3area(1)        ,'c3area(1)           =',c3area(1,1)       
      WRITE(201,15) 's3rtp(1,1)          =',s3rtp(1,1)      ,'c3rtp(1,1)          =',c3rtp(1,1,1)     
      WRITE(201,15) 's3rtp(1,2)          =',s3rtp(1,2)      ,'c3rtp(1,2)          =',c3rtp(1,1,2)
      WRITE(201,15) 's3rtp(1,3)          =',s3rtp(1,3)      ,'c3rtp(1,3)          =',c3rtp(1,1,3)
      WRITE(201,15) 's3rtp(1,4)          =',s3rtp(1,4)      ,'c3rtp(1,4)          =',c3rtp(1,1,4)
      WRITE(201,15) 's3rtp(1,5)          =',s3rtp(1,5)      ,'c3rtp(1,5)          =',c3rtp(1,1,5)
      WRITE(201,15) 's3rtp(1,6)          =',s3rtp(1,6)      ,'c3rtp(1,6)          =',c3rtp(1,1,6)
      WRITE(201,15) 's3rtp(1,7)          =',s3rtp(1,7)      ,'c3rtp(1,7)          =',c3rtp(1,1,7)
      WRITE(201,15) 's3rtp(1,8)          =',s3rtp(1,8)      ,'c3rtp(1,8)          =',c3rtp(1,1,8)
      WRITE(201,15) 's3rtp(1,9)          =',s3rtp(1,9)      ,'c3rtp(1,9)          =',c3rtp(1,1,9)
      WRITE(201,15) 's3rtp(1,10)         =',s3rtp(1,10)      ,'c3rtp(1,10)         =',c3rtp(1,1,10)
      WRITE(201,15) 's3rtp(1,15)         =',s3rtp(1,15)      ,'c3rtp(1,15)         =',c3rtp(1,1,15)
      WRITE(201,15) 's3delp(1)           =',s3delp(1)        ,'c3delp(1)           =',c3delp(1,1)       
WRITE(201,15) 's3brn(1)            =',s3brn(1)         ,'c3brn(1)            =',c3brn(1,1)        
WRITE(201,15) 's3vpgno(1)          =',s3vpgno(1)       ,'c3vpgno(1)          =',c3vpgno(1,1)      
      WRITE(201,15) 's3dpmt(1,6)         =',s3dpmt(1,6)      ,'c3dpmt(1,6)         =',c3dpmt(1,1,6)     
WRITE(201,15) 's3pa(1)             =',s3pa(1)          ,'c3pa(1)             =',c3pa(1,1)         
WRITE(201,15) 's3uf(1)             =',s3uf(1)          ,'c3uf(1)             =',c3uf(1,1)         
WRITE(201,15) 's3ug(1)             =',s3ug(1)          ,'c3ug(1)             =',c3ug(1,1)         
WRITE(201,15) 's3al(1)             =',s3al(1)          ,'c3al(1)             =',c3al(1,1)         
WRITE(201,15) 's3arxq(1)           =',s3arxq(1)        ,'c3arxq(1)           =',c3arxq(1,1)       
      WRITE(201,15) 's3odr(6)            =',s3odr(6)         ,'c3odr(6)            =',c3odr(6)         
      WRITE(201,15) 's3dpm(6)            =',s3dpm(6)         ,'c3dpm(6)            =',c3dpm(6)         
WRITE(201,15) 's3aloold(100,100)   =',s3aloold(100,100),'aloold(100,100)     =',aloold(100,100)     !c3aloold
      WRITE(201,15) 's3packfactor        =',s3packfactor     ,'packfactor          =',packfactor     
WRITE(201,15) 's3rhof(1)           =',s3rhof(1)        ,'c3rhof(1)           =',c3rhof(1,1)       
WRITE(201,15) 's3rhog(1)           =',s3rhog(1)        ,'c3rhog(1)           =',c3rhog(1,1)       

WRITE(201,15) 's3ent(1)            =',s3ent(1)         ,'r3ent(1)            =',r3ent(1,1)        !c3ent
WRITE(201,15) 's3liq(1)            =',s3liq(1)         ,'r3liq(1)            =',r3liq(1,1)        !c3liq 
       
      WRITE(201,15) 'ss3rtp(1,15)        =',ss3rtp(1,15)     ,'cc3rtp(1,15)        =',cc3rtp(1,1,15)    
WRITE(201,15) 's3vpp(1)            =',s3vpp(1)         ,'c3vpp(1)            =',c3vpp(1)        
WRITE(201,15) 's3ngpp(1)           =',s3ngpp(1)        ,'c3ngpp(1)           =',c3ngpp(1)       
WRITE(201,15) 's3vt(1)             =',s3vt(1)          ,'c3vt(1)             =',c3vt(1)         
WRITE(201,15) 's3lt(1)             =',s3lt(1)          ,'c3lt(1)             =',c3lt(1)         
WRITE(201,15) 's3dt(1)             =',s3dt(1)          ,'c3dt(1)             =',c3dt(1)         
WRITE(201,15) 's3ngmf(1,10)        =',s3ngmf(1,10)     ,'c3ngmf(1,10)        =',c3ngmf(1,10)    
	    
      WRITE(201,20)'j1cupid(1)          =',j1cupid(1)          ,'i3cupid(1)          =',i3cupid(1)          
      WRITE(201,20)'j_where             =',j_where             ,'i_where             =',i_where              
WRITE(201,20)'jflag_relap         =',jflag_relap         ,'jflag_relap         =',jflag_relap          
WRITE(201,20)'jflag_cobra         =',jflag_cobra         ,'jflag_cobra         =',jflag_cobra          
WRITE(201,20)'jflag_stop          =',jflag_stop          ,'jflag_stop          =',jflag_stop           
WRITE(201,20)'jflag_bd            =',jflag_bd            ,'jflag_bd            =',jflag_bd             
WRITE(201,20)'jflag_wp            =',jflag_wp            ,'jflag_wp            =',jflag_wp             
WRITE(201,20)'jflag_bd2           =',jflag_bd2           ,'jflag_bd2           =',jflag_bd2            
      WRITE(201,20)'j3bcn(1)            =',j3bcn(1)            ,'i3bcn(1)            =',i3bcn(1,1)            
      WRITE(201,20)'j1Cvoln(1)          =',j1Cvoln(1)          ,'i1Cvoln(1)          =',i1Cvoln(1,1)          
      WRITE(201,20)'j3chan(1)           =',j3chan(1)           ,'i3chan(1)           =',i3chan(1,1)           
      WRITE(201,20)'j3cell(1)           =',j3cell(1)           ,'i3cell(1)           =',i3cell(1,1)           
      WRITE(201,20)'j3mode(1)           =',j3mode(1)           ,'i3mode(1)           =',i3mode(1,1)           
      WRITE(201,20)'j3modet(1,3)        =',j3modet(1,3)        ,'i3modet(1,3)        =',i3modet(1,1,3)        
      WRITE(201,20)'j3nodr              =',j3nodr              ,'i3nodr              =',i3nodr               
      WRITE(201,20)'j3nic(2)            =',j3nic(2)            ,'i3nic(2)            =',i3nic(2)             
      WRITE(201,20)'j3line(30)          =',j3line(30)          ,'i3line(30)          =',i3line(30)           
      WRITE(201,20)'j3dir(1)            =',j3dir(1)            ,'i3dir(1)            =',i3dir(1,1)            
      WRITE(201,20)'j1Cvodn(1)          =',j1Cvodn(1)          ,'i1Cvodn(1)          =',i1Cvodn(1,1)          
      WRITE(201,20)'j1Cvndx(1)          =',j1Cvndx(1)          ,'i1Cvndx(1)          =',i1Cvndx(1,1)          
      WRITE(201,20)'j1jndx(1)           =',j1jndx(1)           ,'i1jndx(1)           =',i1jndx(1,1)           
      WRITE(201,20)'j1Rvodn(1)          =',j1Rvodn(1)          ,'i1Rvodn(1)          =',i1Rvodn(1,1)          
      WRITE(201,20)'j1Rvndx(1)          =',j1Rvndx(1)          ,'i1Rvndx(1)          =',i1Rvndx(1,1)          
WRITE(201,20)'jndxsp(10,1)        =',jndxsp(10,1)        ,'jndxsp(10,1)        =',jndxsp(10,1)        
WRITE(201,20)'j1nic(10,1)         =',j1nic(10,1)         ,'j1nic(10,1)         =',j1nic(10,1)         
WRITE(201,20)'s3sdbvolpack        =',s3sdbvolpack        ,'s3sdbvolpack        =',s3sdbvolpack         
WRITE(201,20)'s3sdbvolpacko(1)    =',s3sdbvolpacko(1)    ,'s3sdbvolpacko(1)    =',s3sdbvolpacko(1)    
WRITE(201,20)'s3c3pack(100,100)   =',s3c3pack(100,100)   ,'s3c3pack(100,100)   =',s3c3pack(100,100)    
      WRITE(201,20)'s3vpackflag(100,100)=',s3vpackflag(100,100),'vpackflag(100,100)  =',vpackflag(100,100) 
      WRITE(201,20)'s3hpackflag(100,100)=',s3hpackflag(100,100),'hpackflag(100,100)  =',hpackflag(100,100) 
WRITE(201,20)'s3overcorrection(1) =',s3overcorrection(1) ,'s3overcorrection(1) =',s3overcorrection(1) 
WRITE(201,20)'j1max(10)           =',j1max(10)           ,'j1max(10)           =',j1max(10)            

    15 format(10x,a,1e12.5,10x,a,1e12.5)
    20 format(10x,a,1i12,10x,a,1i12)  

!++++++++++++++++++++ Writing all variables in c3com_space CWCHOI ++++++++++++++++++++
            
            
      RETURN 
      END SUBROUTINE c3com_print_space

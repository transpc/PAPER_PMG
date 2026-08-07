!comdeck ufiles                                                         
!                                                                       
      INTEGER input,output,rstplt,stripf,plotfl,sth2xt,tty,eoin,jbinfo 
      INTEGER coupfl,inpout,error_i,screen, pltrec
!SUH+ 2018.04.03
	  INTEGER cpdplt
!SUH+ 2018.04.03
      COMMON/ufiles_cupid/input,output,rstplt,stripf,plotfl,sth2xt,tty,eoin,  &
!SUH+ 2018.04.03
!SUH  jbinfo,coupfl,inpout,error_i,screen,pltrec 
      jbinfo,coupfl,inpout,error_i,screen,pltrec,cpdplt
!SUH+ 2018.04.03
      INTEGER max_size,max_c,max_vj,max_hts,max_lv,max_tbls
      COMMON/memsize_cupid/max_size,max_c,max_vj,max_hts,max_lv,max_tbls
!
!  input   input file.
!  output  printed output file.
!  rstplt  restart-plot file.
!  stripf  strip file.
!  plotfl  scratch file for internal plot capability.
!  pltrec  plotrec file for APTPLOT (exclude restart information)
!  sth2xt  used for all thermodynamic property files.
!  tty     online screen file.
!  jbinfo  user file to be copied to job output file.
!  filsch  holds file names used in open statements (except for
!          thermodynamic properties files);  default values set by data
!          statements in blkdta;  user supplied values can be optionally
!          entered in some machine versions.
!  cpdplt  restart-plot file for CUPID. filsch(16)
!
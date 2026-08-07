include ./makefile.in

MYNAME = cupid
MYNAME1 = S_CUPID
VER = 10
EXCNAME = $(MYNAME).x
SRCNAME = $(MYNAME).for
DATE = $(shell date +%Y%m%d_%HH)
BACKUPNAME = $(MYNAME1)_$(VER)_$(DATE).tar
#
D0 = .
include ./makefile_dir.in
#
LFLAGS =
#LFLAGS =-L$(D2) -l$(D2) \
#	-L$(D3) -l$(D3) \
#	-L$(D0_6) -l$(D0_6) \
#	-L$(D0_7) -l$(D0_7)
FFLAGS = $(FOPTFLAGS)
# CFLAGS = $(COPTFLAGS) -I/$(INCDIR)
OFLAGS = $(OOPTFLAGS) 

default : mk_exe

mk_exe: 
	cd $(D0_0); make -j $(p) & \
	cd ../$(D0_d); make -j $(p) & \
	wait;
	cd $(D0_a); make -j $(p) & \
	cd ../../$(D0_b); make -j $(p) & \
	cd ../../$(D0_c); make -j $(p) & \
	cd ../../$(D0_e); make -j $(p) & \
	cd ../../$(D0_g); make -j $(p) & \
	cd ../../$(D0_h); make -j $(p) & \
	cd ../../$(D0_i); make -j $(p) & \
	cd ../../$(D0_j); make -j $(p) & \
	cd ../../$(D0_l); make -j $(p) & \
	cd ../../$(D0_G1); make -j $(p) & \
	wait;
	cd $(D0_k); make -j $(p)
	cd $(D0_1); make -j $(p) & \
	cd ../$(D0_2); make -j $(p) & \
	wait;
	cd $(D0_3); make -j $(p)
	cd $(D0_31); make -j $(p) & \
	cd ../../$(D0_311); make -j $(p) & \
	cd ../../../$(D0_312); make -j $(p) & \
	cd ../../../$(D0_313); make -j $(p) & \
	cd ../../../$(D0_314); make -j $(p) & \
	cd ../../../$(D0_315); make -j $(p) & \
	cd ../../../$(D0_316); make -j $(p) & \
	cd ../../../$(D0_318); make -j $(p) & \
	cd ../../../$(D0_319); make -j $(p) & \
	cd ../../../$(D0_32); make -j $(p) & \
	cd ../../$(D0_321); make -j $(p) & \
	cd ../../../$(D0_322); make -j $(p) & \
	cd ../../../$(D0_323); make -j $(p) & \
	wait;
	cd $(D0_4); make -j $(p)
	cd $(D0_G); make -j $(p)
	cd $(D0_5); make -j $(p)
	cd $(D0_51); make -j $(p) & \
	cd ../../$(D0_52); make -j $(p) & \
	wait;
	cd $(D0_6); make -j $(p) & \
	cd ../$(D0_7); make -j $(p) & \
	wait;
	cd $(D0_71); make -j $(p)
	cd $(D0_8); make -j $(p)
	cd $(D0_8_1); make -j $(p) & \
	cd ../../../$(D0_8_2); make -j $(p) & \
	cd ../../../$(D0_8_3); make -j $(p) & \
	cd ../../../$(D0_8_4); make -j $(p) & \
	cd ../../../$(D0_8_6); make -j $(p) & \
	cd ../../../$(D0_8_8); make -j $(p) & \
	cd ../../../$(D0_8_9); make -j $(p) & \
	cd ../../../$(D0_8_10); make -j $(p) & \
	cd ../../../$(D0_8_11); make -j $(p) & \
	cd ../../../$(D0_8_12); make -j $(p) & \
	cd ../../../$(D0_8_13); make -j $(p) & \
	cd ../../../$(D0_8_14); make -j $(p) & \
	cd ../../../$(D0_8_16); make -j $(p) & \
	cd ../../../$(D0_8_17); make -j $(p) & \
	cd ../../../$(D0_8_18); make -j $(p) & \
	cd ../../../$(D0_8_19); make -j $(p) & \
	cd ../../../$(D0_8_20); make -j $(p) & \
	cd ../../../$(D0_8_21); make -j $(p) & \
	cd ../../../$(D0_8_22); make -j $(p) & \
	cd ../../../$(D0_8_23); make -j $(p) & \
	cd ../../../$(D0_8_24); make -j $(p) & \
	cd ../../../$(D0_8_25); make -j $(p) & \
	cd ../../../$(D0_8_26); make -j $(p) & \
	cd ../../../$(D0_8_27); make -j $(p) & \
	wait;
	cd $(D0_9); make -j $(p) & \
	cd ../$(D0_10); make -j $(p) & \
	cd ../$(D0_11); make -j $(p) & \
	wait;
	cd $(D0_C); make -j $(p) & \
	cd ../$(D0_C2); make -j $(p) & \
	cd ../../$(D0_C1); make -j $(p) & \
	wait;
	cd $(D0_13); make -j $(p)
	cd $(D0_14); make -j $(p)
	$(LD) -o $(EXCNAME) $(OFLAGS) $(OBJS) $(LIBDIR)
clean:
	(cd $(D0_0); make $@)
	(cd $(D0_a); make $@)
	(cd $(D0_b); make $@)
	(cd $(D0_c); make $@)
	(cd $(D0_d); make $@)
	(cd $(D0_e); make $@)
	(cd $(D0_g); make $@)
	(cd $(D0_h); make $@)
	(cd $(D0_j); make $@)
	(cd $(D0_1); make $@)
	(cd $(D0_2); make $@)
	(cd $(D0_3); make $@)
	(cd $(D0_31); make $@)
	(cd $(D0_311); make $@)
	(cd $(D0_312); make $@)
	(cd $(D0_313); make $@)
	(cd $(D0_314); make $@)
	(cd $(D0_315); make $@)
	(cd $(D0_4); make $@)
	(cd $(D0_5); make $@)
	(cd $(D0_6); make $@)
	(cd $(D0_7); make $@)
	(cd $(D0_8); make $@)
	(cd $(D0_8_1); make $@)
	(cd $(D0_8_2); make $@)
	(cd $(D0_8_3); make $@)
	(cd $(D0_8_4); make $@)
	(cd $(D0_8_6); make $@)
	(cd $(D0_8_8); make $@)
	(cd $(D0_8_9); make $@)
	(cd $(D0_8_10); make $@)	
	(cd $(D0_8_11); make $@)	
	(cd $(D0_8_12); make $@)
	(cd $(D0_8_13); make $@)
	(cd $(D0_8_14); make $@)
	(cd $(D0_8_16); make $@)
	(cd $(D0_8_17); make $@)
	(cd $(D0_8_18); make $@)
	(cd $(D0_9); make $@)
	(cd $(D0_10); make $@)
	(cd $(D0_C); make $@)
	rm -f *~
	rm -f *.o 
	rm -f *.mod
	rm -f *__genmod.*
#	rm -f $(EXCNAME)
	rm -f $(LIBNAME)
	rm -f $(SRCNAME)
new : clean default

	(cd $(D0_0); make $@)
	(cd $(D0_a); make $@)
	(cd $(D0_b); make $@)
	(cd $(D0_c); make $@)
	(cd $(D0_d); make $@)
	(cd $(D0_e); make $@)
	(cd $(D0_f); make $@)
	(cd $(D0_g); make $@)
	(cd $(D0_h); make $@)
	(cd $(D0_j); make $@)
	(cd $(D0_1); make $@)
	(cd $(D0_2); make $@)
	(cd $(D0_3); make $@)
	(cd $(D0_31); make $@)
	(cd $(D0_311); make $@)
	(cd $(D0_312); make $@)
	(cd $(D0_313); make $@)
	(cd $(D0_314); make $@)
	(cd $(D0_315); make $@)
	(cd $(D0_4); make $@)
	(cd $(D0_5); make $@)
	(cd $(D0_6); make $@)
	(cd $(D0_7); make $@)
	(cd $(D0_8); make $@)
	(cd $(D0_8_1); make $@)
	(cd $(D0_8_2); make $@)
	(cd $(D0_8_3); make $@)
	(cd $(D0_8_4); make $@)
	(cd $(D0_8_6); make $@)
	(cd $(D0_8_8); make $@)
	(cd $(D0_8_9); make $@)
	(cd $(D0_8_10); make $@)	
	(cd $(D0_8_11); make $@)	
	(cd $(D0_8_12); make $@)
	(cd $(D0_8_13); make $@)			
	(cd $(D0_8_14); make $@)
	(cd $(D0_8_16); make $@)
	(cd $(D0_8_17); make $@)
	(cd $(D0_8_18); make $@)
	(cd $(D0_8_19); make $@)
	(cd $(D0_8_20); make $@)
	(cd $(D0_8_21); make $@)
	(cd $(D0_8_22); make $@)
	(cd $(D0_8_23); make $@)
	(cd $(D0_8_24); make $@)
	(cd $(D0_8_25); make $@)
	(cd $(D0_8_26); make $@)
	(cd $(D0_9); make $@)
	(cd $(D0_10); make $@)
	(cd $(D0_C); make $@)
	cat $(SRC) $(D0_0)/$(D0_0).for $(D0_a)/$(D0_a).for $(D0_b)/$(D0_b).for $(D0_c)/$(D0_c).for $(D0_d)/$(D0_d).for $(D0_e)/$(D0_e).for $(D0_f)/$(D0_f).for $(D0_g)/$(D0_g).for $(D0_h)/$(D0_h).for $(D0_i)/$(D0_h).for $(D0_j)/$(D0_j).for> $(SRCNAME)
	cat $(SRC) $(D0_1)/$(D0_1).for $(D0_2)/$(D0_2).for $(D0_3)/$(D0_3).for $(D0_31)/$(D0_31).for $(D0_311)/$(D3_b).for $(D0_312)/$(D0_312).for $(D0_313)/$(D0_313).for $(D0_314)/$(D0_314).for $(D0_315)/$(D0_315).for $(D0_4)/$(D0_4).for $(D0_5)/$(D0_5).for $(D0_6)/$(D0_6).for $(D0_7)/$(D0_7).for $(D0_8)/$(D0_8).for $(D0_9)/$(D0_9).for $(D0_10)/$(D0_10).for > $(SRCNAME)
	cat $(SRC) $(D0_C)/$(D0_C).for> $(SRCNAME)
	cat $(SRC) $(D0_8_1)/$(D0_8_1).for $(D0_8_2)/$(D0_8_2).for $(D0_8_3)/$(D0_8_3).for $(D0_8_4)/$(D0_8_4).for $(D0_8_6)/$(D0_8_6).for $(D0_8_8)/$(D0_8_8).for $(D0_8_9)/$(D0_8_9).for $(D0_8_10)/$(D0_8_10).for $(D0_8_11)/$(D0_8_11).for $(D0_8_12)/$(D0_8_12).for $(D0_8_13)/$(D0_8_13).for $(D0_8_14).for/$(D0_8_14).for $(D0_8_16).for/$(D0_8_16).for $(D0_8_17).for/$(D0_8_17).for $(D0_8_18).for/$(D0_8_18).for$(D0_8_18).for $(D0_8_19).for/$(D0_8_19).for $(D0_8_20).for/$(D0_8_20).for $(D0_8_21).for/$(D0_8_21).for $(D0_8_22).for/$(D0_8_22).for $(D0_8_23).for/$(D0_8_23).for $(D0_8_24).for/$(D0_8_24).for $(D0_8_25).for/$(D0_8_25).for $(D0_8_26).for/$(D0_8_26).for > $(SRCNAME) 
tag :
	(cd $(D0_0); make $@)
	(cd $(D0_a); make $@)
	(cd $(D0_b); make $@)
	(cd $(D0_c); make $@)
	(cd $(D0_d); make $@)
	(cd $(D0_e); make $@)
	(cd $(D0_f); make $@)
	(cd $(D0_g); make $@)
	(cd $(D0_h); make $@)
	(cd $(D0_j); make $@)
	(cd $(D0_i); make $@)
	(cd $(D0_1); make $@)
	(cd $(D0_2); make $@)
	(cd $(D0_3); make $@)
	(cd $(D0_31); make $@)
	(cd $(D0_311); make $@)
	(cd $(D0_312); make $@)
	(cd $(D0_313); make $@)
	(cd $(D0_314); make $@)
	(cd $(D0_315); make $@)
	(cd $(D0_4); make $@)
	(cd $(D0_5); make $@)
	(cd $(D0_6); make $@)
	(cd $(D0_7); make $@)
	(cd $(D0_8); make $@)
	(cd $(D0_8_1); make $@)
	(cd $(D0_8_2); make $@)
	(cd $(D0_8_3); make $@)
	(cd $(D0_8_4); make $@)
	(cd $(D0_8_6); make $@)
	(cd $(D0_8_8); make $@)
	(cd $(D0_8_9); make $@)
	(cd $(D0_8_10); make $@)
	(cd $(D0_8_11); make $@)
	(cd $(D0_8_12); make $@)
	(cd $(D0_8_13); make $@)
	(cd $(D0_8_14); make $@)
	(cd $(D0_8_16); make $@)
	(cd $(D0_8_17); make $@)
	(cd $(D0_8_18); make $@)
	(cd $(D0_8_19); make $@)
	(cd $(D0_8_20); make $@)
	(cd $(D0_8_21); make $@)
	(cd $(D0_8_22); make $@)
	(cd $(D0_8_23); make $@)
	(cd $(D0_8_24); make $@)
	(cd $(D0_8_25); make $@)
	(cd $(D0_8_26); make $@)
	(cd $(D0_9); make $@)
	(cd $(D0_10); make $@)
	(cd $(D0_C); make $@)
	ctags $(SRC) $(D0_0)/*.f90 $(D0_a)/*.f90 $(D0_b)/*.f90 $(D0_c)/*.f90 $(D0_d)/*.f90 $(D0_e)/*.f90 $(D0_f)/*.f90 $(D0_g)/*.f90 $(D0_h)/*.f90 $(D0_j)/*.f90
	ctags $(SRC) $(D1)/*.f90 $(D2)/*.f90 $(D3)/*.f90 $(D3_a)/*.f90 $(D0_311)/*.f90 $(D0_312)/*.f90 $(D0_313)/*.f90 $(D0_314)/*.f90 $(D0_315)/*.f90 $(D0_4)/*.f90 $(D0_5)/*.f90 $(D0_6)/*.f90 $(D0_7)/*.f90 $(D0_8)/*.f90 $(D0_9)/*.f90 $(D0_10)/*.f90
	ctags $(SRC) $(D0_C)/*.f90
	ctags $(SRC) $(D0_8_1)/*.f90 $(D0_8_2)/*.f90 $(D0_8_3)/*.f90 $(D0_8_4)/*.f90 $(D0_8_6)/*.f90 $(D0_8_8)/*.f90 $(D0_8_9)/*.f90 $(D0_8_10)/*.f90 $(D0_8_11)/*.f90 $(D0_8_12)/*.f90 $(D0_8_13)/*.f90 $(D0_8_14)/*.f90 $(D0_8_16)/*.f90 $(D0_8_17)/*.f90 $(D0_8_18)/*.f90 $(D0_8_19)/*.f90 $(D0_8_20)/*.f90 $(D0_8_21)/*.f90 $(D0_8_22)/*.f90 $(D0_8_23)/*.f90 $(D0_8_24)/*.f90 $(D0_8_25)/*.f90 $(D0_8_26)/*.f90
backup :
	(cd $(D0_0); make $@)
	(cd $(D0_a); make $@)
	(cd $(D0_b); make $@)
	(cd $(D0_c); make $@)
	(cd $(D0_d); make $@)
	(cd $(D0_e); make $@)
	(cd $(D0_f); make $@)
	(cd $(D0_g); make $@)
	(cd $(D0_h); make $@)
	(cd $(D0_j); make $@)
	(cd $(D0_1); make $@)
	(cd $(D0_2); make $@)
	(cd $(D0_3); make $@)
	(cd $(D0_31); make $@)
	(cd $(D0_311); make $@)
	(cd $(D0_312); make $@)
	(cd $(D0_313); make $@)
	(cd $(D0_314); make $@)
	(cd $(D0_315); make $@)
	(cd $(D0_4); make $@)
	(cd $(D0_5); make $@)
	(cd $(D0_6); make $@)
	(cd $(D0_7); make $@)
	(cd $(D0_8); make $@)
	(cd $(D0_8_1); make $@)
	(cd $(D0_8_2); make $@)
	(cd $(D0_8_3); make $@)
	(cd $(D0_8_4); make $@)
	(cd $(D0_8_6); make $@)
	(cd $(D0_8_8); make $@)
	(cd $(D0_8_9); make $@)
	(cd $(D0_8_10); make $@)	
	(cd $(D0_8_11); make $@)	
	(cd $(D0_8_12); make $@)
	(cd $(D0_8_13); make $@)
	(cd $(D0_8_14); make $@)
	(cd $(D0_8_16); make $@)
	(cd $(D0_8_17); make $@)
	(cd $(D0_8_18); make $@)
	(cd $(D0_8_19); make $@)
	(cd $(D0_8_20); make $@)
	(cd $(D0_8_21); make $@)
	(cd $(D0_8_22); make $@)
	(cd $(D0_8_23); make $@)
	(cd $(D0_8_24); make $@)
	(cd $(D0_8_25); make $@)
	(cd $(D0_8_26); make $@)
	(cd $(D0_9); make $@)
	(cd $(D0_10); make $@)
	(cd $(D0_C); make $@)
	tar cvf $(BACKUPNAME) $(BACKUPFILE)
xtar: 
	(cd $(D0_0); tar xvf *.tar)
	(cd $(D0_a); tar xvf *.tar)
	(cd $(D0_b); tar xvf *.tar)
	(cd $(D0_c); tar xvf *.tar)
	(cd $(D0_d); tar xvf *.tar)
	(cd $(D0_e); tar xvf *.tar)
	(cd $(D0_f); tar xvf *.tar)
	(cd $(D0_g); tar xvf *.tar)
	(cd $(D0_h); tar xvf *.tar)
	(cd $(D0_i); tar xvf *.tar)
	(cd $(D0_j); tar xvf *.tar)
	(cd $(D0_1); tar xvf *.tar)
	(cd $(D0_2); tar xvf *.tar)
	(cd $(D0_3); tar xvf *.tar)
	(cd $(D0_31); tar xvf *.tar)
	(cd $(D0_311); tar xvf *.tar)
	(cd $(D0_312); tar xvf *.tar)
	(cd $(D0_313); tar xvf *.tar)
	(cd $(D0_314); tar xvf *.tar)
	(cd $(D0_315); tar xvf *.tar)
	(cd $(D0_4); tar xvf *.tar)
	(cd $(D0_5); tar xvf *.tar)
	(cd $(D0_6); tar xvf *.tar)
	(cd $(D0_7); tar xvf *.tar)
	(cd $(D0_8); tar xvf *.tar)
	(cd $(D0_8_1); tar xvf *.tar)
	(cd $(D0_8_2); tar xvf *.tar)
	(cd $(D0_8_3); tar xvf *.tar)
	(cd $(D0_8_4); tar xvf *.tar)
	(cd $(D0_8_6); tar xvf *.tar)
	(cd $(D0_8_8); tar xvf *.tar)
	(cd $(D0_8_9); tar xvf *.tar)
	(cd $(D0_8_10); tar xvf *.tar)	
	(cd $(D0_8_11); tar xvf *.tar)	
	(cd $(D0_8_12); tar xvf *.tar)
	(cd $(D0_8_13); tar xvf *.tar)
	(cd $(D0_8_14); tar xvf *.tar)
	(cd $(D0_8_16); tar xvf *.tar)
	(cd $(D0_8_17); tar xvf *.tar)
	(cd $(D0_8_18); tar xvf *.tar)
	(cd $(D0_8_19); tar xvf *.tar)
	(cd $(D0_8_20); tar xvf *.tar)
	(cd $(D0_8_21); tar xvf *.tar)
	(cd $(D0_8_22); tar xvf *.tar)
	(cd $(D0_8_23); tar xvf *.tar)
	(cd $(D0_8_24); tar xvf *.tar)
	(cd $(D0_8_25); tar xvf *.tar)
	(cd $(D0_8_26); tar xvf *.tar)
	(cd $(D0_9); tar xvf *.tar)
	(cd $(D0_10); tar xvf *.tar)
	(cd $(D0_C); tar xvf *.tar)

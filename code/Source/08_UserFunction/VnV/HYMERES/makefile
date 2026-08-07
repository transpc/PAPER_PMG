include ../../../makefile.in

MYNAME = cupidUSER
MYNAME1 = s_cupidUSER
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

FFLAGS = $(FOPTFLAGS) -module ../../../Modules
CFLAGS =
#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = udfn_H2P1_initialHeTg.f90\
	udfn_H2P1_inputTg.f90\
	udfn_H2P1_inputTw.f90\
	udfn_H2P1_outputFOV.f90\
	udfn_H2P1_outputHeTg.f90\
	udfn_H2P1_outputTE.f90\
	udfn_H2P1_ouputEF.f90\
	udfn_tw_profile.f90\
	udfn_CUBE_inputBC.f90\
	udfn_CUBE_inputTW.f90\
	udfn_CUBE_outputTg.f90\
	udfn_CUBE_outputTs.f90\
#----------------------------------------------------------------------

FSRC = $(FSRC1)

SRC  = $(FSRC)

OBJS = $(FSRC1:.f90=.o)

BACKUPFILE = *.f90 *.h makefile
.SUFFIXES: .f90 .o
.F.o :
	$(FC) $(FFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<
.f90.o :
	$(FC) $(FFLAGS) -c $<

default : mk_lib

mk_exe: $(OBJS)
	$(FC) -o $(EXCNAME) $(OBJS) $(LFLAGS)

mk_lib: $(OBJS)

clean:
#	rm -f *~
#	rm -f *.o 
#	rm -f *.mod
#	rm -f *__genmod.*
#	rm -f $(EXCNAME)
#	rm -f $(LIBNAME)
#	rm -f $(SRCNAME)
	
new : clean default

src :
	cat $(SRC) > $(SRCNAME)
tag :
	ctags $(SRC)
backup :
	tar cvf $(BACKUPNAME) $(BACKUPFILE)

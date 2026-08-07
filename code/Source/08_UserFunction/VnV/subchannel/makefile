include ../../../makefile.in

MYNAME = cupidUSER
MYNAME1 = s_cupidUSER
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = 1subchannel_geo.f90\
	udfn_mom_source_explicit.f90\
	udfn_mom_source_implicit.f90
FSRCV = udfn_mom_source_evvd.f90\
	udfn_mom_source_em.f90\
	udfn_sij.f90\
	udfn_s_gapij.f90
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../../../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../../../Modules

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

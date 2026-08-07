include ../makefile.in

MYNAME = cupidPRO
MYNAME1 = s_cupidPRO
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = convert_temp2erg.f90\
	init_steamtable.f90\
	mat_prop.f90\
	property_calc.f90\
	read_tpfh2o.f90\
	steamtable.f90\
	sth2x1.f90\
	sth2x3.f90\
	sth2x6.f90\
	stread.f90\
	thcond_cupid.f90\
	viscos_cupid.f90\
	viscos_lw_cupid.f90\
	thcond_lw_cupid.f90\
	convert_temp2erg_gas.f90\
	convert_temp2erg_liq.f90
FSRCV = ncg_cell.f90\
	surftn.f90
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../Modules

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
	$(ARCH) $(ARCHFLAGS) $(LIBNAME) $(OBJS)
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


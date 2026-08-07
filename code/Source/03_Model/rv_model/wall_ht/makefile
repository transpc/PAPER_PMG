include ../../../makefile.in

MYNAME = cupidMODrv7
MYNAME1 = s_cupidMODrv7
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

FFLAGS = $(FOPTFLAGS) -module ../../../Modules
CFLAGS =
#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = CHF_calc.f90\
	condensation_HTC.f90\
	initialize_cond.f90\
	nucl_boiling_HTC.f90\
	psatpd.f90\
	rv_property.f90\
	single_phase_HTC.f90\
	subcooled_boiling.f90\
	trans_film_boiling_HTC.f90\
	trans_film_reflood_HTC.f90\
	condensation_ncg.f90
#----------------------------------------------------------------------

FSRC = $(FSRC1)

SRC  = $(FSRC)

OBJS = $(FSRC1:.f90=.o)

BACKUPFILE = *.f90 makefile
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


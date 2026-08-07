include ../../makefile.in

MYNAME = cupidSOL1
MYNAME1 = s_cupidSOL1
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

# LFLAGS = -L/cm/shared/apps/metis/METIS -lmetis

FFLAGS = $(FOPTFLAGS) -module ../../Modules
# CFLAGS = -I$(COPTFLAGS) \-L$(INCDIR)
CFLAGS = 
#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = grad_conv.f90\
	grad_limiter.f90\
	quality_2nd_conv.f90\
	quality_2nd_conv_imp.f90\
	mom_2nd_conv_imp.f90\
	mom_2nd_conv.f90\
	mult_ncg_2nd_conv_imp.f90\
	mult_ncg_2nd_conv.f90\
	mass_2nd_conv.f90\
	mass_2nd_conv_imp.f90\
	energy_2nd_conv.f90\
	energy_2nd_conv_imp.f90\
	boron_2nd_conv.f90\
	boron_2nd_conv_imp.f90
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
	$(FC) -o $(EXCNAME) $(OBJS) $(CFLAGS) $(LFLAGS)

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

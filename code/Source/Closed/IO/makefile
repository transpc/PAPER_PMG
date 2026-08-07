include ../../makefile.in

MYNAME = cupidCLOSEDio
MYNAME1 = s_cupidCLOSEDio
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

FFLAGS = $(FOPTFLAGS) -module ../../Modules
CFLAGS =
#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = allocate_geo_var.f90\
	allocate_model_var.f90\
	allocate_physical_var.f90\
	allocate_preliminary_var.f90\
	allocate_var_solid.f90\
	find_walls.f90\
	frink_weight.f90\
	restart_read_rv.f90\
	restart_write_rv.f90\
	find_icelltype.f90
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


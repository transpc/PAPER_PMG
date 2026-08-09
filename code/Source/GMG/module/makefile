include ../../makefile.in

MYNAME = cupidCLOSED
MYNAME1 = s_cupidCLOSED
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = MD_matrix.f90\
	MD_parameter.f90\
	MD_MG_index.f90\
	MD_geometry.f90\
	MD_MG_coord.f90\
	MD_MG_matrix.f90\
	MD_connectivity.f90\
	MD_MPI.f90\
	dump_pmg.f90
FSRCV = 
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../../Modules

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
	$(FC) -o $(EXCNAME) $(OBJS)  $(CFLAGS) $(LFLAGS)

dump_pmg.o: MD_matrix.o MD_MPI.o MD_geometry.o MD_parameter.o

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

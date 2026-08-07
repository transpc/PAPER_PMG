include ../makefile.in

MYNAME = MASTER
MYNAME1 = s_MASTER
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

LFLAGS = 

FFLAGS = $(FOPTFLAGS) -module ../Modules
CFLAGS = $(COPTFLAGS) -I. $(INCDIR)

#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = t_masterC.f90\
	t_masterC_debug.f90\
	t_masterC_rod.f90\
	t_masterC_debug.f90\
	CM_Control_Reactivity.f90\
	save_restart_master.f90\
	read_cupid_master.f90\
	check_steady_dtl.f90\
	check_steady_dmass.f90\
	print_ss_status.f90\
	check_steady_dpower.f90\
	rv_power_transient_master.f90\
	rv_power_master_to_cupid.f90\
	core_boron_concentration.f90
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
	$(ARCH) $(ARCHFLAGS) $(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBNAME)

clean:
	rm -f *~
	rm -f *.o 
	rm -f *.mod
	rm -f *__genmod.*
	rm -f $(EXCNAME)
	rm -f $(LIBNAME)
	rm -f $(SRCNAME)
	
new : clean default

src :
	cat $(SRC) > $(SRCNAME)
tag :
	ctags $(SRC)
backup :
	tar cvf $(BACKUPNAME) $(BACKUPFILE)


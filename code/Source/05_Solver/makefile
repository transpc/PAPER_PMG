include ../makefile.in

MYNAME = cupidSOL
MYNAME1 = s_cupidSOL
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = scalar_reset.f90\
	press_work_face.f90
FSRCV = momentum_diffusion.f90\
	momentum_convection.f90\
	pressure_matrix.f90\
	pressure_matrix_smac.f90\
	scalar_matrix.f90\
	scalar_matrix_solid.f90\
	scalar_energy_convection.f90\
	scalar_energy_diffusion.f90\
	scalar_work_convection.f90\
	scalar_xn_convection.f90\
	scalar_mass_diffusion.f90\
	scalar_mass_convection.f90\
	scalar_update.f90\
	scalar_update_smac3.f90\
	pressure_solve.f90\
	ncg_transport.f90\
	boron_convection.f90\
	continuity_for_smac3.f90\
	implicit_scalar.f90\
	implicit_momentum.f90\
	imp_mom_diffusion.f90\
	imp_mom_convection.f90\
	imp_eng_diffusion.f90\
	imp_scalar_convection.f90\
	continuity_for_smac3_imp.f90\
	boron_transport.f90\
	ncg_transport_imp.f90\
	ncg_convection_imp.f90\
	non_orthogonal_src.f90\
	set_dt.f90\
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

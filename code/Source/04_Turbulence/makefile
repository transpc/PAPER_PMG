include ../makefile.in

MYNAME = cupidTURB
MYNAME1 = s_cupidTURB
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = turb_2nd_conv.f90\
	turb_2nd_conv_imp.f90\
	Turb_ke_calc_gas.f90\
	Turb_ke_calc_liq.f90\
	Turb_ke_cond.f90\
	Turb_ke_convection.f90\
	Turb_ke_diffusion.f90\
	Turb_ke_gas.f90\
	Turb_ke_init.f90\
	Turb_ke_liq.f90\
	Turb_ke_product.f90\
	Turb_ke_vis_gas.f90\
	Turb_ke_vis_lg.f90\
	Turb_ke_vis_liq.f90\
	Turb_ke_vt_tauw.f90\
	Turb_ke_vtangent.f90\
	Turb_SST_calc_gas.f90\
	Turb_SST_calc_liq.f90\
	Turb_SST_gas.f90\
	Turb_SST_init.f90\
	Turb_SST_liq.f90\
	Turb_SST_product.f90\
	Turb_SST_vis_gas.f90\
	Turb_SST_vis_liq.f90\
	Turb_SST_vt_tauw.f90\
	Turb_zeq.f90\
	turb_zeq_wallfn.f90\
	Turbulence_MOD.f90\
	Turb_ke_product_real.f90\
	Turb_LES_WALE.f90
FSRCV = imp_diffusion_ke.f90
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../Modules

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


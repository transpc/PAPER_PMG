include ../makefile.in

MYNAME = MODULE
MYNAME1 = s_MODULE
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

FFLAGS = $(FOPTFLAGS) -module ../Modules
CFLAGS = 
#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = Zparam.f90\
	Zare.f90\
	Zbc_index.f90\
	Zb_condition.f90\
	Zbicg.f90\
	Zboron.f90\
	Zconst1.f90 Zconst2.f90\
	Zcoord1.f90 Zcoord2.f90 Zcoord3.f90 Zcoord4.f90\
	Zdaint_ag.f90\
	Zdecoupled.f90\
	Zdel_scalar.f90\
	Zdhda.f90\
	Zenergy_diff.f90\
	Zface.f90\
	Zflowregime.f90\
	Zflux.f90\
	Zgrad_ls_c2d.f90 Zgrad_ls_c3d.f90\
	Zgradoption.f90\
	Zheat_partition.f90\
	Ziat.f90\
	Zimplicit.f90\
	Zinlet_bc.f90\
	Zmars.f90\
	Zmaster.f90\
	Zmass_conv.f90 Zmass_diff.f90\
	Zmodel.f90\
	Zm_src.f90\
	Zncg.f90\
	Zndforce.f90\
	Znode.f90\
	Znormal.f90\
	Znum_cell.f90\
	Zporous.f90\
	Zpress_coeff.f90 Zpress.f90\
	Zrocom_specific.f90 \
	Zqvol.f90\
	Zscalar_coeff.f90 Zscalar_src.f90\
	Ztimecon.f90\
	Zturb.f90 Zturbzeq.f90\
	Zuserdefined.f90\
	Zvector.f90\
	Zvoid.f90\
	Zzone.f90\
	Z2nd_order.f90
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

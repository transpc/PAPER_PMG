include ../makefile.in

MYNAME = cupidCLOSED
MYNAME1 = s_cupidCLOSED
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = md_function.f90\
	    06_solver_pcg_ilu.f90\
	    6_solver_pbcg_mg.f90\
        1_read_input.f90\
        2_Prep_fine_global.f90\
        2_read_mesh_MPI.f90\
        3_Prep_fine_P.f90\
        4_Prep_MG_GarL.f90\
        5_PREP_GMG.f90\
		5_PREP_GMG_global_coarse.f90\
        6_subdomain_infor_mg.f90\
        7_assemble_FVM.f90\
        7_SOLVE_GMG.f90\
		blaslapack_sub_n.f90\
        bubble_sort.f90\
        coarsening_semi.f90\
        connectivity_coarse.f90\
        mt_amux.f90\
        PCG_Dig.f90\
		poly_smooth.f90\
        Domain_infor_FVM_fine.f90\
        Domain_infor_FVM_coarse.f90\
        Domain_infor_FVM_coarsest.f90\
        Relax_GSP.f90\
        send_receive.f90\
        send_receive_csr.f90\
        send_receive_mt.f90\
        SOLVE_GC.f90\
        stiffness_MG.f90\
	    mt_trans.f90\
	    stiffness_GC.f90\
        imapGZ_coarse.f90\
        neighboring.f90\
        write_report.f90\
        interpolation.f90\
        Neighbor_node_ARP.f90\
        linear_algebra.f90

FSRCV = 
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
#


mk_exe: $(OBJS)
	$(FC) -o $(EXCNAME) $(OBJS)  $(CFLAGS) $(LFLAGS)

mk_lib: $(OBJS)
mt_trans.o: md_function.f90
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

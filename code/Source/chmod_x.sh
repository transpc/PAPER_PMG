#!/bin/sh -x
#
cd ../run_vv
ifort -O compare2.f -o compare2
ifort -O compare_profile.f -o compare_profile
ifort -O compare_profile1.f -o compare_profile1
cd ../SourceM
#
chmod +x ../1.rm
chmod +x 1.rm
chmod +x diff0
chmod +x diff_all.sh
chmod +x diff_make.sh
#
chmod +x ../run_vv/diff0
chmod +x ../run_vv/diffsh.sh
chmod +x ../run_vv/diff_all.sh
chmod +x ../run_vv/check_troubles.sh
#
chmod +x ../run_vv/run.sh
chmod +x ../run_vv/run1.sh
chmod +x ../run_vv/run2.sh
chmod +x ../run_vv/run3.sh
chmod +x ../run_vv/run8.sh
chmod +x ../run_vv/run24.sh
chmod +x ../run_vv/run48.sh
chmod +x ../run_vv/runs.sh
#
chmod +x ../run_vv/run1_mpi.sh
chmod +x ../run_vv/run2_mpi.sh
chmod +x ../run_vv/run3_mpi.sh
chmod +x ../run_vv/run4_mpi.sh
chmod +x ../run_vv/run5_mpi.sh
chmod +x ../run_vv/run6_mpi.sh
chmod +x ../run_vv/run7_mpi.sh
chmod +x ../run_vv/run8_mpi.sh
#
chmod +x ../run_vv/run1_mpi_ca.sh
chmod +x ../run_vv/run2_mpi_ca.sh
chmod +x ../run_vv/run3_mpi_ca.sh
chmod +x ../run_vv/run4_mpi_ca.sh
chmod +x ../run_vv/run5_mpi_ca.sh
#
chmod +x ../run_vv/run1_ser.sh
chmod +x ../run_vv/run2_ser.sh
chmod +x ../run_vv/run3_ser.sh
chmod +x ../run_vv/run4_ser.sh
chmod +x ../run_vv/run5_ser.sh
chmod +x ../run_vv/run6_ser.sh
chmod +x ../run_vv/run7_ser.sh
chmod +x ../run_vv/run8_ser.sh
chmod +x ../run_vv/runs_mcc.sh
#
chmod +x ../run_vv/run1_ser_ca.sh
chmod +x ../run_vv/run2_ser_ca.sh
chmod +x ../run_vv/run3_ser_ca.sh
chmod +x ../run_vv/run4_ser_ca.sh
chmod +x ../run_vv/run5_ser_ca.sh
#
chmod +x ../run_vv/copy_comp.sh
chmod +x ../run_vv/copy_comps.sh
chmod +x ../run_vv/copy_comp1.sh
chmod +x ../run_vv/copy_comp2.sh
chmod +x ../run_vv/copy_comp3.sh
chmod +x ../run_vv/copy_comp8.sh
chmod +x ../run_vv/copy_comp24.sh
chmod +x ../run_vv/copy_comp48.sh
chmod +x ../run_vv/copy_comp0.sh
chmod +x ../run_vv/copy_comp0_mcc.sh
#
chmod +x ../run_vv/copy_ref.sh
chmod +x ../run_vv/copy_refs.sh
chmod +x ../run_vv/copy0_ref.sh
chmod +x ../run_vv/compare_ref.sh
chmod +x ../run_vv/compare_refs.sh
#
chmod +x ../run_vv/run-MCC/run1.sh
chmod +x ../run_vv/run-MCC/run2.sh
#
chmod +x ../run_vv/run-MCC/copy_ref.sh
chmod +x ../run_vv/run-MCC/compare_ref.sh
chmod +x ../run_vv/run-MCC/copy0_ref.sh
#
chmod +x ../run_vv/run-MCC/OPR1000_MSLB/run/copy0.sh
chmod +x ../run_vv/run-MCC/OPR1000_MSLB/run/run1.sh
chmod +x ../run_vv/run-MCC/OPR1000_MSLB/run/run2.sh
chmod +x ../run_vv/run-MCC/OPR1000_MSLB/run/run3.sh
chmod +x ../run_vv/run-MCC/OPR1000_MSLB/run/run4.sh
chmod +x ../run_vv/run-MCC/OPR1000_MSLB/run/run.sh
#
chmod +x ../run_vv/run-MCC/APR1400_MSLB/run/copy0.sh
chmod +x ../run_vv/run-MCC/APR1400_MSLB/run/run1.sh
chmod +x ../run_vv/run-MCC/APR1400_MSLB/run/run2.sh
chmod +x ../run_vv/run-MCC/APR1400_MSLB/run/run3.sh
chmod +x ../run_vv/run-MCC/APR1400_MSLB/run/run4.sh
chmod +x ../run_vv/run-MCC/APR1400_MSLB/run/run.sh
#
chmod +x ../run_vv/run-MCC/init_reactor_safety/run/copy0.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety/run/run1.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety/run/run2.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety/run/run3.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety/run/run4.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety/run/run.sh
#
chmod +x ../run_vv/run-MCC/init_reactor_safety_3D/run/copy0.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety_3D/run/run1.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety_3D/run/run2.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety_3D/run/run3.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety_3D/run/run.sh
chmod +x ../run_vv/run-MCC/init_reactor_safety_3D/run/run4.sh

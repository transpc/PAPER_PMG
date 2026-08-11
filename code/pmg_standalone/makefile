# PMG 단독 빌드 (LOOP C006) — 스텁 없이 원본 파일 직접 참조 (소스 이원화 금지, PLAN §4-0)
# 플래그는 ../Source/makefile.in (apptainer 변형) 을 그대로 include — 골든 재현성의 전제
include ../Source/makefile.in

SRCDIR = ../Source
BUILD  = build

# 1) CUPID 본체 리프 Z-모듈 (전부 USE=0 확인, LOG C006) — Zinterface 는 아래 의존들 뒤
ZSRC = $(SRCDIR)/00_Module/Zparam.f90 \
       $(SRCDIR)/00_Module/Zbicg.f90 \
       $(SRCDIR)/00_Module/Zcoord1.f90 \
       $(SRCDIR)/00_Module/Znode.f90 \
       $(SRCDIR)/00_Module/Ztimecon.f90 \
       $(SRCDIR)/00_Module/Zzone.f90 \
       $(SRCDIR)/00_Module/Znum_cell.f90 \
       $(SRCDIR)/00_Module/d_MPI/Zcore.f90 \
       $(SRCDIR)/00_Module/d_MPI/Zmpi.f90 \
       $(SRCDIR)/00_Module/j_rv_model/Zrv_mpi.f90 \
       $(SRCDIR)/00_Module/j_rv_model/Zrv_hts_2d.f90 \
       $(SRCDIR)/00_Module/k_interface/Zinterface.f90

# 2) GMG 전용 모듈 (GMG/module/makefile 의 FSRCS 순서 그대로)
GMOD = $(SRCDIR)/GMG/module/MD_matrix.f90 \
       $(SRCDIR)/GMG/module/MD_parameter.f90 \
       $(SRCDIR)/GMG/module/MD_MG_index.f90 \
       $(SRCDIR)/GMG/module/MD_geometry.f90 \
       $(SRCDIR)/GMG/module/MD_MG_coord.f90 \
       $(SRCDIR)/GMG/module/MD_MG_matrix.f90 \
       $(SRCDIR)/GMG/module/MD_connectivity.f90 \
       $(SRCDIR)/GMG/module/MD_MPI.f90

# 3) GMG 본체 (GMG/makefile 의 FSRCS 순서 그대로 추출 — 6_solver_pbcg_ali/genmod 제외)
GSRC = $(SRCDIR)/GMG/md_function.f90 $(SRCDIR)/GMG/06_solver_pcg_ilu.f90 $(SRCDIR)/GMG/6_solver_pbcg_mg.f90 $(SRCDIR)/GMG/1_read_input.f90 $(SRCDIR)/GMG/2_Prep_fine_global.f90 $(SRCDIR)/GMG/2_read_mesh_MPI.f90 $(SRCDIR)/GMG/3_Prep_fine_P.f90 $(SRCDIR)/GMG/4_Prep_MG_GarL.f90 $(SRCDIR)/GMG/5_PREP_GMG.f90 $(SRCDIR)/GMG/5_PREP_GMG_global_coarse.f90 $(SRCDIR)/GMG/6_subdomain_infor_mg.f90 $(SRCDIR)/GMG/7_assemble_FVM.f90 $(SRCDIR)/GMG/7_SOLVE_GMG.f90 $(SRCDIR)/GMG/blaslapack_sub_n.f90 $(SRCDIR)/GMG/bubble_sort.f90 $(SRCDIR)/GMG/coarsening_semi.f90 $(SRCDIR)/GMG/coarsening_semi_amg.f90 $(SRCDIR)/GMG/connectivity_coarse.f90 $(SRCDIR)/GMG/mt_amux.f90 $(SRCDIR)/GMG/mt_amuxvr1.f90 $(SRCDIR)/GMG/mt_precond.f90 $(SRCDIR)/GMG/PCG_Dig.f90 $(SRCDIR)/GMG/poly_smooth.f90 $(SRCDIR)/GMG/Domain_infor_FVM_fine.f90 $(SRCDIR)/GMG/Domain_infor_FVM_coarse.f90 $(SRCDIR)/GMG/Domain_infor_FVM_coarsest.f90 $(SRCDIR)/GMG/Relax_GSP.f90 $(SRCDIR)/GMG/Relax_GSx.f90 $(SRCDIR)/GMG/send_receive.f90 $(SRCDIR)/GMG/send_receive_csr.f90 $(SRCDIR)/GMG/send_receive_mt.f90 $(SRCDIR)/GMG/solve_CG.f90 $(SRCDIR)/GMG/solve_CG_P.f90 $(SRCDIR)/GMG/SOLVE_GC.f90 $(SRCDIR)/GMG/stiffness_MG.f90 $(SRCDIR)/GMG/mt_trans.f90 $(SRCDIR)/GMG/stiffness_GC.f90 $(SRCDIR)/GMG/imapGZ_coarse.f90 $(SRCDIR)/GMG/neighboring.f90 $(SRCDIR)/GMG/write_report.f90 $(SRCDIR)/GMG/interpolation.f90 $(SRCDIR)/GMG/Neighbor_node_ARP.f90 $(SRCDIR)/GMG/linear_algebra.f90 

# 4) MPI 래퍼: allreduce·communicate 모두 원본 (C010-1 에서 스텁 → 실물 교체, np=1 은 niut=0 으로 항등)
SSRC = $(SRCDIR)/06_MPI/allreduce_fns.f90 \
       $(SRCDIR)/06_MPI/communicate.f90

FSRC = $(ZSRC) $(GMOD) $(SSRC) $(GSRC)
OBJS = $(patsubst %.f90,$(BUILD)/%.o,$(notdir $(FSRC))) $(BUILD)/padiso.o

VPATH = stub:$(SRCDIR)/06_MPI:$(SRCDIR)/00_Module:$(SRCDIR)/00_Module/d_MPI:$(SRCDIR)/00_Module/j_rv_model:$(SRCDIR)/00_Module/k_interface:$(SRCDIR)/GMG/module:$(SRCDIR)/GMG

# 직렬 컴파일 전제 (목록 순서가 모듈 의존 순서) — -j 사용 금지
objs: $(OBJS)
	@echo "== C006 compile OK: $(words $(OBJS)) objects =="

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(FC) $(FOPTFLAGS) -module $(BUILD) -c $< -o $@

$(BUILD)/padiso.o: $(SRCDIR)/GMG/padiso.f | $(BUILD)
	$(FC) $(FOPTFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)

# 링크 클로저 프로브 (C006): 미정의 외부 심볼 인벤토리용
probe: objs $(BUILD)/link_probe.o
	$(FC) -o $(BUILD)/link_probe $(BUILD)/link_probe.o $(OBJS) $(LIBDIR)

$(BUILD)/link_probe.o: driver/link_probe.f90 | $(BUILD)
	$(FC) $(FOPTFLAGS) -module $(BUILD) -c $< -o $@

# 합성 Poisson 테스트 드라이버 (C007)
driver: objs $(BUILD)/driver_pmg.o
	$(FC) -o $(BUILD)/driver_pmg $(BUILD)/driver_pmg.o $(OBJS) $(LIBDIR)

$(BUILD)/driver_pmg.o: driver/driver_pmg.f90 | $(BUILD)
	$(FC) $(FOPTFLAGS) -module $(BUILD) -c $< -o $@

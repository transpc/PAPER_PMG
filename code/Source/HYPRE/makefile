include ../makefile.in

# =====================================================================
#  HYPRE 서브디렉토리 — CUPID 압력 솔버의 hypre(BiCGSTAB+BoomerAMG) 경로
#
#  두 가지를 한다:
#   (1) hypre 라이브러리 부트스트랩 — hypre-<VER>.tar.gz 를 전개해
#       configure/make/install 까지 수행. 산출물이 이미 있으면 건너뛴다.
#       (최초 1회 수 분 소요. 컨테이너/네이티브 어느 쪽에서도 동일 동작)
#   (2) CUPID 래퍼 Fortran 소스 컴파일 — 다른 서브디렉토리와 동일 규약
#       (.o 는 제자리, .mod 는 ../Modules, 링크는 최상위 makefile 이 담당)
#
#  버전 2.33.0 고정: 컨테이너 컴파일러(icc classic 2021.10)와의 호환·재현성.
#  라이브러리 링크는 ../makefile.in 의 LIBDIR 에 -lHYPRE 로 등재되어 있다.
# =====================================================================
HYPRE_VER    = 2.33.0
HYPRE_TOP    = $(CURDIR)/hypre-$(HYPRE_VER)
HYPRE_LIB    = $(HYPRE_TOP)/install/lib/libHYPRE.a
# 이식성 우선 기본값. HPC 실계산용으로 올리려면:
#   make HYPRE_CFLAGS='-O3 -xHost'
HYPRE_CFLAGS ?= -O2
# hypre 자체 빌드의 병렬도 (로그인 노드 배려용으로 분리)
HYPRE_JOBS   ?= $(p)

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = MD_hypre.f90\
        SOLVE_HYPRE.f90\
        HYPRE_FINALIZE.f90

FSRCV =
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS  = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../Modules

.SUFFIXES: .f90 .o
.F.o :
	$(FC) $(FFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<
.f90.o :
	$(FC) $(FFLAGS) -c $<

default : mk_lib

mk_lib : $(HYPRE_LIB) $(OBJS)

# 라이브러리만 따로 세울 때 (Fortran 래퍼와 무관하게 검증/선빌드)
hypre_lib : $(HYPRE_LIB)

SOLVE_HYPRE.o    : MD_hypre.o
HYPRE_FINALIZE.o : MD_hypre.o

# hypre 라이브러리 — 산출물(.a)이 타깃이라 두 번째 make 부터는 no-op
$(HYPRE_LIB) :
	@echo "== hypre $(HYPRE_VER) 최초 빌드 (수 분 소요, 이후 재사용) =="
	@test -d $(HYPRE_TOP) || tar xzf hypre-$(HYPRE_VER).tar.gz
	cd $(HYPRE_TOP)/src && ./configure CC=mpiicc CFLAGS='$(HYPRE_CFLAGS)' \
	    F77=mpiifort FFLAGS='$(HYPRE_CFLAGS)' \
	    --prefix='$(HYPRE_TOP)/install' --without-openmp --disable-shared \
	    > configure.log 2>&1
	cd $(HYPRE_TOP)/src && make -j $(HYPRE_JOBS) > make.log 2>&1 && \
	    make install > install.log 2>&1
	@test -f $@ || { echo "hypre 빌드 실패 — $(HYPRE_TOP)/src/{configure,make,install}.log 확인"; false; }
	@echo "== hypre OK: $@ =="

# 다른 서브디렉토리와 동일하게 clean 은 no-op (수 분짜리 재빌드 실수 유발 방지).
# hypre 트리까지 지우려면 hypre_distclean 을 명시적으로 부를 것.
clean:

hypre_distclean:
	rm -rf $(HYPRE_TOP)

new : clean default

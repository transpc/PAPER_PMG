#!/usr/bin/env bash
# 공통 환경변수 — 모든 스크립트가 source 한다. (LOOP C001, 계약: PLAN §2-2)
# SIF/METIS_LIB 는 외부에서 export 로 재정의 가능 (체계마다 경로가 다를 수 있음)
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_DIR
# 2026-08-11 환경 이주 (root → sjdo, LOOP C013): SIF 는 hpc2023_ubuntu_prc3.3 (동일 oneAPI 2023.2.1,
# ifort 2021.10.0), METIS 는 SIF 내장 /usr/local/lib/libmetis.so 사용 (호스트 빌드본 불필요)
export SIF="${SIF:-$HOME/00_apptainer/hpc2023_ubuntu_prc3.3.sif}"
export METIS_LIB="${METIS_LIB:-/usr/local/lib}"
export CUPID_SRC="$CODE_DIR/Source"
export CUPID_CASE="${CUPID_CASE:-$CODE_DIR/2_iSMR_ECT_res1}"
export CUPID_NP="${CUPID_NP:-4}"        # mpirun rank 수
export MAKE_JOBS="${MAKE_JOBS:-12}"     # 병렬 make 잡 수

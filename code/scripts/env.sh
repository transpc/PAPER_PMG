#!/usr/bin/env bash
# 공통 환경변수 — 모든 스크립트가 source 한다. (LOOP C001, 계약: PLAN §2-2)
# SIF/METIS_LIB 는 외부에서 export 로 재정의 가능 (체계마다 경로가 다를 수 있음)
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_DIR
export SIF="${SIF:-/root/00_apptainer/hpc23.sif}"
export METIS_LIB="${METIS_LIB:-/root/00_apptainer/metis-5.0.2/build/Linux-x86_64/libmetis}"
export CUPID_SRC="$CODE_DIR/Source"
export CUPID_CASE="${CUPID_CASE:-$CODE_DIR/2_iSMR_ECT_res1}"
export CUPID_NP="${CUPID_NP:-4}"        # mpirun rank 수
export MAKE_JOBS="${MAKE_JOBS:-12}"     # 병렬 make 잡 수

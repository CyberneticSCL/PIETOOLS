#!/bin/bash
# regime.sh BLOCK... -- run the 2026-09-26 testing regime blocks in order.   CC, 09/26/2026
#
# One MATLAB process per MATLAB block (bl_regime.m), so a crash or OOM loses one
# block; bl_regime resumes items already 'ok' in the manifest. Bash-only blocks:
#   B0g  GPU hygiene (idle GPU, sm_89 cubins in the libraries that hold the
#        kernels) + clean re-times of three banked runs
#   B6b  clean 1e-4 re-times of the objective-form dumps whose 09-24 timings
#        overlapped other GPU jobs
# Environment: CUADMM_OUT (Windows path), REGIME_DEADLINE (posix s). A block is
# not started past the deadline or if $CUADMM_OUT/STOP exists.
set -u
PKG="$(cd "$(dirname "$0")/.." && pwd)"
OUTW="${CUADMM_OUT:?set CUADMM_OUT}"
OUTU="$(cygpath -u "$OUTW")"
RG="$OUTU/regime"; mkdir -p "$RG"
DMP="$OUTU/baseline/dumps"
EXE=${CUADMM_EXE:-/home/mpeet/solvers/cuADMM/build/cuadmm_exe}
DL="${REGIME_DEADLINE:?set REGIME_DEADLINE}"
PKGW="$(cygpath -w "$PKG")"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$RG/driver.log"; }
wsldir() { local w; w="$(cygpath -w "$1")"; echo "/mnt/$(echo "${w:0:1}" | tr A-Z a-z)$(echo "${w:2}" | tr '\\' '/')"; }

retime() {   # retime <id> <tol> <cap_s> : one clean run, no dump
  local d="$DMP/$1"; [ -f "$d/blk.txt" ] || { log "retime $1: no dump"; return; }
  local wd; wd="$(wsldir "$d")"
  wsl -e bash -c "LD_LIBRARY_PATH=/usr/lib/wsl/lib timeout $3 $EXE '$wd/' $2 1000000 15 0 > '$wd/retime_$2.log' 2>&1"
  local it sv
  it=$(awk '/^ *[0-9]+ \|/{i=$1} END{print i}' "$d/retime_$2.log")
  sv=$(grep -oP '(?<= solve )[0-9.eE+-]+' "$d/retime_$2.log" | tail -1)
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "${it:-NA}" "${sv:-NA}" "$(grep -c 'converged' "$d/retime_$2.log")" >> "$RG/retime.tsv"
  log "retime $1 tol=$2 iters=${it:-NA} solve=${sv:-NA}s"
}

for blk in "$@"; do
  [ -f "$OUTU/STOP" ] && { log "STOP file present; stopping before $blk"; break; }
  now=$(date +%s); [ "$now" -ge "$DL" ] && { log "deadline reached; stopping before $blk"; break; }
  n=$(wsl -e bash -c 'pgrep -c cuadmm_exe' 2>/dev/null | tr -d '\r'); n=${n:-0}
  if [ "$n" != "0" ]; then log "a cuadmm_exe is running ($n); waiting 60 s"; sleep 60; fi
  case "$blk" in
  B0g)
    log "B0g GPU hygiene"
    wsl -e bash -c "nvidia-smi --query-compute-apps=pid,name --format=csv,noheader" >> "$RG/driver.log" 2>&1
    wsl -e bash -c "cd /home/mpeet/solvers/cuADMM/build && sha256sum cuadmm_exe libcuadmm_lib.so psd_projection/libpsd_lib.so && cuobjdump -lelf libcuadmm_lib.so psd_projection/libpsd_lib.so | sort | uniq -c" >> "$RG/driver.log" 2>&1
    printf 'id\ttol\titers\tsolve_s\tconverged\n' > "$RG/retime.tsv"
    retime stab_rd1_hv 1e-6 120        # banked 09-24: 4023 it, 15.8 s
    retime stab2_rd 1e-4 200           # banked: 1316 it, 27.5 s (09-23 same shape: 12.7 ms/it)
    retime nl_fisher_nobnd 1e-4 400    # banked: 4567 it, 169 s under contention; 75 s clean (memory)
    ;;
  B6b)
    log "B6b clean re-times"
    for id in est_rd1 ctrl_rd1 hinfco_rd1 hinfduco_rd1 hinfdu_rd1 nl_fisher_opt; do
      [ "$(date +%s)" -ge "$DL" ] && break
      retime "$id" 1e-4 120
    done
    ;;
  *)
    log "block $blk start"
    matlab -batch "addpath('$PKGW'); bl_regime('$blk')" > "$RG/$blk.log" 2>&1
    log "block $blk exit $? ($(grep -c 'REGIME done' "$RG/$blk.log") items done, $(grep -c 'REGIME ERR' "$RG/$blk.log") ERR)"
    ;;
  esac
done
log "driver finished"

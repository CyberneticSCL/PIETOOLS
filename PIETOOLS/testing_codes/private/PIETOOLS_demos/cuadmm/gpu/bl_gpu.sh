#!/bin/bash
# bl_gpu.sh -- GPU arm. cuADMM over every dump, tightest-first-that-makes-sense.
#
# LD_LIBRARY_PATH=/usr/lib/wsl/lib is MANDATORY: a distro libcuda.so.595.91.07 in
# /usr/lib/x86_64-linux-gnu shadows the WSL driver and every CUDA call then
# returns "no CUDA-capable device". Without it this arm silently produces nothing.
#
# cuadmm_exe EXITS 0 ON TOTAL CUDA FAILURE, so success is judged by the presence
# of a CUADMM_TIMING line, never by $?. timeout's 124 is the one trusted code.
#
# ESCALATION RULE: 1e-6 is attempted ONLY if 1e-4 converged inside its cap.
# Measured justification: ctrl_rd1 (m=181, which Mosek disposes of in 0.05 s)
# burned 155,900 iterations and the full 300 s cap at 1e-4 without converging.
# Spending a second, longer cap at a TIGHTER tolerance on a case that already
# failed the looser one buys nothing and, across ~31 dumps, is the difference
# between a 3-hour arm and a 10-hour one.
#
# Residuals are NOT parsed here beyond the summary block; capped runs have no
# summary. The per-case logs are kept and post-processed, so a capped run still
# yields its last-iteration residuals.
set -u
EXE=${CUADMM_EXE:-/home/mpeet/solvers/cuADMM/build/cuadmm_exe}   # override with CUADMM_EXE
DMP="$1"; OUT="$2"; CAP4=${3:-120}; CAP6=${4:-300}

if [ ! -f "$OUT" ]; then
  printf 'id\ttol\texit\tcapped\tinit\tsolve\ttotal\titers\tpinf\tdinf\tgap\tpsdX\tpsdS\tmem_init\tmem_solve\n' > "$OUT"
fi

run_one() {
  local d="$1" id="$2" tol="$3" cap="$4"
  local log="$d/cu_${tol}.log"
  LD_LIBRARY_PATH=/usr/lib/wsl/lib timeout "$cap" "$EXE" "$d" "$tol" 2000000 15 1 > "$log" 2>&1
  local rc=$? capped=0; [ $rc -eq 124 ] && capped=1
  for f in X_opt y_opt S_opt; do
    [ -f "$d/$f.txt" ] && mv -f "$d/$f.txt" "$d/${f}_${tol}.txt"
  done
  g() { grep -oP "$1" "$log" 2>/dev/null | tail -1; }
  local init solve total pinf dinf gap psdX psdS mi ms its
  init=$(g '(?<=CUADMM_TIMING init )[0-9.eE+-]+'); solve=$(g '(?<=solve )[0-9.eE+-]+')
  total=$(g '(?<=total )[0-9.eE+-]+')
  pinf=$(g '(?<=primal infeasibility = )[0-9.eE+-]+')
  [ -z "${pinf:-}" ] && pinf=$(g '(?<=primal   infeasibility = )[0-9.eE+-]+')
  dinf=$(g '(?<=dual   infeasibility = )[0-9.eE+-]+')
  gap=$(g '(?<=relative gap         = )[0-9.eE+-]+')
  psdX=$(g '(?<=PSD violation X      = )[-0-9.eE+]+'); psdS=$(g '(?<=PSD violation S      = )[-0-9.eE+]+')
  mi=$(g '(?<=after_init_MiB )[0-9.]+'); ms=$(g '(?<=after_solve_MiB )[0-9.]+')
  its=$(awk '/^ *[0-9]+ \|/{n=$1} END{print n+0}' "$log")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$id" "$tol" "$rc" "$capped" "${init:-}" "${solve:-}" "${total:-}" "$its" \
    "${pinf:-}" "${dinf:-}" "${gap:-}" "${psdX:-}" "${psdS:-}" "${mi:-}" "${ms:-}" >> "$OUT"
  echo "GPU $id $tol rc=$rc capped=$capped solve=${solve:-NA} iters=$its gap=${gap:-NA}"
  [ $capped -eq 0 ] && [ -n "${solve:-}" ]
}

for d in "$DMP"/*/; do
  id=$(basename "$d")
  [ -f "$d/blk.txt" ] || continue
  if grep -qP "^${id}\t1e-4\t" "$OUT" 2>/dev/null; then echo "GPU skip $id"; continue; fi
  if run_one "$d" "$id" 1e-4 "$CAP4"; then
    grep -qP "^${id}\t1e-6\t" "$OUT" 2>/dev/null || run_one "$d" "$id" 1e-6 "$CAP6"
  else
    echo "GPU $id 1e-6 SKIPPED (1e-4 did not converge)"
  fi
done
echo "GPUDONE"

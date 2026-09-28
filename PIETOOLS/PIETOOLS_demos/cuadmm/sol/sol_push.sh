#!/bin/bash
# sol_push.sh -- copy this PIETOOLS working tree and chosen dumps to Sol.   CC, 09/27/2026
#
#   bash sol/sol_push.sh <local CUADMM_OUT> <dump id>...
#   bash sol/sol_push.sh 'C:\Users\mpeet\cuadmm_out\regime_0927' stab_tr1 hinf_rd1_hv
#
# The WHOLE tree goes, not just this folder: cuadmm_path asserts that every
# routine resolves once inside the checkout, and the Sol repo copy is an old
# commit.  The exact state is recorded in harness/VERSION (HEAD plus a hash of
# the uncommitted diff), so a Sol result can be traced to the code that made it.
# Dumps are copied, not rebuilt: identical bytes on both machines is what makes
# the desktop/Sol comparison a comparison.  Also copies regime/*.mat (the
# desktop's Mosek brackets, which Sol cannot compute).  Needs the ASU VPN and
# the ssh key (memory: asu-sol-access-and-environment).  No Duo prompt.
set -euo pipefail
OUTW="${1:?local CUADMM_OUT}"; shift
OUTU="$(cygpath -u "$OUTW")"
PKG="$(cd "$(dirname "$0")/.." && pwd)"          # .../PIETOOLS/PIETOOLS_demos/cuadmm
TREE="$(cd "$PKG/../.." && pwd)"                 # .../PIETOOLS (the checkout root)
SOL=/scratch/mpeet/pietools
ssh sol "mkdir -p $SOL/harness $SOL/harness_out/baseline/dumps $SOL/harness_out/regime $SOL/harness_out/slurm"

ver="$(git -C "$TREE" rev-parse HEAD) dirty-diff-sha256=$(git -C "$TREE" diff HEAD | sha256sum | cut -c1-16)"
echo "pushing tree ($ver)"
tar -C "$(dirname "$TREE")" --exclude=.git --exclude=.claude --exclude='*.asv' \
    -czf - "$(basename "$TREE")" | ssh sol "rm -rf $SOL/harness/PIETOOLS && tar xzf - -C $SOL/harness && echo '$ver' > $SOL/harness/VERSION"
# git on Windows may check shell scripts out with CRLF; bash on Sol would then
# fail on the \r (CC, 09/27/2026)
ssh sol "find $SOL/harness/PIETOOLS -name '*.sh' -o -name '*.slurm' | xargs -r sed -i 's/\r\$//'"

for id in "$@"; do
  d="$OUTU/baseline/dumps"
  [ -f "$d/$id.mat" ] || { echo "no dump $id in $d"; exit 1; }
  echo "pushing dump $id ($(du -sh "$d/$id" | cut -f1))"
  tar -C "$d" -czf - "$id.mat" "${id}_meta.mat" "$id" | ssh sol "tar xzf - -C $SOL/harness_out/baseline/dumps"
done
if ls "$OUTU/regime/"*.mat >/dev/null 2>&1; then
  tar -C "$OUTU/regime" -czf - $(cd "$OUTU/regime" && ls *.mat) | ssh sol "tar xzf - -C $SOL/harness_out/regime"
fi
ssh sol "cat $SOL/harness/VERSION; du -sh $SOL/harness $SOL/harness_out"

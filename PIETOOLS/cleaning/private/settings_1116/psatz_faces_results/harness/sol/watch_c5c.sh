#!/bin/bash
# watch_c5c.sh JOBID [start|end]: poll Sol every 120 s for Slurm job JOBID; print its state
# changes and every new marker line (@@, ##, ===, errors) of its log; exit at a terminal
# state, or in 'start' mode as soon as the build line has appeared. SCRATCH, MMP, 11/15/2026.
J=$1; MODE=${2:-end}
OUT=/scratch/mpeet/pietools/psatz_faces/jobs/c5c_box3d_$J.out
prev_state=""; shown=0
while true; do
  resp=$(ssh -o ConnectTimeout=30 -o BatchMode=yes sol "st=\$(squeue -h -j $J -o %T 2>/dev/null); if [ -z \"\$st\" ]; then st=\$(sacct -n -X -j $J -o State%20 2>/dev/null | head -1 | tr -d ' '); fi; echo STATE=\$st; if [ -f $OUT ]; then grep -E '^(@@|##|===)|rror|Killed|OOM|xceeded|Terminated' $OUT; fi" 2>/dev/null | grep -viE 'post-quantum|store now|openssh.com/pq|^[*]')
  state=$(printf '%s\n' "$resp" | grep '^STATE=' | head -1 | cut -d= -f2)
  lines=$(printf '%s\n' "$resp" | grep -v '^STATE=')
  total=$(printf '%s\n' "$lines" | grep -c .)
  if [ "$state" != "$prev_state" ]; then echo "$(date +%H:%M) job $J state ${state:-unknown}"; prev_state=$state; fi
  if [ "$total" -gt "$shown" ]; then printf '%s\n' "$lines" | grep . | tail -n +$((shown+1)); shown=$total; fi
  case "$state" in COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|DEADLINE) echo "$(date +%H:%M) job $J ended: $state"; exit 0;; esac
  if [ "$MODE" = start ] && printf '%s\n' "$lines" | grep -q 'build .* s: m='; then echo "$(date +%H:%M) job $J: build line seen"; exit 0; fi
  sleep 120
done

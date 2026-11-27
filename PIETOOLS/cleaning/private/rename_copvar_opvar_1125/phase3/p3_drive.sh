#!/bin/bash
# p3_drive.sh TREES BATCHES - 11/27/2026: runs the Phase 3 batches with d3run.m, one MATLAB process
# per batch (each stops taking rows after 17 min; a paused batch resumes in a new process, at most
# 3 processes per batch), trees one after another. Log: drive.log.
S=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3
SW='C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3'
ML="/c/Program Files/MATLAB/R2025b/bin/matlab.exe"
TREES=${1:-"old new"}
BATCHES=${2:-"T1 T2 D P"}
for tree in $TREES; do
  out=$S/out_$tree; mkdir -p "$out"
  for b in $BATCHES; do
    want=$(grep -c -P "^$b\t" "$S/q_$tree.tsv")
    for attempt in 1 2 3; do
      have=0; [ -f "$out/results.tsv" ] && have=$(grep -c -P "^$b\t" "$out/results.tsv")
      [ "$have" -ge "$want" ] && break
      [ -f "$out/STOP" ] && break
      echo "$(date '+%F %T') start $tree $b attempt $attempt ($have of $want done)" >> "$S/drive.log"
      "$ML" -batch "addpath('$SW'); d3run('C:/wt1127/$tree','$SW/q_$tree.tsv','$SW/out_$tree','$b')" >> "$out/$b.console.txt" 2>&1
      echo "$(date '+%F %T') end $tree $b exit $?" >> "$S/drive.log"
    done
  done
done
echo "$(date '+%F %T') all done" >> "$S/drive.log"

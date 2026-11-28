#!/bin/bash
# P3S_DRIVE.SH TREE QUEUE OUT BATCHES - 11/28/2026, Phase 3 of the _sop rename (p3_drive.sh of
# ../../rename_copvar_opvar_1125/phase3 with the tree, queue and output folder as arguments):
# the d3run.m batches BATCHES ("T1 D", ...) of QUEUE on the export TREE (holding PIETOOLS/), one
# MATLAB process per batch; a batch takes no new row after 17 min and resumes in a new process,
# at most 3 processes per batch. Log: drive.log beside OUT.
S=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3s
SW='C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/p3s'
ML="/c/Program Files/MATLAB/R2025b/bin/matlab.exe"
tree=$1; q=$2; out=$3; batches=$4
mkdir -p "$S/$out"
for b in $batches; do
  want=$(grep -c -P "^$b\t" "$S/$q")
  for attempt in 1 2 3; do
    have=0; [ -f "$S/$out/results.tsv" ] && have=$(grep -c -P "^$b\t" "$S/$out/results.tsv")
    [ "$have" -ge "$want" ] && break
    [ -f "$S/$out/STOP" ] && break
    echo "$(date '+%F %T') start $out $b attempt $attempt ($have of $want done)" >> "$S/drive.log"
    "$ML" -batch "addpath('$SW'); d3run('$tree','$SW/$q','$SW/$out','$b')" >> "$S/$out/$b.console.txt" 2>&1
    echo "$(date '+%F %T') end $out $b exit $?" >> "$S/drive.log"
  done
done
echo "$(date '+%F %T') done $out $batches" >> "$S/drive.log"

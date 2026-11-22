#!/bin/bash
# rr_run.sh TAG STARTS STATEDIR MAXSESS - runs ex_search (opts.post) on the
# worktree C:/wt1121r in sessions of at most 17 minutes of builds, resuming
# from the row states, until the run prints DONE or MAXSESS sessions end.
# Each session first calls lb_setup(), which puts the worktree tree on the
# path and asserts single resolution; without it MATLAB starts on the saved
# path of the main checkout (corrected 11/21 03:25: the first queues ran there).
# Settings rerun of 11/21/2026 (PLAN.md). MMP, 11/21/2026.
TAG=$1; STARTS=$2; SD=$3; MAX=$4
H=C:/wt1121r/PIETOOLS/testing_codes/private/baselines/2026-11-18_studies/private
L=C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/cleaning/private/settings_rerun_1121/logs
for i in $(seq 1 "$MAX"); do
  echo "=== $TAG session $i $(date '+%m/%d %H:%M:%S')" >> "$L/$TAG.log"
  matlab -batch "cd('$H'); lb_setup(); fprintf('## which: %s | %s | %s\n',which('PIETOOLS_stability_sop'),which('sossolve'),which('post_verify')); ex_search('$STARTS','$TAG',struct('post',true,'stopmin',17,'statedir','$SD'))" >> "$L/$TAG.log" 2>&1
  if grep -q "## EX_SEARCH $TAG DONE" "$L/$TAG.log"; then
    echo "=== $TAG done $(date '+%m/%d %H:%M:%S')" >> "$L/$TAG.log"; break
  fi
done

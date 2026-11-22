#!/bin/bash
# pc_qrun_b.sh (worktree C:/wt1120b at 41c9210c; st_qrun3.sh for the problem-case run of 11/20: worktree C:/wt1120, pc_run)
# Runs st_run(TAG) for the tags of QFILE, taken one at a time from its first
# line (tags may be appended or reordered while it runs; an empty file waits;
# a line STOP ends the runner). One MATLAB -batch process at a time, started
# only when a slot is free: fewer than 2 lock files in the Fable scratchpad
# matlab_slots and fewer than 2 MATLAB -batch processes on the machine. A
# tag stopped at its session stop time is rerun until done (12 sessions at
# most; lb_run resumes). Holds sister_<TAG>.lock while it runs.
#   -s          second queue: start only while no lock other than sister_*
#               exists, no exclusive lock (sister_*_X.lock) and no request
#               for one (*.want) is pending
#   -x          exclusive tags (memory-heavy): announces sister_<TAG>.want
#               while waiting, holds sister_<TAG>_X.lock while running, so
#               second queues stay out
#   -a TASKOUT  begin once "QUEUE DONE" appears in the file TASKOUT
#   -p PID      begin once process PID (Windows) has exited
SLOTS="C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/matlab_slots"
H="C:/wt1120b/PIETOOLS/testing_codes/private/baselines/2026-11-18_studies/private"
LOGS="C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/cleaning/private/problem_cases_1120/logs_41c9210c"
mkdir -p "$LOGS"
MODE=n; AFTER=""; PIDW=""
while [ $# -gt 1 ]; do
  case "$1" in -s) MODE=s; shift;; -x) MODE=x; shift;; -a) AFTER="$2"; shift 2;; -p) PIDW="$2"; shift 2;; *) break;; esac
done
QF="$1"; QN=$(basename "$QF" .txt)
if [ -n "$AFTER" ]; then until grep -q "QUEUE DONE" "$AFTER" 2>/dev/null; do sleep 60; done; fi
if [ -n "$PIDW" ]; then
  while powershell -NoProfile -Command "if (Get-Process -Id $PIDW -ErrorAction SilentlyContinue) { exit 0 } else { exit 1 }"; do sleep 30; done
fi
echo "@@ runner $QN ($MODE) begins $(date +%H:%M:%S)"
nbatch() {
  powershell -NoProfile -Command "@(Get-CimInstance Win32_Process -Filter \"Name='MATLAB.exe'\" | Where-Object { \$_.CommandLine -match '-batch' -and \$_.ExecutablePath -match 'win64' }).Count" | tr -d '\r '
}
waitslot() {   # $1: tag
  local want="$SLOTS/sister_$1.want"
  [ "$MODE" = "x" ] && echo "want $1 $(date +%H:%M)" > "$want"
  while true; do
    nl=$(ls "$SLOTS"/*.lock 2>/dev/null | wc -l)
    nf=$(ls "$SLOTS"/*.lock 2>/dev/null | grep -v '/sister_' | wc -l)
    nx=$(ls "$SLOTS"/sister_*_X.lock 2>/dev/null | wc -l)
    nw=$(ls "$SLOTS"/*.want 2>/dev/null | grep -v "/sister_$1.want" | wc -l)
    nb=$(nbatch); [ -z "$nb" ] && nb=2
    if [ "$nl" -lt 2 ] && [ "$nb" -lt 2 ]; then
      if [ "$MODE" = "x" ]; then break; fi
      if [ "$MODE" = "n" ] && [ "$nw" -eq 0 ]; then break; fi
      if [ "$MODE" = "s" ] && [ "$nf" -eq 0 ] && [ "$nx" -eq 0 ] && [ "$nw" -eq 0 ]; then break; fi
    fi
    sleep 30
  done
  rm -f "$want"
}
while true; do
  TAG=$(head -1 "$QF" 2>/dev/null | tr -d '\r ')
  if [ -z "$TAG" ]; then sleep 60; continue; fi
  sed -i '1d' "$QF"
  if [ "$TAG" = "STOP" ]; then break; fi
  for ATT in 1 2 3 4 5 6 7 8 9 10 11 12; do
    waitslot "$TAG"
    if [ "$MODE" = "x" ]; then LOCK="$SLOTS/sister_${TAG}_X.lock"; else LOCK="$SLOTS/sister_$TAG.lock"; fi
    echo "sister $TAG $(date +%H:%M)" > "$LOCK"
    LOG="$LOGS/${TAG}_$ATT.log"
    echo "@@ $QN $TAG start #$ATT $(date +%H:%M:%S)"
    matlab -batch "cd('$H'); pc_run('$TAG')" > "$LOG" 2>&1
    rc=$?
    rm -f "$LOCK"
    if grep -q "stop time .* reached" "$LOG"; then
      echo "@@ $QN $TAG paused #$ATT $(date +%H:%M:%S) rc=$rc (stop time; resuming)"
      sleep 5; continue
    fi
    echo "@@ $QN $TAG end #$ATT $(date +%H:%M:%S) rc=$rc $(grep -c '^#####' "$LOG") rows"
    break
  done
  sleep 5
done
echo "@@ QUEUE DONE $QN $(date +%H:%M:%S)"

#!/bin/bash
# run_ab.sh TREE HARNESS_DIR [LIMIT_S] [FUNC]: FUNC(TREE) (default ab_dde) from
# HARNESS_DIR in the second MATLAB slot with a watchdog; log
# libpass/logs/FUNC_TREE.log. The MATLAB-side paths use the C:/ form. MMP, 11/17/2026.
H=/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad
HW=C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad
tree=$1; N=$2; lim=${3:-900}; fn=${4:-ab_dde}
log=$H/libpass/logs/${fn}_$tree.log
while [ $(ls $H/matlab_slots | wc -l) -ge 2 ]; do sleep 30; done
echo "${fn}_$tree $(date +%H:%M)" > $H/matlab_slots/${fn}_$tree.lock
start=$(date +%s)
"/c/Program Files/MATLAB/R2025b/bin/matlab.exe" -batch "cd('$N'); addpath('$HW/libpass/ab'); $fn('$tree')" > $log 2>&1 &
pid=$!; winpid=$(cat /proc/$pid/winpid 2>/dev/null); killed=0
while kill -0 $pid 2>/dev/null; do
  sleep 5
  if (( $(date +%s) - start > lim )); then
    echo "## TIMEOUT killed after $(( $(date +%s)-start )) s" >> $log
    powershell.exe -NoProfile -Command "Get-CimInstance Win32_Process | Where-Object { \$_.ParentProcessId -eq $winpid } | ForEach-Object { Stop-Process -Id \$_.ProcessId -Force }; Stop-Process -Id $winpid -Force" >> $log 2>&1
    killed=1; break
  fi
done
wait $pid 2>/dev/null; rc=$?
rm -f $H/matlab_slots/${fn}_$tree.lock
echo "## END rc=$rc killed=$killed wall=$(( $(date +%s)-start ))s" >> $log

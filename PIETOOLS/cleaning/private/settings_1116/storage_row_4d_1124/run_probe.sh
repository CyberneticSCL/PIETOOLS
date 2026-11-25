#!/bin/bash
# run_probe.sh <function> <label> <maxsec> <maxGB>
# Runs <function>('<label>') in a MATLAB batch under a slot lock; polls the
# working set of the MATLAB processes of this run every 10 s; stops the
# run after <maxsec> s or above <maxGB> GB. Diary: <label>_diary.log.
FN="$1"; LABEL="$2"; MAXSEC="$3"; MAXGB="$4"
SCR="C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/poscopvar4d"
SLOT="C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/matlab_slots/af90_${LABEL}.lock"
echo "$(date) $FN $LABEL" > "$SLOT"
cd "$SCR" || exit 1
rm -f "${LABEL}_diary.log" "${LABEL}_stdout.log"
"/c/Program Files/MATLAB/R2025b/bin/matlab.exe" -batch "${FN}('${LABEL}')" > "${LABEL}_stdout.log" 2>&1 &
start=$(date +%s)
sleep 15
peak=0
while true; do
  el=$(( $(date +%s) - start ))
  mem=$(powershell -NoProfile -Command "\$p=Get-CimInstance Win32_Process -Filter \"Name='MATLAB.exe'\" | Where-Object { \$_.CommandLine -like '*${FN}*' }; if (\$p) { [math]::Round((\$p | Measure-Object -Property WorkingSetSize -Sum).Sum/1GB,2) } else { -1 }")
  echo "t=${el}s mem=${mem}GB"
  if [ "$mem" = "-1" ]; then echo "MATLAB exited at t=${el}s"; break; fi
  peak=$(awk "BEGIN{m=$peak; if ($mem>m) m=$mem; print m}")
  if [ "$el" -ge "$MAXSEC" ] || awk "BEGIN{exit !($mem > $MAXGB)}"; then
    sleep 3
    echo "STOP at t=${el}s mem=${mem}GB"
    powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='MATLAB.exe'\" | Where-Object { \$_.CommandLine -like '*${FN}*' } | ForEach-Object { Stop-Process -Id \$_.ProcessId -Force }"
    break
  fi
  sleep 10
done
rm -f "$SLOT"
echo "peak working set ${peak} GB"
echo "===== diary"
cat "${LABEL}_diary.log"

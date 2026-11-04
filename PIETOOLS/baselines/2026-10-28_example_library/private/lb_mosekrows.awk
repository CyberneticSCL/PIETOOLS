# lb_mosekrows.awk, 10/28/2026: per catalog row of an lb_run log (full.log,
# sop.log), the MOSEK runs between the row's start and end markers.
# Usage: awk -f lb_mosekrows.awk sop.log
# Output, one line per finished row: key tier status | n_sdp, MOSEK optimizer
# seconds summed over the row's SDPs | last SDP: problem status, solution
# status, primal obj, primal violations con/var/barvar, dual obj, dual
# violations con/var/barvar (MOSEK "Interior-point solution summary").
/^##### / && /result/ { if (cur!="") { if (csv) printf "%s,%s,%d,%.2f,%s,%s,%s\n", key, tier, n, t, mc, pro, sol; else printf "%s %s | n_sdp %d  mosek_s %.2f  m %s | %s %s | P %s con %s var %s bar %s | D %s con %s var %s bar %s\n", cur, $4, n, t, mc, pro, sol, pobj, pcon, pvar, pbar, dobj, dcon, dvar, dbar }; cur=""; next }
/^##### / { cur=$2" "$3; key=$2; tier=$3; gsub(/[\[\]]/,"",tier); n=0; t=0; mc="-"; pro=sol=pobj=pcon=pvar=pbar=dobj=dcon=dvar=dbar="-"; next }
cur=="" { next }
/^  Constraints  +:/ { mc=$3 }
/^Interior-point solution summary/ { n++ }
/^  Problem status/  { pro=$4 }
/^  Solution status/ { sol=$4 }
/^  Primal\.  obj:/ { pobj=$3; pcon=$8; pvar=$10; pbar=$12 }
/^  Dual\.    obj:/ { dobj=$3; dcon=$8; dvar=$10; dbar=$12 }
/^  Optimizer  +- +time:/ { t+=$4 }

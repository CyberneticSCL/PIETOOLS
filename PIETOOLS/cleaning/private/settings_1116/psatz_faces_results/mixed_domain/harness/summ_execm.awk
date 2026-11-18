# summ_execm.awk: one table row per '@@ EXECM' line of a diagnostics log. SCRATCH, MMP, 11/16/2026.
/^@@ EXECM/ {
  sys=$3; par=""; sto=""; cfg=""; dw=""; st=""; cert=""; ndv=""; m=""; t=""; ne=""; fe=""; eta=""; slv="";
  for (i=4;i<=NF;i++) {
    if ($i ~ /^par=/) par=substr($i,5); else if ($i ~ /^storage=/) sto=substr($i,9);
    else if ($i ~ /^cfg=/) cfg=substr($i,5); else if ($i ~ /^dw=/) dw=substr($i,4);
    else if ($i ~ /^cert=/) {cert=substr($i,6); st=$(i-1)} else if ($i ~ /^ndv=/) ndv=substr($i,5);
    else if ($i ~ /^m=/) m=substr($i,3); else if ($i ~ /^t=/) t=substr($i,3);
    else if ($i ~ /^numerr=/) ne=substr($i,8); else if ($i ~ /^feas=/) fe=substr($i,6);
    else if ($i ~ /^eta=/) {eta=substr($i,5); slv=$(i+1)}
  }
  printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", sys, par, sto, cfg, dw, st, cert, ndv, m, t, ne, fe, eta
}

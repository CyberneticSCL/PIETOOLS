# summ_pin.awk: one line per (N,set,dw) from the '@@ PIN' lines: verdicts per rung, sizes, max solve time
/^@@ PIN/ {
    for (i = 3; i <= NF; i++) {
        n = index($i, "=");
        if (n > 0) { k = substr($i, 1, n-1); v = substr($i, n+1); kv[k] = v; }
    }
    key = kv["N"] " " kv["set"] " o" kv["o"] " pd" kv["pd"] " dw" kv["dw"] " D" kv["D"] " w" kv["w"];
    if (!(key in seen)) { seen[key] = 1; order[++no] = key; }
    s = kv["st"]; if (s == "+1") c = "+"; else if (s == "-1") c = "-"; else if (s == "+0" || s == "0") c = "0"; else c = "?";
    verd[key] = verd[key] c;
    size[key] = "ndec=" kv["ndec"] " m=" kv["m"];
    t = kv["ts"] + 0; if (t > tmax[key]) tmax[key] = t;
    if ($0 ~ /err|rror/) errs[key] = errs[key] "E";
}
END {
    for (i = 1; i <= no; i++) { k = order[i]; printf "%-46s %-7s %-22s tsmax=%4d %s\n", k, verd[k], size[k], tmax[k], errs[k]; }
}

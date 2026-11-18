# summ_exec.awk: compact table of the '@@ EXEC' and '@@ HINF' lines: key fields only
/^@@ EXEC/ || /^@@ HINF/ {
    for (i = 2; i <= NF; i++) {
        n = index($i, "=");
        if (n > 0) { k = substr($i, 1, n-1); v = substr($i, n+1); kv[k] = v; }
    }
    if ($1 == "@@" && $2 == "EXEC") {
        st = $4; for (i = 5; i <= NF; i++) if ($i !~ /=/ && $i != "" ) st = st;
        printf "%-18s %-10s %-9s o%-3s pd%s dw%s  cert=%s ndv=%-8s m=%-6s t=%s\n", $3, kv["storage"], kv["cfg"], kv["o"], kv["pd"], kv["dw"], kv["cert"], kv["ndv"], kv["m"], kv["t"];
    } else {
        printf "HINF %-7s %-9s o%-3s pd%s dw%s  cert=%s gam=%-10s gap=%-10s ndv=%-8s m=%-6s t=%s\n", kv["storage"], kv["cfg"], kv["o"], kv["pd"], kv["dw"], kv["cert"], kv["gam"], kv["gap"], kv["ndv"], kv["m"], kv["t"];
    }
}

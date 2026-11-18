# summ_heat.awk: compact table of the '@@ HEAT' lines
/^@@ HEAT/ {
    for (i = 3; i <= NF; i++) {
        n = index($i, "=");
        if (n > 0) { k = substr($i, 1, n-1); v = substr($i, n+1); kv[k] = v; }
    }
    printf "%-9s o%-4s pd%s QD=%-6s lo=%-11s hi=%-11s reach=%-8s res=%s nsolve=%-3s t=%-5s nx=%-8s m=%-5s nnz=%s\n", kv["cfg"], kv["o"], kv["pd"], kv["QD"], kv["lo"], kv["hi"], kv["reach"], kv["res"], kv["nsolve"], kv["t"], kv["nx"], kv["m"], kv["nnz"];
}

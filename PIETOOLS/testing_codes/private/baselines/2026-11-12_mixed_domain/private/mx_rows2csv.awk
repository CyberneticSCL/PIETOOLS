# MX_ROWS2CSV.AWK turns the result lines of MX_EXEC, MX_EXEC2 and MX_EXEC3
# in the run logs into CSV rows: one row per executive call.
# Usage: awk -v run=N8 -v tree=t -v tier=light -f mx_rows2csv.awk N8.log
# A result line starts with the system name (S.., C.., D..) and has five
# fields separated by " | ": system, parameters, executive and form; status
# and certified flag; gam, expected value and gap; verdict; time, rows m,
# decision variables, Gram block sizes Ks, degree raise dw and a message.
BEGIN { FS = " [|] " }
/^(S[0-9]+|C[0-9]+|D[0-9]+)_[A-Za-z0-9_]+ +\{/ && NF >= 5 {
    split($1, a, " +");  sys = a[1];  par = a[2];  exe = a[3];  form = a[4];
    split($2, b, " +");  status = b[1];  cert = b[3];
    split($3, c, " +");  gam = c[2];  expect = c[4];  gap = c[6];
    verdict = $4;  gsub(/^ +| +$/, "", verdict);
    rest = $5;
    t = rest;  sub(/s .*/, "", t);  gsub(/^ +/, "", t);
    m = rest;  sub(/.* m +/, "", m);  sub(/ .*/, "", m);
    ndv = rest;  sub(/.* ndv +/, "", ndv);  sub(/ .*/, "", ndv);
    ks = rest;  sub(/.* Ks /, "", ks);  sub(/ dw .*/, "", ks);
    dw = rest;  sub(/.* dw /, "", dw);  msg = dw;  sub(/ .*/, "", dw);
    if (msg ~ / /) { sub(/^[^ ]+ /, "", msg) } else { msg = "" }
    gsub(/"/, "\"\"", par);  gsub(/"/, "\"\"", msg);
    printf "%s,%s,%s,%s,\"%s\",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\"%s\",\"%s\"\n", \
        run, tree, tier, sys, par, exe, form, status, cert, dw, ndv, m, t, gam, expect, gap, verdict, ks, msg
}

# RUNS.AWK: prints FILE:START-END:LENGTH for each run of at least N (default
# 15; set with -v N=...) consecutive code lines that carry no comment. Blank
# lines neither extend nor break a run. A comment line, a %{ %} block, a code
# line with a % comment outside strings, or text after a ... continuation
# breaks the run (with -v B=1 a blank line breaks it too). Strings and
# transposes are told apart as in CODEIDENT.AWK.
BEGIN { if (N == "") N = 15; blk = 0; run = 0; st = 0 }
FNR == 1 && NR > 1 { flush(prevfile); blk = 0 }
{
    prevfile = FILENAME
    s = $0; sub(/\r$/, "", s)
    t = s; gsub(/^[ \t]+|[ \t]+$/, "", t)
    if (t == "%{") { blk++; flush(FILENAME); next }
    if (blk > 0) { if (t == "%}") blk--; next }
    if (t == "") { if (B == 1) flush(FILENAME); next }
    if (substr(t, 1, 1) == "%") { flush(FILENAME); next }
    com = 0; n = length(s); i = 1; prev = ""
    while (i <= n) {
        c = substr(s, i, 1)
        if (c == "%") { com = 1; break }
        if (c == "." && substr(s, i, 3) == "...") {
            r = substr(s, i + 3); gsub(/[ \t]/, "", r); if (r != "") com = 1; break
        }
        if (c == "\"" || (c == "'" && prev !~ /[A-Za-z0-9_)\]}.']/)) {
            j = i + 1
            while (j <= n) {
                d = substr(s, j, 1)
                if (d == c) { if (substr(s, j + 1, 1) == c) { j += 2; continue } break }
                j++
            }
            prev = "\""; i = j + 1; continue
        }
        if (c != " " && c != "\t") prev = c; else prev = " "
        i++
    }
    if (com) { flush(FILENAME); next }
    if (run == 0) st = FNR
    run++; last = FNR
}
END { flush(prevfile) }
function flush(f) {
    if (run >= N) printf "%s:%d-%d:%d\n", f, st, last, run
    run = 0
}

# MX_CHECKS2CSV.AWK turns the check lines of TEST_MIXED_DOMAIN_PDE in a run
# log into CSV rows: run, part, name, status, detail.
# Usage: awk -v run=N7 -f mx_checks2csv.awk N7.log
# A check line is: part (convert, ops, exact, exec, piesim), the name (which
# may hold spaces: system, executive, parameters), the status (pass, FAIL,
# KNOWN, rejected) and the detail.
/^(convert|ops|exact|exec|piesim) / {
    part = $1
    rest = substr($0, length($1) + 1);  sub(/^ +/, "", rest);  sub(/\r$/, "", rest)
    if (!match(rest, / +(pass|FAIL|KNOWN|rejected) +/)) next
    name = substr(rest, 1, RSTART - 1)
    st = substr(rest, RSTART, RLENGTH);  gsub(/ /, "", st)
    detail = substr(rest, RSTART + RLENGTH)
    gsub(/"/, "\"\"", name);  gsub(/"/, "\"\"", detail)
    printf "%s,%s,\"%s\",%s,\"%s\"\n", run, part, name, st, detail
}

# CODEIDENT.AWK: prints the code of a MATLAB file with every comment removed,
# one line per non-empty code line, trailing blanks dropped, indentation kept.
# Comments: text after % outside strings; %{ ... %} blocks; text after a ...
# continuation. A ' is a transpose when the character before it is a letter,
# digit, _, ), ], }, . or ', and starts a string otherwise; "..." strings and
# doubled quotes inside strings are handled.
BEGIN { blk = 0 }
{
    s = $0; sub(/\r$/, "", s)
    t = s; gsub(/^[ \t]+|[ \t]+$/, "", t)
    if (t == "%{") { blk++; next }
    if (blk > 0) { if (t == "%}") blk--; next }
    out = ""; n = length(s); i = 1; prev = ""
    while (i <= n) {
        c = substr(s, i, 1)
        if (c == "%") break
        if (c == "." && substr(s, i, 3) == "...") { out = out "..."; break }
        if (c == "\"") {
            j = i + 1
            while (j <= n) {
                d = substr(s, j, 1)
                if (d == "\"") { if (substr(s, j + 1, 1) == "\"") { j += 2; continue } break }
                j++
            }
            out = out substr(s, i, j - i + 1); prev = "\""; i = j + 1; continue
        }
        if (c == "'") {
            if (prev ~ /[A-Za-z0-9_)\]}.']/) { out = out c; prev = c; i++; continue }
            j = i + 1
            while (j <= n) {
                d = substr(s, j, 1)
                if (d == "'") { if (substr(s, j + 1, 1) == "'") { j += 2; continue } break }
                j++
            }
            out = out substr(s, i, j - i + 1); prev = "\""; i = j + 1; continue
        }
        out = out c
        if (c != " " && c != "\t") prev = c; else prev = " "
        i++
    }
    sub(/[ \t]+$/, "", out)
    if (out != "") print out
}

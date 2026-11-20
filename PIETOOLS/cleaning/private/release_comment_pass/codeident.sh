#!/bin/bash
# CODEIDENT.SH [FILE...]: for each .m file (paths relative to the PIETOOLS
# folder, or absolute), compares its code with comments removed against the
# same file at HEAD. Prints OK, DIFF (with the first differing lines) or NEW
# (not in HEAD). Without arguments: every modified .m file in the working tree.
# Exit status 1 if any file is DIFF.
REPO=/c/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS
P=$REPO/PIETOOLS
HERE=$(cd "$(dirname "$0")" && pwd)
if [ $# -eq 0 ]; then
    set -- $(git -C "$REPO" diff --name-only HEAD -- '*.m' | sed 's#^PIETOOLS/##')
fi
bad=0
for f in "$@"; do
    f=${f#$P/}
    if ! git -C "$REPO" cat-file -e "HEAD:PIETOOLS/$f" 2>/dev/null; then echo "NEW   $f"; continue; fi
    a=$(git -C "$REPO" show "HEAD:PIETOOLS/$f" | awk -f "$HERE/codeident.awk")
    b=$(awk -f "$HERE/codeident.awk" "$P/$f")
    if [ "$a" == "$b" ]; then echo "OK    $f"
    else
        bad=1; echo "DIFF  $f"
        diff <(printf '%s\n' "$a") <(printf '%s\n' "$b") | head -6 | sed 's/^/      /'
    fi
done
exit $bad

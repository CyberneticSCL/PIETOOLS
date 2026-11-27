# Phase 1b – 2025 call forms (D6), 11/27/2026

Maintainer, 11/27/2026: "commit them and continue to phase 1b", then "run autonomously for the
next 12 hours to finish up this copvar/opvar switch" (each phase committed locally on passing its
checks). Committed abfb0e83.

## Change

One test at the top of `poslpivar`, `lpivar` and `mat2opvar`, on `nargin` and the types of the
second and third inputs; it raises `<name>:call2025` with a message that gives the current call.

- `poslpivar`, `lpivar`: a non-empty numeric count array, then nothing (`nargin` 2), or a third
  input that is a non-empty number, empty with two or more counts, a `struct` without field
  `out` (the 2-D degrees dx, dy, d2), or a cell holding a non-empty number (the 1-D degree cell).
- `mat2opvar`: fewer than four inputs; or two or more numeric counts with a third input that is a
  `polynomial`, char, cellstr or non-empty number (the 2025 VARS or DOM).

Without the test each such call fails further down (SPACES missing, or `parse_opvar_spaces`), so
no current call changes. Not caught: `mat2opvar` with a scalar PDIM (2025: R^n → R^n), a pair of
names not ending in `_dum` and a 1 x 2 DOM – also a valid current call, a multiplier on
L2^n[s,t]. One line of help in each header. The three files are new in this release: no log line,
no markers (D1).

## Measured

- Before (export of 6135b775 at `C:/wt1127/pre1b`): 31 of 31 calls in the 2025 forms error –
  11 with "Not enough input arguments." (MATLAB:minrhs), 20 with a message of the space parser
  (`pre1b_6135b775.log`).
- After (the checkout, = abfb0e83): `test_call2025_sop` passes its 54 checks – 31 identifiers,
  and class and size for 23 calls in the current forms; the outputs of the 23 equal those of
  6135b775 under `isequaln` (`post_checkout.log`).
- `checkcode`: no message in the three files on either tree; one ISCLSTR suggestion in
  `mat2opvar` is suppressed with its reason (a 2025 VARS is never a string).
- Cost (`d6_time.m`, `time_1127.txt`): the test alone takes 3.6, 5.0, 8.3 and 14.7 µs for 2, 4,
  8 and 16 spaces and reads no decision variable; the end-to-end medians of the three functions
  in 1-D and 2-D lie within the run-to-run spread of the two trees (`poslpivar` 2-D: 26.6 and
  29.2 ms before, 32.4 and 26.8 ms after).

## Files

`d6_run.m` (path to one tree, resolution check, the test, comparison of the current-form
outputs, `checkcode`), `d6_time.m`, the two logs, `time_1127.txt`.

# Storage row on a space in 4 or more variables – 11/24/2026

Report of the cleanup session (11/24, forwarded on the maintainer's word "send this issue
back to New Posopvar agent to consider"): with every shipped tier, `poscopvar` cannot declare
storage on a space in 4 or more variables. Every executive declares its storage through this
path (`pos_row` of `exec_tools_sop` calls `poscopvar` without DEG).

## Reproduced

`probe4d_mmp.m` (this folder): `poscopvar(lpiprogram(vars,dom,tier),1,{vars})` on
L2[s1,...,sn], n = 3, 4, 5, each of the six tiers. R2025b, `poscopvar.m` at HEAD 4d06fa46.

| n | stripped, light, custom, heavy, veryheavy | extreme |
|---|---|---|
| 3 | OK | OK |
| 4, 5 | `A length row of lpi.pos has 3 entries for a space in n variables` (fit_len) | `Unable to perform assignment ... left side is 1-by-n ... right side is 1-by-3` (the word 'aaa') |

Full output: `before_diary.log`.

## Cause

Entry 3 of each 1 x 3 field of `lpi.pos` serves every space in 3 or more variables (help of
`poscopvar`, DEG OMITTED), but the tier files write it with one value per direction for
exactly 3 (`len_inc [0 0 0]`, `len_com [2 2 2]`, `len_sos` 3 x 3, `cap [Inf Inf Inf]`,
the word 'aaa' of extreme), and `fit_len` accepted a scalar or exactly n entries. `len_sos`
took column 1 on any width mismatch, silently. The rule for more than 3 variables was left
unspecified when `row_storage` was written (11/17–18); nothing else in the executive path has
this limit: the per-dimension rows of `st.sop` and `lpi.ineq.psatz` are read at `min(nv,end)`.

## Fix (`lpi_programming_sopvar/poscopvar.m`, uncommitted)

Rule, now in the help (DEG OMITTED): a row, word or row of codes shorter than the number of
variables extends to the directions beyond it when its entries agree; unequal entries raise,
since nothing assigns a length to the extra directions.

- `fit_len`: one branch, `numel(v) < n && all(v == v(1))` extends; the message reads
  "A row of lpi.pos has k entries for a space in n variables" (was "A length row", since
  words and `len_sos` rows now pass through it).
- Words: `A(i,:) = fit_len(codes of the word, n)`. Rows of codes: each row through
  `fit_len` (was `A = pt`, unchecked: a 3-column row on 4 variables indexed out of bounds).
- `len_sos`, both sites: each of the rows a, b, c through `fit_len` (was: column 1 repeated
  on any width mismatch, unequal columns included).

Annotation: `poscopvar.m` is new since PIETOOLS 2025; by the release decision of 11/19–20 new
files carry no log line and no markers (the bug-fix pass of 11/21 edited this file the same
way). The maintainer may ask for markers.

Cost: comparisons on rows of at most n entries per space; independent of the number of
decision variables.

## Measured after the fix

Same probe (`after_diary.log`); ndv = decision variables of the program after the call.

| n | stripped | light | custom | heavy | veryheavy | extreme |
|---|---|---|---|---|---|---|
| 3 | 4 dv, 0.57 s | 4,225 dv, 0.22 s | 4,225, 0.08 s | 4,225, 0.03 s | 4,225, 0.16 s | 1 dv, 0.10 s |
| 4 | 4 dv, 0.04 s | 66,049 dv, 0.58 s | 66,049, 0.05 s | 66,049, 0.05 s | 66,049, 0.05 s | 1 dv, 0.02 s |
| 5 | 4 dv, 0.03 s | 1,050,625 dv, 44.3 s | 1,050,625, 1.08 s | 1,050,625, 1.16 s | 1,050,625, 1.37 s | 1 dv, 0.14 s |

Rows at n = 4: `[1 1 1 1 1 1; 1 4 4 4 4 256]` (a, then b as the full integral with Gram 256)
for light, custom, heavy, veryheavy; stripped `[...; 1 4 4 4 4 1]`; extreme the multiplier
alone. The 44 s of light at n = 5 is the position-map cache of the 1024-wide shape, built
once per MATLAB session; custom, heavy and veryheavy reuse it.

Tests (`tests_after.log`): the cleanup session's uncommitted
`testing_codes/private/lpi_programming_sopvar/tests/test_settings_prog_sop.m` passes 13 of 13
in 7.7 s, part (13) included, with its file untouched; `test_executives_sop(false)` passes 20
of 20 in 23.8 s.

## The three questions

(a) Equal-entry extension is now the specified meaning of "the last entry for every higher
e" within a row. The alternatives were not taken: scalars in entry 3 of the six tier files
would change six files and still leave the word 'aaa'; a help text alone would leave every
tier unusable above 3 variables. A user writing a row for a 4-variable problem writes 4
entries, which `fit_len` accepts as before.

(b) Yes. `len_sos` rows pass through `fit_len`, so an unequal row of the wrong width raises
instead of taking column 1. On every shipped tier and on the synthetic row of part (11) the
result is unchanged (widths equal to the entry number, or width 1).

(c) Untested. The record holds class-level runs at 4 variables only (`test_posmopvar`, nv 3
and 4, 09/21/2026; an `int_semisep` memory wall at nv = 4 for overlapping spaces). The storage
alone costs 66,049 decision variables at light on 4 variables and 1,050,625 on 5; an
executive would add the converter, the slack of `lpi_ineq_sop` and the equality rows in 4
variables. No 4-D executive was attempted.

## Held for the maintainer

- Commit of `poscopvar.m` (this fix) – nothing committed.
- Whether the file should carry a log line and markers despite the new-file rule.
- Part (13) of the cleanup session's test covers the row of (11) and light; the cleanup
  session intended to extend it to all six tiers once fixed (its file, uncommitted).

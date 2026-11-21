# dS ledger – PIESOS items for the report to Declan, 11/21/2026

Agent dS. Items: C62, C65, C72 (credited DJ); C63, C64, C66–C71, C73, C111 (credited SS, AT, CR).
HEAD 905b0d26 (fa0deda4 plus one MMP commit outside PIESOS). `git diff 78641839 HEAD -- PIESOS` is
empty, so the CANDIDATES.md line numbers hold at HEAD; no uncommitted change in PIESOS/.
Scripts and outputs: `/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/dS/`
(`<item>_repro.m` writes `<item>_repro.txt`; `run_all_dS.m` runs all in one process, 30 s; run
02:37 and 02:38, 11/21/2026, rc 0; `snippets_dS.m` runs the reproduction lines of the report
entries as written – all behave as stated, 02:45). Authors: `git blame -w -C -C -M` and `git log -S` over the file
history (@quadPoly: opvar/ -> ndopvar/ (1c3157f5, SS, 01/20/2026) -> PIESOS/polyopvar/ (6757b6d3, DJ,
11/08/2026)). 6fd03c15 (DJ, 10/15/2026, "Update PIESOS codes.") is the first commit on any ref that
contains each DJ line below (`git log --all -S`); the file headers read "DJ, 10/15/2026: Initial coding".

| Item | Status | Author (true) | Class |
|---|---|---|---|
| C62 | CONFIRMED (run) | DJ | silent wrong result (rare path) |
| C63 | CONFIRMED (run) | SS | error, reachable path |
| C64 | CONFIRMED (run) | CR | error, reachable path |
| C65 | CONFIRMED (run) | DJ | error, rare path |
| C66 | CONFIRMED (run) | AT | error, rare path |
| C67 | CONFIRMED (run) | AT | error, rare path |
| C68 | CONFIRMED (run) | SS; DJ lines 38-53 | silent wrong result (P + c) |
| C69 | CONFIRMED (run) | SS; help MMP | dead code (misleading help) |
| C70 | CONFIRMED (run) | SS | dead code / no effect |
| C71 | CONFIRMED (reading) | SS | dead code / no effect |
| C72 | CONFIRMED (run) | DJ | misleading message |
| C73 | CONFIRMED (reading) | CR | misleading text (license) |
| C111 | CONFIRMED (reading) | CR | misleading message |

## C62 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@intop/isvalid.m:291`, local is_consistent (unchanged since 78641839).
- Author: DJ, 6fd03c15, 10/15/2026 (as PIESOS/@intop/isvalid.m); `git log --all -S'safe_vars = [safe_vars,find(safe_idcs)]'`: 6fd03c15 only (plus two ledger commits of MMP).
- Mechanism (read): from pass 2 on, find(safe_idcs) gives positions in rem_vars, not variable numbers. The truly resolved variable is then missing from safe_vars (resolution delayed by one pass) and a lower-numbered variable is wrongly marked resolved, which removes it from every limit list. Wrong entries only remove limits, so verdicts can only change from false to true: acyclic rows are never rejected (inferred). A cycle is missed only if every node has its walk round the cycle cut by a wrong entry at distance >= 3 added in time: two wrong entries from pass 2 on a cycle of length >= 6, plus a full-domain variable below them and two variables resolved in pass 2 → d >= 9 (inferred).
- Script: C62_repro.m. Measured: row `[0,10, 0,3, 0,5, 0,1, 0,6, 0,7, 0,9, 0,1, 0,2]` (6-cycle s2->s3->s5->s6->s7->s9->s2): topological peel 0, copy with rem_vars(safe_idcs) 0, verbatim copy of HEAD 1, `intop({polynomial(1)},{row},{'s'},[0,1])` ACCEPTED, isvalid(K) = 1. Control 2-cycle `[0,4, 0,3, 0,2]`: rejected, "The intervals of integration are incompatible.". Sweep (HEAD copy vs truth, fixed copy vs truth): 0/4, 0/81, 0/4096 wrong for d = 1, 2, 3 (all rows); 0/20000 for each d = 4..9 (random, rng(1)). canonicalize_intmat does not catch the row (no variable has two variable limits).
- Consequence: silent wrong result on a rare path – the constructor accepts an intmat whose integration limits depend on themselves; a user reaches it only by declaring such an intmat of degree >= 9.
- Note: `testing_codes/private/PIESOS/Testfolder/intop/remove_badints.m:110` repeats the line; Test_remove_badints_fuzz uses di <= 4, below the smallest failing degree.
- Entry: see Report entries.

## C63 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@quadPoly/horzcat.m:51`, `vertcat.m:49`.
- Author: SS – `git log -S'isnumeric(X)'`: 23edd836 (SS, 01/29/2026, "Modified quadPoly class") in both files; re-indented in 0992dfd0 (SS, 02/11/2026).
- Script: C63_repro.m. Measured: `which X` empty; [P,P] dim [1 2], [P;P] dim [2 1]; [P,3], [3,P] → "Unrecognized function or variable 'X'." at horzcat:51; [P;3], [3;P] → same at vertcat:49.
- Consequence: error on a reachable path (any concatenation with a numeric array; the help, horzcat.m:16, supports it).

## C64 – CONFIRMED (run)
- HEAD: `PIESOS/custom/filter_Pcell.m:115`.
- Author: CR, df0172cc, 10/07/2026 (created as polyopvar/Examples/LocalStability/Interpret_LF/filter_Pcell.m; `git log -S'entry.coefficient(:) ~= 0'`).
- Script: C64_repro.m. Measured: polynomial block → positions [1 1;2 2]; numeric block `{[1,1e-9;1e-9,2]}`, threshold 1e-6, with two outputs and with one output → "Dot indexing is not supported for variables of this type." at filter_Pcell:115 (the loop runs regardless of nargout).
- Consequence: error on a reachable path (help line 13: "Numeric blocks are also supported"). No caller in the tree outside cleaning/.

## C65 – CONFIRMED (run)
- HEAD: `PIESOS/internal/otimes_params.m:656`, local otimes_quadPolyCell.
- Author: DJ, 6fd03c15, 10/15/2026 (as PIESOS/misc/otimes_params.m); `git log --all -S`: 6fd03c15 only.
- Script: C65_repro.m. Measured, S the struct form of the constant quadPoly 3 (fields C = {sparse(3)}, Zs, Zt, ns, nt, dim): otimes_params(S,2) → C = 6, dim [1 1]; otimes_params(2,S) → "Both parameters must be specified as objects of type 'quadPoly'." at line 657 (A converted to a quadPoly at line 650, then ~isa(A,'struct') is true); otimes_params(S,'abc') and (S,true) → "Dot indexing is not supported for variables of this type." at line 661; otimes_params(S,1+x) → "Property ns is not gettable".
- Consequence: error on a rare path – (a) a documented pair (help lines 13-14: B any type listed for A) is rejected; (b) invalid B with struct A fails later with a less informative error. PIESOS caller `@intop/mtimes_tensopvar.m:791` passes two structs (read), so only direct calls reach it.

## C66 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@quadPoly/subs.m:73`.
- Author: AT (Talitsky), 17bb1f7d, 02/25/2026 (`git log -S'if size(vars) ~= size(vals)'`).
- Script: C66_repro.m. P = 1 + 2*s2 + 3*s1 + 4*s1*s2, P(1,2) = 16. Measured: {'s1','s2'},[1 2] → 16 (control); {'s1','s2'},[1 2 3] → 16, no error (third value ignored); {'s1'},[1;2] → quadPoly in s2, no error; {'s1','s2'},1 → "Index exceeds the number of array elements. Index must not exceed 1." at subs:305; {'s1','s2'},[1;2] (equal counts) → "number of vars and vals should be the same" at subs:74.
- Consequence: error on a rare path; a count mismatch in the same orientation passes silently; a column of values is rejected with a false message.

## C67 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@quadPoly/reshape.m:1, 57`.
- Author: AT, 914b3721, 02/18/2026 (signature and line 57, `git log -S`).
- Script: C67_repro.m. Measured: reshape(P,4,1) → dim [4 1]; reshape(P,[4,1]) → "Not enough input arguments." at reshape:57; reshape([1 3;2 4],[4 1]) works for a double.
- Consequence: error on a rare path. No fix suggested (the help documents only RESHAPE(A,M,N)).

## C68 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@quadPoly/plus.m:51-88`.
- Author: lines 54-88 SS (5bba3db0 01/16/2026, 23edd836 01/29/2026, e67d118e 02/15/2026, 0992dfd0 02/11/2026); lines 38-53 – the 0 shortcut, the conversion of a double operand (45-50) and the comment "should add error checking (e.g. matching matrix dimensions)" (51) – DJ, 9c8ac083, 06/07/2026 ("support simple but non-standard cases (e.g. addition with 0)"). liftIndex (private) SS 23edd836; DJ f169b4fe (06/26/2026) changed only orientation (i(:), v(:)).
- Script: C68_repro.m. F = [s1 0; 0 s1]. Measured: F + 5 → [5+s1 0; 0 s1] (polynomial sum [s1+5 5; 5 s1+5]; from H.C at s1 = 0.5: [5.5 0;0 0.5], expected [5.5 5;5 0.5]); F + t1 → t1 in entry (1,1) only; 5 + F → "Arrays have incompatible sizes for this operation." at plus:80; F + quadPoly(ones(3)) → "Index exceeds array bounds." at liftIndex:63; quadPoly([1;2]) + quadPoly([3 4]) (fast path) → dim [2 1], size(C) [2 2], C = [4 5;5 6], quadPoly2polynomial → [4;5] (double: [4 5;5 6]).
- Consequence: silent wrong result – P + c for a non-scalar quadPoly P with variables and a nonzero double c (admitted by DJ lines 45-50, lifted with F.dim by SS line 72); before 9c8ac083 a double operand errored at G.ns (inferred from the diff). Whether a PIESOS call reaches P + c was not determined.

## C69 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@quadPoly/subsasgn_1.m:1`.
- Author: SS – created as subsasgn.m in 0992dfd0 (02/11/2026) with `error("...").` (trailing period, a parse error); renamed to subsasgn_1.m and the period removed in dc7378f0 (SS, 02/18/2026, "Rudimentary mtimes for sopvar", R067). Help lines 2-10: MMP, 8e335783 (11/20/2026), stating that assignment throws an error.
- Script: C69_repro.m. Measured: methods('quadPoly') has subsasgn_1, not subsasgn; `which subsasgn(P)` → built-in; `P.dim = [5,5]` → no error, dim [5 5] with a 2 x 1 C; subsasgn_1(P,S) → the intended error.
- Consequence: dead code; the help is misleading. Code outside the class assigns quadPoly properties: `sopvar/constructors/int2sopvar.m:163` (Kmat.nt = ...), `sopvar/@sopvar/inv_1D.m:498-501`, so a working overload named subsasgn would break int2sopvar.

## C70 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@quadPoly/quadPoly.m:98-101` (help lines 65-66 and comment line 97, MMP 8e335783, already say "unused").
- Author: SS – line 98 dc7378f0 (02/18/2026), lines 99-100 e67d118e (02/15/2026).
- Script: C70_repro.m. Measured: Zs = {[2;0;1]} stored as given for presorted 0, 1 and 6 inputs (isequal 1); methods('quadPoly'): clean 0, combine 1; `which clean(P0)` → C:\YALMIP\...\extras\clean.m, and clean(P0) → "Invalid data type. Argument must be numeric, char, or logical.".
- Consequence: dead code; an unsorted basis is accepted although the property comments (lines 40-43) require sorted exponent vectors.

## C71 – CONFIRMED (reading)
- HEAD: `PIESOS/polyopvar/@quadPoly/var_swap.m:125-126, 232`.
- Author: SS, 0992dfd0, 02/11/2026.
- Evidence: PosMap/posMap occur only at lines 125, 126, 219, 229, 232; posMap = cell(1,k) is never filled (help line 229 says so).
- Consequence: dead code or no effect.

## C72 – CONFIRMED (run)
- HEAD: `PIESOS/polyopvar/@tensopvar/isvalid.m:127-128`; the same text at `initialize.m:65-66`.
- Author: DJ, 6fd03c15, 10/15/2026 (both files; `git log --all -S` on the test and the message: 6fd03c15 only).
- Script: C72_repro.m (T = tensopvar(mat2sopvar(1,...)), s in [0,1]). Measured: tensopvar({T}) accepted (factor class tensopvar); isvalid with ops{1,1} = T → 1; tensopvar({5}) → "Factors in tensor-PI operator must be of type 'sopvar'." at initialize:66; isvalid with ops{1,1} = 5 → "... Factors in the tensor-PI operator must be specified as objects of type 'sopvar'.".
- Consequence: misleading message; reached by tensopvar(cell) with an invalid factor, by assignment to ops, and through @tensopmat/isvalid.m:160.

## C73 – CONFIRMED (reading)
- HEAD: `PIESOS/custom/kernel_degree.m:28`.
- Author: CR, ca05d251, 09/07/2026 (block copied from polyopvar/Examples/LocalStability/innerprod_v2.m).
- Evidence: outside cleaning/, 1065 .m files read "version 2 of the License"; the only .m file with a GPL block that does not is kernel_degree.m ("version 3"); the repository LICENSE is GPL "Version 2, June 1991".
- Consequence: misleading text (license).

## C111 – CONFIRMED (reading)
- HEAD: `PIESOS/custom/DMB.m:110`.
- Author: CR, ca05d251, 09/07/2026 (created as polyopvar/Examples/LocalStability/polyopvar_product.m).
- Evidence: "Def. 4 product has an inconsistent output component count." names a definition of a source that neither DMB.m nor PIESOS/documentation/polyopvar.pdf (no numbered definitions, pdftotext) names; all error identifiers read 'polyopvar_product:*', the former file name. The guard follows otimes(repelem(TA,kB,1),repmat(TB,kA,1),[1,2]) of two kA*kB-row operators, so it fires only if otimes returns a wrong row count (read). DMB has no caller outside cleaning/private/dead_code.
- Consequence: misleading message (guard path).

## Incidental (by reading, not items)
- `PIESOS/polyopvar/@tensopvar/initialize.m:62-64` (DJ, 6fd03c15) converts factors of classes 'opvar'/'opvar2d', which are not in the repository; isa(...) is false for every object, so the branch is unreachable.
- `sopvar/@sopvar/inv_1D.m:497` (SS, 0992dfd0), local zeroKernelLike, calls quadPoly() with no input, which quadPoly.m:84-85 rejects – measured: quadPoly() → "quadPoly must have 1 or 6 inputs, nothing else." (snippets_dS.txt) – on the branch of zero kernels (lines 86-101); inv_1D has no caller outside itself and a test; the branch itself was not run.
- The help of subsasgn_1.m (MMP, 8e335783) describes the overload as active (C69).

## Report entries

**C62** `PIESOS/polyopvar/@intop/isvalid.m:291`, local is_consistent (DJ, 6fd03c15, 2026-10-15) –
from the second pass on, `safe_vars = [safe_vars,find(safe_idcs)]` records positions in the reduced list rem_vars, not variable numbers: a resolved variable is missed and another is marked resolved. For the degree-9 row below (s2 -> s3 -> s5 -> s6 -> s7 -> s9 -> s2 is a cycle; s4 and s8, resolved in pass 2, are recorded as s3 and s7) isvalid returns true and the constructor accepts the row; the 2-cycle `[0,4, 0,3, 0,2]` is rejected ("The intervals of integration are incompatible."), and with rem_vars(safe_idcs) so is the row. Silent wrong result on a rare path: a declared intmat whose limits depend on themselves passes the check. Wrong entries only remove limits, so acyclic rows are never rejected, and a missed cycle needs at least 6 variables, degree 9 or more (both inferred); measured: no wrong verdict in all rows for degree <= 3 and in 20000 random rows for each degree 4 to 9. Fix: `rem_vars(safe_idcs)` for `find(safe_idcs)` at line 291 (rem_vars is reduced at line 292); `testing_codes/private/PIESOS/Testfolder/intop/remove_badints.m:110` repeats the line.
Reproduction:
```matlab
row = [0,10, 0,3, 0,5, 0,1, 0,6, 0,7, 0,9, 0,1, 0,2];
K = intop({polynomial(1)},{row},{'s'},[0,1]);    % accepted; expected: error
isvalid(K)                                        % 1
K2 = intop({polynomial(1)},{[0,4, 0,3, 0,2]},{'s'},[0,1]);   % error, as intended
```

**C65** `PIESOS/internal/otimes_params.m:656`, local otimes_quadPolyCell (DJ, 6fd03c15, 2026-10-15) –
the check of B tests `~isa(A,'struct')`. otimes_params(2,S), S the struct form of a quadPoly (help, lines 10-14), is rejected with "Both parameters must be specified as objects of type 'quadPoly'." because A was converted to a quadPoly at line 650, while otimes_params(S,2) returns 2*S; an invalid B with a struct A passes and fails at line 661 ("Dot indexing is not supported for variables of this type." for a char, "Property ns is not gettable" for a polynomial). Error on a rare path: direct calls only; @intop/mtimes_tensopvar.m:791 passes two structs. Fix: `isa(B,'struct')`.
Reproduction:
```matlab
S = struct('C',{{sparse(3)}},'Zs',{{}},'Zt',{{}},'ns',{{}},'nt',{{}},'dim',[1,1]);
C = otimes_params(S,2);  full(C.C{1})    % 6
otimes_params(2,S)                       % error: Both parameters must be ... 'quadPoly'.
otimes_params(S,'abc')                   % error at line 661: Dot indexing ...
```

**C72** `PIESOS/polyopvar/@tensopvar/isvalid.m:127-128`, `initialize.m:65-66` (DJ, 6fd03c15, 2026-10-15) –
both tests accept 'sopvar', 'tensopvar' and 'intop' factors, the messages name only 'sopvar': tensopvar({T}) with a tensopvar T is accepted, tensopvar({5}) errors with "Factors in tensor-PI operator must be of type 'sopvar'." (initialize:66), and isvalid with ops{1,1} = 5 returns "... must be specified as objects of type 'sopvar'.". Misleading message; reached by tensopvar(cell) with an invalid factor, by assignment to ops, and through @tensopmat/isvalid.m:160. Fix: name the three classes in both messages.
Reproduction:
```matlab
T = tensopvar(mat2sopvar(1,struct('in',{{'s'}},'out',{{'s'}}),struct('in',[0,1],'out',[0,1])));
Tn = tensopvar({T});     % accepted
Tb = tensopvar({5});     % error: Factors in tensor-PI operator must be of type 'sopvar'.
```

### PIESOS items of other authors (Declan maintains the folder)

**C63** `PIESOS/polyopvar/@quadPoly/horzcat.m:51`, `vertcat.m:49` (SS, 23edd836, 2026-01-29) –
`elseif isnumeric(X)`: X is never assigned, so [P,3], [3,P], [P;3] and [3;P] error with "Unrecognized function or variable 'X'." Error on a reachable path: concatenation with a numeric array, which the help supports (horzcat.m:16). Fix: `isnumeric(varargin{k})`.
Reproduction:
```matlab
P = quadPoly([1;2],{[0;1]},{},[1,1],{'s1'},{});
[P,3]     % Unrecognized function or variable 'X'. (horzcat:51)
[P;3]     % the same at vertcat:49
```

**C64** `PIESOS/custom/filter_Pcell.m:115` (CR, df0172cc, 2026-10-07) –
the loop for the second output reads entry.coefficient, so a numeric Pcell{1,1}, which the help supports and the first loop filters (lines 89-91), errors with "Dot indexing is not supported for variables of this type.", with one output as well (the loop runs regardless of nargout). Error on a reachable path. Fix: for a numeric block, test `entry ~= 0`.
Reproduction:
```matlab
F = filter_Pcell({[1,1e-9;1e-9,2]},1e-6)   % error at filter_Pcell:115
```

**C66** `PIESOS/polyopvar/@quadPoly/subs.m:73` (AT, 17bb1f7d, 2026-02-25) –
`if size(vars) ~= size(vals)` is true only if both dimensions differ. For P = 1 + 2*s2 + 3*s1 + 4*s1*s2: subs(P,{'s1','s2'},[1 2 3]) returns 16 (the third value ignored), subs(P,{'s1'},[1;2]) returns P(1,s2), subs(P,{'s1','s2'},1) errors at line 305 (index exceeds), and subs(P,{'s1','s2'},[1;2]), equal counts, errors with "number of vars and vals should be the same". Error on a rare path; a count mismatch in the same orientation passes silently. Fix: `numel(vars) ~= numel(vals)`.
Reproduction:
```matlab
P = quadPoly([1;2;3;4],{[0;1],[0;1]},{},[1,1],{'s1','s2'},{});
R = subs(P,{'s1','s2'},[1 2 3]);  full(R.C)   % 16, no error
subs(P,{'s1','s2'},[1;2])                     % error: number of vars and vals should be the same
```

**C67** `PIESOS/polyopvar/@quadPoly/reshape.m:1, 57` (AT, 914b3721, 2026-02-18) –
only RESHAPE(A,M,N) is accepted: reshape(P,[4,1]) errors at line 57 with "Not enough input arguments.", where reshape(P,4,1) works and reshape(A,[4 1]) works for a double array. Error on a rare path.
Reproduction:
```matlab
P = quadPoly(sparse([1 3;2 4]),{},{},[2,2],{},{});
reshape(P,4,1);      % works
reshape(P,[4,1])     % Not enough input arguments. (reshape:57)
```

**C68** `PIESOS/polyopvar/@quadPoly/plus.m:51-88` (SS, lines 54-88, 5bba3db0 to e67d118e, 2026-01-16 to 2026-02-15; DJ, lines 38-53, 9c8ac083, 2026-06-07: a double operand converted to a 1 x 1 quadPoly, and the comment of line 51) –
no dimension check; the sum takes F.dim, and liftIndex lifts both operands with F.dim. For F = [s1 0; 0 s1]: F + 5 returns [5+s1 0; 0 s1] (expected [s1+5 5; 5 s1+5]); F + t1 adds t1 to entry (1,1) only; 5 + F errors at line 80 ("Arrays have incompatible sizes for this operation."); a 3 x 3 G errors in liftIndex ("Index exceeds array bounds."); quadPoly([1;2]) + quadPoly([3 4]) returns dim [2 1] with a 2 x 2 coefficient matrix, read back by quadPoly2polynomial as [4;5]. Silent wrong result: P + c for a non-scalar quadPoly P with variables and a nonzero double c. Fix: the dimension check that line 51 asks for, with a 1 x 1 operand expanded first.
Reproduction:
```matlab
C = sparse(4,2);  C(2,1) = 1;  C(4,2) = 1;
F = quadPoly(C,{[0;1]},{},[2,2],{'s1'},{});      % [s1 0; 0 s1]
quadPoly.quadPoly2polynomial(F + 5)             % [5+s1 0; 0 s1]; expected [s1+5 5; 5 s1+5]
5 + F                                           % error: Arrays have incompatible sizes
```

**C69** `PIESOS/polyopvar/@quadPoly/subsasgn_1.m:1` (SS: subsasgn.m in 0992dfd0, 2026-02-11, renamed in dc7378f0, 2026-02-18) –
the file defines `function P = subsasgn(P,S)`, and MATLAB names a method by its file, so the class has subsasgn_1 and no subsasgn: `P.dim = [5,5]` succeeds and leaves dim [5 5] with a 2 x 1 coefficient matrix, while the help (lines 2-10, added 11/20/2026) says assignment throws an error. Dead code, misleading help. sopvar/constructors/int2sopvar.m:163 and sopvar/@sopvar/inv_1D.m:498-501 assign quadPoly properties from outside the class, so an overload named subsasgn would break them. Fix: delete the file.
Reproduction:
```matlab
P = quadPoly([1;2],{[0;1]},{},[1,1],{'s1'},{});
any(strcmp(methods('quadPoly'),'subsasgn'))   % 0 (subsasgn_1: 1)
P.dim = [5,5];                                % no error
```

**C70** `PIESOS/polyopvar/@quadPoly/quadPoly.m:98-101` (SS, dc7378f0, 2026-02-18; lines 99-100 e67d118e, 2026-02-15) –
`if ~presorted && 0` never runs, so presorted has no effect (the help says so): an unsorted basis Zs = {[2;0;1]} is stored as given for presorted 0 and 1, although the property comments (lines 40-43) require sorted exponent vectors. The branch would call clean(obj); quadPoly has no clean method, and `clean` resolves to YALMIP extras/clean.m on the test machine, which errors for a quadPoly. Dead code. Fix: delete lines 97-101.
Reproduction: by reading: `~presorted && 0` is false for every presorted; methods('quadPoly') lists no clean.

**C71** `PIESOS/polyopvar/@quadPoly/var_swap.m:125-126, 232` (SS, 0992dfd0, 2026-02-11) –
sPosMap and tPosMap are assigned and never read; rebuildBasisFromDegrees returns posMap = cell(1,k), never filled. Dead code. Fix: drop the second outputs.
Reproduction: by reading: PosMap/posMap occur only at lines 125, 126, 219, 229 and 232.

**C73** `PIESOS/custom/kernel_degree.m:28` (CR, ca05d251, 2026-09-07) –
the license block says "version 3 of the License"; the other 1065 .m files outside cleaning/ say version 2, as does the repository LICENSE (GPL Version 2). Misleading text (license). Fix: "version 2".
Reproduction: by reading: the only .m file outside cleaning/ with "version 3 of the License".

**C111** `PIESOS/custom/DMB.m:110` (CR, ca05d251, 2026-09-07) –
"Def. 4 product has an inconsistent output component count." cites a definition of a source that neither the file nor PIESOS/documentation/polyopvar.pdf names; the identifiers 'polyopvar_product:*' carry the former file name. The guard fires only if otimes returns a wrong row count. Misleading message. Fix: name the product as the help does (the tensor product of A and B).
Reproduction: by reading: DMB.m:108-111.

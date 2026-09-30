# sopvar/private/dead_code

Dead code of the sopvar classes that other developers wrote, kept here rather than deleted
(MMP, 09/29/2026). A transparency audit of 09/29/2026 found that nothing live calls these files. The
check covered a whole-repository grep of code and strings, `feval`/`str2func`/`eval` that could build
the names, and the MATLAB profiles of 1-D, 2-D and 3-D program builds. Dead files that MMP wrote
(`@sopvar/randsopvar.m`, `@sopvar/test_MP.m`) were deleted; git history holds them.

**Not on the MATLAB path, by construction.** `genpath` (and so `pietools_path_update`) skips folders
named `private` and everything below them. Private-function resolution sees only files directly
inside a `private` folder, never its subfolders, so nothing here is callable from `sopvar/`. The
subfolders avoid `@` and `private` in their names, so adding this folder to the path by hand would
not create a second definition of a class or a private function. To use a file, copy it to where it
came from.

Each file is byte-identical to its version at commit 55e8d42c. "Created" is the author and date of
the first commit that added the file, traced with `git log --follow`. All files were last changed in
2135d87f (08/30/2026, the move of sopvar to a top-level folder) unless stated.

| file (here) | original location | created | why dead |
|---|---|---|---|
| sopvar_class/sopvar_old.m | @sopvar/ | Sachin Shivakumar 2026-01-19 (header "SS,MP, 01/15/2026") | quadPoly-era constructor; listed by `methods` but cannot be called |
| sopvar_class/ctranspose_old.m | @sopvar/ | Sachin Shivakumar 2026-01-19 (header "MMP, SS - 1_16_2026") | quadPoly-era adjoint; `ctranspose_old(P)` failed on a sopvar (field `vars_in`) |
| sopvar_class/randOpvar.m | @sopvar/ | Sachin Shivakumar 2026-01-20 | no caller |
| sopvar_class/mtimestest.m | @sopvar/ | Sachin Shivakumar 2026-05-10 | script in the class folder, no assertions |
| sopvar_class/test_script.m | @sopvar/ | Sachin Shivakumar 2026-01-20 | script in the class folder, no assertions |
| sopvar_class_private/int_2b.m | @sopvar/private/ | Sachin Shivakumar 2026-03-18 | 9-input variant; the local `int_2b` in `@sopvar/mtimes.m` takes precedence; only `termCompose` called it |
| sopvar_class_private/leftPermuteVec.m | @sopvar/private/ | Sachin Shivakumar 2026-02-11 | no caller |
| sopvar_class_private/leftShiftMonomials_old.m | @sopvar/private/ | Sachin Shivakumar 2026-04-15 | superseded by `leftShiftMonomials_SS`; its errors carry that file's label |
| sopvar_class_private/leftshiftMonomoials_AT.m | @sopvar/private/ | Talitsky 2026-03-17 | reads quadPoly fields a sparse parameter lacks |
| sopvar_class_private/rightshiftMonomials_AT.m | @sopvar/private/ | Talitsky 2026-03-17 | reads quadPoly fields a sparse parameter lacks |
| sopvar_class_private/mapAlphaBetaToGamma.m | @sopvar/private/ | Sachin Shivakumar 2026-02-11 | no caller |
| sopvar_class_private/monomial_outer.m | @sopvar/private/ | Sachin Shivakumar 2026-01-29 | no caller |
| sopvar_class_private/monomial_shift_left.m | @sopvar/private/ | Sachin Shivakumar 2026-06-26 | no caller |
| sopvar_class_private/monomial_shift_right.m | @sopvar/private/ | Sachin Shivakumar 2026-06-26 | no caller |
| sopvar_class_private/termCompose.m | @sopvar/private/ | Sachin Shivakumar 2026-02-11 | reads quadPoly fields; no caller |
| sopvar_class_private/termCompose_old.m | @sopvar/private/ | Sachin Shivakumar 2026-02-18 | reads quadPoly fields; no caller |
| sdopvar_class/mtimes_AT.m | @sdopvar/ | DanBraghini 2026-08-17 (header entries AT 09/07, MMP 09/09-09/10) | uncalled variant of `mtimes`, 574 of its 8-line windows duplicated there; last changed 6a5021f3 (09/11/2026) |
| sdopvar_class/plus_decparam.m | @sdopvar/ | Sachin Shivakumar 2026-08-10 | no caller |
| sdopvar_class/plus_decparam_batch.m | @sdopvar/ | Sachin Shivakumar 2026-08-17 | no caller; last changed 0b6c2268 (09/07/2026) |
| sdopvar_class/minus_decparam.m | @sdopvar/ | Sachin Shivakumar 2026-08-10 | calls `plus_dpvar`, which exists nowhere |
| sdopvar_class/lrmultiply_batch.m | @sdopvar/ | Sachin Shivakumar 2026-08-17 | no caller; failed on an sdopvar (`cellfun:NotACell`) |
| sdopvar_class/rand_sdopvar.m | @sdopvar/ | Declan 2026-06-07 | stale copy: no caller passes an sdopvar, so every call reaches `Testfolder/sdopvar/rand_sdopvar.m` |
| sdopvar_class/test_script.m | @sdopvar/ | DanBraghini 2026-08-10 | script in the class folder, no assertions |
| sdopvar_class_private/MatrixMultiply.m | @sdopvar/private/ | DanBraghini 2026-08-10 | its two calls in `@sdopvar/mtimes.m` are commented out |
| sdopvar_class_private/CombineDecisionBasis.m | @sdopvar/private/ | DanBraghini 2026-08-10 | callers moved to `sync_basis` on 09/07/2026; last changed 1b353e15 (08/30/2026) |
| sdopvar_class_private/UnionBasisMonomials.m | @sdopvar/private/ | Declan 2026-06-08 (header "MMP, SS, AT, DJ, 06/08/2026") | last caller, `sync_basis`, now reads positions directly (09/29/2026); `@sopvar/private` keeps its own live copy |
| Testfolder/test_leftshift_monomials_sopvar.m | sopvar/Testfolder/ | Talitsky 2026-07-31 | calls class-private functions, which a test cannot see |
| Testfolder/example_of_leftshiftmonomials_AT.m | sopvar/Testfolder/ | Talitsky 2026-04-01 | calls class-private functions, which a test cannot see |
| Testfolder/leftshiftMonomoials_ATv2.m | sopvar/Testfolder/ | Talitsky 2026-03-17 | called only by the example above |
| Testfolder/rightshiftMonomials_ATv2.m | sopvar/Testfolder/ | Talitsky 2026-04-14 | called only by the example above |

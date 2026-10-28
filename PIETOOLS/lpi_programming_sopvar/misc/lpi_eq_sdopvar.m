function prog = lpi_eq_sdopvar(prog,P,opts)                                 % MMP, 09/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_EQ_SDOPVAR(PROG,P,OPTS) adds to the LPI program PROG the
% equality constraints P == 0 for the 'sdopvar' decision operator
%
%   (P x)(s) = sum_gamma int_{s'} K_gamma(s,s') I_gamma(s-s') x(s') ds',
%   K_gamma(s,s') = kron(I_m,ZL(s))'*C_gamma(d)*kron(I_n,ZR(s')),
%   vec(C_gamma(d)) = A_gamma + B_gamma'*d,
%
% gamma in {1,2,3}^n3 over the n3 variables common to the input and output
% spaces (1 a multiplier, a factor delta(s_k-s_k'); 2 a lower and 3 an
% upper integral), m x n = P.dims, d the decision variables P.Zd: for every
% parameter cell gamma, the equations A_gamma(c) + B_gamma(:,c)'*d = 0 at
% each position c of C_gamma that the canonical multiplier form allows
% (the columns whose ZR monomial has degree 0 in every multiplier
% direction of gamma), all-zero equations omitted. In canonical multiplier
% form, which every 'sdopvar' object has, these equations hold exactly
% when P is the zero operator.
%
% INPUT
% - prog:   'struct', an LPI program (lpiprogram or lpiprogram_sop) whose
%           prog.decvartable holds every name of P.Zd;
% - P:      'sdopvar' object; an empty or all-zero parameter adds no
%           equation;
% - opts:   (optional) 'char', 'symmetric' for a P that is self-adjoint for
%           every value of d: of each pair of parameter cells mapped to each
%           other by exchanging lower (2) and upper (3) integrals in every
%           direction, only the one of lower linear index is constrained
%           (default []: every cell).
%
% OUTPUT
% - prog:   'struct', PROG with the equations appended to prog.expr as
%           'eq' entries, in the order of the parameter cells and, within a
%           cell, of the positions c.
%
% NOTES
% P of a class other than 'sdopvar' (a 'sopvar', 'polynomial' or 'double'
% included), P not in canonical multiplier form, a name of P.Zd absent from
% prog.decvartable, and OPTS other than 'symmetric' are errors. To impose
% P == Q, call LPI_EQ_SDOPVAR(PROG,P-Q).
%
% See also LPI_EQ, LPI_EQ_NDOPVAR, POSSOPVAR, SOSEQ, SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_eq_sdopvar
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% MMP, 08/28/2026: Initial coding
% MMP, 08/29/2026: Drop the multiplier contraction. The canonical multiplier
%                  form is now a class invariant, enforced by the 'sdopvar'
%                  constructor, so the stored coefficients are determined by
%                  the operator and can be constrained directly. What remains
%                  is a selection of the positions the form allows to be
%                  nonzero, which is pure indexing.
% MMP, 09/26/2026: Name handling in O(q) per call removed from the hot path.
%                  (1) Use P.Zd and prog.decvartable as they are when they are
%                  already cellstr, instead of a cellstr(string(.)) round trip
%                  over q names each; (2) optional 4th input lets
%                  'lpi_eq_cdopvar' skip the membership check for blocks whose
%                  list it has already verified; (3) hand soseq only the
%                  decision variables a parameter uses, so sosconstr's
%                  combine/compress no longer sort all q names per call.
%                  Measured, 2-D container Hinf build (q = 7.6e5, 6 calls,
%                  14 soseq): the check and conversions took 3.7 s and
%                  combine/compress 1.6 s; same program, bit for bit.
% MMP, 09/26/2026: Membership check by one unique over [table; names]
%                  instead of ismember(names,table). The reverse scan
%                  getequation took today (ismember(table,names)) does not
%                  help here: the names are half the table (3.8e5 of 7.6e5),
%                  so both sides are large; measured 261 vs 275 ms. The
%                  unique route, same process vs HEAD: 2-D Hinf 411 -> 212
%                  ms, q = 1e6 over 3 variables 461 -> 191 ms; peak +8 MB
%                  at 1.1e6 names (6.5 -> 14.5 MB). A single name keeps
%                  ismember (a strcmp there). Same verdict, same error.
% MMP, 09/26/2026: One soseq per ~q nonzeros instead of one per parameter.
%                  Every soseq scans all q names of prog.decvartable in
%                  getequation (ismember(decvartable,dvarname)), and each
%                  call hashes its names again. Measured, 2-D container Hinf
%                  (q = 7.6e5): 14 calls, 28 ms fixed each (1-entry soseq),
%                  738k names hashed for a 3.8e5 list; soseq 1.74 of the
%                  2.29 s of S6b. The parameters' columns are now collected
%                  in the old call order and imposed as row dpvars joined up
%                  to q nonzeros each (impose_rows). A second output returns
%                  them unimposed and P = that output imposes them, so
%                  'lpi_eq_cdopvar' joins all its blocks. The cap, not one
%                  soseq for everything: soseq holds ~9 copies of its input,
%                  so one uncapped call peaked at 3.4x the per-parameter
%                  memory. Measured 2-D Hinf, q = 7.6e5 / 2.3e6, old -> new:
%                  2.2 -> 1.1 s / 9.5 -> 5.0 s; peak 90 -> 138 MB / 291 ->
%                  447 MB (stock lpi_eq_2d 1.5 s, 174 MB / 6.3 s, 540 MB).
%                  Item (3) of the first 09/26 entry now applies per soseq:
%                  the names used by any of its parameters. prog.expr has
%                  fewer, longer entries; the At, b, c, K sossolve builds
%                  are bit-identical, row order included. Cost: soseq calls
%                  ~ nnz/q + name-list changes, instead of the parameter
%                  count (3^n3 per block, multiplicative in dimension);
%                  memory ~9 x 16 B x min(q + largest parameter, nnz) per
%                  soseq plus one 16 B/nnz copy of the collected blocks.
%                  Measured in 1-D and 2-D only.
% MMP, 09/30/2026: The gamma of parameter k and the cell of its adjoint
%                  (direction 1 fastest; Sec. 1 and 4 of
%                  sopvar_implementation_notes.pdf) come from the
%                  shared 'gamma_of_cell' and 'cell_of_gamma'
%                  (sopvar/misc/conventions) instead of ind2sub/sub2ind over
%                  sz_C, a map several files wrote locally. Same integers;
%                  3.5 us less per parameter (4.2 -> 0.75 us per call,
%                  measured standalone).
% MMP, 09/30/2026: One routine per mode. The signature is the public one
%                  again, prog = lpi_eq_sdopvar(prog,P,opts): collect the
%                  rows ('collect_eq_rows'), impose them ('impose_eq_rows'),
%                  both moved verbatim to lpis_sopvar/private. A second
%                  output used to mean "collect, do not impose" and a struct
%                  P "impose only", a switch invisible at the call site in
%                  'lpi_eq_cdopvar', which now calls the two directly. So
%                  the 09/26/2026 entries' 4th input, second output and
%                  batch form (and their help, 13 lines stamped 09/26/2026)
%                  are no longer here; the check, the cap and the row order
%                  they describe hold in the private routines. Same program
%                  bit for bit (tr_ab, 91 SDP leaves over w1-w3).
% MMP, 10/02/2026: NOTES: the joined rows are written by 'lpi_soseq', not
%                  soseq (see 'impose_eq_rows'); no code change here. The
%                  second 09/26/2026 entry's cost model (each soseq scans
%                  all q names; ~9 copies per soseq) and the 09/30/2026
%                  "soseq calls are unchanged" now describe only the soseq
%                  fallback; by default each batch costs O(nnz) and ~80 B
%                  per nonzero (lpi_soseq), the batches unchanged.
% MMP, 10/16/2026: moved, with private/, from sopvar/lpis_sopvar/ to
%                  lpi_programming_sopvar/; 'lpis_sopvar' in the 09/30/2026
%                  entry and BEGIN block now reads lpi_programming_sopvar.
%                  The help text names the new folder. Code unchanged.
% MMP, 10/21/2026: moved to lpi_programming_sopvar/misc/ and the private
%                  helpers to misc/private/: the folder root keeps the user
%                  functions only. The help text and the locations above
%                  read accordingly. Code unchanged.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), step comments in the body; code
%                  unchanged. The deleted header text was the comparison
%                  with 'lpi_eq', the argument that constraining the
%                  canonical coefficients is exact, and notes, with their
%                  (was) lines, on the joined prog.expr entries and on the
%                  private routines that do the work. The (was) lines of
%                  the earlier signatures moved from under the function
%                  line to here, so that help starts with the description:
% function [prog,Eq] = lpi_eq_sdopvar(prog,P,opts,dvars_checked)            % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% function prog = lpi_eq_sdopvar(prog,P,opts,dvars_checked)                 % MMP, 09/26/2026 (was)
% function prog = lpi_eq_sdopvar(prog,P,opts)                               % MMP, 09/26/2026 (was)


% BEGIN MMP, 09/30/2026: the three modes split. Deleted here, moved verbatim
% to lpis_sopvar/private: the input checks and the collection loop (lines
% stamped MMP 08/28, 08/29, 09/26 and 09/30/2026) with the local function
% canonical_positions, to 'collect_eq_rows'; the local function impose_rows
% (MMP, 09/26/2026), to 'impose_eq_rows'. Deleted outright: the batch-form
% branch for a struct P and the nargout switch (MMP, 09/26/2026), through
% which 'lpi_eq_cdopvar' reached the two halves; it now calls them itself.
% Collect, then impose: the rows and the soseq calls are unchanged.
if nargin<3,    opts = [];      end                                         % MMP, 09/30/2026
% The equations of P, every name of P.Zd checked against prog
Eq = collect_eq_rows(prog,P,opts,false);                                    % MMP, 09/30/2026
% Appended to prog.expr
prog = impose_eq_rows(prog,Eq);                                             % MMP, 09/30/2026
% END MMP, 09/30/2026

end

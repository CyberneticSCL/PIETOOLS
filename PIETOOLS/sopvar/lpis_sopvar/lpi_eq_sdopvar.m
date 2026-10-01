function prog = lpi_eq_sdopvar(prog,P,opts)                                 % MMP, 09/30/2026
% function [prog,Eq] = lpi_eq_sdopvar(prog,P,opts,dvars_checked)            % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% function prog = lpi_eq_sdopvar(prog,P,opts,dvars_checked)                 % MMP, 09/26/2026 (was)
% function prog = lpi_eq_sdopvar(prog,P,opts)                               % MMP, 09/26/2026 (was)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_EQ_SDOPVAR(PROG,P) takes an LPI optimization program structure
% 'prog' and an 'sdopvar' decision variable P, and adds equality constraints
% enforcing P==0. It is the 'sdopvar' counterpart of 'lpi_eq', and applies in
% particular to the operators returned by 'possopvar'.
%
% An 'sdopvar' object represents the operator
%
%   (P x)(s) = sum_gamma int_{s'} K_gamma(s,s') I_gamma(s-s') x(s') ds'
%
%   K_gamma(s,s') = (I_m o ZL(s))^T unvec(A_gamma + B_gamma'*d) (I_n o ZR(s'))
%
% for gamma in {1,2,3}^n3, where index 1 denotes a multiplier (a factor
% delta(s_k-s_k')), 2 a lower integral and 3 an upper integral.
%
% Since the 3^n3 indicator patterns I_gamma act on disjoint regions, the
% operator vanishes exactly when every term vanishes, and a term vanishes
% exactly when its coefficients do. The second half of that relies on the
% canonical multiplier form, which every 'sdopvar' satisfies: in a direction
% where gamma_k=1 the factor delta(s_k-s_k') identifies s_k with s_k', so
% degrees could otherwise be moved freely between ZL and ZR without changing
% the operator, and constraining the stored coefficients would be strictly
% stronger than constraining the operator. The canonical form pins that
% freedom by keeping all of the degree on the left, so the coefficients are
% determined by the operator and constraining them is exactly right. See
% 'canonicalize_multiplier' and the class documentation.
%
% This routine therefore constrains the stored coefficients directly, after
% dropping the positions the canonical form requires to be zero.
%
% INPUT
% - prog:   'struct' specifying an LPI program structure to modify (see also
%           'lpiprogram'). All decision variables of P must already appear in
%           'prog.decvartable';
% - P:      'sdopvar' object of which to enforce P==0. Parameters that are
%           empty or identically zero generate no constraints;
% - opts:   (optional) set opts = 'symmetric' if P is known to be
%           self-adjoint, to constrain only one parameter per adjoint pair.
%           Taking the adjoint exchanges lower and upper integrals in every
%           direction, so the parameters pair up under the index map
%           1->1, 2->3, 3->2;
%
% OUTPUT
% - prog:   Same LPI program structure as the input, but with the added
%           constraints enforcing P==0.
%
% NOTES
% To enforce P==Q, call LPI_EQ_SDOPVAR(PROG,P-Q); '@sdopvar/plus' aligns the
% monomial and decision variable bases of the two operators.
%
% The constraints of one call are joined into as few 'soseq' calls as a cap % MMP, 09/26/2026
% of numel(prog.decvartable) nonzeros each allows, so prog.expr gets a few  % MMP, 09/26/2026
% entries per call, not one per parameter. The rows sossolve assembles from % MMP, 09/26/2026
% prog.expr (At and b, in order) are the same either way.                   % MMP, 09/26/2026
%
% The work is two private routines of lpis_sopvar: 'collect_eq_rows' builds % MMP, 09/30/2026
% the rows without imposing them and 'impose_eq_rows' imposes them.         % MMP, 09/30/2026
% 'lpi_eq_cdopvar' calls the two directly, to join the rows of all blocks   % MMP, 09/30/2026
% that share a decision variable list.                                      % MMP, 09/30/2026
%
% See also LPI_EQ, LPI_EQ_NDOPVAR, POSSOPVAR, SOSEQ, SDOPVAR.

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


% BEGIN MMP, 09/30/2026: the three modes split. Deleted here, moved verbatim
% to lpis_sopvar/private: the input checks and the collection loop (lines
% stamped MMP 08/28, 08/29, 09/26 and 09/30/2026) with the local function
% canonical_positions, to 'collect_eq_rows'; the local function impose_rows
% (MMP, 09/26/2026), to 'impose_eq_rows'. Deleted outright: the batch-form
% branch for a struct P and the nargout switch (MMP, 09/26/2026), through
% which 'lpi_eq_cdopvar' reached the two halves; it now calls them itself.
% Collect, then impose: the rows and the soseq calls are unchanged.
if nargin<3,    opts = [];      end                                         % MMP, 09/30/2026
Eq = collect_eq_rows(prog,P,opts,false);                                    % MMP, 09/30/2026
prog = impose_eq_rows(prog,Eq);                                             % MMP, 09/30/2026
% END MMP, 09/30/2026

end

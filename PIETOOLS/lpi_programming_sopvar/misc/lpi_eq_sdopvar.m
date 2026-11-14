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
% when P is the zero operator. Of each set of these equations that are
% multiples of one another (to 1e-12 of the largest coefficient), only
% one is written; equations without a decision variable are all written.
%
% INPUT
% - prog:   'struct', an LPI program (lpiprogram) whose
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
% See also LPI_EQ_SOP, LPI_EQ_CDOPVAR, LPI_EQ_NDOPVAR, POSSOPVAR, SOSEQ,
%          SDOPVAR.
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
% MMP, 08/29/2026: The multiplier contraction dropped: the canonical
%                  multiplier form is a class invariant (the 'sdopvar'
%                  constructor), so the stored coefficients are determined
%                  by the operator and are constrained directly, at the
%                  positions the form allows.
% MMP, 09/26/2026: Name handling off the hot path: P.Zd and prog.decvartable
%                  used as cellstr without conversion, membership by one
%                  unique over [table; names], soseq handed only the names a
%                  parameter uses, and one soseq per ~q nonzeros instead of
%                  one per parameter (the parameters' columns collected in
%                  the old order and imposed as joined row dpvars; the cap
%                  bounds the transient of soseq, which holds ~9 copies of
%                  its input). Same program bit for bit; measured in 1-D and
%                  2-D (2-D H-infinity, q = 7.6e5: 2.2 -> 1.1 s).
% MMP, 09/30/2026: The gamma of a parameter and the cell of its adjoint from
%                  the shared 'gamma_of_cell' and 'cell_of_gamma'
%                  (sopvar/misc/conventions). One routine per mode: the
%                  collection ('collect_eq_rows') and the imposition
%                  ('impose_eq_rows') moved verbatim to private helpers, and
%                  the 4th input, second output and struct-P batch form of
%                  09/26/2026, through which 'lpi_eq_cdopvar' reached them,
%                  are gone; it calls the two directly. Same program bit for
%                  bit. Earlier signatures:
% function [prog,Eq] = lpi_eq_sdopvar(prog,P,opts,dvars_checked)            % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% function prog = lpi_eq_sdopvar(prog,P,opts,dvars_checked)                 % MMP, 09/26/2026 (was)
% function prog = lpi_eq_sdopvar(prog,P,opts)                               % MMP, 09/26/2026 (was)
% MMP, 10/02/2026: Help only: the joined rows are written by 'lpi_soseq',
%                  not soseq (see 'impose_eq_rows'); each batch costs O(nnz)
%                  and ~80 B per nonzero.
% MMP, 10/16/2026: Moved, with private/, from sopvar/lpis_sopvar/ to
%                  lpi_programming_sopvar/; 10/21/2026 to its misc/ (helpers
%                  to misc/private/). Code unchanged.
% MMP, 10/28/2026: Header to the header protocol, step comments in the body;
%                  code unchanged. Deleted: the comparison with 'lpi_eq',
%                  the argument that constraining the canonical coefficients
%                  is exact, and notes on the joined prog.expr entries and
%                  the private routines.
% MMP, 11/01/2026: Help only: 'collect_eq_rows' drops the equations that
%                  are multiples of another, and the description says so.
% MMP, 11/04/2026: Help only: PROG comes from LPIPROGRAM; LPIPROGRAM_SOP
%                  is retired.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed
%                  (measurements to one figure), See also without the
%                  retired lpi_eq, the BEGIN note names misc/private; code
%                  unchanged.
% MMP, 11/10/2026: Help only: 'collect_eq_rows' writes every equation
%                  without a decision variable.


% BEGIN MMP, 09/30/2026: the three modes split. Moved verbatim to
% misc/private: the input checks and the collection loop (lines stamped MMP
% 08/28, 08/29, 09/26 and 09/30/2026) with the local function
% canonical_positions, to 'collect_eq_rows'; the local function impose_rows
% (MMP, 09/26/2026), to 'impose_eq_rows'. Deleted outright: the batch-form
% branch for a struct P and the nargout switch (MMP, 09/26/2026), through
% which 'lpi_eq_cdopvar' reached the two halves; it now calls them itself.
% Collect, then impose: the rows are unchanged.
if nargin<3,    opts = [];      end                                         % MMP, 09/30/2026
% The equations of P, every name of P.Zd checked against prog
Eq = collect_eq_rows(prog,P,opts,false);                                    % MMP, 09/30/2026
% Appended to prog.expr
prog = impose_eq_rows(prog,Eq);                                             % MMP, 09/30/2026
% END MMP, 09/30/2026

end

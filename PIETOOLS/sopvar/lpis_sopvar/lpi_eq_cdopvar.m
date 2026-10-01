function prog = lpi_eq_cdopvar(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_EQ_CDOPVAR(PROG,P) takes an LPI optimization program structure
% 'prog' and a 'cdopvar' decision operator P on a concatenation of mixed
% spaces, and adds equality constraints enforcing P==0. It is the container
% counterpart of 'lpi_eq_sdopvar', and applies in particular to the operators
% returned by 'poscopvar' and to sums and compositions of them.
%
% A container acts blockwise,
%
%   (P x)_i = sum_j P.C{i,j} x_j,
%
% on x = [x_1;...;x_N] with the x_j independent, so P==0 holds exactly when
% every block is the zero operator. Each block is an ordinary 'sdopvar' and
% is constrained by 'lpi_eq_sdopvar'; there is nothing further to do at the
% container level. A block stored as [] is structurally zero and generates no
% constraints.
%
% INPUT
% - prog:   'struct' specifying an LPI program structure to modify (see also
%           'lpiprogram'). All decision variables of P must already appear in
%           'prog.decvartable';
% - P:      'cdopvar' object of which to enforce P==0;
% - opts:   (optional) set opts = 'symmetric' if P is known to be
%           self-adjoint, to constrain only one block per adjoint pair.
%           Block (j,i) of a self-adjoint container is the adjoint of block
%           (i,j), so the off-diagonal blocks are constrained for i<j only,
%           and each diagonal block is passed 'symmetric' in turn so that its
%           own parameters are paired the same way. As in 'lpi_eq_sdopvar'
%           this is taken on the caller's word: if P is not self-adjoint the
%           resulting constraint is strictly weaker than P==0;
%
% OUTPUT
% - prog:   Same LPI program structure as the input, but with the added
%           constraints enforcing P==0.
%
% NOTES
% To enforce P==Q, call LPI_EQ_CDOPVAR(PROG,P-Q); '@cdopvar/plus' aligns    % MMP, 09/26/2026
% the blocks, monomial bases and decision variable lists of the two
% containers, and a zero block on one side passes through untouched.        % MMP, 09/26/2026
% (was) To enforce P==Q, call LPI_EQ_CDOPVAR(PROG,P+(-1)*Q); '@cdopvar/plus' aligns % MMP, 09/26/2026 (was)
% (was) containers, and a zero block on one side passes through untouched. Note % MMP, 09/26/2026 (was)
% (was) that neither 'cdopvar' nor 'copvar' defines 'minus' or 'uminus', so P-Q is % MMP, 09/26/2026 (was)
% (was) not executable; the scalar branch of '@cdopvar/mtimes' is what negates a % MMP, 09/26/2026 (was)
% (was) container.                                                          % MMP, 09/26/2026 (was)
%
% Where one block of an adjoint pair is structurally zero and the other is
% not, the nonzero one is constrained, so 'symmetric' does not depend on
% which of the two a sum happened to populate. For a genuinely self-adjoint
% container the two cases cannot differ in substance.
%
% All blocks sharing a decision variable list are imposed together, joined  % MMP, 09/26/2026
% into as few 'soseq' calls as a cap of numel(prog.decvartable) nonzeros    % MMP, 09/26/2026
% each allows (see 'lpi_eq_sdopvar'), so prog.expr gains a few entries per  % MMP, 09/26/2026
% call. The rows sossolve assembles (At and b, in order) are those one      % MMP, 09/26/2026
% 'soseq' per parameter gave.                                               % MMP, 09/26/2026
%
% See also LPI_EQ_SDOPVAR, LPI_EQ, POSCOPVAR, COPQUADVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_eq_cdopvar
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
% Initial coding MMP, 09/21/2026
% MMP, 09/21/2026: Fix the 'symmetric' skip condition, which tested the
%                   mirror block for emptiness rather than for carrying
%                   decision variables, so an adjoint pair made of a stored
%                   all-zero 'sopvar' and a nonzero 'sdopvar' had NEITHER
%                   block constrained -- the opposite of what the NOTES
%                   above promise. Reproduced by counting 'prog.expr.num'.
%                   Also corrected the P-Q recipe in the NOTES: no container
%                   class defines 'minus'.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  lpi_eq_mdopvar -> lpi_eq_cdopvar, posmopvar -> poscopvar.
%                  File was 'lpi_eq_mdopvar.m'.
% MMP, 09/26/2026: NOTES recipe is P-Q again: minus/uminus exist for both
%                  containers since 09/25/2026. The 09/21 note that no
%                  container class defines 'minus' is obsolete. Doc only.
% MMP, 09/26/2026: Verify each block's decision variables against the
%                  program once per distinct Zd list, not once per block.
%                  lpi_eq_sdopvar's check hashes all q names of
%                  prog.decvartable; blocks of a container sum share one list
%                  (@cdopvar/plus), and soseq adds no variables, so one check
%                  covers every block with that list. Measured, 2-D container
%                  Hinf build (q = 7.6e5, 6 blocks): 2.3 s of ismember became
%                  one call plus an O(q) isequal (~10 ms) per block. Same
%                  program, bit for bit.
% MMP, 09/26/2026: Blocks sharing a Zd list are imposed together, joined up
%                  to q nonzeros per soseq, instead of one soseq per nonzero
%                  parameter of every block: each soseq scans all q names of
%                  prog.decvartable in getequation. Blocks are collected with
%                  lpi_eq_sdopvar's second output and imposed by it, columns
%                  in the old order. 2-D container Hinf: 14 soseq -> 4, S6b
%                  1.9x faster at +50% stage peak, below stock's in both
%                  (numbers: lpi_eq_sdopvar, MMP 09/26/2026). prog.expr has
%                  fewer entries; the
%                  At, b, c, K sossolve builds are bit-identical, row order
%                  included. The 09/21 entry's check by counting
%                  prog.expr.num no longer discriminates (blocks now share
%                  entries); count equality rows (columns of At) instead.
% MMP, 09/30/2026: Collect with 'collect_eq_rows' and impose with
%                  'impose_eq_rows' (lpis_sopvar/private), the two halves
%                  of 'lpi_eq_sdopvar', called by name instead of through
%                  its second output, 4th input and struct-P modes, which
%                  are gone. The previous entry's "collected with
%                  lpi_eq_sdopvar's second output and imposed by it" now
%                  reads: collected by collect_eq_rows, imposed by
%                  impose_eq_rows. Same calls in the same order, same
%                  program bit for bit.

if isa(P,'copvar')
    error("Input of type 'copvar' carries no decision variables; use 'eq' "...
          +"on its blocks to test whether it is zero.")
end
if ~isa(P,'cdopvar')
    error("Input must be of type 'cdopvar'.")
end
if ~isa(prog,'struct') || ~isfield(prog,'decvartable')
    error("First input must be an LPI program structure; see 'lpiprogram'.")
end

symm = false;
if nargin>2 && ~isempty(opts)
    if (ischar(opts) || isstring(opts)) && strcmpi(char(opts),'symmetric')
        symm = true;
    else
        error("Third input is not recognized; the only supported option is "...
              +"'symmetric'.")
    end
end

[M,N] = size(P.C);
if symm && M~=N
    error("Only a square container can be self-adjoint; this one is "...
          +num2str(M)+" x "+num2str(N)+".")
end

% Last Zd list lpi_eq_sdopvar verified against prog.decvartable (NaN:       % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% Last Zd list collect_eq_rows verified against prog.decvartable (NaN:      % MMP, 09/30/2026
% none). prog.decvartable does not change in this loop, so the verdict      % MMP, 09/26/2026
% carries over to every later block with an equal list.                     % MMP, 09/26/2026
Zd_ok = NaN;                                                                % MMP, 09/26/2026
% Constraints collected, not yet imposed, all over the list Zd_ok. Every    % MMP, 09/26/2026
% soseq scans all q names of prog.decvartable, so all blocks sharing a list % MMP, 09/26/2026
% (the common case, see Zd_ok) are imposed together, ~q nonzeros per soseq. % MMP, 09/26/2026
batch = [];                                                                 % MMP, 09/26/2026
for i = 1:M
    for j = 1:N
        if symm && j<i
            % Skip only when the mirror will ACTUALLY be constrained, which % MMP, 09/21/2026
            % means it carries decision variables. Testing isempty alone     % MMP, 09/21/2026
            % also skipped this block when the mirror was a stored all-zero  % MMP, 09/21/2026
            % 'sopvar', and the mirror's own visit then returns early at the  % MMP, 09/21/2026
            % fixed-block branch below, so NEITHER was constrained -- the     % MMP, 09/21/2026
            % behaviour this file's own NOTES rule out. Measured on a 2 x 2   % MMP, 09/21/2026
            % container with block (1,2) replaced by a zero 'sopvar':         % MMP, 09/21/2026
            % 'symmetric' emitted 3 constraint expressions, exactly what the  % MMP, 09/21/2026
            % diagonal alone emits, against 4 for the all-decision case.      % MMP, 09/21/2026
            % Only reachable when the container is not in fact self-adjoint,  % MMP, 09/21/2026
            % since a zero block of a self-adjoint container has a zero       % MMP, 09/21/2026
            % mirror -- but a guard whose only job is to be safe under a      % MMP, 09/21/2026
            % violated precondition should not fail silently there.           % MMP, 09/21/2026
            if isa(P.C{j,i},'sdopvar')                                       % MMP, 09/21/2026
%           if ~isempty(P.C{j,i}) || isempty(P.C{i,j})                       % MMP, 09/21/2026 (was)
                continue
            end
        end
        Bij = P.C{i,j};
        if isempty(Bij)
            continue                        % structurally zero
        end
        if isa(Bij,'sopvar')
            % A fixed block among decision ones is legal in this container
            % (see the class documentation), but it carries no decision
            % variables, so it either is zero already or cannot be made so.
            if ~all(cellfun(@(c) isempty(c) || ~any(c(:)),Bij.params))
                error("Block ("+num2str(i)+","+num2str(j)+") is a fixed "...
                      +"'sopvar' and is not zero, so no choice of decision "...
                      +"variables can make the container vanish.")
            end
            continue
        end
        % isequal is O(q) with a small constant (~10 ms at q = 3.8e5); the  % MMP, 09/26/2026
        % check it replaces hashes the whole program table.                 % MMP, 09/26/2026
        checked = isequal(Bij.Zd,Zd_ok);                                    % MMP, 09/26/2026
        if ~checked && ~isempty(batch)                                      % MMP, 09/26/2026
            % New list: the pending rows index the old one, so impose them. % MMP, 09/26/2026
%           prog = lpi_eq_sdopvar(prog,batch);                              % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
            prog = impose_eq_rows(prog,batch);                              % MMP, 09/30/2026
            batch = [];                                                     % MMP, 09/26/2026
        end                                                                 % MMP, 09/26/2026
        % Second output: collect this block's constraints, do not impose.   % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
        % Collect this block's constraints, do not impose (prog unchanged). % MMP, 09/30/2026
        if symm && i==j
%           prog = lpi_eq_sdopvar(prog,Bij,'symmetric');                    % MMP, 09/26/2026 (was)
%           prog = lpi_eq_sdopvar(prog,Bij,'symmetric',checked);            % MMP, 09/26/2026 (was)
%           [prog,Eq] = lpi_eq_sdopvar(prog,Bij,'symmetric',checked);       % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
            Eq = collect_eq_rows(prog,Bij,'symmetric',checked);             % MMP, 09/30/2026
        else
%           prog = lpi_eq_sdopvar(prog,Bij);                                % MMP, 09/26/2026 (was)
%           prog = lpi_eq_sdopvar(prog,Bij,[],checked);                     % MMP, 09/26/2026 (was)
%           [prog,Eq] = lpi_eq_sdopvar(prog,Bij,[],checked);                % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
            Eq = collect_eq_rows(prog,Bij,[],checked);                      % MMP, 09/30/2026
        end
        % Appended after the earlier blocks' columns, the order one soseq   % MMP, 09/26/2026
        % per block imposed them in, so sossolve's rows keep their order.   % MMP, 09/26/2026
        if isempty(batch),  batch = Eq;                                     % MMP, 09/26/2026
        else,               batch.Cs = [batch.Cs, Eq.Cs];                   % MMP, 09/26/2026
        end                                                                 % MMP, 09/26/2026
        % Verified now: lpi_eq_sdopvar errors on an unknown variable.       % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
        % Verified now: collect_eq_rows errors on an unknown variable.      % MMP, 09/30/2026
        Zd_ok = Bij.Zd;                                                     % MMP, 09/26/2026
    end
end
% The rest: all blocks if they share one list. lpi_eq_sdopvar caps each     % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% soseq at q nonzeros, bounding its transient; 'batch' holds one copy of    % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% the collected coefficients, O(nnz), until then (see impose_rows there).   % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% The rest: all blocks if they share one list. impose_eq_rows caps each     % MMP, 09/30/2026
% soseq at q nonzeros, bounding its transient; 'batch' holds one copy of    % MMP, 09/30/2026
% the collected coefficients, O(nnz), until then (see its help).            % MMP, 09/30/2026
if ~isempty(batch)                                                          % MMP, 09/26/2026
%   prog = lpi_eq_sdopvar(prog,batch);                                      % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
    prog = impose_eq_rows(prog,batch);                                      % MMP, 09/30/2026
end                                                                         % MMP, 09/26/2026

end

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
% To enforce P==Q, call LPI_EQ_CDOPVAR(PROG,P+(-1)*Q); '@cdopvar/plus' aligns
% the blocks, monomial bases and decision variable lists of the two
% containers, and a zero block on one side passes through untouched. Note
% that neither 'cdopvar' nor 'copvar' defines 'minus' or 'uminus', so P-Q is
% not executable; the scalar branch of '@cdopvar/mtimes' is what negates a
% container.
%
% Where one block of an adjoint pair is structurally zero and the other is
% not, the nonzero one is constrained, so 'symmetric' does not depend on
% which of the two a sum happened to populate. For a genuinely self-adjoint
% container the two cases cannot differ in substance.
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
        if symm && i==j
            prog = lpi_eq_sdopvar(prog,Bij,'symmetric');
        else
            prog = lpi_eq_sdopvar(prog,Bij);
        end
    end
end

end

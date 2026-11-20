function prog = lpi_eq_cdopvar(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_EQ_CDOPVAR(PROG,P,OPTS) adds to the LPI program PROG the
% equality constraints P == 0 for the 'cdopvar' decision container P,
%
%   (P x)_i = sum_j P.C{i,j} x_j,
%
% block by block. An 'sdopvar' block, with parameter cells gamma of
% coefficients vec(C_gamma(d)) = A_gamma + B_gamma'*d over the decision
% variables d = P.Zd, gets the equations A_gamma(c) + B_gamma(:,c)'*d = 0
% at each position c of C_gamma that its canonical multiplier form allows
% (the columns whose ZR monomial has degree 0 in every multiplier
% direction of gamma), all-zero equations omitted; these hold exactly
% when the block is the zero operator. Of each set of equations of a block
% that are multiples of one another, only one is written (see entries.dup
% in OUTPUT); equations without a decision variable are all written. A
% block [] or an all-zero 'sopvar' block gets no equation.
%
% INPUT
% - prog:   'struct', an LPI program whose prog.decvartable holds every
%           name of P.Zd; with the field sopeq (a cell, e.g. {}), the rows
%           of this call are recorded there (see OUTPUT);
% - P:      'cdopvar' object, M x N blocks P.C{i,j}, each an 'sdopvar'
%           object, a 'sopvar' object or [];
% - opts:   (optional) 'char', 'symmetric' for a P that is self-adjoint for
%           every value of d: of the blocks (i,j) and (j,i), i < j, only
%           (i,j) is constrained, or (j,i) when (i,j) is not an 'sdopvar';
%           in a diagonal block, of each pair of parameter cells mapped to
%           each other by exchanging lower (2) and upper (3) integrals in
%           every direction, only the one of lower linear index (default
%           []: every block and cell).
%
% OUTPUT
% - prog:   'struct', PROG with the equations appended to prog.expr as
%           'eq' entries, in the block order (1,1), (1,2), ..., (M,N) and,
%           within a block, in the order of the parameter cells and
%           positions; when PROG has the field sopeq, a 'struct' appended to
%           prog.sopeq with fields
%             rows       1 x 2 double, the first and last row of the
%                        equations among the rows of prog.expr;
%             entries    struct array, one entry per 'sdopvar' block except
%                        those 'symmetric' skips, in row order, with fields
%                        i, j (double, the block), cells (double row, its
%                        parameter cells with equations), cols (cell, per
%                        cell the linear positions c of its equations, in
%                        row order), blk ('struct' with the ZL, ZR, dims,
%                        vars and dom of the block, NL and NR (double, its
%                        numbers of left and right monomials), psize (the
%                        size of its parameter array) and symmetric) and
%                        dup ('struct' with fields cell, col, rcell, rcol
%                        and factor, double columns: per equation of the
%                        block not written because it is a multiple, its
%                        parameter cell and position, those of the written
%                        equation, and f with (equation not written) =
%                        f*(equation written), to 1e-12 of its largest
%                        coefficient);
%             M, N       double, the block grid;
%             symmetric  logical, OPTS was 'symmetric';
%             space_out, space_in, vars, dom
%                        those of P.
%
% NOTES
% P a 'copvar' or of a class other than 'cdopvar', PROG without the field
% decvartable, OPTS other than 'symmetric', 'symmetric' with M ~= N, a
% nonzero 'sopvar' block, a block not in canonical multiplier form, and a
% name of P.Zd absent from prog.decvartable are errors. To impose P == Q,
% call LPI_EQ_CDOPVAR(PROG,P-Q).
%
% See also LPI_EQ_SDOPVAR, LPI_EQ_SOP, POSCOPVAR, COPQUADVAR, CDOPVAR.
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

% Inputs: a decision container and an LPI program
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

% Option 'symmetric'
symm = false;
if nargin>2 && ~isempty(opts)
    if (ischar(opts) || isstring(opts)) && strcmpi(char(opts),'symmetric')
        symm = true;
    else
        error("Third input is not recognized; the only supported option is "...
              +"'symmetric'.")
    end
end

% Block grid; 'symmetric' needs a square one
[M,N] = size(P.C);
if symm && M~=N
    error("Only a square container can be self-adjoint; this one is "...
          +num2str(M)+" x "+num2str(N)+".")
end

% Last Zd list verified against prog.decvartable (NaN: none).
% prog.decvartable does not change in this loop, so the verdict carries
% over to every later block with an equal list.
Zd_ok = NaN;
% Constraints collected, not yet imposed, all over the list Zd_ok. The
% check of the list against prog.decvartable runs once per list, so all
% blocks sharing one (the common case, see Zd_ok) are imposed together.
batch = [];
% Row bookkeeping for the dual: the row count of the program before this
% call and, per block written, the parameters and positions in the order
% the rows are produced. Opt-in: only a program that already carries the
% field 'sopeq' is recorded; any other program keeps the stock fields.
record = isfield(prog,'sopeq');
r0 = 0;
if record && isfield(prog,'expr') && isfield(prog.expr,'At')
    for e = 1:numel(prog.expr.At),  r0 = r0 + size(prog.expr.At{e},2);  end
end
% dup: the equations dropped as multiples
ents = struct('i',{},'j',{},'cells',{},'cols',{},'blk',{},'dup',{});
% Blocks in the order (1,1), (1,2), ..., (M,N)
for i = 1:M
    for j = 1:N
        % Under 'symmetric', a block below the diagonal is skipped when its
        % mirror (j,i) is an 'sdopvar', which is constrained instead
        if symm && j<i
            % Skip only when the mirror carries decision variables: a
            % stored all-zero 'sopvar' mirror returns early below, so an
            % isempty test left NEITHER block constrained. Reachable only
            % when the container is not self-adjoint; no silent failure.
            if isa(P.C{j,i},'sdopvar')
                continue
            end
        end
        % Block (i,j): [] and an all-zero 'sopvar' add nothing, a nonzero
        % 'sopvar' is an error
        Bij = P.C{i,j};
        if isempty(Bij)
            continue                        % structurally zero
        end
        if isa(Bij,'sopvar')
            % A fixed block among decision ones is legal in this container,
            % but it carries no decision variables, so it either is zero
            % already or cannot be made so.
            if ~all(cellfun(@(c) isempty(c) || ~any(c(:)),Bij.params))
                error("Block ("+num2str(i)+","+num2str(j)+") is a fixed "...
                      +"'sopvar' and is not zero, so no choice of decision "...
                      +"variables can make the container vanish.")
            end
            continue
        end
        % A list equal to the last verified one needs no check; a new list
        % first imposes the rows collected over the old one
        % isequal is O(q) with a small constant (~10 ms at q = 3.8e5); a
        % full check hashes the whole program table.
        checked = isequal(Bij.Zd,Zd_ok);
        if ~checked && ~isempty(batch)
            % New list: the pending rows index the old one, so impose them.
            prog = impose_eq_rows(prog,batch);
            batch = [];
        end
        % Collect the constraints of this block, do not impose (prog unchanged).
        if symm && i==j
            Eq = collect_eq_rows(prog,Bij,'symmetric',checked);
        else
            Eq = collect_eq_rows(prog,Bij,[],checked);
        end
        % Appended after the columns of the earlier blocks, so the rows of
        % the SDP keep the block order. batch keeps the Eq.pos (table rows)
        % of its first block: that one is collected unchecked (a batch
        % starts only at Zd_ok = NaN or after a flush, i.e. ~checked), every
        % later one has an equal list, and nothing here changes
        % prog.decvartable. Otherwise the rows, written at the trusted pos,
        % land on the wrong variables.
        if isempty(batch),  batch = Eq;
        else,               batch.Cs = [batch.Cs, Eq.Cs];
        end
        % Row record of the block
        if record
            ents(end+1) = struct('i',i,'j',j,'cells',Eq.cells,'cols',{Eq.cols},'blk',Eq.blk,'dup',Eq.dup); %#ok<AGROW>
        end
        % Verified now: the collection errors on an unknown variable.
        Zd_ok = Bij.Zd;
    end
end
% The rest: all blocks if they share one list. Each prog.expr entry is
% capped at q nonzeros, bounding its transient; 'batch' holds one copy of
% the collected coefficients, O(nnz), until then.
if ~isempty(batch)
    prog = impose_eq_rows(prog,batch);
end
% Record the rows of this call (prog.sopeq{end})
if record
    % Rows written: one per recorded position
    nrow = 0;
    for e = 1:numel(ents),  nrow = nrow + sum(cellfun(@numel,ents(e).cols));    end
    prog.sopeq{end+1} = struct('rows',[r0+1,r0+nrow],'entries',ents,'M',M,'N',N, ...
                               'symmetric',symm,'space_out',P.space_out, ...
                               'space_in',P.space_in,'vars',{P.vars},'dom',P.dom);
end

end

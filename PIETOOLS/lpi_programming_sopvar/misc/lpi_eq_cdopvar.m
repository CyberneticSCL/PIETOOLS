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
% when the block is the zero operator. A block [] or an all-zero 'sopvar'
% block gets no equation.
%
% INPUT
% - prog:   'struct', an LPI program (lpiprogram or lpiprogram_sop) whose
%           prog.decvartable holds every name of P.Zd; with the field sopeq
%           (a cell, e.g. {}), the rows of this call are recorded there
%           (see OUTPUT);
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
%                        row order) and blk ('struct' with the block's ZL,
%                        ZR, dims, vars and dom, NL and NR (double, its
%                        numbers of left and right monomials), psize (the
%                        size of its parameter array) and symmetric);
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
% MMP, 10/08/2026: Row bookkeeping for the dual: prog.sopeq{end} records the
%                  rows this call produced (first and last row of the
%                  program, the blocks in writing order with their
%                  parameters, positions and bases, from collect_eq_rows),
%                  so 'lpigetdual_sop' can assemble prog.solinfo.y into the
%                  signed dual kernel of this equality. Rows unchanged.
%                  Opt-in (same day): recorded only when the program already
%                  carries the field 'sopeq', so that no other caller's
%                  program gains a field. Checked: b'y equals the primal
%                  objective with sos_opts.simplify 0 and 1 (sossolve puts
%                  the multipliers back in the original row coordinates).
% MMP, 10/02/2026: Comments only. impose_eq_rows now writes each joined
%                  batch with 'lpi_soseq' (the rows from collect_eq_rows'
%                  check), so no soseq scans the program names; 'soseq
%                  calls' below read 'prog.expr entries'. The batch keeps
%                  the Eq.pos of its first block, valid for every later
%                  block with an equal list.
% MMP, 10/16/2026: moved, with private/, from sopvar/lpis_sopvar/ to
%                  lpi_programming_sopvar/; 'lpis_sopvar' in the entries
%                  above now reads lpi_programming_sopvar. Code unchanged.
% MMP, 10/21/2026: moved to lpi_programming_sopvar/misc/ and the private
%                  helpers to misc/private/: the folder root keeps the user
%                  functions only. The locations in the entries above read
%                  accordingly. Code unchanged.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), step comments in the body; code
%                  unchanged. The deleted header text was the comparison
%                  with 'lpi_eq_sdopvar', the argument that blockwise
%                  equations are exact, the P-Q recipe's (was) lines, the
%                  reason for the rule on adjoint pairs with one zero block
%                  and notes on joining blocks into prog.expr entries. The
%                  rule that the 09/21/2026 entry and the comment at the
%                  'symmetric' skip call the NOTES is now in the OPTS entry.

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

% Last Zd list lpi_eq_sdopvar verified against prog.decvartable (NaN:       % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% Last Zd list collect_eq_rows verified against prog.decvartable (NaN:      % MMP, 09/30/2026
% none). prog.decvartable does not change in this loop, so the verdict      % MMP, 09/26/2026
% carries over to every later block with an equal list.                     % MMP, 09/26/2026
Zd_ok = NaN;                                                                % MMP, 09/26/2026
% Constraints collected, not yet imposed, all over the list Zd_ok. Every    % MMP, 09/26/2026
% soseq scans all q names of prog.decvartable, so all blocks sharing a list % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
% (the common case, see Zd_ok) are imposed together, ~q nonzeros per soseq. % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
% check of the list against prog.decvartable runs once per list, so all     % MMP, 10/02/2026
% blocks sharing one (the common case, see Zd_ok) are imposed together.     % MMP, 10/02/2026
batch = [];                                                                 % MMP, 09/26/2026
% Row bookkeeping for the dual read-back ('lpigetdual_sop'): the program's % MMP, 10/08/2026
% row count before this call and, per block written, the parameters and    % MMP, 10/08/2026
% positions in the order the rows are produced. Opt-in: only a program     % MMP, 10/08/2026
% that already carries the field 'sopeq' (the container executives set     % MMP, 10/08/2026
% prog.sopeq = {} after lpiprogram) is recorded; every other caller's       % MMP, 10/08/2026
% program keeps the stock fields.                                           % MMP, 10/08/2026
record = isfield(prog,'sopeq');                                             % MMP, 10/08/2026
r0 = 0;                                                                     % MMP, 10/08/2026
if record && isfield(prog,'expr') && isfield(prog.expr,'At')                % MMP, 10/08/2026
    for e = 1:numel(prog.expr.At),  r0 = r0 + size(prog.expr.At{e},2);  end % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026
ents = struct('i',{},'j',{},'cells',{},'cols',{},'blk',{});                 % MMP, 10/08/2026
% Blocks in the order (1,1), (1,2), ..., (M,N)
for i = 1:M
    for j = 1:N
        % Under 'symmetric', a block below the diagonal is skipped when its
        % mirror (j,i) is an 'sdopvar', which is constrained instead
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
        % Block (i,j): [] and an all-zero 'sopvar' add nothing, a nonzero
        % 'sopvar' is an error
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
        % A list equal to the last verified one needs no check; a new list
        % first imposes the rows collected over the old one
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
        % batch keeps the Eq.pos (table rows) of its first block: that one  % MMP, 10/02/2026
        % is collected unchecked (a batch starts only at Zd_ok = NaN or     % MMP, 10/02/2026
        % after a flush, i.e. ~checked), every later one has an equal list, % MMP, 10/02/2026
        % and nothing here changes prog.decvartable. Break that, and        % MMP, 10/02/2026
        % lpi_soseq, which trusts pos, puts rows on the wrong variables.    % MMP, 10/02/2026
        if isempty(batch),  batch = Eq;                                     % MMP, 09/26/2026
        else,               batch.Cs = [batch.Cs, Eq.Cs];                   % MMP, 09/26/2026
        end                                                                 % MMP, 09/26/2026
        % Row record of the block
        if record                                                           % MMP, 10/08/2026
            ents(end+1) = struct('i',i,'j',j,'cells',Eq.cells,'cols',{Eq.cols},'blk',Eq.blk); %#ok<AGROW> % MMP, 10/08/2026
        end                                                                 % MMP, 10/08/2026
        % Verified now: lpi_eq_sdopvar errors on an unknown variable.       % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
        % Verified now: collect_eq_rows errors on an unknown variable.      % MMP, 09/30/2026
        Zd_ok = Bij.Zd;                                                     % MMP, 09/26/2026
    end
end
% The rest: all blocks if they share one list. lpi_eq_sdopvar caps each     % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% soseq at q nonzeros, bounding its transient; 'batch' holds one copy of    % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% the collected coefficients, O(nnz), until then (see impose_rows there).   % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% The rest: all blocks if they share one list. impose_eq_rows caps each     % MMP, 09/30/2026
% soseq at q nonzeros, bounding its transient; 'batch' holds one copy of    % MMP, 09/30/2026 % MMP, 10/02/2026 (was)
% prog.expr entry at q nonzeros, bounding its transient; 'batch' holds      % MMP, 10/02/2026
% one copy of                                                               % MMP, 10/02/2026
% the collected coefficients, O(nnz), until then (see its help).            % MMP, 09/30/2026
if ~isempty(batch)                                                          % MMP, 09/26/2026
%   prog = lpi_eq_sdopvar(prog,batch);                                      % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
    prog = impose_eq_rows(prog,batch);                                      % MMP, 09/30/2026
end                                                                         % MMP, 09/26/2026
% Record the rows of this call (prog.sopeq{end}) for 'lpigetdual_sop'.      % MMP, 10/08/2026
if record                                                                   % MMP, 10/08/2026
    % Rows written: one per recorded position
    nrow = 0;                                                               % MMP, 10/08/2026
    for e = 1:numel(ents),  nrow = nrow + sum(cellfun(@numel,ents(e).cols));    end % MMP, 10/08/2026
    prog.sopeq{end+1} = struct('rows',[r0+1,r0+nrow],'entries',ents,'M',M,'N',N, ...% MMP, 10/08/2026
                               'symmetric',symm,'space_out',P.space_out, ...% MMP, 10/08/2026
                               'space_in',P.space_in,'vars',{P.vars},'dom',P.dom); % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026

end

function Pmop = opvar2copvar(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMOP = OPVAR2COPVAR(POP) takes an 'opvar' object POP and returns the
% equivalent 2 x 2 'copvar' container.
%
% INPUTS
% - Pop:    'opvar' object, P: R^m2 x L2^n2[s] -> R^m1 x L2^n1[s], with     % MMP, 09/25/2026
%           dim = [m1 m2; n1 n2]; a zero-dimensional space is dropped and   % MMP, 09/25/2026
%           a zero component becomes a [] block. Unlike                     % MMP, 09/25/2026
% - (was)   'opvar' object, P: R^m x L2^n[s] -> R^m x L2^n[s]. Unlike       % MMP, 09/25/2026 (was)
%           'opvar2sopvar', ALL FOUR components may be nonempty - handling
%           the full 4-PI operator is the point of this routine;
%
% OUTPUTS
% - Pmop:   'copvar' from (R^m2, L2^n2[s]) to (R^m1, L2^n1[s]), with        % MMP, 09/25/2026
% - (was)   'copvar' over the spaces (R^m, L2^n[s]), with                   % MMP, 09/25/2026 (was)
%
%               Pmop.C{1,1} = P    (R^m2  -> R^m1)                          % MMP, 09/25/2026
%               Pmop.C{1,2} = Q1   (L2^n2 -> R^m1)                          % MMP, 09/25/2026
%               Pmop.C{2,1} = Q2   (R^m2  -> L2^n1)                         % MMP, 09/25/2026
%               Pmop.C{2,2} = R    (L2^n2 -> L2^n1)                         % MMP, 09/25/2026
%     (was)     Pmop.C{1,1} = P    (R^m  -> R^m)                            % MMP, 09/25/2026 (was)
%     (was)     Pmop.C{1,2} = Q1   (L2^n -> R^m)                            % MMP, 09/25/2026 (was)
%     (was)     Pmop.C{2,1} = Q2   (R^m  -> L2^n)                           % MMP, 09/25/2026 (was)
%     (was)     Pmop.C{2,2} = R    (L2^n -> L2^n)                           % MMP, 09/25/2026 (was)
%
%           Row is output, column is input, which is the container's
%           convention and matches opvar's own component naming;
%
% NOTES
% The per-block conversion is delegated to 'opvar2sopvar' rather than
% reimplemented. That routine accepts only a single-component 'opvar', so
% each of the four components is first placed in an otherwise empty 'opvar'
% whose dim matrix selects just that block. opvar's dim is
% [out_R in_R; out_L2 in_L2], so the four selectors are
%
%   P  [m1 m2;0 0]   Q1 [m1 0;0 n2]   Q2 [0 m2;n1 0]   R  [0 0;n1 n2]       % MMP, 09/25/2026
% (was) P  [m m;0 0]     Q1 [m 0;0 n]     Q2 [0 m;n 0]     R  [0 0;n n]     % MMP, 09/25/2026 (was)
%
% A component that is empty or identically zero is left as a structurally
% zero block ([]) rather than converted, so a sparse operator stays sparse
% and the container does not carry bases for blocks that hold nothing.
% Exception: a row or column that would have no block at all gets one       % MMP, 09/25/2026
% explicit zero block, since 'verify' and copvar(C) read its space off it.  % MMP, 09/25/2026
%
% A zero-dimensional space is skipped: for m1 = 0 there is no R row, for    % MMP, 09/25/2026
% m2 = 0 no R column, likewise n1, n2 for L2, so e.g. dim = [0 0; n n]      % MMP, 09/25/2026
% gives a 1 x 1 container. Including such a space would give a row or       % MMP, 09/25/2026
% column of zero dimension.                                                 % MMP, 09/25/2026
% (was) A zero-dimensional side is skipped the same way: for m = 0 the whole first % MMP, 09/25/2026 (was)
% (was) row and column are absent, and the result is a 1 x 1 container. Rows and % MMP, 09/25/2026 (was)
% (was) columns of a 'copvar' cannot be entirely empty, since the container reads % MMP, 09/25/2026 (was)
% (was) each row's space and dimension off a populated block, so an operator with % MMP, 09/25/2026 (was)
% (was) m = 0 or n = 0 must not produce a 2 x 2 grid.                       % MMP, 09/25/2026 (was)
%
% See also COPVAR2OPVAR, OPVAR2SOPVAR, OPVAR2D2COPVAR, RAND_COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar2copvar
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
% Initial coding MMP, 09/18/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  opvar2mopvar -> opvar2copvar. File was 'opvar2mopvar.m'.
% MMP, 09/25/2026: Accept RECTANGULAR opvars, dim = [m1 m2; n1 n2], and
%                  operators with an all-zero row or column. The PIE's input
%                  and output operators are of that kind - B1: R^nw ->
%                  R^n0 x L2^n1, C1, D11 - and the H-infinity KYP operator
%                  cannot be assembled in a container without them. Two
%                  changes: the four component selectors take the general
%                  dim, and the metadata is stated from the opvar's own dim,
%                  var1 and I instead of derived from the blocks, since the
%                  derivation needs a populated block in every row and column
%                  and a zero component (D11 = 0, or a B1 whose ODE part is
%                  zero) leaves one empty. For dim = [m m; n n] the selectors
%                  are the four used before and the result is unchanged. The
%                  'opvar2copvar:nonSquare' error of 09/18/2026 is gone.
% MMP, 09/25/2026: A row or column left with no block by zero components
%                  gets one explicit zero block. With only the stated
%                  metadata such a container failed 'verify' ("Row 1 contains
%                  no populated block") and could not be rebuilt by copvar(C)
%                  ('copvar:emptyRow'), e.g. D11 = 0. A zero block costs one
%                  degree-0 basis.

if ~isa(Pop,'opvar')
    error('opvar2copvar:badInput','Input must be an opvar object.')
end
d = Pop.dim;
% % % BEGIN change MMP, 09/25/2026: general dim. Deleted: the square check  % MMP, 09/25/2026
% % % and its 'nonSquare' error, and the square selectors (09/18/2026).     % MMP, 09/25/2026
% if d(1,1)~=d(1,2) || d(2,1)~=d(2,2)                                       % MMP, 09/25/2026 (was)
%     error('opvar2copvar:nonSquare',...
%         ['This routine expects an opvar with equal input and output '...
%          'dimensions on each space, i.e. dim = [m m; n n]; got %s.'],mat2str(d)) % MMP, 09/25/2026 (was)
% end                                                                       % MMP, 09/25/2026 (was)
% m = d(1,1);     n = d(2,1);                                               % MMP, 09/25/2026 (was)

% Which of the four spaces are present at all. A zero-dimensional space is
% not a space: including it would give the container an all-empty row.
% keep_R  = m>0;                                                            % MMP, 09/25/2026 (was)
% keep_L2 = n>0;                                                            % MMP, 09/25/2026 (was)
% if ~keep_R && ~keep_L2                                                    % MMP, 09/25/2026 (was)
%     error('opvar2copvar:empty','Operator has no dimensions to convert.')  % MMP, 09/25/2026 (was)
% end                                                                       % MMP, 09/25/2026 (was)
rows = [d(1,1)>0, d(2,1)>0];         % output spaces R, L2 that exist       % MMP, 09/25/2026
cols = [d(1,2)>0, d(2,2)>0];         % input spaces R, L2 that exist        % MMP, 09/25/2026
if ~any(rows) || ~any(cols)                                                 % MMP, 09/25/2026
    error('opvar2copvar:empty','Operator has no dimensions to convert.')    % MMP, 09/25/2026
end                                                                         % MMP, 09/25/2026

% One selector per component, in the [out_R in_R; out_L2 in_L2] layout.
% sel = {[m m;0 0], [m 0;0 n];                                              % MMP, 09/25/2026 (was)
%        [0 m;n 0], [0 0;n n]};                                             % MMP, 09/25/2026 (was)
% For d = [m m;n n] these are the four selectors above.                     % MMP, 09/25/2026
sel = {[d(1,1) d(1,2);0 0], [d(1,1) 0;0 d(2,2)];                            % MMP, 09/25/2026
       [0 d(1,2);d(2,1) 0], [0 0;d(2,1) d(2,2)]};                           % MMP, 09/25/2026
nm  = {'P' ,'Q1';
       'Q2','R'};
% rows = [keep_R, keep_L2];       cols = [keep_R, keep_L2];                 % MMP, 09/25/2026 (was)
% % % END change MMP, 09/25/2026                                            % MMP, 09/25/2026

% % % BEGIN change MMP, 09/25/2026: no row or column entirely [].           % MMP, 09/25/2026
% Blocks to convert: the nonzero components, plus one explicit zero block   % MMP, 09/25/2026
% in each row, then each column, they leave empty - 'verify' and            % MMP, 09/25/2026
% copvar(C) read a row's or column's space off a block. The opvar holds     % MMP, 09/25/2026
% such a component as zeros of the right size, which converts as is.        % MMP, 09/25/2026
keep = false(2,2);                                                          % MMP, 09/25/2026
for i = find(rows),     for j = find(cols)                                  % MMP, 09/25/2026
    keep(i,j) = ~is_zero_component(Pop.(nm{i,j}));                          % MMP, 09/25/2026
end,                    end                                                 % MMP, 09/25/2026
ir = find(rows);        jc = find(cols);                                    % MMP, 09/25/2026
for i = ir,     if ~any(keep(i,:)),     keep(i,jc(1)) = true;   end,    end % MMP, 09/25/2026
for j = jc,     if ~any(keep(:,j)),     keep(ir(1),j) = true;   end,    end % MMP, 09/25/2026
% % % END change MMP, 09/25/2026                                            % MMP, 09/25/2026

C = cell(sum(rows),sum(cols));
ri = 0;
for i = 1:2
    if ~rows(i),    continue,   end
    ri = ri+1;      ci = 0;
    for j = 1:2
        if ~cols(j),    continue,   end
        ci = ci+1;
        comp = Pop.(nm{i,j});
%       if is_zero_component(comp)                                          % MMP, 09/25/2026 (was)
        if ~keep(i,j)                                                       % MMP, 09/25/2026
            continue                            % structurally zero block
        end
        Pblk = opvar();
        Pblk.I = Pop.I;     Pblk.var1 = Pop.var1;   Pblk.var2 = Pop.var2;
        Pblk.dim = sel{i,j};
        Pblk.(nm{i,j}) = comp;
        C{ri,ci} = opvar2sopvar(Pblk);
    end
end
% Pmop = copvar(C);                                                         % MMP, 09/25/2026 (was)
% % % BEGIN change MMP, 09/25/2026: metadata from the opvar, not blocks.    % MMP, 09/25/2026
% A zero component leaves a row or column with no block, which the block    % MMP, 09/25/2026
% derivation rejects; the opvar states the spaces itself. The registry      % MMP, 09/25/2026
% holds var1 only if an L2 space survives.                                  % MMP, 09/25/2026
isL2 = [false, true];                                                       % MMP, 09/25/2026
if any(rows & isL2) || any(cols & isL2)                                     % MMP, 09/25/2026
    vars = {char(Pop.var1.varname{1})};     dom = double(Pop.I(1,:));       % MMP, 09/25/2026
    so = isL2(rows)';                       si = isL2(cols)';               % MMP, 09/25/2026
else                                                                        % MMP, 09/25/2026
    vars = cell(1,0);                       dom = zeros(0,2);               % MMP, 09/25/2026
    so = false(nnz(rows),0);                si = false(nnz(cols),0);        % MMP, 09/25/2026
end                                                                         % MMP, 09/25/2026
meta = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',si,...
    'dim_out',d(rows,1),'dim_in',d(cols,2));                                % MMP, 09/25/2026
Pmop = copvar(C,meta);                                                      % MMP, 09/25/2026
% % % END change MMP, 09/25/2026                                            % MMP, 09/25/2026

end

% ------------------------------------------------------------------------
function tf = is_zero_component(comp)
% True when a component carries nothing. R is a struct of R0, R1, R2, so a
% scalar emptiness test does not reach it.
if isstruct(comp)
    tf = true;
    f = fieldnames(comp);
    for k = 1:numel(f)
        tf = tf && is_zero_component(comp.(f{k}));
    end
    return
end
if isempty(comp)
    tf = true;      return
end
comp = polynomial(comp);
tf = isempty(comp.coefficient) || ~any(comp.coefficient(:));
end

function Pmop = opvar2copvar(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMOP = OPVAR2COPVAR(POP) takes an 'opvar' object POP and returns the
% equivalent 2 x 2 'copvar' container.
%
% INPUTS
% - Pop:    'opvar' object, P: R^m2 x L2^n2[s] -> R^m1 x L2^n1[s], with
%           dim = [m1 m2; n1 n2]; a zero-dimensional space is dropped and
%           a zero component becomes a [] block. All four components may
%           be nonempty;
%
% OUTPUTS
% - Pmop:   'copvar' from (R^m2, L2^n2[s]) to (R^m1, L2^n1[s]), with
%
%               Pmop.C{1,1} = P    (R^m2  -> R^m1)
%               Pmop.C{1,2} = Q1   (L2^n2 -> R^m1)
%               Pmop.C{2,1} = Q2   (R^m2  -> L2^n1)
%               Pmop.C{2,2} = R    (L2^n2 -> L2^n1)
%
%           Row is output, column is input, as in the component naming of
%           the 'opvar';
%
% NOTES
% Each of the four components is placed in an otherwise empty 'opvar'
% whose dim matrix selects just that block, and converted to an 'sopvar'
% block. The dim of an 'opvar' is [out_R in_R; out_L2 in_L2], so the four
% selectors are
%
%   P  [m1 m2;0 0]   Q1 [m1 0;0 n2]   Q2 [0 m2;n1 0]   R  [0 0;n1 n2]
%
% A component that is empty or identically zero is left as a structurally
% zero block ([]), so a sparse operator stays sparse and the container
% does not carry bases for blocks that hold nothing.
% Exception: a row or column that would have no block at all gets one
% explicit zero block, since the container reads the space of a row or
% column off one of its blocks.
%
% A zero-dimensional space is skipped: for m1 = 0 there is no R row, for
% m2 = 0 no R column, likewise n1, n2 for L2, so e.g. dim = [0 0; n n]
% gives a 1 x 1 container. Including such a space would give a row or
% column of zero dimension.
% With no space at all on one side (m1 = n1 = 0, or m2 = n2 = 0) the
% container has no space on that side: an empty channel of a PIE.
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

% Check the input
if ~isa(Pop,'opvar')
    error('opvar2copvar:badInput','Input must be an opvar object.')
end
d = Pop.dim;

% Which of the four spaces are present at all. A zero-dimensional space is
% not a space: including it would give the container an all-empty row.
rows = [d(1,1)>0, d(2,1)>0];         % output spaces R, L2 that exist
cols = [d(1,2)>0, d(2,2)>0];         % input spaces R, L2 that exist
% No space on a side (an empty channel of a PIE, as C1 without a regulated
% output): the container has no space on that side, as an opvar has a
% zero dimension; its products and concatenations contribute nothing.

% One selector per component, in the [out_R in_R; out_L2 in_L2] layout.
sel = {[d(1,1) d(1,2);0 0], [d(1,1) 0;0 d(2,2)];
       [0 d(1,2);d(2,1) 0], [0 0;d(2,1) d(2,2)]};
nm  = {'P' ,'Q1';
       'Q2','R'};

% Blocks to convert: the nonzero components, plus one explicit zero block
% in each row, then each column, they leave empty - the container reads
% the space of a row or column off a block. The opvar holds such a
% component as zeros of the right size, which converts as is.
keep = false(2,2);
for i = find(rows),     for j = find(cols)
    keep(i,j) = ~is_zero_component(Pop.(nm{i,j}));
end,                    end
ir = find(rows);        jc = find(cols);
% With no space on the other side there is no block to place: the row or
% column stays empty and the stated metadata carries its space
for i = ir,     if ~isempty(jc) && ~any(keep(i,:)),     keep(i,jc(1)) = true;   end,    end
for j = jc,     if ~isempty(ir) && ~any(keep(:,j)),     keep(ir(1),j) = true;   end,    end

% Convert each kept component, placed alone in an 'opvar' with the
% selector of its block, to an 'sopvar' block
C = cell(sum(rows),sum(cols));
ri = 0;
for i = 1:2
    if ~rows(i),    continue,   end
    ri = ri+1;      ci = 0;
    for j = 1:2
        if ~cols(j),    continue,   end
        ci = ci+1;
        comp = Pop.(nm{i,j});
        if ~keep(i,j)
            continue                            % structurally zero block
        end
        Pblk = opvar();
        Pblk.I = Pop.I;     Pblk.var1 = Pop.var1;   Pblk.var2 = Pop.var2;
        Pblk.dim = sel{i,j};
        Pblk.(nm{i,j}) = comp;
        C{ri,ci} = opvar2sopvar(Pblk);
    end
end
% Metadata of the container, stated from the opvar (a row or column may
% hold no block). The registry holds var1 only if an L2 space survives.
isL2 = [false, true];
if any(rows & isL2) || any(cols & isL2)
    vars = {char(Pop.var1.varname{1})};     dom = double(Pop.I(1,:));
    so = isL2(rows)';                       si = isL2(cols)';
else
    vars = cell(1,0);                       dom = zeros(0,2);
    so = false(nnz(rows),0);                si = false(nnz(cols),0);
end
meta = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',si,...
    'dim_out',d(rows,1),'dim_in',d(cols,2));
Pmop = copvar(C,meta);

end

% ------------------------------------------------------------------------
function tf = is_zero_component(comp)
% True when a component carries nothing. R is a struct of R0, R1, R2, so a
% scalar emptiness test does not reach it.
% INPUT:  comp - component of an 'opvar': 'double', 'polynomial', or a
%                'struct' of these (checked field by field)
% OUTPUT: tf   - logical scalar
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

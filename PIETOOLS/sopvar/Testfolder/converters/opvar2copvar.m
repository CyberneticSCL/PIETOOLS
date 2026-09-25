function Pmop = opvar2copvar(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMOP = OPVAR2COPVAR(POP) takes an 'opvar' object POP and returns the
% equivalent 2 x 2 'copvar' container.
%
% INPUTS
% - Pop:    'opvar' object, P: R^m x L2^n[s] -> R^m x L2^n[s]. Unlike
%           'opvar2sopvar', ALL FOUR components may be nonempty - handling
%           the full 4-PI operator is the point of this routine;
%
% OUTPUTS
% - Pmop:   'copvar' over the spaces (R^m, L2^n[s]), with
%
%               Pmop.C{1,1} = P    (R^m  -> R^m)
%               Pmop.C{1,2} = Q1   (L2^n -> R^m)
%               Pmop.C{2,1} = Q2   (R^m  -> L2^n)
%               Pmop.C{2,2} = R    (L2^n -> L2^n)
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
%   P  [m m;0 0]     Q1 [m 0;0 n]     Q2 [0 m;n 0]     R  [0 0;n n]
%
% A component that is empty or identically zero is left as a structurally
% zero block ([]) rather than converted, so a sparse operator stays sparse
% and the container does not carry bases for blocks that hold nothing.
%
% A zero-dimensional side is skipped the same way: for m = 0 the whole first
% row and column are absent, and the result is a 1 x 1 container. Rows and
% columns of a 'copvar' cannot be entirely empty, since the container reads
% each row's space and dimension off a populated block, so an operator with
% m = 0 or n = 0 must not produce a 2 x 2 grid.
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

if ~isa(Pop,'opvar')
    error('opvar2copvar:badInput','Input must be an opvar object.')
end
d = Pop.dim;
if d(1,1)~=d(1,2) || d(2,1)~=d(2,2)
    error('opvar2copvar:nonSquare',...
        ['This routine expects an opvar with equal input and output '...
         'dimensions on each space, i.e. dim = [m m; n n]; got %s.'],mat2str(d))
end
m = d(1,1);     n = d(2,1);

% Which of the four spaces are present at all. A zero-dimensional space is
% not a space: including it would give the container an all-empty row.
keep_R  = m>0;
keep_L2 = n>0;
if ~keep_R && ~keep_L2
    error('opvar2copvar:empty','Operator has no dimensions to convert.')
end

% One selector per component, in the [out_R in_R; out_L2 in_L2] layout.
sel = {[m m;0 0], [m 0;0 n];
       [0 m;n 0], [0 0;n n]};
nm  = {'P' ,'Q1';
       'Q2','R'};
rows = [keep_R, keep_L2];       cols = [keep_R, keep_L2];

C = cell(sum(rows),sum(cols));
ri = 0;
for i = 1:2
    if ~rows(i),    continue,   end
    ri = ri+1;      ci = 0;
    for j = 1:2
        if ~cols(j),    continue,   end
        ci = ci+1;
        comp = Pop.(nm{i,j});
        if is_zero_component(comp)
            continue                            % structurally zero block
        end
        Pblk = opvar();
        Pblk.I = Pop.I;     Pblk.var1 = Pop.var1;   Pblk.var2 = Pop.var2;
        Pblk.dim = sel{i,j};
        Pblk.(nm{i,j}) = comp;
        C{ri,ci} = opvar2sopvar(Pblk);
    end
end
Pmop = copvar(C);

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

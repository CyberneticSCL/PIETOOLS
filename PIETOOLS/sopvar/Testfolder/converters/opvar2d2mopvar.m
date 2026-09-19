function Pmop = opvar2d2mopvar(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMOP = OPVAR2D2MOPVAR(POP) takes an 'opvar2d' object and returns the
% equivalent 4 x 4 'mopvar' container.
%
% INPUTS
% - Pop:    'opvar2d' object over R x L2[x] x L2[y] x L2[x,y]. Unlike
%           'opvar2d2sopvar', all sixteen components may be nonempty;
%
% OUTPUTS
% - Pmop:   'mopvar' over the four spaces in opvar2d's own order,
%           (R, L2[x], L2[y], L2[x,y]), so that
%
%               C{i,j}  <->  Pop.R<out><in>,  out,in in {0,x,y,2}
%
%           Row is output, column is input. That is opvar2d's component
%           naming read as a grid, and it is already the convention of
%           'opvar2d2sopvar', whose internal name table is exactly
%
%               R00 R0x R0y R02
%               Rx0 Rxx Rxy Rx2
%               Ry0 Ryx Ryy Ry2
%               R20 R2x R2y R22
%
% NOTES
% Each component is delegated to 'opvar2d2sopvar' rather than converted by
% hand. That routine takes a single-component 'opvar2d' and identifies which
% component is live from the dim matrix, so each of the sixteen is placed in
% an otherwise empty 'opvar2d' whose 4 x 2 dim selects just that block.
%
% Components that are empty or identically zero are left as structurally
% zero blocks. A PIE's 2D operator is mostly zero, so this matters: the
% container then carries no basis for them and 'mtimes' skips their
% products entirely.
%
% Spaces with dimension zero are dropped, along with their row and column,
% because a 'mopvar' row must have a populated block to define its space
% and dimension. The surviving spaces keep opvar2d's relative order.
%
% The cell-valued components (Rxx and R2x are 3 x 1, Ryy and R2y are 1 x 3,
% R22 is 3 x 3) are copied whole; the 3-way split is the alpha index of the
% PI operator in that direction and 'opvar2d2sopvar' maps it into the
% sopvar parameter cell.
%
% See also MOPVAR2OPVAR2D, OPVAR2D2SOPVAR, OPVAR2MOPVAR, RAND_MOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar2d2mopvar
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

if ~isa(Pop,'opvar2d')
    error('opvar2d2mopvar:badInput','Input must be an opvar2d object.')
end
nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
d = Pop.dim;
if ~isequal(size(d),[4,2])
    error('opvar2d2mopvar:badDim','opvar2d dim must be 4x2; got %s.',mat2str(size(d)))
end
rows = d(:,1)>0;        cols = d(:,2)>0;
if ~any(rows) || ~any(cols)
    error('opvar2d2mopvar:empty','Operator has no dimensions to convert.')
end

C = cell(sum(rows),sum(cols));
ri = 0;
for i = 1:4
    if ~rows(i),    continue,   end
    ri = ri+1;      ci = 0;
    for j = 1:4
        if ~cols(j),    continue,   end
        ci = ci+1;
        comp = Pop.(nm{i,j});
        if is_zero_component(comp)
            continue                        % structurally zero block
        end
        Pblk = opvar2d();
        Pblk.I = Pop.I;     Pblk.var1 = Pop.var1;   Pblk.var2 = Pop.var2;
        sel = zeros(4,2);
        sel(i,1) = d(i,1);  sel(j,2) = d(j,2);
        Pblk.dim = sel;
        Pblk.(nm{i,j}) = comp;
        C{ri,ci} = opvar2d2sopvar(Pblk);
    end
end
Pmop = mopvar(C);

end

% ------------------------------------------------------------------------
function tf = is_zero_component(comp)
% True when a component carries nothing. Seven of the sixteen are CELLS of
% polynomials (the alpha index of a PI direction), so a scalar emptiness
% test does not reach them.
if iscell(comp)
    tf = true;
    for k = 1:numel(comp),  tf = tf && is_zero_component(comp{k});  end
    return
end
if isempty(comp)
    tf = true;      return
end
comp = polynomial(comp);
tf = isempty(comp.coefficient) || ~any(comp.coefficient(:));
end

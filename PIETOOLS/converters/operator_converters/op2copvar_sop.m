function Pm = op2copvar_sop(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM = OP2COPVAR_SOP(POP) converts a PI operator to a 'copvar' container:
%   opvar, dopvar -> the 2 x 2 container over (R^m, L2[s]);
%   opvar2d    -> a copvar over the present spaces of POP (R^n, L2[x],
%                 L2[y], L2[x,y]), one block per nonzero component, with a
%                 zero block kept in every all-zero row and column so that
%                 no space is dropped; with no space at all on one side,
%                 a container with no space on that side;
%   copvar/cdopvar -> returned as is.
%
% INPUT
% - Pop:    'opvar', 'dopvar', 'opvar2d', 'copvar' or 'cdopvar' object;
%           other classes throw an error.
%
% OUTPUT
% - Pm:     'copvar' (or 'cdopvar') object representing the same operator,
%           one space per nonzero dimension of POP.dim, in the order R^n,
%           L2[x], L2[y], L2[x,y].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - op2copvar_sop
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MMP, 10/08/2026

% Containers are returned as is; 1-D operators convert block by block
if isa(Pop,'copvar') || isa(Pop,'cdopvar'),     Pm = Pop;   return,     end
if isa(Pop,'opvar') || isa(Pop,'dopvar')
    Pm = opvar2copvar(Pop);     return
end
if ~isa(Pop,'opvar2d')
    error('op2copvar_sop:class','POP should be an opvar, opvar2d or a container.');
end
% Names of the 16 components, row = output space, column = input space
nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
% Present output and input spaces, and the nonzero components between them
d = Pop.dim;
rows = find(d(:,1)>0);  cols = find(d(:,2)>0);
keep = false(4,4);
for i = rows',  for j = cols',  keep(i,j) = ~is_zero(Pop.(nm{i,j}));  end,  end
% One zero block in each row and column left empty; a side with no space
% has no block to place
for i = rows',  if ~isempty(cols) && ~any(keep(i,cols)),  keep(i,cols(1)) = true;   end,    end
for j = cols',  if ~isempty(rows) && ~any(keep(rows,j)),  keep(rows(1),j) = true;   end,    end
% Convert each kept component, placed alone in an 'opvar2d' with the
% dimensions of its row and column, to an 'sopvar' block
C = cell(numel(rows),numel(cols));
for a = 1:numel(rows)
    for b = 1:numel(cols)
        i = rows(a);    j = cols(b);
        if ~keep(i,j),  continue,   end
        Pblk = opvar2d();
        Pblk.I = Pop.I;     Pblk.var1 = Pop.var1;   Pblk.var2 = Pop.var2;
        sel = zeros(4,2);   sel(i,1) = d(i,1);  sel(j,2) = d(j,2);
        Pblk.dim = sel;
        Pblk.(nm{i,j}) = Pop.(nm{i,j});
        C{a,b} = opvar2d2sopvar(Pblk);
    end
end
% Assemble the container
if ~isempty(rows) && ~isempty(cols)
    Pm = copvar(C);
else
    % No space on one side (an empty channel of a PIE): no block to read
    % the spaces off, so they are stated from POP.dim. Space k of R, L2[x],
    % L2[y], L2[x,y] holds the variables of row k of MASK; the registry is
    % those the present spaces hold, sorted, with their domains
    nmv = reshape(Pop.var1.varname,1,[]);
    mask = logical([0 0; 1 0; 0 1; 1 1]);
    used = any(mask(unique([rows;cols]),:),1);
    [vs,ix] = sort(nmv(used));
    dI = Pop.I(used,:);         mk = mask(:,used);      mk = mk(:,ix);
    meta = struct('vars',{vs},'dom',dI(ix,:),'space_out',mk(rows,:),...
                  'space_in',mk(cols,:),'dim_out',d(rows,1),'dim_in',d(cols,2));
    Pm = copvar(C,meta);
end
end


function tf = is_zero(comp)
% TF = IS_ZERO(COMP) is true if the component COMP carries no nonzero
% coefficient.
% INPUT:  comp - component of an 'opvar2d': 'double', 'polynomial', or a
%                cell of these (checked element by element)
% OUTPUT: tf   - logical scalar
if iscell(comp)
    tf = true;
    for k = 1:numel(comp),  tf = tf && is_zero(comp{k});  end
    return
end
if isempty(comp),   tf = true;  return,     end
comp = polynomial(comp);
tf = isempty(comp.coefficient) || ~any(comp.coefficient(:));
end

function val = getfield_opvar2d(P,fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% val = getfield_opvar2d(P,fname) returns the value that the field fname
% would have if the PI operator
%   P: R^m0 x L2^mx[x] x L2^my[y] x L2^m2[x,y]
%                           -> R^n0 x L2^nx[x] x L2^ny[y] x L2^n2[x,y]
% were a 2025 'dopvar2d' object, with components
%   [ R00  R0x  R0y  R02 ]
%   [ Rx0  Rxx  Rxy  Rx2 ]
%   [ Ry0  Ryx  Ryy  Ry2 ]
%   [ R20  R2x  R2y  R22 ],
% Rij mapping the space j to the space i, e.g.
%   (Rxy v)(x) = int_c^d Rxy(x,y) v(y) dy,
%   (Rxx v)(x) = Rxx{1}(x) v(x) + int_a^x Rxx{2}(x,t) v(t) dt
%                               + int_x^b Rxx{3}(x,t) v(t) dt,
% and R22{i,j} the kernel with the multiplier (i,j=1), the integral over
% [a,x] or [c,y] (i,j=2), or over [x,b] or [y,d] (i,j=3) in x and in y.
%
% INPUT
% P:      dopvar class object with the 2-D field layout: at most
%         two spatial variables, and the block rows (columns) ordered by
%         space, first R, then L2[x], then L2[y], then L2[x,y]
% fname:  char or string, a field name of the 2-D layout: 'R00', 'R0x',
%         'R0y', 'R02', 'Rx0', 'Rxx', 'Rxy', 'Rx2', 'Ry0', 'Ryx', 'Ryy',
%         'Ry2', 'R20', 'R2x', 'R2y', 'R22', 'I', 'var1', 'var2' or 'dim'
%
% OUTPUT
% val:    the value of the field:
%         'Rij':    ni x mj 'dpvar' if neither x nor y is in both spaces i
%                   and j; else a cell of these, 3 x 1 if x is in both
%                   ('Rxx','Rx2','R2x'), 1 x 3 if y is ('Ryy','Ry2','R2y'),
%                   and 3 x 3 if both are ('R22'). It is in var1 (x,y) and,
%                   in a direction with an integral over part of the
%                   domain, in var2, and affine in the decision variables
%                   of P.Zd that appear in it;
%         'I':      2 x 2 array [a,b; c,d], the domain P.dom of x and y;
%         'var1':   2 x 1 'polynomial', the spatial variables [x;y];
%         'var2':   2 x 1 'polynomial', the dummy variables, named as var1
%                   with '_dum' appended;
%         'dim':    4 x 2 array [n0,m0; nx,mx; ny,my; n2,m2];
%
% NOTES:
% - x is P.vars{1} and y is P.vars{2}. With fewer
%   variables, 'I', 'var1', 'var2' have fewer rows, and the components on
%   an absent space have zero rows or columns.
% - The block rows of each space are stacked in order into the rows of
%   that space, likewise for the block columns. The structure of P ensures
%   this keeps the order of the components.
% - A component is assembled from the kernels of the blocks, so it is in
%   the output variables and the dummy variables of the input variables.
%   The parameters of a block form a 3 x ... x 3 cell, one direction per
%   variable in both spaces, sorted, so x before y as in the 2025 class. A variable only in
%   the input space is integrated over its whole domain, and the 2025 class
%   names it as var1, so its dummy is renamed (e.g. Rxy in x and y, not
%   y_dum).
% - A block may be 'sdopvar', 'sopvar' (no decision variables, kernels
%   'polynomial') or []. Each component is returned as 'dpvar', as the
%   components of the 2025 class are, also if no block of it has decision
%   variables.
% - A structurally zero block ([]) gives zeros.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getfield_opvar2d (dopvar)
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
%
% DJ, 11/15/2026: Initial coding;

% Check the structure of P, giving the space of each block row and column,
% and the field name
[is_2d,idx_out,idx_in] = is_opvar2d(P);
if ~is_2d
    error('dopvar:notOpvar2d',['The operator does not have the 2-D field layout (R00 to R22) of an operator '...
        'R x L2[x] x L2[y] x L2[x,y] -> R x L2[x] x L2[y] x L2[x,y]; see ''is_opvar2d''.'])
end
if ~(ischar(fname) || (isstring(fname) && isscalar(fname)))
    error('dopvar:badField','The field name must be a char or string.')
end
fname = char(fname);

% Names of the variables and their dummies, as 'sdopvar' names them in kernels
vname = reshape(P.vars,[],1);           dname = strcat(vname,'_dum');
% Component (i,j) maps space j to space i, the spaces 1 to 4 being R,
% L2[x], L2[y], L2[x,y]; has_x(i) (has_y(i)) is true if space i has x (y).
nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
has_x = [false,true,false,true];        has_y = [false,false,true,true];

% Extract the desired field. Use the assemble_kernels helper to extract
% the kernels associated with a particular component.
switch fname
    case 'I'
        val = P.dom;
    case 'var1'
        val = polynomial(vname);
    case 'var2'
        val = polynomial(dname);
    case 'dim'
        val = zeros(4,2);
        for t=1:4
            val(t,:) = [sum(P.dim_out(idx_out==t)), sum(P.dim_in(idx_in==t))];
        end
    otherwise
        [i,j] = find(strcmp(nm,fname));
        if isempty(i)
            error('dopvar:badField',['''%s'' is not in the 2-D field layout; use one of '...
                'R00, R0x, R0y, R02, Rx0, Rxx, Rxy, Rx2, Ry0, Ryx, Ryy, Ry2, R20, '...
                'R2x, R2y, R22, I, var1, var2, dim.'],fname)
        end
        % One parameter per multiplier/integral index of x (rows of the
        % cell) and y (columns) if the variable is in both spaces
        sz = [1+2*(has_x(i) && has_x(j)), 1+2*(has_y(i) && has_y(j))];
        K = assemble_kernels(P,idx_out==i,idx_in==j,sz);
        % A variable only in the input space keeps its own name in the 2025 class
        rn = find([has_x(j) && ~has_x(i), has_y(j) && ~has_y(i)]);
        for k=1:numel(K)
            vnames = K{k}.varname;
            if isempty(rn) || isempty(vnames),  continue,   end
            for t=rn
                vnames(strcmp(vnames,dname{t})) = vname(t);
            end
            K{k} = dpvar(K{k}.C,K{k}.degmat,vnames,K{k}.dvarname,K{k}.matdim);
        end
        if isscalar(K)
            val = K{1};
        else
            val = K;
        end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function K = assemble_kernels(P,rows,cols,sz)
% K = assemble_kernels(P,rows,cols,sz) stacks the kernels of the blocks
% P.C(rows,cols) into one 'dpvar' per parameter.
%
% INPUT
% P:    dopvar class object with the 2-D field layout
% rows: M x 1 logical, the block rows to stack
% cols: N x 1 logical, the block columns to stack
% sz:   1 x 2 array, the size of the cell of parameters of each block:
%       [1,1], [3,1], [1,3] or [3,3]
%
% OUTPUT
% K:    sz cell of sum(P.dim_out(rows)) x sum(P.dim_in(cols)) 'dpvar',
%       K{k} the kernel of parameter k of every block, with zeros for a
%       structurally zero block
%
% NOTES
% The kernels are computed using the '@sdopvar/kernels' method ('dpvar'),
% or '@sopvar/kernels' for a block without decision variables
% ('polynomial'), converted to 'dpvar' with 'poly2dpvar', so that every
% kernel is 'dpvar'. Their cell has the size of the parameters, 3 x 1 for
% one shared variable, and is reshaped to sz. The 'dpvar' concatenation
% merges the decision variables of the blocks.

ir = find(rows);    jc = find(cols);
K = cell(sz);
% No block: zeros of the size of the absent (or empty) space
if isempty(ir) || isempty(jc)
    K(:) = {dpvar(zeros(sum(P.dim_out(rows)),sum(P.dim_in(cols))))};
    return
end
Kblk = cell(numel(ir),numel(jc),prod(sz));
for a=1:numel(ir)
    for b=1:numel(jc)
        blk = P.C{ir(a),jc(b)};
        if isempty(blk)
            Kblk(a,b,:) = {dpvar(zeros(P.dim_out(ir(a)),P.dim_in(jc(b))))};
        else
            Kb = blk.kernels;
            if isa(blk,'sopvar')
                % 'polynomial' kernels; {}: no variable is a decision
                % variable (poly2dpvar(K) would take any 'coeff*' name as one)
                Kb = cellfun(@(Kk) poly2dpvar(Kk,{}),Kb,'UniformOutput',false);
            end
            Kblk(a,b,:) = Kb(:);
        end
    end
end
for k=1:numel(K)
    Krows = cell(numel(ir),1);
    for a=1:numel(ir)
        Krows{a} = horzcat(Kblk{a,:,k});
    end
    K{k} = vertcat(Krows{:});
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

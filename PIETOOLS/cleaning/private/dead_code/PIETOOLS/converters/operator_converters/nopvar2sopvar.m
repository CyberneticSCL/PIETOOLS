function obj = nopvar2sopvar(obj_ndopvar)
% OBJ = SOPVAR2OPVAR(OBJSOPVAR) takes a nopvar object representing a 4-PI
% operator component and returns an sopvar object representing the same
% operator.
%
% INPUTS
% - objndopvar:  'nopvar' object representing a nd PI operator. It can not
%               map between different function space, 
%               (i.e. L2^n to L2^n);
%
% OUTPUTS
% - obj:        'sopvar' object representing the same operator as the input;
%                        in the dsopvar form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ndopvar2sdopvar
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
% AT, 08/05/2026: Initial coding 


% Extract the dimension, variables and domain of the operator

if ~isa(obj_ndopvar, 'nopvar') 
    error('The input should be ndopvar')
end
P = obj_ndopvar;
dims = [P.dim(1), P.dim(2)];
dom = struct('in', [P.dom], 'out', [P.dom]); % maps space to itself
% ndopvar include dummy variables
% P.vars(:, 2) -- only dummy variables
vars = struct('in', {P.vars(:, 1).varname}, 'out', {P.vars(:, 1).varname}); 
% now construct monomial basis
% ZL = ZR;
ZL = cell(1, length(P.deg));
ZR = cell(1, length(P.deg));
for deg_idx = 1:length(P.deg)
    ZL{deg_idx} = (0:P.deg(deg_idx))';
end
ZR = ZL; 

% now transform coefficient matrix 
% in ndopvar it is (Ik o [1;d])^T Cj in R^{k times }
% Cj has the size dim(1)*len(Zl)*(len(dvarnames)+1) times (len(Zr)*dim(2))
% Cj(1:dim(1)*len(Zl), :) are our A

sz_C = size(P.C);
N = size(P.dom,1);
left_size  = dims(1)*prod(cellfun(@(x) length(x), ZL));
right_size = dims(2)*prod(cellfun(@(x) length(x), ZR));
params = cell(size(P.C));

% In nopvar P.C have different dimensions
% In sopvar .params have the same dimensions for all {ii}
for ii=1:numel(P.C) 
    % new C is dim(1)*len(Zl)*(len(dvarnames)+1) times (len(Zr)*dim(2))
    % need to determinated which columns to choose
    % Determine the index of element ii along each dimension of the cell C
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs); % array of indeces 
    % If element ii corresponds to an integral, we need to account for the
    % monomial basis in the associated dummy variable
    is_int = logical(idcs-1);
    column_idx = 1;% column_idx indicates 
    for dim_idx = 1:length(P.deg)
        if is_int(dim_idx)
            column_idx = kron(column_idx, ones(1, P.deg(dim_idx) + 1));
        else
            var_temp= zeros(1, P.deg(dim_idx) + 1);
            var_temp(1) = 1;
            column_idx = kron(column_idx, var_temp); 
        end
    end
    column_idx = kron(ones(1, dims(2)), column_idx);
    [~, cols_mon, ~] = find(column_idx);
    [rows, cols, vs] = find(P.C{ii});
    new_cols = cols_mon(cols);
    C_new = sparse(rows, new_cols, vs, left_size, right_size);
    params{ii} = C_new;
end

obj = sopvar(params, vars,  ZL, ZR, dom, dims);

end
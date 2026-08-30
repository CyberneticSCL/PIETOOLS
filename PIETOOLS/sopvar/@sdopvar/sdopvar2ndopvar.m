function obj = sdopvar2ndopvar(objsdopvar)
% OBJ = SOPVAR2OPVAR(OBJSOPVAR) takes a sdopvar object representing a 4-PI
% operator component and returns an ndopvar object representing the same
% operator.
%
% INPUTS
% - objSopvar:  'sdopvar' object representing a 1D PI operator. It can not
%               map between different function space, 
%               (i.e. it maps L2^n to L2^n);
%
% OUTPUTS
% - obj:        'ndopvar' object representing the same operator as the input;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sopvar2opvar
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


if ~isa(objsdopvar, 'sdopvar') 
    error('The input should be sdopvar')
end



P = objsdopvar;
dims = [P.dims(1), P.dims(2)];
dvarname = P.Zd; %
if ~isequal(P.dom.in, P.dom.out)
    error('Input/Output domains dismatch')
end
if any(~strcmp(char(P.vars.in),char(P.vars.out)))
    error('Input/Output vars dismatch')
end
% ndopvar include dummy variables
% P.vars(:, 2) -- only dummy variables
vars_old = string(char(P.vars.in));

for idx = 1:length(vars_old)
    [var1, var1_dummy] = pvar(vars_old{idx}, [vars_old{idx}, '_dum']);
    if idx == 1
        vars_new = [var1, var1_dummy];
    else
        vars_new = [vars_new; [var1, var1_dummy]];
    end
end
% now construct monomial basis
% ZL = ZR;
ZL = P.ZL;
ZR = P.ZR;
degL = cellfun(@(x) max(x), ZL);
degR = cellfun(@(x) max(x), ZL);
if ~isequal(degR, degL)
    error('left and right monomials have different degrees')
end
deg = degL;
% now transform coefficient matrix 
% in ndopvar it is (Ik o [1;d])^T Cj in R^{k times }
% Cj has the size dim(1)*len(Zl)*(len(dvarnames)+1) times (len(Zr)*dim(2))
% Cj(1:dim(1)*len(Zl), :) are our A

N = size(P.dom.in,1);
left_size_block = dims(1)*prod(cellfun(@(x) length(x), ZL));
left_size_full = left_size_block*(length(dvarname) + 1);
right_size = dims(2)*prod(cellfun(@(x) length(x), ZR));

A = P.params.A;
B= P.params.B;
sz_C = size(A);
n_dvarnames = length(dvarname);
C_new = cell(size(A));
for ii=1:numel(A) 
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
    for dim_idx = 1:length(deg)
        if is_int(dim_idx)
            column_idx = kron(column_idx, ones(1, deg(dim_idx) + 1));
        else
            var_temp= zeros(1, deg(dim_idx) + 1);
            var_temp(1) = 1;
            column_idx = kron(column_idx, var_temp); 
        end
    end
    column_idx = kron(ones(1, dims(2)), column_idx);
    [~, cols_mon, ~] = find(column_idx);

    A_sub = reshape(A{ii}, left_size_block, right_size);
    B_sub = reshape(B{ii}.', left_size_block*n_dvarnames, right_size);

    C_sub = [A_sub; B_sub];
    C_new{ii} = C_sub(:, cols_mon); 
    % cdim = n*prod(deg(is_int)+1);
    % % Set sparse coefficients of dimension rdim x cdim
    % rho = (q+10)/(rdim*cdim);
    % Pop.C{ii} = sprand(rdim,cdim,rho);
end


obj = ndopvar(); % empty nopvar
obj.dom =  P.dom.in;
obj.deg =  deg;
obj.vars = vars_new;
obj.dvarname = dvarname;
obj.C = C_new;

end
function C = otimes(A,B)
% C = OTIMES(A,B) returns the tensor product of coefficient operators A and
% B, so that (A*Zd(x)).*(B*Zd(y)) = C*(Zd(x) o Zd(y))
%
% INPUTS
% - A:  m x n1 'intop' object
% - B:  m x n2 'intop' object
%
% OUTPUTS
% - C:  m x n1*n2 'intop' object representing the tensor product of A
%       and B
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - otimes
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
% DJ, 03/01/2026: Initial coding


% Make sure the row dimensions of the functionals match
[m1,n1] = size(A);
[m2,n2] = size(B);
if m1~=m2
    if m1==1 
        A.params = ones(m2,1)*A;
    elseif m2~=1
        B.params = ones(m1,1)*B;
    else
        error("Row dimension of the operators must match")
    end 
end

% Only support scalar-valued monomials for now
if n1>1 || n2>1
    error("Multiplication of functionals on vector-valued monomials is currently not supported.")
end

% For multiplication with degree-0 monomial, we can just multiply the
% coefficients
if isa(A,'double') || isa(A,'polynomial') || isa(A,'dpvar')
    C = B;
    C.params = A.*B.params;
    return
elseif ~isa(A,'intop')
    error("Functionals must be specified as 'intop' objects.")
end
if isa(B,'double') || isa(B,'polynomial') || isa(B,'dpvar')
    C = A;
    C.params = A.params.*B;
    return
elseif ~isa(B,'intop')
    error("Functionals must be specified as 'intop' objects.")
end

% Make sure the domains of the variables match
if any(A.dom~=B.dom)
    error("Domains of integration of the functionals must match.")
end

% Extract parameters
Aparams = A.params;     omat1 = A.omat;
Bparams = B.params;     omat2 = B.omat;

% Declare a full set of dummy variables for integration 
var1 = A.pvarname;      Apar_varname = Aparams.varname;
var2 = B.pvarname;      Bpar_varname = Bparams.varname;
var_new = cell(1,numel(var1)+numel(var2));
for j=1:numel(var1)
    % Generate a new name for the jth variable
    var_new{j} = ['s_',num2str(j)];
    % Replace variable in the parameters defining A
    var_idx = ismember(Apar_varname,var1(j));
    if any(var_idx)
        Apar_varname(var_idx) = var_new(j);
    end
end
Aparams.varname = Apar_varname;
for j=1:numel(var2)
    % Generate a new name for the (j+k)th variable
    var_new{j+numel(var1)} = ['s_',num2str(j+numel(var1))];
    % Replace variable in the parameters defining B
    var_idx = ismember(Bpar_varname,var2(j));
    if any(var_idx)
        Bpar_varname(var_idx) = var_new(j+numel(var1));
    end
end
Bparams.varname = Bpar_varname;

% Update variable indices
omat2 = omat2 + numel(var1);

omat_new = zeros(0,2*numel(var_new)+1);
params_new = zeros(size(Bparams,1),0);
for i=1:size(omat1,1)
    % Establish all orders of integration in the new variable list based on
    % the order of integration in term i of A
    omat_i = omat2;
    prev_idx = -ones(size(omat2,1),1);
    for j=1:numel(var1)
        % Place variable j in every possible position
        omat_j = omat_i;
        omat_i = zeros(0,size(omat_j,2)+1);
        prev_idx_new = zeros(0,1);
        for k=0:size(omat_j,2)
            % Check for which rows in omat_j we can place var j in position
            % k without violating the ordering imposed by A
            is_geq = k>prev_idx;
            nr = sum(is_geq);
            omat_tmp = [omat_j(is_geq,1:k),omat1(i,j)*ones(nr,1),omat_j(is_geq,k+1:end)];
            omat_i = [omat_i; omat_tmp];
            prev_idx_new = [prev_idx_new; k*ones(nr,1)];
        end
        prev_idx = prev_idx_new;
    end
    % Add the order to the full list of terms
    omat_new = [omat_new; omat_i];

    % Declare the parameters associated with each integral
    params_i = Aparams(:,i).*Bparams;
    params_new = [params_new,repmat(params_i,[1,size(omat_i,1)/size(omat2,1)])];
end

% Declare the operator representing the tensor product
C = intop(params_new,omat_new,var_new,A.dom);

end
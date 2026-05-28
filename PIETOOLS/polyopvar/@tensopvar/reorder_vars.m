function Cop_out = reorder_vars(Cop,var_order)
% COP_OUT = REORDER_VARS(COP,VAR_ORDER) takes a tensor-PI operator
% COP and reorders the integrals defining this functional to match the
% order VAR_ORDER = [i1,...,id] of the state variables
%   COP_OUT*(x_{i1} o ... o x_{id}) = COP*(x1 o ... o xd)
% where o denotes the tensor product.
%
% INPUTS
% - Cop:        m x n 'tensopvar' object representing a tensor-PI operator
%               on a degree-d distributed monomial
% - var_order:  1 x d array spcifying the new order of the state variables
%               in the distributed monomial;
%
% OUTPUTS
% - Cop_out:    m x n 'tensopvar' object representing the same tensor-PI
%               operator as the input but now acting on the reordered
%               distributed monomial;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - reorder_vars
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
% DJ, 04/21/2026: Initial coding

% % Check the inputs
if nargin<=1
    error("Not enough inputs.")
end
if ~isa(Cop,'tensopvar')
    error("Tensor-PI operator must be specified as 'tensopvar' object.")
end
nvars = size(Cop.ops,2);
if ~isa(var_order,'double')
    error("Order of variables must be specified as  d x 1 array of integers.")
elseif numel(var_order)~=nvars || ~isempty(setdiff(1:nvars,var_order))
    error("A new index must be specified for each of the variables in the functional.")
end

% % Reorder the factors in the tensor-PI operator
Cop_out = Cop;
if Cop.type(2) && size(Cop,2)>1
    % For Kronecker products, we cannot simply change the order of
    % the operators
    Cop_out.order = var_order;
else
    Cop_out.ops = Cop.ops(:,var_order);
end

end
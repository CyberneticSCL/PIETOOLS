function Kop_out = reorder_vars(Kop,var_order)
% KOP_OUT = REORDER_VARS(KOP,VAR_ORDER) reorders the integrals defining the
% PI functional KOP to match the specified order VAR_ORDER, so that
%   KOP_OUT*(x_{i1} o ... o x_{id}) = KOP*(x1 o ... o xd)
% where o denotes the tensor product.
%
% INPUTS
% - Kop:        m x n 'intop' object representing a functional on
%               L2[[a1,b1]^{d1} x ... x [aN,bN]^{dN}];
% - var_order:  1 x N cell, with each element a 1 x di array specifying the
%               new order of the different degrees of each spatial
%               variable. If N=1, may also be specified as a 1 x d1 array;
%
% OUTPUTS
% - Kop_out:    m x n 'intop' object representing the same functional as
%               the input, acting on the reordered monomial;
%
% NOTES
% - For univariate, scalar-valued states, x1(s1),...,xd(s1), reordering the
%   spatial variables in Kop is equivalent to reordering the state
%   variables,
%      KOP_OUT*(x_{i1} o ... o x_{id}) = KOP*(x1 o ... o xd)
%   for VAR_ORDER = [i1,...,id]. However, for multivariate or vector-valued
%   state variables, such a reordering of the spatial variables is not
%   sufficient.
% - The dummy variables of a variable are the columns of the monomial
%   degrees Z.in{i}, so reordering them only permutes these columns. The
%   coefficient matrices Kop.params are not changed.

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
% DJ, 10/15/2026: Initial coding;

% Check that the functional is properly specified
if nargin<=1
    error("Not enough inputs.")
end
if ~isa(Kop,'intop')
    error("Functional must be specified as 'intop' object.")
end
if ~Kop.chkval && ~isvalid(Kop)
    error("The functional is not properly declared.")
end
N = numel(Kop.vars.in);
degs = Kop.depmat.in;

% Check that the variable order is properly specified
if isnumeric(var_order) && N==1
    var_order = {var_order};
elseif ~isa(var_order,'cell') || numel(var_order)~=N
    error("Order of variables must be specified as 1 x N cell for N spatial variables.")
end
is_id = false(1,N);
for i=1:N
    if ~isnumeric(var_order{i}) || numel(var_order{i})~=degs(i)
        error("Element 'var_order{i}' must be a 1 x di array of integers, for degree di along variable i.")
    elseif any(~ismember(var_order{i},1:degs(i))) || numel(unique(var_order{i}))~=degs(i)
        error("Element 'var_order{i}' must be a permutation of the vector (1:di), for degree di along variable i.")
    end
    var_order{i} = reshape(var_order{i},1,[]);
    if isequal(var_order{i},1:degs(i))
        is_id(i) = true;
    end
end
if all(is_id)
    % The order of the variables does not actually change
    Kop_out = Kop;
    return
end

Kop_out = Kop;
for i=find(~is_id)
    % The dummy variable that was called si_dumk is now called si_dumj, with
    % j the position of k in the new order: permute the monomial degrees
    Kop_out.Z.in{i} = Kop.Z.in{i}(:,var_order{i});
    % Update the domains of integration to match the new variable names
    % Extract the lower and upper bounds of the domain
    intmat_i = Kop.intmat{i};
    Amat_i = intmat_i(:,1:2:end);
    Bmat_i = intmat_i(:,2:2:end);
    % Account for updated order
    Amat_i = Amat_i(:,var_order{i});
    Bmat_i = Bmat_i(:,var_order{i});
    [~,old_order] = sort(var_order{i});
    Amat_i(Amat_i>0) = old_order(Amat_i(Amat_i>0));
    Bmat_i(Bmat_i<=degs(i)) = old_order(Bmat_i(Bmat_i<=degs(i)));
    intmat_i(:,1:2:end) = Amat_i;
    intmat_i(:,2:2:end) = Bmat_i;
    Kop_out.intmat{i} = intmat_i;
end
Kop_out.chkval = true;

end

function g = weighted_sobolev_norm(u,alpha)
% G = WEIGHTED_SOBOLEV_NORM(U,ALPHA) returns the squared weighted Sobolev
% norm
%   g(x) = sum_i sum_k alpha{i}(k_1+1,...,k_N+1) * ||D^k u_i||_{L2}^2
% of the PDE state u = T*x, expressed in terms of the fundamental state x.
% Here, u_i is the ith component of u, a function of N spatial variables
% s_1,...,s_N (N=0 if u_i takes values in R), the sum is over the
% multi-indices k = (k_1,...,k_N) of derivative orders, and
%   D^k u_i = d^{k_1}/ds_1^{k_1} ... d^{k_N}/ds_N^{k_N} u_i.
% The components of u may depend on different numbers of spatial variables.
%
% INPUTS
% - u:      'polyopvar' representing the state, with m components u_i: the
%           blocks of rows of u.C (m = numel(u.dims.out), which equals
%           size(u,1) if all components are scalar-valued). Component u_i
%           is a function of the N spatial variables in u.vars.out for
%           which u.depmat.out(i,:) is nonzero (in the order of u.vars.out),
%           e.g. N = 1, 2, 0 for u = [u1(s1); u2(s1,s2); u3]. All entries
%           of a vector-valued component get the same weights;
% - alpha:  m x 1 cell, with alpha{i} specifying the weights of the
%           derivatives of u_i: alpha{i}(k_1+1,...,k_N+1) is the weight of
%           ||D^k u_i||^2, where k_j is the order of the derivative in the
%           jth variable of u_i, as ordered above. For derivatives up to order
%           (d_1,...,d_N), alpha{i} is a (d_1+1) x ... x (d_N+1) array,
%           so a scalar if N=0, a vector if N=1, and a matrix if N=2.
%           Derivatives with zero weight are not taken and need not exist:
%           for an L2 state u_i, use alpha{i} = 1 (or any nonnegative
%           scalar) for the weight of ||u_i||^2.
%           Alternatively, alpha can be a single array, used for all
%           components: a scalar then weights ||u_i||^2 for each
%           component, whereas an array with N>0 dimensions requires each
%           component to depend on N variables;
%
% OUTPUTS
% - g:      scalar 'polyopvar' object representing the squared weighted
%           Sobolev norm (the number 0 if all weights are zero);
%
% EXAMPLE
% For u = [u1(s1); u2(s1,s2); u3], with u3 in R, and
%   alpha = {[1;0;2]; [1,3;0,0]; 5},
% the output is
%   g = ||u1||^2 + 2*||d^2 u1/ds1^2||^2 + ||u2||^2 + 3*||d u2/ds2||^2 + 5*u3^2.
%
% NOTES
% - The derivatives are taken of the output of u, see 'diff' for
%   'polyopvar', so they must again be operators acting on x: an error is
%   thrown for a derivative of nonzero weight if u_i involves a multiplier
%   along the differentiated variable (a derivative of x would be needed);
% - g is a norm only if all weights are nonnegative, which is not checked;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - weighted_sobolev_norm
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
% If you modify this code, document all changes carefully and
% include date authorship, and a brief description of modifications
%
% CR, DJ, 11/09/2026: Initial coding;

if ~isa(u,'polyopvar')
    error("weighted_sobolev_norm: the state must be specified as 'polyopvar' object.")
end

% A single array of weights is used for all components
m = numel(u.dims.out);
if ~iscell(alpha)
    alpha = repmat({alpha},m,1);
elseif numel(alpha)~=m
    error("weighted_sobolev_norm: alpha must specify the weights for each of the %d components of u.",m)
end

g = 0;
for i=1:m
    A = alpha{i};
    if ~any(A(:))
        continue                % no terms: also allows L2 states, for any alpha{i}=0
    end

    % Component i only, as the other components may depend on other
    % variables, or be L2 states that cannot be differentiated. Its
    % variables are those of u_i.vars.out with nonzero depmat.out, as
    % u_i.vars.out lists the variables of all operators in u.
    u_i = u{i};
    if any(u_i.depmat.out>1)
        error("weighted_sobolev_norm: component %d of u is a tensor product of functions of the same variable, which is not supported.",i)
    end
    svars = u_i.vars.out(u_i.depmat.out>0);
    N = numel(svars);
    if ~isnumeric(A) || ndims(A)>max(N,2) || (N<2 && ~isvector(A)) || (N==0 && ~isscalar(A))
        error("weighted_sobolev_norm: the weights of component %d of u, depending on N=%d spatial variables, must be a scalar for N=0, a vector for N=1, or an array with at most N dimensions for N>1.",i,N)
    end
    if N==1
        A = A(:);               % column, so that the index is the order+1
    end

    % Loop over the derivatives with nonzero weight, k(j) the order in svars{j}
    for idx=find(A(:))'
        sz = [size(A,1:N),1,1]; % ind2sub needs at least two dimensions
        k = cell(1,numel(sz));
        [k{:}] = ind2sub(sz,idx);
        k = [k{1:N}] - 1;
        du_i = u_i;
        if any(k)
            try
                du_i = diff(u_i,svars(k>0),k(k>0));
            catch ME
                error("weighted_sobolev_norm: derivative of order %s of component %d of u is not defined (%s)",mat2str(k),i,ME.message)
            end
        end
        g = g + A(idx)*innerprod(du_i,du_i);
    end
end

end

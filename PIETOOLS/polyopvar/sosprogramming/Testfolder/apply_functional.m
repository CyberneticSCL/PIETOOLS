function fval = apply_functional(Kop,xvals,degmat)
% FVAL = APPLY_FUNCTIONAL(KCELL,XVALS,IDX_MAT,VARS,DOM) computes the value
% of the integral
% f(x) = sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           KCELL{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1
% where kj = idx_ii(i,j) for j=1,...,d and i=1,...,m for m, for the
% given values of the polynomial functions xi;
%
% INPUTS
% - Kop:    struct with fields
%               params: m x q*n 'polynomial' or 'dpvar' object,
%                       with elements (1:m,(i-1)*n+1:i*n) corresponding to 
%                       the kernels in the ith term of the functional
%                       in linear format
%               omat:   q x p array of integers, where p is the number of
%                       dummy variables for integration, with row i
%                       specifying the order of the variables in the
%                       integral associated with the ith kernel.
%                       Specifically, if omat(i,:) = [k1,k2,...,kd], 
%                       then term i is defined by the integral
%                           int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b};
%               matdim: 1 x 2 array specifying the dimensions of the
%                       kernels K{i};
%               vars:   p x 1 'polynomial' array specifying the names of
%                       the dummy variables used in definition of the
%                       kernels in Kop.params;
%               dom:    interval [a,b] over which to integrate;
% - xvals:  d x 1 array of type 'polynomial' 
% - degmat: 1 x d array of integers, specifying the degrees of the
%           state variables appearing in the distributed monomials;
%
% OUTPUTS
% - fval:   The value of the function f(x) for the specified K and x;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_functional
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
% DJ, 01/26/2026: Initial coding

d = sum(degmat);
nstates = size(degmat,2);
if size(xvals,2)~=nstates
    error("Distributed monomial must be specified as 1 x d array")
end
if size(xvals,1)~=1
    error("Vector-valued state variables are currently not supported.")
end
if size(degmat,1)~=1
    error("Polynomials involving multiple monomials are not supported.")
end

Kparams = Kop.params;
idx_mat = Kop.omat;
vars = Kop.vars;
dom = Kop.dom;
mdim = Kop.matdim(1);
ndim = Kop.matdim(2);

% Establish for each factor in the monomial which state variable is
% considered
state_idcs = [];
for i=1:nstates
    state_idcs = [state_idcs; i*ones(degmat(i),1)];
end
xvals_full = polynomial(zeros(ndim,d));
for i=1:d
    state_num = state_idcs(i);
    var_ii = xvals(:,state_num).varname;
    if isscalar(var_ii)
        xvals_full(:,i) = xvals(:,state_num);
        xvals_full(:,i).varname = vars(i).varname;
    elseif isempty(var_ii)
        xvals_full(:,i) = xvals(state_num);
    elseif numel(var_ii)>1
        error("Each state variable can depend on at most one independent variable")
    end
end

if isempty(idx_mat)
    % List all possible orders of the variables t1 through td
    %   idx_mat(l,:) = [i,j,k] means a <= ti <= tj <= tk <= b in term l
    idx_mat = 1;
    for i=2:d
        n_ords = size(idx_mat,1);
        idx_mat_new = zeros(i*n_ords,i);
        for j=1:i
            % Place variable t_ii in position jj
            idx_mat_new(j:i:end,:) = [idx_mat(:,1:j-1),i*ones(n_ords,1),idx_mat(:,j:end)];
        end
        idx_mat = idx_mat_new;
    end
end


fval = zeros(mdim,ndim);
n_terms = size(Kparams,2)/ndim;
if d>2 && n_terms>factorial(d)
    error("Number of operators should be d! for monomial degree d.")
end
%Kf = Kparams*prod(xvals_full(1,:),2);   % requires state variables to be scalar!
if all(idx_mat(1,:)==0)
    % Account for multiplier term
    has_multiplier = true;
    fval_i = Kparams(:,1:ndim);
    for j=2:d
        fval_i = fval_i.*subs(xvals_full(:,j),vars(j),vars(1));
    end
    fval_i = int_simple(fval_i.*xvals_full(:,1),vars(1),dom(1),dom(2));
    fval = fval + fval_i;
else
    has_multiplier = false;
end
for i=has_multiplier+1:size(idx_mat,1)
    % Check the order of the variables:
    %   idx_ii = [i,j,k] implies a <= ti <= tj <= tk <= b
    idx_ii = idx_mat(i,:);
    
    % Perform the appropriate integrals
    %   int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
    %           Kcell{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1 
    % where kj = idx_ii(i,j) for j=1,...,d and i=1,...,m for m=d!.
    fval_i = Kparams(:,(i-1)*ndim+1:i*ndim);
    for j=d:-1:2
        var_num = idx_ii(j);   % integrating over the jth biggest variable
        L = vars(idx_ii(j-1)); % integrating from the j-1th biggest variable up to b
        fval_i = fval_i*xvals_full(:,var_num);
        fval_i = int_simple(fval_i,vars(var_num),L,dom(2));
    end
    % Finally, perform the integral from a to b in the smallest variable
    var_num = idx_ii(1); 
    fval_i = fval_i*xvals_full(:,var_num);
    fval_i = int_simple(fval_i,vars(var_num),dom(1),dom(2));
    fval = fval + fval_i;
end

if d==2 && n_terms==size(idx_mat,1)+1
    % Third term corresponds to the multiplier
    %   int_{a}^{b} P(s)*x(s)^2 ds
    fval = fval + int(Kparams(:,end-ndim+1:end)*xvals_full(:,1).*subs(xvals_full(:,2),vars(2),vars(1)),vars(1),dom(1),dom(2));
end


end
function Kop = combine_terms(Kop,state_idcs)
% KOP = COMBINE_TERMS(KOP,STATE_IDCS) takes a functional defined by KOP,
% acting on a distributed monomial, and combines different terms in the 
% functional to account for the fact that x(t1)*x(t2)=x(t2)*x(t1), so that
% e.g.
%   int_{a}^{b} int_{a}^{t1} K1(t1,t2)*x(t1)*x(t2) dt2 dt1
%       + int_{a}^{b} int_{a}^{t2} K2(t1,t2)*x(t1)*x(t2) dt1 dt2
%   = int_{a}^{b} int_{a}^{t1} [K1(t1,t2)+K2(t2,t1)]*x(t1)*x(t2) dt2 dt1
%
% INPUTS
%  OPTION 1
% - Kop:        'polyopvar' object representing a distributed polynomial
%               functional defined by a single monomial, so that
%               size(Kop.degmat,1) = 1. In this case, the degrees of the
%               different state variables in the monomial will be checked
%               to determine which terms in the integral can be combined;
%
%  OPTION 2
% - Kop:        struct representing a functional acting on a distributed
%               monomial with d factors. This struct must have fields:
%               params: 1 x d 'polynomial' or 'dpvar' object,
%                       with each element i corresponding to 
%                       the kernels in the ith term of the functional;
%               omat:   q x d array of integers, with row i
%                       specifying the order of the variables in the
%                       integral associated with the ith kernel.
%                       Specifically, if omat(i,:) = [k1,k2,...,kd], 
%                       then term i is defined by the integral
%                           int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b};
%               vars:   d x 1 'polynomial' array specifying the names of
%                       the dummy variables used in definition of the
%                       kernels in Kop.params;
%               dom:    interval [a,b] over which to integrate;
% - state_idcs: 1 x d array of integers specifying for each of the 
%               independent variables over which to integrate what state
%               variable in the monomial it corresponds, so that if
%               state_idcs(i) = j_{i} for each i, then
%                   Kop*Z(x) 
%                       = int K(t1,...,td)*x_{j_{i}}(t1)*...*x_{j_{i}}(td);
%
% OUTPUTS
%  OPTION 1
% - Kop:        'polyopvar' object representing the same polynomial
%               functional as the input, but now accounting for 
%               non-uniqueness of the monomials;
%
%  OPTION 2
% - Kop:        struct with same elements as the input, representing the
%               same functional operator, but now accounting for
%               non-uniqueness of the monomials. In particular, different
%               rows of Kop.omat corresponding to the same integral (for
%               the given monomial) will have been combined, and the
%               associated elements of Kop.params will have been added.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - combine_terms
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


% Allow functional to be specified as 'polyopvar' object.
if nargin==1 && isa(Kop,'polyopvar')
    if size(Kop.degmat,1)>1
        error("Polynomial must be defined in terms of a single monomial.")
    end
    Kfun = Kop;
    state_idcs = 1:size(Kop.degmat,2);
    state_idcs = repelem(state_idcs,Kop.degmat);
    Kop = Kop.C.ops{1};
end

% Extract pararmeters defining the functional operator
omat = Kop.omat;
vars = Kop.vars;
Kparams = Kop.params;
ndim = Kop.matdim(2);
% Account for possible multiplier term in case of monomial of degree 2
Kparams1 = Kparams(:,1:size(omat,1)*ndim);
Kparams2 = Kparams(:,size(omat,1)*ndim+1:end);

% Establish the order of the independent variables in each term
state_nums = unique(state_idcs);
[~,var_orders] = sort(omat,2);

for ii=1:numel(state_nums)
    % For each state variable, check if it is integrated with respect to
    % multiple variables
    var_idcs_ii = find(state_idcs==state_nums(ii));
    if isscalar(var_idcs_ii)
        continue
    end
    for jj=1:numel(var_idcs_ii)-1
        % For each pair of variables, (var1,var2), with respect to which x
        % is integrated, retain only the integrals with var2 <= var1
        var1_idx = var_idcs_ii(jj);
        var1 = vars(var1_idx);
        for kk=jj+1:numel(var_idcs_ii)
            var2_idx = var_idcs_ii(kk);
            var2 = vars(var2_idx);
            % Find all terms with var2 >= var1
            var1_order = var_orders(:,var1_idx);
            var2_order = var_orders(:,var2_idx);
            is_greater = var2_order>var1_order;
            % Reorder such that var2 <= var1
            lindcs1 = sub2ind(size(omat),find(is_greater),var1_order(is_greater));
            lindcs2 = sub2ind(size(omat),find(is_greater),var2_order(is_greater));
            omat(lindcs1) = var2_idx;
            omat(lindcs2) = var1_idx;
            var_orders(is_greater,[var1_idx,var2_idx]) = var_orders(is_greater,[var2_idx,var1_idx]);
            is_greater_params = find(logical(kron(is_greater,ones(1,ndim))));
            Kparams1(:,is_greater_params) = var_swap(Kparams1(:,is_greater_params),var1,var2);
            % Combine terms involving the same integral
            [P,omat] = uniquerows_integerTable(omat);
            Kparams1 = Kparams1*kron(P,speye(ndim));
            P = P.*(1./sum(P,1));
            var_orders = P'*var_orders;
        end
    end
end

% Return the functional operator with combined terms
Kop.omat = omat;
Kop.params = [Kparams1, Kparams2];
ztol  = 1e-12;
if isa(Kop.params,'double')
    Kop.params(abs(Kop.params)<ztol) = 0;
else
    Kop.params.C(abs(Kop.params.C)<ztol) = 0;
end
if nargin==1
    % If the input was a 'polyopvar' object, we return a 'polyopvar' object
    % representing the same polynomial.
    Kfun.C.ops{1} = Kop;
    Kop = Kfun;
end

end
function dP = diff(P,vars,d)
% DP = DIFF(P,VARS,D) returns the derivative of the output of the
% distributed polynomial P with respect to the output variables VARS,
%   DP(x) = d^D/(d vars^D) P(x).
% VARS and D can be row vectors to differentiate with respect to several
% variables at once, and D a matrix to take several derivatives, stacked
% vertically.
%
% INPUTS
% - P:      m x 1 'polyopvar' object, not an adjoint;
% - vars:   1 x r 'cellstr', char, or 'polynomial' (pvar) array, the
%           output variables of P with respect to which to differentiate;
% - d:      (optional) q x r array of nonnegative integers, row k the
%           orders of differentiation in the variables vars for the k-th
%           derivative. A scalar is used for all variables. Defaults to 1;
%
% OUTPUTS
% - dP:     q*m x 1 'polyopvar' object, the q derivatives of P stacked
%           vertically, dP(x) = [D1 P(x); ...; Dq P(x)], for the
%           differential operators D1,...,Dq defined by the rows of d;
%
% NOTES
% - The output of P is differentiated, dP(x) = (D C)*Z(x), with the same
%   monomials Z. Differentiation with respect to the states (the
%   derivative of the polynomial in x) is not supported.
% - The derivatives of the operators C must be tensor-PI operators: an
%   error is thrown e.g. if the output of P involves a multiplier along a
%   variable in VARS, as the derivative then involves derivatives of the
%   states.
% - The derivative with respect to a variable that is not an output
%   variable of P is zero.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - diff
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
% DJ, 10/20/2026: Initial coding;

% Check that the polynomial is properly specified
if ~isa(P,'polyopvar')
    error("Distributed polynomial must be specified as object of type 'polyopvar'.")
elseif ~P.chkval
    P = initialize(P);
end
if P.is_adjoint
    error("polyopvar:adjoint_not_supported","'diff' is not supported for adjoint 'polyopvar' objects (u'), see 'ctranspose'.")
end

% Check that the variables and orders of differentiation are properly
% specified
if nargin<=2
    [vars,d] = parse_diff_inputs(vars);
else
    [vars,d] = parse_diff_inputs(vars,d);
end
if any(ismember(vars,P.varname))
    error("Differentiation with respect to state variables is not supported.")
end
p = size(d,1);

% Differentiation with no variables returns the original polynomial
if isempty(vars)
    dP = repmat(P,p,1);
    return
end

% Differentiate the operators acting on the monomials, which stacks the
% outputs of the p derivatives. The monomials remain the same.
dP = P;
dP.C = diff(P.C,vars,d);
% The output space is that of the new operators: let initialize infer the
% fields describing it, and check the rest.
fdims = P.dims;         fdims.out = zeros(0,1);
fvars = P.vars;         fvars.out = cell(1,0);
fdom = P.dom;           fdom.out = zeros(0,2);
fdepmat = P.depmat;     fdepmat.out = zeros(0,0);
dP.dims = fdims;        dP.vars = fvars;
dP.dom = fdom;          dP.depmat = fdepmat;
dP = initialize(dP);

end

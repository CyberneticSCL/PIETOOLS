function dP = diff(P,vars,d)
% DP = DIFF(P,VARS,D) takes the composition of the differential operator
%   d^D/(d vars^D)
% in the output variables VARS with the distributed polynomial P(x),
%   DP(x) = d^D/(d vars^D) P(x)
% VARS and D can be row vectors to simultaneously differentiate with
% respect to multiple variables. D can also be a column vector or matrix to
% take multiple derivatives, concatenated vertically.
%
% INPUTS
% - P:      m x 1 'polyopvar' object representing a distributed polynomial
%           P(x) = C*Z(x) in the states x;
% - vars:   1 x r 'cellstr' or 'polynomial' (pvar) array specifying the
%           names of the spatial variables, in terms of which the output of
%           P is defined, with respect to which to take the derivative;
% - d:      (optional) p x r array of nonnegative integers specifying the
%           order of the derivatives to be taken with respect to each of
%           the spatial variables. Defaults to 1;
%
% OUTPUTS
% - dP:     m*p x 1 'polyopvar' representing the derivatives of P, stacked
%           vertically, so that dP(x) = [D1 P(x); ...; Dp P(x)] for the
%           differential operators D1,...,Dp defined by the rows of d;
%
% NOTES
% - The derivative is taken of the output of P, by composing the
%   tensor-PI operators C acting on the monomials Z(x) with the
%   differential operator, see 'diff' for 'tensopmat': dP(x) = (D C)*Z(x)
%   with the same monomials Z. Differentiation with respect to the state
%   variables of P, that is, the derivative of the polynomial in x, is not
%   supported;
% - Derivatives of the tensor-PI operators must again be tensor-PI
%   operators, see 'diff' for 'tensopvar' and 'sopvar'. For instance,
%   an error is thrown if the output of P involves a multiplier along a
%   variable in VARS, as then the derivative involves derivatives of the
%   states;
% - A variable that is not an output variable of P is allowed, returning
%   a zero polynomial;
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
% DJ, 10/20/2026: Initial coding;

% Check that the polynomial is properly specified
if ~isa(P,'polyopvar')
    error("Distributed polynomial must be specified as object of type 'polyopvar'.")
elseif ~P.chkval
    P = initialize(P);
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

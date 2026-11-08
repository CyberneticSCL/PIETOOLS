function P_out = lift_output(P,vars,dom,dims,deps)
% P_OUT = LIFT_OUTPUT(P,VARS,DOM,DIMS,DEPS) returns the 'polyopvar' object
% P_OUT representing the distributed polynomial P with its output lifted
% to the function space defined by VARS, DOM and DEPS: the output is
% extended with variables it was not yet defined in, along which it is
% constant, so that
%   P_out(x)(s,t) = P(x)(s),
% for s the output variables of P, and t the new output variables. This
% allows e.g. the addition of a polynomial in an ODE state, with output in
% R, to a polynomial in a PDE state, with output in L2[a,b].
%
% INPUTS
% - P:      m x 1 'polyopvar' object;
% - vars:   1 x M cellstr specifying the names of the output variables of
%           the lifted polynomial. Must contain all output variables of P.
%           Can also be specified as struct with field 'out' containing
%           this cellstr (other fields, e.g. 'in', are ignored);
% - dom:    M x 2 'double' array specifying the intervals on which these
%           variables are defined. Must match the domains of the output
%           variables of P. Can also be specified as struct with field
%           'out' (other fields are ignored);
% - dims:   (optional) p x 1 array of integers (or struct with field
%           'out') specifying the dimension of each row-block of the output
%           (see P.dims.out), or a scalar for the same dimension of all
%           row-blocks. Defaults to P.dims.out. A scalar-valued row-block
%           can be lifted to a vector-valued one, repeating the scalar
%           output in each row. If P has a single (scalar-valued)
%           row-block, dims can be a q x 1 array, to repeat the row-block q
%           times, with dimension dims(i) for the ith copy;
% - deps:   (optional) p x M array of nonnegative integers (or struct with
%           field 'out') specifying the degree of each row-block of the
%           output in each of the variables vars (see P.depmat.out), or a
%           1 x M array for the same degrees of all row-blocks. Must be at
%           least the degrees of P, and defaults to these degrees, or to 1
%           for variables that a row-block is not yet defined in;
%
% OUTPUTS
% - P_out:  m x 1 'polyopvar' object representing the lifted polynomial,
%           with the same state variables and monomials as P, and with
%               P_out.vars.out = sort(vars);
%               P_out.depmat.out = deps (sorted accordingly);
%
% NOTES
% - The coefficient operator P.C is lifted by 'lift_output' for 'tensopmat'
%   objects, see there for the lifting of the different types of blocks.
%   The input function spaces, the states, and the monomials are
%   unchanged.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lift_output
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
% DJ, 11/03/2026: Initial coding;
% DJ, 11/06/2026: Output space of P_out follows from P_out.C (see set.C);
% DJ, 11/07/2026: Reject adjoints, see 'ctranspose';

% Check that the polynomial is properly declared
if ~isa(P,'polyopvar')
    error("First argument must be of type 'polyopvar'.")
end
if P.is_adjoint
    error("polyopvar:adjoint_not_supported","'lift_output' is not supported for adjoint 'polyopvar' objects (u'), see 'ctranspose'.")
end
if ~P.chkval
    P = initialize(P);
end
if nargin<=3
    dims = [];
end
if nargin<=4
    deps = [];
end

% Lift the coefficient operator, which determines the output space of the
% polynomial (all inputs are checked there)
C = lift_output(P.C,vars,dom,dims,deps);

% Set the output space of the polynomial to that of the new operator, by the
% assignment of P_out.C, see set.C. The input space, and the monomials, are
% unchanged
P_out = P;
P_out.C = C;
P_out.chkval = false;
P_out = initialize(P_out);

end

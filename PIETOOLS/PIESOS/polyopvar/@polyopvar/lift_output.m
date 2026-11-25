function P_out = lift_output(P,vars,dom,dims,deps)
% P_OUT = LIFT_OUTPUT(P,VARS,DOM,DIMS,DEPS) lifts the output of the
% distributed polynomial P to the function space defined by VARS, DOM and
% DEPS, extending it to the new variables, along which it is constant,
%   P_out(x)(s,t) = P(x)(s),
% for s the output variables of P, and t the new output variables. This
% allows e.g. adding a polynomial in an ODE state, with output in R, to a
% polynomial in a PDE state, with output in L2[a,b].
%
% INPUTS
% - P:      m x 1 'polyopvar' object, not an adjoint, with nr row-blocks;
% - vars:   1 x M 'cellstr', or 'struct' with field 'out', the output
%           variables of the lifted polynomial, including all output
%           variables of P;
% - dom:    M x 2 'double' array, or 'struct' with field 'out', the
%           domains of these variables, matching those of P;
% - dims:   (optional) nr x 1 array of integers (or 'struct' with field
%           'out'), the dimension of each row-block of the output, or a
%           scalar for all row-blocks. Defaults to P.dims.out. A
%           scalar-valued row-block can be lifted to a vector-valued one,
%           repeating the scalar in each row. If P has a single
%           scalar-valued row-block, dims can be a q x 1 array, repeating
%           the row-block q times, with dimension dims(i) for copy i;
% - deps:   (optional) nr x M array of nonnegative integers (or 'struct'
%           with field 'out'), the degree of each row-block in each of the
%           variables vars, or a 1 x M array for all row-blocks. At least
%           the degrees of P. Defaults to these degrees, and to 1 for
%           variables that a row-block is not yet defined in;
%
% OUTPUTS
% - P_out:  'polyopvar' object representing the lifted polynomial, with
%           the same states and monomials as P, and with
%               P_out.vars.out = sort(vars);
%               P_out.depmat.out = deps (sorted accordingly);
%
% NOTES
% - A row-block containing a functional ('intop') cannot be lifted along
%   any variable, and one consisting only of blocks of type 'double',
%   'polynomial' or 'dpvar' only to degree at most 1 in each variable.
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/03/2026: Initial coding;

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
% assignment of P_out.C. The input space, and the monomials, are
% unchanged
P_out = P;
P_out.C = C;
P_out.chkval = false;
P_out = initialize(P_out);

end

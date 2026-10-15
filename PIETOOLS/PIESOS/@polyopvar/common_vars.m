function [A_out,B_out] = common_vars(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_VARS(A_IN,B_IN) takes two distributed polynomials
% and expresses them in terms of the same state variables and the same
% output variables.
%
% INPUTS
% - A_in, B_in:     m x n 'polyopvar' objects;
%
% OUTPUTS
% - A_out, B_out:   m x n 'polyopvar' objects representing the same
%                   distributed polynomials as A_in, B_in, respectively,
%                   but expressed in terms of the same variables, so that
%                   A_out.varname = B_out.varname;
%                   A_out.vars = B_out.vars;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_vars
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

% Unify the state variables (varname/dims.in/vars.in/dom.in/depmat.in/Z)
[A_out,B_out] = common_states(A_in,B_in);

% Also unify the output space: output variables will follow from
% coefficients
[C_A,C_B] = common_vars(A_out.C,B_out.C);
A_out.C = C_A;                       B_out.C = C_B;
A_out.vars.out = C_A.vars.out;       B_out.vars.out = C_B.vars.out;
A_out.dom.out = C_A.dom.out;         B_out.dom.out = C_B.dom.out;
A_out.depmat.out = C_A.depmat.out;   B_out.depmat.out = C_B.depmat.out;

end
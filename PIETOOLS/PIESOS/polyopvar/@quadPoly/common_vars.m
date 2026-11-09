function [A_out,B_out] = common_vars(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_VARS(A_IN,B_IN) takes two 'quadPoly' objects and
% expresses them in terms of the same left- and right-variables.
%
% INPUTS
% - A_in, B_in:     'quadPoly' objects representing polynomials in
%                   (potentially) different variables
%
% OUTPUTS
% - A_out, B_out:   'quadPoly' objects representing the same polynomials as
%                   the inputs, but now expressed in terms of the same 
%                   variables, so that
%                       A_out.ns = B_out.ns;    A_out.nt = B_out.nt;
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
% DJ, 06/19/2026: Initial coding

% Extract the variable names
vars1_A = A_in.ns;      vars1_B = B_in.ns;
vars2_A = A_in.nt;      vars2_B = B_in.nt;

% Combine into unique list
vars1 = unique([vars1_A,vars1_B]);
vars2 = unique([vars2_A,vars2_B]);

% Set the new variables of each object
A_out = set_vars(A_in,vars1,vars2);
B_out = set_vars(B_in,vars1,vars2);

end
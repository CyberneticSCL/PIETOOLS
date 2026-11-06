function [dP,is_zero] = diff(P,vars,d)
% [DP,IS_ZERO] = DIFF(P,VARS,D) takes the composition of the PI operator P
% (a 'copvar' object) with a Dth order differential operator in variables
% VARS,
%   DP = d^D/(d vars^D) o P,
% for each of the output spaces of P. VARS and D can be row vectors to
% simultaneously differentiate with respect to multiple variables. D can
% also be a column vector or matrix to take multiple derivatives,
% concatenated vertically. See the 'diff' method of 'sopvar' for the
% operators that are supported.
%
% INPUTS
% - P:      M x N 'copvar' object representing a PI operator;
% - vars:   1 x r 'cellstr' or 'polynomial' (pvar) array specifying the
%           names of the spatial variables with respect to which to take
%           the derivative;
% - d:      (optional) p x r array of nonnegative integers specifying the
%           order of the derivatives to be taken with respect to each of
%           the spatial variables. Defaults to 1;
%
% OUTPUTS
% - dP:     M x N 'copvar' object representing the composition of the
%           requested differential operator with P, with the same input
%           spaces and output spaces as P, and with p*q_i components in the
%           i-th output space, for q_i those of P;
% - is_zero:    boolean value set to true if dP is known to be the zero
%               operator (false does not imply it is not zero);
%
% NOTES
% - Block row i of P acts into L_2^{q_i}[s^i], and (P*u)_i depends only on
%   the variables s^i. A variable that is not in s^i therefore does not
%   affect (P*u)_i, and block row i of dP is zero for any derivative in that
%   variable. In particular, a variable that is in none of the output spaces
%   gives the zero operator;
% - Each nonzero block P.C{i,j} is differentiated by the 'diff' method of
%   'sopvar', so the derivative of block (i,j) is the composition of the
%   differential operator with that block alone. Since (P*u)_i is the sum of
%   P.C{i,j}*u_j, this is the derivative of the full operator;
% - A multiplier along the variable of differentiation cannot be
%   differentiated, in which case an error is thrown by the 'diff' method of
%   'sopvar';
% - Zero blocks remain []. A block that differentiates to zero stays an
%   explicit zero 'sopvar', since a container may not have a row or column
%   without any block;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - diff (copvar)
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
% DJ, 11/06/2026: Initial coding;

if ~isa(P,'copvar')
    error("PI operator must be specified as object of type 'copvar'.")
end

% Check and normalize the variables and orders once, so that the 'sopvar'
% blocks can skip the check
if nargin<=2
    [vars,d] = parse_diff_inputs(vars);
else
    [vars,d] = parse_diff_inputs(vars,d);
end
p = size(d,1);

% Differentiate each nonzero block; 'sopvar' returns a zero block for the
% blocks that do not depend on the variables.
C = P.C;
is_zero = true;
for ii=find(~cellfun(@isempty,C(:)))'
    [C{ii},is_zero_ii] = diff(C{ii},vars,d,true);
    is_zero = is_zero && is_zero_ii;
end

% The p derivatives are stacked in the output of each block, q_i --> p*q_i.
% Spaces and domains do not change, so the metadata of P is reused
meta = metadata(P);
meta.dim_out = p*P.dim_out;
dP = copvar(C,meta);

end

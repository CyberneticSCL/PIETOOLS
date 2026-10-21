function dKop = diff(Kop,vars,d)
% DKOP = DIFF(KOP,VARS,D) takes the composition of the functional KOP with
% a Dth order differential operator in variables VARS
%   DKOP = d^D/(d vars^D) o KOP
% VARS and D can be row vectors to simultaneously differentiate with
% respect to multiple variables. D can also be a column vector or matrix to
% take multiple derivatives, concatenated vertically.
%
% INPUTS
% - Kop:    m x n 'intop' object representing a functional;
% - vars:   1 x r 'cellstr' or 'polynomial' (pvar) array specifying the
%           names of the spatial variables with respect to which to take
%           the derivative;
% - d:      (optional) p x r array of nonnegative integers specifying the
%           order of the derivatives to be taken with respect to each of
%           the spatial variables. Defaults to 1;
%
% OUTPUTS
% - dKop:   m*p x n 'intop' object representing the composition of the
%           requested differential operator with the functional Kop;
%
% NOTES
% - A functional maps to R^m, which does not depend on any spatial
%   variable, so for any derivative of positive order the composition is
%   the zero operator, with the same input space as Kop. For orders of
%   differentiation all zero, Kop is returned (repeated p times);
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

% Check that the operator is properly specified
if ~isa(Kop,'intop')
    error("Functional must be specified as object of type 'intop'.")
end

% Check that the variables and orders of differentiation are properly
% specified
if nargin<=2
    [vars,d] = parse_diff_inputs(vars);
else
    [vars,d] = parse_diff_inputs(vars,d);
end
p = size(d,1);

dKop = Kop;
if isempty(vars)
    % Differentiation with no variables returns the original functional:
    % stack p copies, as the rows of params are ordered (matrix row,
    % decision variable), so a vertical stack of the functionals is the
    % vertical stack of the parameters.
    dKop.params = cellfun(@(P) repmat(P,p,1),Kop.params,'UniformOutput',false);
else
    % The output of a functional does not depend on the spatial variables:
    % the composition is zero, with the rows of each parameter stacked p
    % times.
    dKop.params = cellfun(@(P) sparse(p*size(P,1),size(P,2)),Kop.params,'UniformOutput',false);
end

end

function vars = native_vars(op)
% VARS = NATIVE_VARS(OP) returns the spatial variables that appear in the
% matrix OP of type 'double', 'polynomial', or 'dpvar', as a sorted 1 x M
% cellstr array. Only variables appearing in a term with a nonzero
% coefficient are counted, and decision variables are not spatial
% variables.
%
% INPUTS
% - op:     m x n matrix of type 'double', 'polynomial' or 'dpvar';
%
% OUTPUTS
% - vars:   1 x M cellstr array with the names of the spatial variables
%           appearing in op, sorted, empty (1 x 0) for a 'double' op;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - native_vars
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

vars = cell(1,0);
if isa(op,'double')
    return
end
deg = full(op.degmat);
if isa(op,'dpvar')
    % Coefficients C: rows (decision variable inner, matrix row outer), columns (monomial inner,
    % matrix column outer), so that column k belongs to monomial mod(k-1,nmon)+1
    nmon = size(deg,1);
    mons = unique(mod(find(any(op.C,1))-1,nmon)+1);
else
    % Coefficients: row k holds the coefficients of monomial k
    mons = find(any(op.coefficient,2));
end
used = any(deg(mons,:)>0,1);
vars = reshape(sort(op.varname(used)),1,[]);

end

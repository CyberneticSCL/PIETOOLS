function H = jacobian(F, vars)
%JACOBIAN Compute the Jacobian of a quadPoly object F.
%
% INPUTS
% - F:      n x 1 or 1 x n 'quadPoly' object;
% - vars:   (optional) 'char', 'string' array or 'cell' of k variable
%           names; defaults to all variables of F, left (F.ns) then
%           right (F.nt);
%
% OUTPUTS
% - H:      'quadPoly' object, the derivatives of F in the k variables,
%           stacked vertically (k x n) for a row vector F, horizontally
%           (n x k) for a column vector F (k x 1 for a scalar F).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - jacobian
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
% AT, DJ, 11/24/2026: Initial coding;

% Only 'quadPoly' objects; by default all variables of F
if ~isa(F, 'quadPoly')
    error('quadPoly:jacobian:badType', 'jacobian supports only quadPoly ');
end

if nargin < 2 || isempty(vars)
    vars = union(F.ns, F.nt, 'stable');
end


% Only vector-valued F
dim = F.dim;

if (dim(1) > 1) && (dim(2) > 1)
    error('quadPoly:jacobian:dimMismatch', 'jacobian supports only 1xn or nx1 quadpoly ')
end

% Variables as a cell array of names
if isa(vars, 'char')
    vars = string(vars);
end

if ~isa(vars, 'cell')
    vars = num2cell(vars);
end

% Concatenate the derivatives in the variables, in order
for var_ind = 1:length(vars)
    var_name = vars{var_ind};

    diff_obj = diff(F, var_name);

    if var_ind == 1
        H = diff_obj;  
    else
        if dim(1) == 1
            H = [H; diff_obj];
        else
            H = [H, diff_obj];
        end
    end
end

end

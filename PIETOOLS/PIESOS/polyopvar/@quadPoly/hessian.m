function H = hessian(F, vars)
%HESSIAN compute hessian of quadpoly object F.
% 
% INPUTS
% - F:    1x1 'quadpoly' class object.
% - vars:  (optional) cell array of names of vars, k of them; defaults to
%          all variables of F, left (ns) then right (nt)
% 
% OUTPUTS
% - H:      'quadpoly' object differentiated with respect all vars: the
%           k x k matrix H(i,j) = d^2 F/(d vars{i} d vars{j})
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% AT, 02/18/2026: Initial coding;


% Only 'quadPoly' objects; by default all variables of F
if ~isa(F, 'quadPoly')
    error('quadPoly:hessian:badType', 'hessian supports only quadPoly ');
end

if nargin < 2 || isempty(vars)
    vars = union(F.ns, F.nt, 'stable');
end


% Only scalar F
dim = F.dim;

if (dim(1) > 1) || (dim(2) > 1)
    error('quadPoly:hessian:dimMismatch', 'hessian supports only 1x1 quadpoly ')
end

% Variables as a cell array of names
if isa(vars, 'char')
    vars = string(vars);
end

if ~isa(vars, 'cell')
    vars = num2cell(vars);
end

% Row i of H: the derivatives of dF/d(vars{i}) in each variable
for var_ind1 = 1:length(vars)
    var_name1 = vars{var_ind1};
    for var_ind2 = 1:length(vars)
        var_name2 = vars{var_ind2};
        diff_obj = diff(F, var_name1);
        diff_obj = diff(diff_obj, var_name2);
        if var_ind2 == 1
            Hrow = diff_obj;  
        else
            Hrow = [Hrow, diff_obj];
        end
    end

    % Stack the rows
    if var_ind1 == 1
        H = Hrow;
    else
        H = [H; Hrow];
    end

end

end

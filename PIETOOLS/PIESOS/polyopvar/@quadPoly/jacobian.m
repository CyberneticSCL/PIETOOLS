function H = jacobian(F, vars)
%JACOBIAN compute jacobian of quadpoly object F.
%
% INPUTS
% - P1:   nx1 or 1xn 'quadpoly' class object.
% - vars:  cell array of names of vars
% 
% OUTPUTS
% - H:      'quadpoly' object differentiated with respect all vars
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

if ~isa(F, 'quadPoly')
    error('quadPoly:jacobian:badType', 'jacobian supports only quadPoly ');
end

if nargin < 2 || isempty(vars)
    vars = union(F.ns, F.nt, 'stable');
end


dim = F.dim;

if (dim(1) > 1) && (dim(2) > 1)
    error('quadPoly:jacobian:dimMismatch', 'jacobian supports only 1xn o nx1 quadpoly ')
end

if isa(vars, 'char')
    vars = string(vars);
end

if ~isa(vars, 'cell')
    vars = num2cell(vars);
end

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

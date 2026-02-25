function [P_diff] = diff(P1, vars)
%DIFF differentiate variables in a quadPoly
% 
% INPUTS
% - P1:   'quadpoly' class object.
% - vars:  cell array of names of vars
% 
% OUTPUTS
% - P_diff:      'quadpoly' object differentiated with respect all vars
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


if isa(P1, 'double')
    P_diff = 0;
    return;
end
if ~isa(P1, 'quadPoly')
    error('quadPoly:diff:badType', 'diff supports only quadPoly or double');
end

if isa(vars, 'char')
    vars = string(vars);
end

if ~isa(vars, 'cell')
    vars = num2cell(vars);
end


% extract monomials and coefficients.
Zt = P1.Zt;
Zs = P1.Zs;
nt = P1.nt;
nt = string(nt);
ns = P1.ns;
ns = string(ns);
C = P1.C;
dim = P1.dim;

% locate var names in ns or nt
% subs_in_ns = {};
% subs_in_nt = {};

vals_in_ns = cell(size(ns));
vals_in_nt = cell(size(nt));
for vars_ind = 1:length(vars)
    ind_in_ns = strcmp(ns, vars{vars_ind});
    if sum(ind_in_ns) > 0
        location_vars_ind = find(ind_in_ns);
        vals_in_ns{location_vars_ind} = 1;
    else
        ind_in_nt = strcmp(nt, vars{vars_ind});

        if sum(ind_in_nt) > 0
            location_vars_ind = find(ind_in_nt);
            vals_in_nt{location_vars_ind} = 1;
        else
            warning(['quadPoly:diff:NotFound ', char(vars{vars_ind}), ' is not in quadPoly. Returned 0'])
            P_diff = 0*P1;
            return
        end
    end
end


left_multiplier = 1;
hatZs = {};
for left_index = 1:length(ns)
    if isempty(vals_in_ns{left_index})
        n_mon = length(Zs{left_index});
        mon = speye(n_mon, n_mon);

        hatZs{end + 1} = Zs{left_index};
    else 
        mon = diag(Zs{left_index});
        mon = mon(:, ~all(mon == 0, 1));

        hatZs{end + 1} = Zs{left_index} -1;
        hatZs{end} = hatZs{end}(hatZs{end} >= 0);
        
    end
    left_multiplier = kron(left_multiplier, mon);
end
left_multiplier = kron(speye(dim(1)), left_multiplier);


hatZt = {};
right_multiplier = 1;
for right_index = 1:length(nt)
    if isempty(vals_in_nt{right_index})
        n_mon =length(Zt{right_index});
        mon = speye(n_mon, n_mon);

        hatZt{end + 1} = Zt{right_index};
    else
        mon = diag(Zt{right_index});
        mon = mon(:, ~all(mon == 0, 1));

        hatZt{end + 1} = Zt{right_index} -1;
        hatZt{end} = hatZt{end}(hatZt{end} >= 0);
    end
    right_multiplier = kron(right_multiplier, mon);
end
right_multiplier = kron(speye(dim(2)), right_multiplier);




hatC = left_multiplier'*C*right_multiplier;

P_diff = quadPoly(hatC, hatZs, hatZt, dim, ns, nt);
% we need to call clean quadpoly after diff
end

function [P_int] = int(P1, vars, low_limit, high_limit)
%int integrate variables in a quadPoly
%
% 
% INPUTS
% - P1:   'quadpoly' class object.
% - vars:  cell array of names of vars
% - low_limit (optional) : 'double' array
% - high_limit (optional): 'double' array
%
% OUTPUTS
% - P_int:      'quadpoly' object integrated with respect all vars
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
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



if ~isa(P1, 'quadPoly')
    error('quadPoly:int:badType', 'int supports only quadPoly');
end

if isa(vars, 'char')
    vars = string(vars);
end

if ~isa(vars, 'cell')
    vars = num2cell(vars);
end


if nargin > 2
    if nargin < 4
        high_limit = low_limit;
        low_limit = 0*high_limit;
    end

    P_int = int(P1, vars);
    P_int = subs(P_int, vars, high_limit) - subs(P_int, vars, low_limit);
    return
end

  

% extract monomials and coefficients.
Zt = P1.Zt;
Zs = P1.Zs;
nt = P1.nt;
nt_str = string(nt);
ns = P1.ns;
ns_str = string(ns);
C = P1.C;
dim = P1.dim;

% locate var names in ns or nt
% subs_in_ns = {};
% subs_in_nt = {};

vals_in_ns = cell(size(ns_str));
vals_in_nt = cell(size(nt_str));
for vars_ind = 1:length(vars)
    ind_in_ns = strcmp(ns_str, vars{vars_ind});
    if sum(ind_in_ns) > 0
        location_vars_ind = find(ind_in_ns);
        vals_in_ns{location_vars_ind} = 1;
    else
        ind_in_nt = strcmp(nt_str, vars{vars_ind});
        
        if sum(ind_in_nt) > 0
            location_vars_ind = find(ind_in_nt);
            vals_in_nt{location_vars_ind} = 1;
        else
            warning(['quadPoly:diff:NotFound ', char(vars{vars_ind}), ' is not in quadPoly'])
        end
    end
end

% Update left-monomials to account for antiderivative
left_multiplier = 1;
hatZs = cell(1,numel(ns));
for left_index = 1:numel(ns)
    if isempty(vals_in_ns{left_index})
        % If variable is not integrated, monomial remains the same
        n_mon = length(Zs{left_index});
        mon = speye(n_mon, n_mon);
        hatZs{left_index} = Zs{left_index};
    else 
        % Set the factor with which the multiply monomial
        mon = diag(1./(Zs{left_index} + 1)); 
        % Increase the monomial degree
        hatZs{left_index} = Zs{left_index} + 1;        
    end
    left_multiplier = kron(left_multiplier, mon);
end
left_multiplier = kron(speye(dim(1)), left_multiplier);


% Update left-monomials to account for antiderivative
hatZt = cell(1,numel(nt));
right_multiplier = 1;
for right_index = 1:numel(nt)
    if isempty(vals_in_nt{right_index})
        % If variable is not integrated, monomial remains the same
        n_mon =length(Zt{right_index});
        mon = speye(n_mon, n_mon);
        hatZt{right_index} = Zt{right_index};
    else
        % Set the factor with which the multiply monomial
        mon = diag(1./(Zt{right_index} + 1));
        % Increase the monomial degree
        hatZt{right_index} = Zt{right_index} + 1; 
    end
    right_multiplier = kron(right_multiplier, mon);
end
right_multiplier = kron(speye(dim(2)), right_multiplier);

% Update coefficients
hatC = left_multiplier'*C*right_multiplier;

% Set the new quadPoly
P_int = quadPoly(hatC, hatZs, hatZt, dim, ns, nt);
% we need to call clean quadpoly after diff
end

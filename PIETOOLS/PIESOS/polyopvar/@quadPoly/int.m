function [P_int] = int(P1, vars, low_limit, high_limit)
%INT Integrate a quadPoly with respect to variables.
%
% INPUTS
% - P1:         m x n 'quadPoly' object;
% - vars:       'char', 'string' array or 'cell' of names of variables; a
%               variable not in P1 is added, on the right (P1.nt) if its
%               name starts with 't', else on the left (P1.ns);
% - low_limit:  (optional) lower limits, one per variable, passed to subs
%               as vals: a 'double' array, or a 'cell' of scalars and
%               variable names; with three inputs, the upper limits, the
%               lower limits then being 0 (numeric limits only);
% - high_limit: (optional) upper limits, as for low_limit;
%
% OUTPUTS
% - P_int:      m x n 'quadPoly' object: without limits, the
%               antiderivative in vars that vanishes at 0; else the
%               iterated integral over vars, in the order given (a limit
%               may be a later variable of vars).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int
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



% Only 'quadPoly' objects are integrated
if ~isa(P1, 'quadPoly')
    error('quadPoly:int:badType', 'int supports only quadPoly');
end

% Variables as a cell array of names
if isa(vars, 'char')
    vars = string(vars);
end

if ~isa(vars, 'cell')
    vars = num2cell(vars);
end


% Definite integral: the antiderivative at the upper minus the lower limits
if nargin > 2
    if nargin < 4
        high_limit = low_limit;
        low_limit = 0*high_limit;
    end

    % One variable at a time, in the order of vars (an iterated integral):
    % G(b)-G(a) in all variables at once would drop the terms in which some
    % variables are at their lower and others at their upper limit
    P_int = P1;
    for k=1:numel(vars)
        G = int(P_int, vars(k));
        P_int = subs(G, vars(k), high_limit(k)) - subs(G, vars(k), low_limit(k));
    end
    return
end


% A variable on which P1 does not depend is added with degree 0, so that
% its integral is var*P1 rather than skipped: on the right if its name
% starts with 't', else on the left
var_names = reshape(cellfun(@char,vars,'UniformOutput',false),1,[]);
new_names = var_names(~ismember(var_names,[P1.ns(:)',P1.nt(:)']));
if ~isempty(new_names)
    is_t = startsWith(new_names,'t');
    P1 = set_vars(P1,sort([P1.ns(:)',new_names(~is_t)]),sort([P1.nt(:)',new_names(is_t)]));
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


% Update right-monomials to account for antiderivative
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
end

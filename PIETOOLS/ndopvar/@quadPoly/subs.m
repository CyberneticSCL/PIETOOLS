function [P_subs] = subs(P1, vars, vals)
%SUBS Substitute variables in a quadPoly with tensor-decomposed exponent bases.
%
% % INPUTS
% - P1:   nx1 or 1xn 'quadpoly' class object.
% - vars:  cell array of names of vars
% - vals:  cell array of 'double'
% % OUTPUTS
% - P_subs:      'quadpoly' object differentiated with respect all vars
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
 

if ~isa(P1, 'quadPoly')
    error('quadPoly:subs:badType', 'subs supports only quadPoly ');
end

if isa(vars, 'char')
    vars = string(vars);
end

if ~isa(vars, 'cell')
    vars = num2cell(vars);
end

if ~isa(vals, 'cell')

    if isa(vals, 'double')
        vals = num2cell(vals);
    elseif isa(vals, 'string')
        vals = {vals};
    elseif isa(vals, 'char')
        vals = {string(vals)};
    else
        error('quadPoly:subs:badType', 'subs supports cell array or double as vals input')
    end
end

if size(vars) ~= size(vals)
    error('quadPoly:subs:badType', 'number of vars and vals should be the same');
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


% check vals for symbolic variables 
isa_symb_var = cellfun(@isa_char_or_string, vals);
if sum(isa_symb_var) > 0
    sym_to_sym_vars = vars(isa_symb_var);
    sym_to_sym_vals = vals(isa_symb_var);
    sym_to_num_vars = vars(~isa_symb_var);
    sym_to_num_vals = vals(~isa_symb_var);
    
    sym_to_sym_vals = cellfun(@char, sym_to_sym_vals, 'UniformOutput',false);

    is_a_ns = cellfun(@(name) find(strcmp(ns, name)), sym_to_sym_vars, 'UniformOutput', false);
    is_a_nt = cellfun(@(name) find(strcmp(nt, name)), sym_to_sym_vars, 'UniformOutput', false);

    chosen_vals_ns = ~cellfun(@isempty, is_a_ns);
    chosen_vals_nt = ~cellfun(@isempty, is_a_nt);

    is_a_ns = cell2mat(is_a_ns);
    is_a_nt = cell2mat(is_a_nt);

    ns_new = ns;
    ns_new(is_a_ns) = sym_to_sym_vals(chosen_vals_ns);
    nt_new = nt;
    nt_new(is_a_nt) = sym_to_sym_vals(chosen_vals_nt);


    P_sub1 = quadPoly(C, Zs, Zt, dim, ns_new, nt_new, 0);
    
    if sum(isa_symb_var) < length(vars)
        P_subs = subs(P_sub1, sym_to_num_vars, sym_to_num_vals);
    else
        P_subs = P_sub1;
    end
    return
end

% locate var names in ns or nt
% subs_in_ns = {};
% subs_in_nt = {};

vals_in_ns = cell(size(ns));
vals_in_nt = cell(size(nt));
for vars_ind = 1:length(vars)
    ind_in_ns = strcmp(ns, vars{vars_ind});
    if sum(ind_in_ns) > 0
        location_vars_ind = find(ind_in_ns);
        vals_in_ns{location_vars_ind} = vals{vars_ind};
    else
        ind_in_nt = strcmp(nt, vars{vars_ind});

        if sum(ind_in_nt) > 0
            location_vars_ind = find(ind_in_nt);
            vals_in_nt{location_vars_ind} = vals{vars_ind};
        else
            warning(['quadPoly:subs:NotFound ', char(vars{vars_ind}), ' is not in quadPoly'])
        end
    end
    
end

left_multiplier = 1;
hatZs = {};
hatns = {};
for left_index = 1:length(ns)
    if isempty(vals_in_ns{left_index})
        n_mon = length(Zs{left_index});
        mon = speye(n_mon, n_mon);

        hatZs{end + 1} = Zs{left_index};
        hatns{end + 1} = ns{left_index};
    else
        value_to_subs = vals_in_ns{left_index};
        mon = value_to_subs.^(Zs{left_index});
    end
    left_multiplier = kron(left_multiplier, mon);
end
left_multiplier = kron(speye(dim(1)), left_multiplier);


hatZt = {};
hatnt = {};
right_multiplier = 1;
for right_index = 1:length(nt)
    if isempty(vals_in_nt{right_index})
        n_mon =length(Zt{right_index});
        mon = speye(n_mon, n_mon);

        hatZt{end + 1} = Zt{right_index};
        hatnt{end + 1} = nt{right_index};
    else
        value_to_subs = vals_in_nt{right_index};
        mon = value_to_subs.^(Zt{right_index});
    end
    right_multiplier = kron(right_multiplier, mon);
end
right_multiplier = kron(speye(dim(2)), right_multiplier);

hatC = left_multiplier'*C*right_multiplier;



P_subs = quadPoly(hatC, hatZs, hatZt, dim, hatns, hatnt);
end


function [ans] = isa_char_or_string(val)
    if isa(val, 'char') || isa(val, 'string')
        ans = true;
    else
        ans = false;
    end
end
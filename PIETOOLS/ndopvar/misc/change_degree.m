function Cop = change_degree(Aop,degree)
% COP = change_degree(Aop,degree) returns the 'nopvar/ndopvar' object Cop with basis of degree d
%
% INPUTS
% - Aop:     'nopvar/ndopvar' object
% --degree:  new degree
%
% OUTPUS
% - Cop:    'nopvar/ndopvar' copy of Aop with new monomial basis
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - change_degree
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
% AT, 01/21/2026: Initial coding



% Check that the operators can indeed be added
if ~isa(Aop, 'nopvar') && ~isa(Aop, 'ndopvar')
    error('The first argument must be "nopvar" or "ndopvar"')
end

if isscalar(degree)
    degree = degree*ones(size(Aop.deg));
end

dec_var = [];
if isa(Aop, 'ndopvar')
    dec_var = Aop.dvarname;
end
% Exclude more complicated addition operations
% --> will need to be included later!
if Aop.deg==degree
    Cop = Aop;
    return
    % error("Addition of operators with different monomial degrees is currently not supported.")
end
if Aop.deg > degree
    error('Decreasing degree is currently unsupported')
end



% Finally, let's actually concatenate
if isa(Aop, 'nopvar')
    Cop = nopvar(); % empty nopvar
elseif isa(Aop, 'ndopvar')
    Cop = ndopvar();
    Cop.dvarname = Aop.dvarname;
else
    error('The first argument must be "nopvar" or "ndopvar"')
end
Cop.dom = Aop.dom;
Cop.deg = degree;
Cop.vars = Aop.vars; 
N = size(Aop.dom,1);
Cop.C = cell([3*ones(1,N),1]);
binStr = dec2base(0:(3^N-1), 3);
% construct degmat for common basis of monomials
for dim = 1:N
    Zd_old_out{dim} = (0:Aop.deg(dim))';
    Zd_old_in{dim}  = (0:Aop.deg(dim))';
    Zd_new_out{dim} = (0:degree(dim))';
    Zd_new_in{dim}  = (0:degree(dim))';
end

% % % consturct monomial basis for nopvar with new degree
degmat_old_out = Zd_old_out{1};
degmat_new_out = Zd_new_out{1};

for dim=2:N
    degmat_old_out = [kron(degmat_old_out, ones(length(Zd_old_out{dim}), 1)), kron(ones(length(degmat_old_out), 1), Zd_old_out{dim})];
    degmat_new_out = [kron(degmat_new_out, ones(length(Zd_new_out{dim}), 1)), kron(ones(length(degmat_new_out), 1), Zd_new_out{dim})];
end

% find location of new monomials in the old monomials
size_left = length(degmat_new_out);
[~, Loc_temp] =ismember(degmat_old_out,  degmat_new_out, 'rows');
for iter_idx = 1:Aop.dim(1)
    % shift m times if (I_m otimes Z(s))^T
    if iter_idx == 1
        Loc_l = Loc_temp;
    else
        Loc_l = [Loc_l; Loc_temp+size_left*(iter_idx - 1)];
    end
end
% if ndopvar Loc_l shows rows of (I otimes xi)^T
% To choose rows of C we need to 
if isa(Aop, 'ndopvar')
    Loc_L_C = (Loc_l -1)*(1+ length(dec_var)) + 1; % starting location of C rows
    Loc_L_C = [kron(Loc_L_C, ones((1+ length(dec_var)), 1)), kron(ones(length(Loc_L_C), 1), (0:length(dec_var))')];
    Loc_L_C = sum(Loc_L_C, 2);
else
    Loc_L_C = Loc_l;
end

for iter_idx = 1:size(binStr, 1)

    substr = binStr(iter_idx, :);
    Qindx = str2num(substr')'+1; % array of .C indeces
    Qindx_cell = num2cell(Qindx); % cell array of .C indeces

    index_for_monom_in_theta = Qindx == 1; % include monomials for theta

    subdegree1 = degree;
    subdegree1(index_for_monom_in_theta) = 0;
    subdegree2 = Aop.deg;
    subdegree2(index_for_monom_in_theta) = 0;

    % construct degmat for Z(theta)
    degmat_new_in = 0;
    degmat_old_in = 0;
    if subdegree1(1) ~= 0
        degmat_new_in = Zd_new_in{1};
    end
    if subdegree2(1) ~= 0
        degmat_old_in = Zd_old_in{1};
    end

    for dim=2:N
        if subdegree1(dim) == 0 % do not perform kron product if no monomials
            continue
        end
        degmat_old_in = [kron(degmat_old_in, ones(length(Zd_old_in{dim}), 1)),kron(ones(length(degmat_old_in), 1), Zd_old_in{dim})];
        degmat_new_in = [kron(degmat_new_in, ones(length(Zd_new_in{dim}), 1)),kron(ones(length(degmat_new_in), 1), Zd_new_in{dim})];
    end

    % construct monomial basis for right multiplier (depends on iter_idx)

    [~, Loc_temp] =ismember(degmat_old_in,  degmat_new_in, 'rows');
    size_right = length(degmat_new_in);
    for iter_idx = 1:Aop.dim(2)  % shift n times if I_n otimes Z(theta)
        if iter_idx == 1
            Loc_r = Loc_temp;
        else
            Loc_r = [Loc_r; Loc_temp+size_right*(iter_idx - 1)];
        end
    end
    % change Cop values
    full_matrix = sparse(size_left*Aop.dim(1)*(length(dec_var) + 1), size_right*Aop.dim(2));
    full_matrix(Loc_L_C, Loc_r) = Aop.C{Qindx_cell{:}};
    Cop.C{Qindx_cell{:}} = full_matrix;
end


end
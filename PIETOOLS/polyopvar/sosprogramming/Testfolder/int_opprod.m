function [Pvec,idx_mat,var2] = int_opprod(Pops)
% [PVEC,IDX_MAT,VAR2] = INT_OPPROD(POPS) takes a cell of 'nopvar' objects
% POPS representing 2-PI operators POPS{i} and computes a cell PVEC of
% kernels such that
%   int_{a}^{b} (POPS{1}*x1)(s).*...*(POPS{d}*xd)(s) ds
%    = sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           Pvec{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1 
% where kj = IDX_MAT(i,j) for j=1,...,d and i=1,...,m for m=d!.
%
% INPUTS
% - Pops:   1 x d cell of 'nopvar' objects, representing a tensor product
%           of 2-PI operators, i.e. (POPS{1}*x1)(s)*...*(POPS{d}*xd)(s);
%
% OUTPUTS
% - Pvec:   d! x 1 cell of 'quadpoly' objects, with the ith element
%           representing the kernel in the ith term in the integral
%           operator acting on x1 o ... o xd (for o the tensor product);
% - idx_mat:    d! x d array, with each row specifying the order of the
%               variables in the associated term in the integral operator.
%               In particular, if idx_mat(i,:) = [k1,k2,...,kd], then term
%               i is defined by the integral
%                   int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b}
%               so that a <= t_{k1} <= t_{k2} <= ... <= t_{kd} <= b;
% - var2:   d x 1 'polynomial' array of 'pvar' objects, specifying the
%           names of the variables t_{j} used in the definition of the
%           kernels in Pvec;
%
% NOTES
% The products of the operators are assumed to be elementwise, so that
%   (Pops{1}*x1)(s).*(Pops{2}*x2)(s).*...*(Pops{d}*xd)(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int_opprod
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
% DJ, 01/21/2026: Initial coding


if isa(Pops,'nopvar')
    Pops = {Pops};
end
dtot = numel(Pops);

% Convert the coefficients defining each parameter to a 'nopvar' object
Pop_params = cell(size(Pops));
var2 = polynomial(zeros(dtot,1));
has_multiplier = zeros(dtot,1);
for kk=1:dtot
    % Extract the kth operator
    Pop_kk = Pops{kk};
    if ~isa(Pop_kk,'nopvar')
        error("Operators must be specified as a cell of 'nopvar' objects.")
    end
    N = size(Pop_kk.vars,1);
    if N>1
        error("Multivariate operators are currently not supported.")
    end
    % Extract the dimension of the operator
    dim_kk = Pop_kk.dim;
    if kk==1
        % Extract primary and dummy variable
        var1 = Pop_kk.vars(1);
        var2_name = Pop_kk.vars(2).varname{1};
        % Extract the domain
        dom = Pop_kk.dom;
    end
    % Check that the domains match
    if any(Pop_kk.dom~=dom)
        error("Spatial domains of the operators must match.")
    end
    % Set the kth dummy variable
    var2_kk_name = [var2_name,'_',num2str(kk)];
    var2_kk = polynomial({var2_kk_name});
    var2(kk) = var2_kk;
    % % Build monomial vectors in primary and dummy variables
    deg_kk = Pop_kk.deg;
    % Z1 = var1.^(0:Pop_kk.deg(1))';
    % Z2 = var2_kk.^(0:Pop_kk.deg(1))';
    % Check that we have no multiplier term
    Pop_kk_params = Pop_kk.C;
    if ~isempty(Pop_kk.C{1}) && any(any(Pop_kk.C{1}))
        has_multiplier(kk) = 1;
        Pop_kk_params{1} = coeff2poly(Pop_kk.C{1},dim_kk,[deg_kk,0],[var2_kk,var1]);  % <-- substitute R0(s=t_k)
    else
        Pop_kk_params{1} = zeros(dim_kk);
    end
    % Convert integral terms to polynomial
    for ii=2:3
        %Pop_kk_params.C{ii} = quadPoly(Pop_kk.C{ii},Z1,Z2,dim_kk,1,1);
        Pop_kk_params{ii} = coeff2poly(Pop_kk.C{ii},dim_kk,deg_kk,[var1,var2_kk]);
    end
    Pop_params{kk} = Pop_kk_params;
end
if sum(has_multiplier)>1
    error("At most one of the operators may include a multiplier component.")
else
    m_num = find(has_multiplier);
    var_m = var2(m_num);
end


% List all possible orders of the variables t1 through td
%   idx_mat(l,:) = [i,j,k] means a <= ti <= tj <= tk <= b in term l
idx_mat = 1;
for ii=2:dtot
    n_ords = size(idx_mat,1);
    idx_mat_new = zeros(ii*n_ords,ii);
    for jj=1:ii
        % Place variable t_ii in position jj
        idx_mat_new(jj:ii:end,:) = [idx_mat(:,1:jj-1),ii*ones(n_ords,1),idx_mat(:,jj:end)];
    end
    idx_mat = idx_mat_new;
end

% For each possible ordering of the variables ti, compute the associated
% kernel
n_ords = size(idx_mat,1);
% ord_mat = zeros(n_ords,dtot);
Pvec = cell(n_ords,1);
for ii=1:n_ords
    Pvec_ii = 0;
    idx_ii = idx_mat(ii,:);
    for jj=1:dtot+1
        % For the considered order of the variables ti, place variable s in
        % position jj
        idx_jj = [idx_ii(1,1:jj-1),dtot+1,idx_ii(1,jj:end)];
        % For each variable ti, determine whether s<= ti or s>=ti
        [~,ord] = sort(idx_jj);
        is_geq_s = ord(1:dtot)>=ord(end);
        % Extract the associated parameter from each operator Pop_kk
        Rjj = 1;
        for kk=1:dtot
            Rjj = Rjj.*Pop_params{kk}{is_geq_s(kk)+2};
        end
        % Integrate the product of the parameters over the interval
        if jj==1
            L = dom(1);
        else
            L = var2(idx_ii(1,jj-1));
        end
        if jj==dtot+1
            U = dom(2);
        else
            U = var2(idx_ii(1,jj));
        end
        Pvec_ii = Pvec_ii + int(Rjj,var1,L,U);
    end
    % Also account for possible multiplier term: tj = s
    if any(has_multiplier)
        R_tmp = 1;
        m_idx = find(idx_ii==m_num,1,'first');
        for kk=idx_ii(1:m_idx-1)
            % Loop over all variables k for which t_k <= t_m = s
            R_tmp = R_tmp.*subs(Pop_params{kk}{2},var1,var_m);
        end
        % Deal with the operator which has a multiplier
        R_tmp = R_tmp.*Pop_params{m_num}{1};
        for kk=idx_ii(m_idx+1:end)
            % Loop over all variables k for which t_k >= t_m = s
            R_tmp = R_tmp.*subs(Pop_params{kk}{3},var1,var_m);
        end
        Pvec_ii = Pvec_ii + R_tmp;
    end
    Pvec{ii} = Pvec_ii;
end

end
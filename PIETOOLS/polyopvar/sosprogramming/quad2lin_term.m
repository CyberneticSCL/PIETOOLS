function [Kcell,idx_mat,var2] = quad2lin_term(Lop,Pmat,Rop)
% [KCELL,IDX_MAT,VAR2] = QUAD2LIN_TERM(ZOPL,PMAT,ZOPR) takes a distributed
% polynomial in quadratic form
%   p(x) = <Lop*Z1(x) , Pmat*Rop*Z2(x)>_{L2}
% and converts it to the linear form, computing a call 
% takes a cell of 'nopvar' objects
% POPS representing 2-PI operators POPS{i} and computes a cell KCELL of
% kernels such that
%  p(x) = sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           Kcell{i}(t1,...,td)*Z1(x)(t_{1},...,t_{d1})*Z2(x)(t_{d1+1},...,t_{d1+d2}) dt_kd ... dt_k1 
% where kj = IDX_MAT(i,j) for j=1,...,d and i=1,...,m for m=d!, with 
% d=d1+d2.
%
% INPUTS
% - Lop:    m x 1 'tensopvar' object representing an operator acting on a
%           scalar distributed monomial Z1(x) of degree d1;
% - Pmat:   m x n 'double' or 'dpvar' object parameterizing the inner
%           product;
% - Rop:    n x 1 'tensopvar' object representing an operator acting on a
%           scalar distributed monomial Z2(x) of degree d2;
%
% OUTPUTS
% - Kcell:  d! x 1 cell of 'quadpoly' objects, with the ith element
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
%   (Lop{1}*x1)(s).*(Lop{2}*x2)(s).*...*(Lop{d1}*xd1)(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - quad2lin_term
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

if ~isa(Lop,'tensopvar') || ~isa(Rop,'tensopvar')
    error("Operators parameterizing the inner product should be specified as 'tensopvar' objects.")
end
if isa(Lop.ops{1},'nopvar')
    d1 = 1;
    Lops = {Lop.ops{1}};
else
    Lops = Lop.ops{1};
    d1 = numel(Lops);
end
if isa(Rop.ops{1},'nopvar')
    d2 = 1;
    Rops = {Rop.ops{1}};
else
    Rops = Rop.ops{1};
    d2 = numel(Rops);
end
dtot = d1+d2;

% Convert the coefficients defining each parameter to a 'nopvar' object
Pop_params = cell(1,d1+d2);
var2 = polynomial(zeros(dtot,1));
has_multiplier = zeros(dtot,1);
for kk=1:dtot
    % Extract the kth operator
    if kk<=d1
        Pop_kk = Lops{kk};
    else
        Pop_kk = Rops{kk-d1};
    end
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
        if kk<=d1
            Pop_kk_params{1} = Pop_kk_params{1}';
        end
    else
        Pop_kk_params{1} = zeros(dim_kk);
    end
    % Convert integral terms to polynomial
    for ii=2:3
        %Pop_kk_params.C{ii} = quadPoly(Pop_kk.C{ii},Z1,Z2,dim_kk,1,1);
        Pop_kk_params{ii} = coeff2poly(Pop_kk.C{ii},dim_kk,deg_kk,[var1,var2_kk]);
        if kk<=d1
            Pop_kk_params{ii} = Pop_kk_params{ii}';     % Lop operators get transposed in inner product
        end
    end
    Pop_params{kk} = Pop_kk_params;
end
if sum(has_multiplier)>1 && (d1~=1 || d2~=1)
    error("At most one of the operators may include a multiplier component.")
else
    m_num = find(has_multiplier,1,'first');
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
Kcell = cell(n_ords,1);
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
        Ljj = 1;
        for kk=1:d1
            Ljj = Ljj.*Pop_params{kk}{is_geq_s(kk)+2};
        end
        Rjj = 1;
        for kk=d1+(1:d2)
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
        Pvec_ii = Pvec_ii + int(Ljj*Pmat*Rjj,var1,L,U);
    end
    % Also account for possible multiplier term: tj = s
    if any(has_multiplier)
        [~,ord_ii] = sort(idx_ii);
        ord_m = ord_ii(m_num);
        L_tmp = 1;
        for kk=1:d1
            if ord_ii(kk)<ord_m
                % t_k <= t_m = s
                L_tmp = L_tmp.*subs(Pop_params{kk}{2},var1,var_m);
            elseif ord_ii(kk)==ord_m
                % t_k == t_m = s
                L_tmp = L_tmp.*Pop_params{m_num}{1};
            else
                % t_k >= t_m = s
                L_tmp = L_tmp.*subs(Pop_params{kk}{3},var1,var_m);
            end
        end
        R_tmp = 1;
        for kk=d1+(1:d2)
            if ord_ii(kk)<ord_m
                % t_k <= t_m = s
                R_tmp = R_tmp.*subs(Pop_params{kk}{2},var1,var_m);
            elseif ord_ii(kk)==ord_m
                % t_k == t_m = s
                R_tmp = R_tmp.*Pop_params{m_num}{1};
            else
                % t_k >= t_m = s
                R_tmp = R_tmp.*subs(Pop_params{kk}{3},var1,var_m);
            end
        end
        Pvec_ii = Pvec_ii + L_tmp*Pmat*R_tmp;

        % m_idx = find(idx_ii==m_num,1,'first');
        % for kk=idx_ii(1:m_idx-1)
        %     % Loop over all variables k for which t_k <= t_m = s
        %     R_tmp = R_tmp.*subs(Pop_params{kk}{2},var1,var_m);
        % end
        % % Deal with the operator which has a multiplier
        % R_tmp = R_tmp.*Pop_params{m_num}{1};
        % for kk=idx_ii(m_idx+1:end)
        %     % Loop over all variables k for which t_k >= t_m = s
        %     R_tmp = R_tmp.*subs(Pop_params{kk}{3},var1,var_m);
        % end
        % Pvec_ii = Pvec_ii + R_tmp;
    end
    Kcell{ii} = Pvec_ii;
end
% In the case of a quadratic function
%   < Zop*x , P*Zop*x>
% we do allow both terms Zop to have multipliers, 
%   < Zop_0*x, P*Zop_0*x> = int_{a}^{b} Z0(s)^T*Pmat*Z0(s)*x(s)*x(s) ds
% We store the parameter Z0(s)^T*Pmat*Z0(s) as an extra element of the
% cell.
if d1==1 && d2==1 && sum(has_multiplier)==2
    Kcell = [Kcell; {Pop_params{1}{1}*Pmat*subs(Pop_params{2}{1},var2(2),var2(1))}];
end

end
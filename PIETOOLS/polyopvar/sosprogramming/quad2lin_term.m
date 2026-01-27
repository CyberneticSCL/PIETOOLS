function [Kop] = quad2lin_term(Pmat,Lmon,Rmon,dom,var1,var2)
% [KOP] = QUAD2LIN_TERM(PMAT,LMON,RMON) takes a distributed
% polynomial in quadratic form
%   V(x) = <LMON(x) , PMAT*RMON(x))>_{L2},
% for 
%   LMON(x) = C1*(x1^d1*...*xm^dm)(s)  and  RMON = C2*(x1^p1*...*xm^pm)
% and converts it to the linear form, computing a cell KCELL of
% kernels such that
%  V(x) = sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           K_{k}(t1,...,td)*
%           x1(t1)*...*x1(t_{d1+p1})*...*xm(t_{d-dm-pm+1})*...*xm(t_{d})
%                                                       dt_kd ... dt_k1 
% where kj = IDX_MAT(i,j) for j=1,...,d and i=1,...,m for m=d!, with 
% d=d1...+dm+p1+...+pm.
%
% INPUTS
% - Pmat:   m1 x n1 'double' or 'dpvar' object parameterizing the inner
%           product;
% - Lmon:   m1 x m 'polyopvar' object representing an distributed monomial
%           involving a single (scalar) distributed monomial Z1(x);
% - Rmon:   n1 x n 'polyopvar' object representing an distributed monomial
%           involving a single (scalar) distributed monomial Z1(x);
% - dom:    (optional) 1x2 array specifing the interval, [a,b], over which
%           integration is performed.
% - var1:   (optional) 1x1 'pvar' object specifying the variable used in
%           the integral defining the inner product;
% - var2:   (optional) dx1 'pvar' object specifying the variables used in
%           the multi-integral defining the inner product in linear format;
%
% OUTPUTS
% - Kop:    struct with fields
%               params: m x q*n 'polynomial' or 'dpvar' object,
%                       with elements (1:m,(i-1)*n+1:i*n) corresponding to 
%                       the kernels in the ith term of the functional
%                       in linear format
%               omat:   q x p array of integers, where p is the number of
%                       dummy variables for integration, with row i
%                       specifying the order of the variables in the
%                       integral associated with the ith kernel.
%                       Specifically, if omat(i,:) = [k1,k2,...,kd], 
%                       then term i is defined by the integral
%                           int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b};
%               matdim: 1 x 2 array specifying the dimensions of the
%                       kernels K{i};
%               vars:   p x 1 'polynomial' array specifying the names of
%                       the dummy variables used in definition of the
%                       kernels in Kop.params;
%               dom:    interval [a,b] over which to integrate;
%

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

% Process the left and right factors
if isa(Lmon,'polyopvar') && isa(Rmon,'polyopvar')
    % Make sure the two factors are expressed in terms of the same variables
    [Lmon,Rmon] = common_vars(Lmon,Rmon);
    degL = Lmon.degmat;
    degR = Rmon.degmat;
elseif (isa(Lmon,'polynomial') || isa(Lmon,'dpvar')) || isempty(Lmon) && isa(Rmon,'polyopvar')
    degL = 0;
    degR = Rmon.degmat;
else
    error("Monomials terms in the inner product should be specified as 'polyopvar' objects.")
end
if size(degL,1)~=1 || size(degR,1)~=1
    error("The left and right polynomial in the inner product must be defined by exactly one monomial.")
end

% Check that the other arguments are properly specified
if isempty(Pmat)
    Pmat = 1;
end
if nargin<=3
    dom = [];
end
if nargin<=4
    var1 = [];
end



% Add the monomial degrees of the left and right monomials:
%   [x1^d1*...*xm^dm] * [x1^p1*...*xm^pm] = x1^(d1+p1)*...*xm^(dm+pm)
deg_new = degL+degR;
d1 = sum(degL);
d2 = sum(degR);
dtot = d1+d2;
nvars = size(deg_new,2);
if nargin<=5
    var2 = polynomial(zeros(dtot,1));
elseif numel(var2)~=dtot
    error("Number of dummy variables must match the cumulative degree of the monomials.")
end
if any(d1>0)
% % Establish the new order of the independent variables:
% If we have e.g. the product of a monomial in 5 factors with a monomial in 
% 4 factors,
%   [x1(t1)*x1(t2)*x2(t3)*x2(t4)*x2(t5)] * [x1(t6)*x1(t7)*x2(t8)*x2(t9)]
% we need to re-index the independent variables (t1,...,t9) to get
%    x1(t1)*x1(t2)*x1(t3)*x1(t4) * x2(t5)*x2(t6)*x2(t7)*x2(t8)*x2(t9);
% First, assign each factor k in Lmon(x) an index j=p1_idcs(k)=k, meaning 
% that factor k depends on variable t_{j}=t_{k}
p1_idcs = [[0,cumsum(degL(1:end-1))];cumsum(degL)];
p1_idcs = mat2cell(p1_idcs,2,ones(1,nvars));
p1_idcs = cellfun(@(a)(a(1)+1:a(2))',p1_idcs,'UniformOutput',false);
% Similarly, assign each factor k in Rmon(x) an index j=p2_idcs(k)=k+d1,
% meaning that factor k depends on variable t_{j}=t_{k+d1}
p2_idcs = d1+[[0,cumsum(degR(1:end-1))];cumsum(degR)];
p2_idcs = mat2cell(p2_idcs,2,ones(1,nvars));
p2_idcs = cellfun(@(a)(a(1)+1:a(2))',p2_idcs,'UniformOutput',false);
% Assign each factor in the product Lmon(x)*Rmon(x) an index
p_idcs = [p1_idcs;p2_idcs];
quad_var_order = cell2mat(p_idcs(:))';     % quad_var_order(i) = k means that factor i in the linear representation corresponds to factor k in the quadratic one
% Reorder so that factor k in the product again depends on variable t_{k}
[~,var_order] = sort(quad_var_order);      % var_order(k) = i means that factor k in the old monomials is assigned variable t_{i}
else
    var_order = 1:dtot;
end


% Extract the operators acting on the left and right monomials
if isa(Lmon,'polyopvar')
    if isa(Lmon.C.ops{1},'nopvar')
        Lops = Lmon.C.ops;
    else
        Lops = Lmon.C.ops{1};
    end
    mdim = Lops{1}.dim(2);
else
    if isempty(Lmon)
        Lfctr = 1;
    else
        Lfctr = Lmon;
    end
    Lops = {};
    mdim = size(Lfctr,1);
end
if isa(Rmon.C.ops{1},'nopvar')
    Rops = Rmon.C.ops;
else
    Rops = Rmon.C.ops{1};
end
ndim = Rops{1}.dim(2);
if size(Lops,1)>1 || size(Rops,1)>1
    error("Conversion of polynomials with multiple terms is not supported.")
end
if numel(Lops)~=d1 || numel(Rops)~=d2
    error("Number of factors in coefficient operator should match the power of the monomial.")
end

% Convert the coefficients defining each parameter to a 'nopvar' object
Pop_params = cell(1,dtot);
has_multiplier = zeros(1,dtot);
for kk=1:numel(Pop_params)
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
        if isempty(var1)
            % Set the variable of integration
            var1 = Pop_kk.vars(1);
        end
        if isempty(dom)
            % Extract the domain
            dom = Pop_kk.dom;
        end        
    end
    % % Check that the domains match
    if any(~isequal(Pop_kk.dom,dom))
        error("Spatial domains of the operators must match.")
    end
    % Set the kth dummy variable
    if isequal(var2(var_order(kk)),0)
        var2_kk_name = [var1.varname{1},'_',num2str(var_order(kk))];
        var2_kk = polynomial({var2_kk_name});
        var2(var_order(kk)) = var2_kk;
    else
        var2_kk = var2(var_order(kk));
    end
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
    m_nums = find(has_multiplier);
    vars_m = var2(var_order(m_nums));
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
if isa(Pmat,'dpvar') || isa(Lmon,'dpvar')
    Kparams = dpvar(zeros(mdim,ndim*n_ords));
else
    Kparams = polynomial(zeros(mdim,ndim*n_ords));
end
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
            Ljj = Ljj.*Pop_params{kk}{is_geq_s(var_order(kk))+2};
        end
        Rjj = 1;
        for kk=d1+(1:d2)
            Rjj = Rjj.*Pop_params{1,kk}{is_geq_s(var_order(kk))+2};
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
        if d1==0
            Pvec_ii = Pvec_ii + int(Lfctr*Pmat*Rjj,var1,L,U);
        else
            Pvec_ii = Pvec_ii + int(Ljj*Pmat*Rjj,var1,L,U);
        end
    end
    % Also account for possible multiplier term: tj = s
    for jj=1:numel(m_nums)
        [~,ord_ii] = sort(idx_ii);
        ord_m = ord_ii(var_order(m_nums(jj)));
        L_tmp = 1;
        for kk=1:d1
            if ord_ii(var_order(kk))<ord_m
                % t_k <= t_m = s
                L_tmp = L_tmp.*subs(Pop_params{kk}{2},var1,vars_m(jj));
            elseif ord_ii(var_order(kk))==ord_m
                % t_k == t_m = s
                L_tmp = L_tmp.*Pop_params{m_nums(jj)}{1};
            else
                % t_k >= t_m = s
                L_tmp = L_tmp.*subs(Pop_params{kk}{3},var1,vars_m(jj));
            end
        end
        if d1==0
            L_tmp = subs(Lfctr,var1,vars_m(jj));
        end
        R_tmp = 1;
        for kk=d1+(1:d2)
            if ord_ii(var_order(kk))<ord_m
                % t_k <= t_m = s
                R_tmp = R_tmp.*subs(Pop_params{1,kk}{2},var1,vars_m(jj));
            elseif ord_ii(var_order(kk))==ord_m
                % t_k == t_m = s
                R_tmp = R_tmp.*Pop_params{1,m_nums(jj)}{1};
            else
                % t_k >= t_m = s
                R_tmp = R_tmp.*subs(Pop_params{1,kk}{3},var1,vars_m(jj));
            end
        end
        Pvec_ii = Pvec_ii + subs(L_tmp*Pmat*R_tmp,var1,vars_m(jj));
    end
    if ~isempty(Pvec_ii)
        Kparams(:,(ii-1)*ndim+1:ii*ndim) = Pvec_ii;
    end
end
% In the case of a quadratic function
%   < Zop*x , P*Zop*x>
% we do allow both terms Zop to have multipliers, 
%   < Zop_0*x, P*Zop_0*x> = int_{a}^{b} Z0(s)^T*Pmat*Z0(s)*x(s)*x(s) ds
% We store the parameter Z0(s)^T*Pmat*Z0(s) as an extra element of the
% cell.
if d1==1 && d2==1 && sum(has_multiplier)==2
    Kparams = [Kparams, Pop_params{1}{1}*Pmat*subs(Pop_params{2}{1},var2(2),var2(1))];
end

Kop = struct();
Kop.params = Kparams;
Kop.omat = idx_mat;
Kop.vars = var2;
Kop.dom = dom;
Kop.matdim = [mdim,ndim];

end
function [KCfun,omat,var2,dom] = int(Cop,Kfun,var1,dom,var2)
% [KCFUN,OMAT,VAR2,DOM] = int(COP,KFUN,DOM,VAR1,VAR2) takes a tensor-PI
% operator defined by Cops,
%   (Pop*x)(s) := sum_{i=1}^{m} (COP.ops{i,1}*x1)(s)*...*(COP.ops{i,d}*xd)(s)
% where COPS{i,j} are all 1D PI operators, and computes the integral
%   V(x) = int_{DOM(1)}^{DOM(2)} KFUN(s)*(Pop*x)(s) ds,
% returning the kernel function KCFUN such that
%  V(x) = sum_{i=1}^{q} int_{a}^{b} int_{t_k1(i)}^{b} ... int_{t_kd(i)}^{b} 
%          KC_{i}(t1,...,td)*x1(t1)*...*xd(t_{d}) dt_{kd(i)} ... dt_{k1(i)} 
% where kj(i) = OMAT(i,j) for j=1,...,d, DOM=[a,b], and VAR2=(t1,...,td);
%
% INPUTS
% - Cop:    p x n 'tensopvar' object representing a tensor-PI operator
%           acting on a single degree-d distributed monomial;
% - Kfun:   m x p array of type 'double', 'polynomial', or 'dpvar',
%           representing the kernel function in the integral;
% - var1:   (optional) 1 x 1 cellstr object specifying the name of the 
%           variable (s) used in the integral defining the inner product. 
%           Defaults to Cops{1}.var1.varname;
% - dom:    (optional) 1 x 2 array specifing the interval, [a,b], over
%           which integration is performed. The limits can be variable, 
%           dom = [r1,r2], in which case we assume r1,r2 in [a,b], where
%           [a,b] = Cops{1}.dom. This is also the default domain;
% - var2:   (optional) 1 x d cellstr specifying the desired names of the
%           dummy variables (t1,...,td) used in the final multi-integral;
%
% OUTPUTS
% - KCfun:  m x n*q array of the same type as 'Kfun', specifying the
%           kernels defining the integral of the specified tensor-PI
%           operator, where q is the number of terms in the integral. For
%           each i in {1,...,q}, elements KCfun(1:m,(i-1)*n+1:i*n)
%           correspond to the kernel in the ith integral term;
% - omat:   q x d array of integers, specifying for each of the q terms in
%           the integral the order of the dummy variables in this integral.
%           Specifically, if omat(i,:) = [k1,...,kd], then the dummy
%           variables are ordered as t_{k1} <= t_{k2} <= ... <= t_{kd}, so
%           that the limits of the integral will be
%               int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b};
% - var2:   1 x p cellstr specifying the actual names of the dummy
%           variables (r1,t1,..,td,r2) used in the final multi-integral. 
%           Here, r1=dom(1) if this lower boundary is variable, and 
%           r2=dom(2) if this upper boundary is variable. Otherwise, p=d,
%           and we have var2 = (t1,...,td);
% - dom:    1 x 2 array specifing the interval, [a,b], over the output
%           variables var2 should be integrated. If the input domain
%           involves any variables, the output domain will be different;
%

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
% DJ, 02/03/2026: Initial coding

% Allow for just integration with kernel 1:
%   int(Cop,var1,dom)
% if nargin>=2 && (isa(Kfun,'polynomial') || iscellstr(Kfun))
%     if nargin==4
%         var2 = dom;
%     end
%     if nargin>=3
%         dom = var1;
%     end
%     if nargin>=2
%         var1 = Kfun;
%     end
%     Kfun = 1;
% end

% Extract the PI operators defining the tensor-PI operator
if ~isa(Cop,'tensopvar')
    error("Tensor-PI operator to integrate must be specified as 'tensopvar' object.")
elseif numel(Cop.ops)>1
    error("Tensor-PI operator must act on a single monomial.")
end
Cops = Cop.ops{1};
if isa(Cops,'nopvar')
    Cops = {Cops};
end

% Check that the dimensions of the kernel and tensor-PI operator match
if size(Kfun,2)~=size(Cops{1},1) && ~all(size(Kfun)==1)
    error("Column dimension of the kernel must match row dimension of the tensor-PI operator.")
end
mdim = size(Kfun,1);
ndim = size(Cops{1},2);
d = size(Cops,2);
ntrms = size(Cops,1);

% Check that the other arguments are properly specified
if isempty(Kfun)
    Kfun = 1;
end
if nargin<=2
    % If not specified, use the spatial variable of the tensor-PI operator
    % as dummy variable for integration.
    pvar1 = Cops{1}.var1;
elseif isa(var1,'polynomial')
    pvar1 = var1;
elseif iscellstr(var1)
    pvar1 = polynomial(var1);
end
if numel(pvar1)>1
    error("Only a single dummy variable for integration can be specified.")
end
if nargin<=3
    % If not specified, assume we integrate over the entire domain on which
    % the tensor-PI operator is defined.
    dom = Cops{1}.dom;
end
if nargin<=4
    % If no variables are specifeid
    var2 = polynomial(zeros(d,1));
elseif isa(var2,'polynomial') && ispvar(var2)
    var2name = cell(numel(var2),1);
    for i=1:numel(var2name)
        var2name(i) = var2(i).varname;
    end
    var2 = unique(var2name);
elseif ~iscellstr(var2)
    error("Variables must be specified as 'cellstr' or array of 'pvar' objects")
end
if numel(var2)~=d
    error("Number of dummy variables must match the cumulative degree of the monomials.")
end


% Extract the parameters defining each of the PI operators in the tensor-PI
% operator
[Cparams,has_multiplier,var2,opdom] = compute_params(Cops,pvar1,var2);
pvars2 = polynomial(var2)';
% Allow only multipliers acting on linear monomials
if sum(has_multiplier)>1 %|| (sum(has_multiplier)==1 && d>1)
    error("Multiplier terms in tensor producs of PI operators are not supported.")
else
    m_nums = find(has_multiplier);
    vars_m = pvars2(m_nums);
end

% If the domain of integration is variable, dom=[r1,r2], add the limits r1
% and/or r2 to the output list of variables
if ~isequal(dom(1),opdom(1))
    if ~isa(dom(1),'polynomial') || ~ispvar(dom(1))
        error("For fixed domains of integration, domain must match that of the operator.")
    end
    is_vardom1 = true;
    var2 = [dom(1).varname, var2];
else
    is_vardom1 = false;
end
if ~isequal(dom(2),opdom(2))
    if ~isa(dom(2),'polynomial') || ~ispvar(dom(2))
        error("For fixed domains of integration, domain must match that of the operator.")
    end
    is_vardom2 = true;
    var2 = [var2, dom(2).varname];
else
    is_vardom2 = false;
end
% Add the (variable) boundaries of the domain to the list of variables.
pvars2 = [dom(1); pvars2; dom(2)];
df = d+2;
dom = opdom;


% List all possible orders of the variables (r1,t1,...,td,r2)
%   idx_mat(l,:) = [i,j,k] means a <= ti <= tj <= tk <= b in term l
omat = [ones(1,is_vardom1),(d+1+is_vardom1)*ones(1,is_vardom2)];   % r1<=r2 must always hold
for i=1:d
    n_ords = size(omat,1);
    n_pos = size(omat,2)+1;
    omat_new = zeros(n_pos*n_ords,n_pos);
    for j=1:n_pos
        % Place variable t_i in position j
        omat_new(j:n_pos:end,:) = [omat(:,1:j-1),(i+is_vardom1)*ones(n_ords,1),omat(:,j:end)];
    end
    omat = omat_new;
end
% In case of fixed domain limits, also add those to the list of variables
omat_f = [ones(size(omat,1),1-is_vardom1), omat+(1-is_vardom1), df*ones(size(omat,1),1-is_vardom2)];


% For each possible ordering of the variables (r1,t1,...,td,r2), compute 
% the kernel in the associated term in the functional
n_ords = size(omat_f,1);
if isa(Kfun,'dpvar')
    KCfun = dpvar(zeros(mdim,ndim*n_ords));
else
    KCfun = polynomial(zeros(mdim,ndim*n_ords));
end
for i=1:n_ords
    param_i = 0;
    % Establish the order of the variables in the integral in term i
    ord_f = omat_f(i,:);
    % Check which variables ti are actually in the domain dom = [r1,r2]
    r1_pos = find(ord_f==1,1,'first');
    r2_pos = find(ord_f==df,1,'last');
    rtn_idcs = setdiff(1:df,[r1_pos,r2_pos]);
    ord_i = ord_f(rtn_idcs)-1;
    % For all variables ti<=r1, we have ti <= r1 <= s
    % --> extract the parameter associated with int_{a}^{s} dti
    Cparam1 = num2cell(ones(ntrms,1));
    %Cparam1 = ones(ntrms_R*kdim,mdim);
    for k=ord_f(1:r1_pos-1)
        for trm=1:ntrms
            % ridcs = (trm-1)*kdim+(1:kdim);
            % Cparam1(ridcs,:) = Cparam1(ridcs,:).*Cparams{trm,k-1}{2};
            Cparam1{trm} = Cparam1{trm}.*Cparams{trm,k-1}{2};
        end
    end
    % For all variables ti>=r2, we have ti >= r2 >= s
    % --> extract the parameter associated with int_{s}^{b} dti
    for k=ord_f(r2_pos+1:end)
        for trm=1:ntrms
            % ridcs = (trm-1)*kdim+(1:kdim);
            % Cparam1(ridcs,:) = Cparam1(ridcs,:).*Cparams{trm,k-1}{3};
            Cparam1{trm} = Cparam1{trm}.*Cparams{trm,k-1}{3};
        end
    end
    % For the variables within the domain [r1,r2], the parameter we need
    % depends on where s is placed in [t1,t2], e.g. s>=ti or s<=ti
    ord_j = ord_f(r1_pos:r2_pos);
    for j=2:r2_pos-r1_pos+1
        % Place variable s in position j between r1 and r2
        ord_tmp = [ord_j(1,2:j-1),df+1,ord_j(1,j:end-1)];
        % For each variable ti, determine whether s<= ti or s>=ti
        [ord_sort,pos] = sort(ord_tmp);
        is_geq_s = pos(1:end-1)>=pos(end);
        % Extract the associated parameter from each operator Rop_k
        Cj = 0;
        for trm=1:ntrms
            %ridcs = (trm-1)*kdim+(1:kdim);
            %Ctmp = Cparam1(ridcs,:);
            Ctmp = Cparam1{trm};
            for k=1:numel(ord_sort)-1
                varnum = ord_sort(k)-1;
                Ctmp = Ctmp.*Cparams{trm,varnum}{is_geq_s(k)+2};
            end
            Cj = Cj + Ctmp;
        end
        % Integrate the product of the parameters over the interval
        % s in [L,U]
        L = pvars2(ord_j(j-1));     % will be dom(1) if j==2
        U = pvars2(ord_j(j));       % will be dom(2) if j==numel(indom_idcs)
        param_i = param_i + int_simple(Kfun*Cj,pvar1,L,U);
    end

    % % Also account for possible multiplier term, where t_m = s
    % Since s in [r1,r2], this also requires t_m in [r1,r2]
    is_indom = ismember(m_nums,ord_f(r1_pos+1:r2_pos-1)-1);
    m_nums_rtn = m_nums(is_indom);
    vars_m_rtn = vars_m(is_indom);
    for j=1:numel(m_nums_rtn)
        % Establish in which position variable tm is, ignoring r1 and r2
        pos_m = find(ord_i==m_nums_rtn(j),1,'first');
        % Determine which variables are smaller and greater than t_m = s
        leq_vars = ord_i(1:pos_m-1);
        geq_vars = ord_i(pos_m+1:end);
        Cj = 0;
        for trm=1:ntrms
            % Multiply the appropriate parameters, starting with t_m = s
            Ctmp = Cparams{trm,m_nums_rtn(j)}{1};
            for k=leq_vars
                % t_k <= t_m = s
                Ctmp = Ctmp.*Cparams{trm,k}{2};
            end
            for k=geq_vars
                % t_k >= t_m = s
                Ctmp = Ctmp.*Cparams{trm,k}{3};
            end
            Cj = Cj+Ctmp;
        end
        param_i = param_i + subs(Kfun*Cj,pvar1,vars_m_rtn(j));
    end
    if ~isempty(param_i)
        KCfun(:,(i-1)*ndim+1:i*ndim) = param_i;
    end
end

end




%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function [Pparams,has_multiplier,var2name,dom] = compute_params(Pops,var1,var2name)
% For a m x d cell of 1D nopvar objects Pops, return the parameters 
% defining the associated operators as an m x d cell, where each element is
% a 1x3 cell specifying the parameters for the multiplier term, lower
% integral, and upper integral;

Pparams = cell(size(Pops));
has_multiplier = false(1,size(Pops,2));
for k=1:numel(Pparams)
    % Determine the row and column index in Pop
    [ridx,cidx] = ind2sub(size(Pops),k);
    % Extract the kth operator
    Pop_k = Pops{k};
    if ~isa(Pop_k,'nopvar')
        error("Operators must be specified as a cell of 'nopvar' objects.")
    end
    N = size(Pop_k.vars,1);
    if N>1
        error("Multivariate operators are currently not supported.")
    end
    % Extract the dimension of the operator
    dim_k = Pop_k.dim;
    % Extract the domain of the operator
    if k==1
        dom = Pop_k.dom;
    elseif any(Pop_k.dom~=dom)
        error("Spatial domains on which the PI operators are defined should match.")
    end
    % Set the kth dummy variable
    if ridx==1 && isempty(var2name{cidx})
        var2name{cidx} = [var1.varname{1},'_',num2str(cidx)];        
    end
    pvar2_k = polynomial(var2name(cidx));
    % % Build monomial vectors in primary and dummy variables
    deg_k = Pop_k.deg;
    % Check that we have no multiplier term
    Pop_k_params = Pop_k.C;
    if ~isempty(Pop_k.C{1}) && any(any(Pop_k.C{1}))
        has_multiplier(cidx) = true;
        %Pop_k_params{1} = coeff2poly(Pop_k.C{1},dim_k,[deg_k,0],[pvar2_k,var1]);  % <-- substitute R0(s=t_k)
        Pop_k_params{1} = coeff2poly(Pop_k.C{1},dim_k,[deg_k,0],[var1,pvar2_k]);
    else
        Pop_k_params{1} = zeros(dim_k);
    end
    % Convert integral terms to polynomial
    for i=2:3
        %Pop_kk_params.C{ii} = quadPoly(Pop_kk.C{ii},Z1,Z2,dim_kk,1,1);
        Pop_k_params{i} = coeff2poly(Pop_k.C{i},dim_k,deg_k,[var1,pvar2_k]);
    end
    Pparams{k} = Pop_k_params;
end


end
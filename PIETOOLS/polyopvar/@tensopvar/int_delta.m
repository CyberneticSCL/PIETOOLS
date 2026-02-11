function [KCfun,omat,var2,dom] = int_delta(Cop,Kfun,var1,dom,var2)
% [KCFUN,OMAT,VAR2,DOM] = int_delta(COP,KFUN,DOM,VAR1,VAR2) takes a
% tensor-PI operator defined by Cop,
%   (Pop*x)(s) := sum_{i=1}^{m} (COP.ops{i,1}*x1)(s)*...*(COP.ops{i,d}*xd)(s)
% where COP.ops{i,j} are all 1D PI operators, and computes the integral
%   V(x,y) = int_{DOM(1)}^{DOM(2)} KFUN(s)*y(s)*(Pop*x)(s) ds,
% returning the kernel function KCFUN such that
%  V(x,y) = sum_{i=1}^{q} int_{a}^{b} int_{t_k1(i)}^{b} ... int_{t_kd(i)}^{b} 
%       KC_{i}(s,t1,...,td)*y(s)*x1(t1)*...*xd(t_{d}) dt_{kd(i)} ... dt_{k1(i)} ds 
% where kj(i) = OMAT(i,j) for j=1,...,d+1, DOM=[a,b], and 
% VAR2=(s,t1,...,td);
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
%           which integration is performed. The limits must be fixed;
% - var2:   (optional) 1 x d cellstr specifying the desired names of the
%           dummy variables (t1,...,td) used in the final multi-integral;
%
% OUTPUTS
% - KCfun:  m x n*q array of the same type as 'Kfun', specifying the
%           kernels defining the integral of the specified tensor-PI
%           operator, where q is the number of terms in the integral. For
%           each i in {1,...,q}, elements KCfun(1:m,(i-1)*n+1:i*n)
%           correspond to the kernel in the ith integral term;
% - omat:   q x (d+1) array of integers, specifying for each of the q terms
%           in the integral the order of the dummy variables in this 
%           integral. Specifically, if omat(i,:) = [k1,...,kd,k(d+1)], then 
%           the dummy variables are ordered as 
%               t_{k1} <= t_{k2} <= ... <= t_{kd},
%           so that the limits of the integral will be
%               int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b};
% - var2:   1 x (d+1) cellstr specifying the actual names of the dummy
%           variables (s,t1,..,td) used in the final multi-integral;
% - dom:    1 x 2 array specifing the interval, [a,b], over the output
%           variables var2 should be integrated. Will be the same as the 
%           input dom;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int_delta
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
% DJ, 02/10/2026: Initial coding


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
[~,pdim] = size(Kfun);
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
    var1 = pvar1.varname;
elseif isa(var1,'polynomial')
    pvar1 = var1;
    var1 = pvar1.varname;
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
if any(opdom~=dom)
    error("Domain of integration must match spatial domain of the tensor-PI operator.")
end
% Allow only multipliers acting on linear monomials
if sum(has_multiplier)>1 || (sum(has_multiplier)==1 && d>1)
    error("Multiplier terms in tensor producs of PI operators are not supported.")
end
var2 = [var1,var2];

% List all possible orders of the variables (r1,t1,...,td,r2)
%   idx_mat(l,:) = [i,j,k] means a <= ti <= tj <= tk <= b in term l
omat = 0;
for i=1:d
    n_ords = size(omat,1);
    n_pos = size(omat,2)+1;
    omat_new = zeros(n_pos*n_ords,n_pos);
    for j=1:n_pos
        % Place variable t_i in position j
        omat_new(j:n_pos:end,:) = [omat(:,1:j-1),i*ones(n_ords,1),omat(:,j:end)];
    end
    omat = omat_new;
end

% For each possible ordering of the variables (r1,t1,...,td,r2), compute 
% the kernel in the associated term in the functional
n_ords = size(omat,1);
if isa(Kfun,'dpvar')
    Cfun = dpvar(zeros(pdim,ndim*n_ords));
else
    Cfun = polynomial(zeros(pdim,ndim*n_ords));
end
for i=1:n_ords
    param_i = 0;
    % Establish in which position the variable s is placed
    [~,pos_i] = find(omat(i,:)==0,1);
    leq_vars = omat(i,1:pos_i-1);
    geq_vars = omat(i,pos_i+1:end);
    for trm=1:ntrms
        % Multiply the appropriate parameters
        Ctmp = 1;
        for k=leq_vars
            % t_k <= s
            Ctmp = Ctmp.*Cparams{trm,k}{2};
        end
        for k=geq_vars
            % t_k >= s
            Ctmp = Ctmp.*Cparams{trm,k}{3};
        end
        param_i = param_i+Ctmp;
    end
    if ~isempty(param_i)
        Cfun(:,(i-1)*ndim+1:i*ndim) = param_i;
    end
end
omat = omat+1;
if any(has_multiplier)
    % Also add multiplier term
    omat = [zeros(1,d+1); omat];
    Cfun = [Cparams{1}{1},Cfun];
end
% Account for the kernel K
KCfun = Kfun*Cfun;


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
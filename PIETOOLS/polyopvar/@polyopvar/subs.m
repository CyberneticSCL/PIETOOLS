function fnew = subs(f,var,val)
% FNEW = SUBS(F,VAR,VAL) takes a distributed polynomial functional F, 
% and replaces factor VAR in the monomial defining F by the distributed
% polynomial VAL
%
% INPUTS
% - f:      'polyopvar' object representing a distributed polynomial
%           functional defined by a single distributed monomial, 
%               Z(x)=x1^{\otimes d1} \otimes ... \otimes xp^{\otimes dp},
%           taking the form f(x) = Kop*Z(z) for functional (intop) Kop;
% - var:    integer value between 1 and d=d1+...+dp, specyfing which of the
%           factors in the monomial defining f to substitute;
% - val:    'polyopvar' object representing the polynomial function,
%               p(x)(s) = (Aop*x)(s) + (Bop*x^{o 2})(s) + ....
%           with which to replace the specified factor in the monomial 
%           defining f;
%
% OUTPUTS
% - fnew:   'polyopvar' object representing the distributed polynomial
%           functional defined by f, but with the factor "var" replaced by
%           the specified polynomial, e.g. if var=i<d1,
%           fnew(x):=Kop(x1^{o(i-1)} o p(x) o x1^{o(d1-i)} o...o xp^{o dp};
%
% NOTES
% - The functional f can only be defined by a single monomial,
%   size(f.degmat,1) = 1;
% - f must be a functional, not a function, meaning that f.C.ops{1} must be
%   an 'intop' object
% - val  must be a polynomial function
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subs
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


if isa(var,'polyopvar')
    error("Substitution of state variables is not currently supported.")
end

% Make sure we are only substituting a single monomial
if size(f.degmat,1)==0
    return
elseif size(f.degmat,1)>1
    error("Substitution of individual factor in multiple monomials is not supported.")
end
if var>sum(f.degmat)
    error("Index of factor to substitute cannot exceed cumulative degree of monomial.")
end
% Only support substitution of polynomial functionals for now
if ~all(all(f.C.depmat1==0))
    error("Only substitution of polynomial functionals is currently suppoted.")
elseif size(f.C.ops,1)>1
    error("Substitution for multiple rows of polynomial functionals is currently not suppoted.")
end
% Only support substitution with polynomial functions for now
if ~isa(val,'polyopvar')
    error("Substitution of variables with non-'polyopvar' objects is currently not supported.")
end
if isa(val.C.ops{1},'intop')
    error("The new value cannot be a polynomial functional")
end
% Express the functional and new value in terms of the same variables
if ~isequal(f.varname,val.varname)
    [f,val] = common_vars(f,val);
end
% Extract the monomial defining f
degs1 = f.degmat;
xvar = find(var<=cumsum(degs1),1,'first');
degs1(xvar) = degs1(xvar)-1;
d1 = sum(degs1);
% Extract the operator acting on the monomial defining f
Kop = f.C.ops{1};
omat = Kop.omat;
var2 = Kop.pvarname;
var1 = var2(var);
Kdom = Kop.dom;
pvar2 = polynomial(var2);
ndim = Kop.matdim(2);
% Separate potential multiplier term
has_mult = false;
if all(omat(1,:)==0,2)
    if d1>1
        error("Integration with delta functions is supported only for quadratic functionals.")
    end
    has_mult = true;
end

fnew = f;
fnew.degmat = degs1+val.degmat;
fnew.C.ops = cell(1,size(val.degmat,1));
fnew.C.depmat2 = zeros(size(val.degmat,1),size(fnew.C.vars,1));
for i=1:size(val.degmat,1)
    % Extract the ith monomial defining val
    degs2 = val.degmat(i,:);
    d2 = sum(degs2);
    % Extract the tensor-PI operator acting on the ith monomial
    Cop = val.C.ops{1,i};
    %Cop.ops = Cop.ops(:,i);
    
    % Declare a full list of dummy variables used for integration of the 
    % monomial after substitution
    var2_full = cell(1,d1+d2);
    var2_full(1:d1) = var2(setdiff((1:d1+1),var));
    for k=1:d2
        var2_full{d1+k} = [var1{1},'_',num2str(k)];
    end
    % Establish which dummy variables are used in integrating each state
    % variable
    varnums1 = cumsum([0,degs1(1:end)]);
    varnums1 = mat2cell([varnums1(1:end-1);varnums1(2:end)],2,ones(1,size(degs1,2)));
    varnums1 = cellfun(@(a) (a(1)+1:a(2))',varnums1,'UniformOutput',false);
    varnums2 = d1+cumsum([0,degs2(1:end)]);
    varnums2 = mat2cell([varnums2(1:end-1);varnums2(2:end)],2,ones(1,size(degs1,2)));
    varnums2 = cellfun(@(a) (a(1)+1:a(2))',varnums2,'UniformOutput',false);
    varnums = [varnums1;varnums2];        % variables in column i correspond to state variable i
    varnums = cell2mat(varnums(:));
    % Reorder the dummy variables to account for the order of the state
    % variables in the new monomial
    [~,new2old_nums] = sort(varnums(:));
    var2_new = var2_full(varnums(:)');
    % Declare an empty polynomial in the new monomial
    fi = f;
    fi.degmat = zeros(0,size(fi.degmat,2));
    fi.C.ops = {};

    % Compute the composition of the delta-function integral with the
    % tensor-PI operator Cop
    if has_mult
        % int_{a}^{b} Kfun(s)*x(s)*(Cop{d}*x^{o d})(s) ds
        % Extract the parameters associated with the multiplier terms
        Kfun = Kop.params(:,1:ndim);
        % Impose the delta: var1 = var2
        var_tmp = Kfun.varname;
        var_tmp(strcmp(var_tmp,var1)) = var2_full(1);
        Kfun.varname = var_tmp;
        % Compute the kernels representing the integral
        [KCfun,omat_tmp] = int_delta(Cop,Kfun,var2_full(1),Kdom,var2_full(d1+1:end));
        % Account for the reordering of the variables
        iszero_trm = all(omat_tmp==0,2);
        omat_new = [omat_tmp(iszero_trm,:); new2old_nums(omat_tmp(~iszero_trm,:))];
        Kparam0 = KCfun(1,find(iszero_trm));
        Kparam_rem = KCfun(1,find(~iszero_trm));
        if isa(Kparam0,'dpvar') || isa(Kparam0,'polynomial')
            if isscalar(Kparam0.varname)
                Kparam0.varname = var2_new(1);
            elseif numel(Kparam0.varname)>1
                error("Something is wrong, the delta term depends on multiple variables.")
            end
        end
        KCfun = [Kparam0,Kparam_rem];
        % Declare a functional defined by the kernels KCfun acting on the
        % monomial i
        fnew_i = fi;
        fnew_i.C = tensopmat(intop(KCfun,omat_new,var2_new,Kdom));
        fnew_i.degmat = degs1+degs2;
        fnew_i = combine_terms(fnew_i);
    else
        fnew_i = 0;
    end

    % Compute the composition of the remaining integrals with the tensor-PI
    % operator Cop
    for j=has_mult+1:size(omat,1)
        % Establish the kernel of the jth term in the functional f
        cidcs = (j-1)*ndim+1:j*ndim;
        Kfun = Kop.params(:,cidcs);
        % Establish the interval over which factor specified by "var" is
        % integrated
        pos = find(omat(j,:)==var);
        dom = polynomial(Kdom);
        if pos>1
            dom(1) = pvar2(omat(j,pos-1));
        end
        if pos<d1+1
            dom(2) = pvar2(omat(j,pos+1));
        end
        % Compute the kernel defining the composition of the Kfun and the
        % tensor-PI operator Cop
        [KCfun,omat_tmp] = int(Cop,Kfun,var1,dom,var2_full(d1+1:end));
        % Change index of variables to match that in var2_new
        omat_j = omat(j,:);
        omat_j(omat_j>var) = omat_j(omat_j>var)-1;
        omat_new = omat_tmp+d1-(pos>1);
        if pos>1
            % Change index of variable pos-1 to match that in var2_new
            omat_new(omat_tmp==1) = omat_j(pos-1);
        end
        if pos<d1+1
            % Change index of variable pos+1 to match that in var2_new
            omat_new(omat_tmp==size(omat_tmp,2)) = omat_j(pos+1);
        end

        % Add variables to omat that are smaller than dom(1)
        if pos>1
            max_idcs = (omat_new==omat_j(pos-1))*(1:size(omat_new,2))';
            [omat_new,idcs] = combine_omats(omat_j(1:pos-2),omat_new,max_idcs);
            idcs_cell = mat2cell([idcs'-1*ndim;idcs'*ndim],2,ones(1,numel(idcs)));
            idcs_cell = cellfun(@(a) (a(1)+1:a(end)),idcs_cell,'UniformOutput',false);
            KCfun = KCfun(:,cell2mat(idcs_cell));
        end
        % Add variables to omat that are greater than dom(2)
        if pos<d1+1
            min_idcs = (omat_new==omat_j(pos+1))*(1:size(omat_new,2))';
            [omat_new,idcs] = combine_omats(omat_j(pos+2:end),omat_new,[],min_idcs);
            idcs_cell = mat2cell([idcs'-1*ndim;idcs'*ndim],2,ones(1,numel(idcs)));
            idcs_cell = cellfun(@(a) (a(1)+1:a(end)),idcs_cell,'UniformOutput',false);
            KCfun = KCfun(:,cell2mat(idcs_cell));
        end

        % Account for the reordering of the variables
        omat_new = new2old_nums(omat_new);
        % Declare a functional defined by the kernels KCfun acting on the
        % monomial i
        fj = fi;
        fj.C = tensopmat(intop(KCfun,omat_new,var2_new,Kdom));
        fj.degmat = degs1+degs2;
        fj = combine_terms(fj);
        % Add to the other polynomial
        fnew_i = fnew_i+fj;
    end
    % Set the functional acting on the ith monomial
    fnew.C.ops{1,i} = fnew_i.C.ops{1};
    fnew.C.depmat2(i) = fnew_i.C.depmat2(1,:);
end

end


function [omat,idcs] = combine_omats(omat1,omat2,max_idcs,min_idcs)
% Returns all possible orderings of the variables appearing in omat1 and 
% omat2, that preserves the ordering specified by omat1 and omat2
%
% INPUTS
% - omat1:  1 x d1 array of integers specifying some order of variables;
% - omat2:  m x d2 array of integers specifying the order of variables
%           distinct from the variables listed in omat1;
% - max_idcs:   m x 1 array specifying for each of the rows of omat2 the
%               maximal position in which the variables of omat1 can be
%               placed, enforcing an ordering
%               omat1(i)<omat2(j,max_idcs(j)) for each element i and j;
% - min_idcs:   m x 1 array specifying for each of the rows of omat2 the
%               minimal position in which the variables of omat1 can be
%               placed, enforcing an ordering
%               omat1(i)>omat2(j,min_idcs(j)) for each element i and j;
%
% OUTPUTS
% - omat:   n x d1+d2 array specifying all possible orderings of the
%           variables that preserve the orderings from omat1 and omat2.
%           That is, if omat(i,j1)=omat2(l,k1) and omat(i,j2)=omat2(l,k2) 
%           for some row i and l, columns j1, j2, k1, and k2, then k1<=k2 
%           implies j1<=j2, and similarly ofr omat1. 
% - idcs:   n x 1 array speciyfing for each of the rows of omat which row
%           of omat2 it is based on.

% Check the first input
if size(omat1,1)~=1
    error("First order matrix must contain at most one element.")
end
nvars1 = size(omat1,2);

% Initialize the new order matrix as omat2
omat = omat2;
nvars = size(omat,2);

% Keep track of which row in omat2 each row in omat corresponds to
idcs = (1:size(omat2,1))';

% Keep track of the maximal position each variable in omat1 can be placed
if nargin<=2 || isempty(max_idcs)
    max_idcs = (nvars+1)*ones(size(omat,1),1);
elseif numel(max_idcs)~=size(omat,1)
    error("Number of indices must match number of rows of omat2.")
end
% Keep track of the lowest position each variable in omat2 can be placed
if nargin<=3 || isempty(min_idcs)
    min_idcs = zeros(size(omat,1),1);
elseif numel(min_idcs)~=size(omat,1)
    error("Number of indices must match number of rows of omat2.")
end

for j=1:nvars1
    % Initialize empty omat including variables omat1(1:j)
    nvars = nvars+1;
    omat_full = zeros(0,nvars);
    min_idcs_full = [];
    max_idcs_full = [];
    idcs_full = [];
    % Extract variable from omat1 to add to omat
    var_num = omat1(j);
    for i=1:size(omat,1)
        % Extract the ith row of the order matrix
        omat_i = omat(i,:);
        min_idx = min_idcs(i);
        max_idx = max_idcs(i);
        if min_idx>=max_idx
            continue
        end

        % Place variable var_num in each position in ith row of omat below 
        % the max_idx
        Dmat = var_num*diag(ones(1,nvars));
        Dmat = Dmat(min_idx+1:max_idx,:);
        omat_i = repmat(omat_i,max_idx-min_idx,1);
        omat_i = tril([omat_i,zeros(max_idx-min_idx,1)],-1+min_idx)+Dmat+triu([zeros(max_idx-min_idx,1),omat_i],1+min_idx);
        omat_full = [omat_full; omat_i];

        % Keep track of position in which var_num is placed
        max_idcs_full = [max_idcs_full; max_idx+ones(max_idx-min_idx,1)];
        min_idcs_full = [min_idcs_full; min_idx+(1:max_idx-min_idx)'];
        % Keep track of which row of omat2 each row of omat corresponds to
        idcs_full = [idcs_full; idcs(i)*ones(max_idx-min_idx,1)];
    end
    omat = omat_full;
    min_idcs = min_idcs_full;
    max_idcs = max_idcs_full;
    idcs = idcs_full;
end

end
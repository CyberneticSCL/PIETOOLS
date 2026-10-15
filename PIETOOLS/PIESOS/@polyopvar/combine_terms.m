function f = combine_terms(f)
% F = COMBINE_TERMS(F) takes a distributed polynomial and adds the
% coefficients acting on the same monomial, to minimize the size of the
% distributed monomial basis.
%
% INPUTS
% - f:  'polyopvar' object representing a distributed polynomial;
%
% OUTPUTS
% - f:  'polyopvar' object representing the same polynomial as the input,
%       with only unique monomials, so that each monomial (row of f.Z)
%       appears once (see NOTES);
%
% NOTES
% - Monomials with the same Z but a different number of factors (factors
%   that are multipliers from R, see 'polyopvar') are different monomials.
% - If f has no output spatial variables (f.C.depmat.out is zero), the
%   monomials of the same degree in each state, in any order, are combined
%   as well, by reordering the integrals in the operators acting on them,
%   e.g.   int_{a}^{b} int_{c}^{d} y(s1) x(s2) ds2 ds1
%        = int_{c}^{d} int_{a}^{b} x(s1) y(s2) ds1 ds2.
%   The states in the resulting monomials are in order of increasing
%   index. Monomials with operators of type 'tensopvar' are not combined.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - combine_terms
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
% DJ, 10/15/2026: Initial coding;

if ~f.chkval
    f = initialize(f);
end

% Extract the coefficients
Cop = f.C;
ops = Cop.ops;
depmat = Cop.depmat.in;
nvars = numel(f.varname);

% Identify a unique set of monomials  in terms of which f is defined
Zmat = f.Z;
% Monomials with the same Z but different numbers of factors are different monomials: e.g.
% x o c, with c a multiplier from R, Z=[1,0], and x, Z=[1], which are the same row once padded
% to the same width, see `polyopvar`
[~,nopsZ,is_empty] = num_factors(f);
nopsZ = factor_key(Zmat,nopsZ,is_empty);        % empty monomials merge with those of the same Z
[Pmat,Zmat_new] = uniquerows_integerTable([sum(Zmat>0,2),nopsZ,Zmat]);
Zmat_new = Zmat_new(:,3:end);
nZ_new = size(Zmat_new,1);

% Add operators acting on the same monomial
ops_new = cell(size(Cop.ops,1),size(Zmat_new,1));
depmat_new = zeros(nZ_new,size(depmat,2));
dims_new = zeros(nZ_new,1);
degmat = zeros(nZ_new,nvars);
for j=1:nZ_new
    % Establish which operators act on the jth monomial
    param_idcs = find(Pmat(:,j));
    depmat_new(j,:) = depmat(param_idcs(1),:);
    dims_new(j) = Cop.dims.in(param_idcs(1));
    % Add the different operators acting on this monomial
    for i=1:size(ops,1)
        Ctmp = ops{i,param_idcs(1)};
        for k=2:numel(param_idcs)
            Ctmp = Ctmp + ops{i,param_idcs(k)};
        end
        % Set the operator acting on the jth monomial
        ops_new{i,j} = Ctmp;
    end
    % Also determine the degree of this monomial in each variable
    for var=1:nvars
        degmat(j,var) = sum(Zmat_new(j,:)==var);
    end
end
Cop.ops = ops_new;
Cop.depmat.in = depmat_new;
Cop.dims.in = dims_new;

if any(any(Cop.depmat.out))
    f.C = Cop;
    f.Z = Zmat_new;
    f.chkval = true;
    return
end

% For polynomial functionals (mapping to R), we also combine monomials
% that are of the same degree, by re-ordering the variables, e.g.
%   int_{a}^{b} int_{c}^{d} y(s1) x(s2) ds2 ds1
%       = int_{c}^{d} int_{a}^{b} x(s1) y(s2) ds1 ds2

% Determine the function spaces on which the states are defined
xdims = f.dims.in;
xdepmat = f.depmat.in;

% Determine the degrees of each state in the monomials
% e.g. x o y o x o z is of degree [2,1,1] in [x,y,z]
[degmat2,idcs1,idcs2] = unique(degmat,"rows","stable");

% Collect monomials that share the same degrees
% e.g. x o y o x o z and x o x o y o z
Zcell = cell(1,size(degmat2,1));
opcell = cell(1,size(degmat2,1));
old_deg_nums = num2cell(idcs1(:)');
dmax = size(Zmat,2);
for j=1:size(degmat2,1)
    % Find all monomials of the considered degrees
    is_degj = idcs2==j;
    Zmat_j = Zmat_new(is_degj,:);        
    ops_j = Cop.ops(:,is_degj);
    degs_j = degmat2(j,:);
    if isa(ops_j{1},'tensopvar')
        old_deg_nums{j} = idcs1(j)*ones(1,nnz(is_degj));
        Zcell{j} = Zmat_j';
        opcell{j} = ops_j;
        continue
    end
    % Declare a canonical monomial of the specified degrees
    %   s1^p1 o ... o sd^pd
    Ztrue = repelem(1:nvars,degs_j);
    Ztrue = [Ztrue,zeros(1,dmax-numel(Ztrue))];
    Zcell{j} = Ztrue';
    opcell{j} = num2cell(zeros(size(ops_j,1),1));
    % Reorder integrals to match new ordering of variables
    for i=1:size(ops_j,1)
        for k=1:size(Zmat_j,1)
            if isempty(ops_j{i,k}) || isequal(ops_j{i,k},zeros(size(ops_j{i,k})))
                continue
            else
                state_vals = Zmat_j(k,Zmat_j(k,:)>0);
                [~,fctr_order] = sort(state_vals);                    
                op_tmp = reorder_states(ops_j{i,k},fctr_order,xdims(state_vals),xdepmat(state_vals,:));
                opcell{j}{i} = opcell{j}{i} + op_tmp;
            end
        end
        if isa(opcell{j}{i},'double')
            % All operators acting on the monomials of these degrees are zero
            opcell{j}{i} = [];
        elseif ~isempty(opcell{j}{i})
            % Combine terms in the functional that represent the same
            % integral on the specified monomial
            opcell{j}{i} = combine_terms(opcell{j}{i},degs_j,xdims,xdepmat);
        end
    end
end
% Update coefficient operator and monomials to remove duplicates
old_deg_nums = [old_deg_nums{:}];
Cop.ops = [opcell{:}];
Cop.depmat.in = Cop.depmat.in(old_deg_nums,:);
Cop.dims.in = Cop.dims.in(old_deg_nums);
Zmat_new = sparse([Zcell{:}]');

% Set the new polynomial
f.C = Cop;
f.Z = Zmat_new;
f.chkval = true;

end
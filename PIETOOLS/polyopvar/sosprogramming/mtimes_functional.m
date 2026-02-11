function [Cfun] = mtimes_functional(Kop,Fx,cntr)
% CFUN = MTIMES_FUNCTIONAL(KOP,FX) computes the composition of a functional 
% operator,
%   KOP*w = sum_{i=1}^{q} int_{a}^{b} ... int_{t_{ord(i,d-1)}}^{b} 
%               K_{i}(t1,...,td)*w(t1,...,td) dt_{ord(i,d)} ... dt_{ord(i,1)}
% with a tensor-PI operator
%   (F(x))(s_{1},...,s_{d}) = F{1}(x)(s_{1}) * ... * F{d}(x)(s_{d})
% where 
%   F{i}(x)(s_{i}) = Fop{i}*Z(x) = (T1*x1)(s_{i}) * ... * (Tn*xn)(s_{i})
% for a single monomial term Z(x)(s1,...,sn) = x1(s1)*...*xn(sn), and 3-PI
% operators Tj.
%
% INPUTS
% - Kop:    'struct' representing a functional operator, mapping a degree d
%           distributed monomial to a real scalar. Must have fields:
%               params: 1 x q 'polynomial' or 'dpvar' object,
%                       with each element i corresponding to 
%                       the kernels in the ith term of the functional;
%               omat:   q x d array of integers, with row i
%                       specifying the order of the variables in the
%                       integral associated with the ith kernel.
%                       Specifically, if omat(i,:) = [k1,k2,...,kd], 
%                       then term i is defined by the integral
%                           int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b};
%               vars:   d x 1 'polynomial' array specifying the names of
%                       the dummy variables used in definition of the
%                       kernels in Kop.params;
%               dom:    interval [a,b] over which to integrate;
% - Fx:     1 x d cell, with each element a 'polyopvar' object representing
%           a tensor PI operator applied to a single monomial. 
% - cntr:   internal variable to keep track of which call to the function
%           this is. Should not be specified by user!
%
% OUTPUT
% - Cfun:   'polyopvar' object representing the composition of the
%           functional defined by Kop with the product of the polynomials
%           specified by Fx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes_functional
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
% DJ, 01/25/2026: Initial coding


% Extract the parameters defining the functional Kop
Kparams = Kop.params;
idx_mat = Kop.omat;
Kvars = Kop.vars;
Kdom = Kop.dom;
if any(Kop.matdim~=1)
    error("Functional must be defined by scalar-valued kernels.")
end

if ~isa(Fx,'cell')
    error("Factors to integrate must be specified as a cell of 'polyopvar' objects.")
elseif numel(Fx)>numel(Kvars)
    error("Number of factors should not exceed the number of variables of integration.")
end
if nargin<=2
    cntr = numel(Fx);
end


% Loop over all the terms in the functional defined by Kop
for ii=1:size(idx_mat,1)
    % Extract the order of the dummy variables for the considered term
    idx_ii = idx_mat(ii,:);

    % Find the largest dummy variable for which we have a non-identity
    % factor
    o_idx = find(idx_ii==cntr);
    midx = idx_ii(o_idx);
    var_max = Kvars(midx);
    Fii = Fx{midx};

    % Declare a set of dummy variables associated with the operator midx
    n_fctrs = sum(Fii.degmat(1,:));
    var2 = polynomial(zeros(n_fctrs,1));
    for jj=1:n_fctrs
        var2_jj = [var_max.varname{1},'_',num2str(jj)];
        var2(jj) = polynomial({var2_jj});
    end


    % Establish the limits of integration: 
    %       var(idx_ii(m-1)) <= var(m) <= var(idx_ii(m+1))
    if o_idx==numel(idx_ii)
        U = Kdom(2);
    else
        U = Kvars(idx_ii(o_idx+1));
    end
    if o_idx==1
        L = Kdom(1);
    else
        L = Kvars(idx_ii(o_idx-1));
    end
    dom_ii = [L,U];

    % Perform the composition of the integral operator with the factor
    % along variable var_max
    %Cii = quad2lin_term(1,Kcell{ii},Fii,dom_ii,var_max);
    Cii = mtimes_functional_fctr(Kparams(ii),Fii,dom_ii,var_max,var2);

    % Include the new variables into the total list of variables
    vars_full = [Kvars(1:midx-1,:); Cii.vars; Kvars(midx+1:end,:)];   % add the variables of Cii to Kop
    nvars_ii = numel(Cii.vars);    
    idx_ii(idx_ii>midx) = idx_ii(idx_ii>midx) + nvars_ii-1; % relabel some of the old variables
    nvars_full = numel(vars_full);

    % Also augment the order matrix, interpreting
    % omat(i,j)=-k   as   var2(k) <= dom_ii(1) = Kvars(idx_ii(cntr-1))
    % NOTE: you can still have e.g. var2(k) >= Kvars(idx_ii(cntr-2))
    omat_ii = Cii.omat + (midx-1).*sign(Cii.omat);  % relabel the newly introduced variables
    omat_full = zeros(0,nvars_full);
    param_idcs = [];
    for jj=1:size(omat_ii,1)
        neg_idcs = find(omat_ii(jj,:)<0);
        neg_vals = abs(omat_ii(jj,neg_idcs));
        neg_idcs = neg_idcs + o_idx-1;
        im_idcs = find(real(omat_ii(jj,:))~=omat_ii(jj,:));
        im_vals = abs(omat_ii(jj,im_idcs));
        im_idcs = im_idcs + o_idx-1;
        nrows = cumprod([1,neg_idcs-1,nvars_full-im_idcs(end:-1:1)]);
        omat_jj = zeros(nrows(end),nvars_full);
        omat_jj(1,:) = [idx_ii(1:o_idx-1), omat_ii(jj,:), idx_ii(o_idx+1:end)];
        for kk=1:numel(neg_idcs)
            % Place the variable in all positions below Kvars(idx(cntr-1))
            neg_idx = neg_idcs(kk);
            Dmat_kk = neg_vals(kk)*eye(neg_idx-1,nvars_full);
            omat_kk = zeros(nrows(kk+1),nvars_full);
            for ll=1:nrows(kk)
                omat_ll = repmat(omat_jj(ll,:),[neg_idx-1,1]);
                omat_rem = omat_ll(:,[1:neg_idx-1,neg_idx+1:end]);
                Lmat_kk = [tril(omat_rem,-1),zeros(neg_idx-1,1)];
                Umat_kk = [zeros(neg_idx-1,1),triu(omat_rem)];
                omat_ll = Lmat_kk + Dmat_kk + Umat_kk;
                omat_kk((ll-1)*(neg_idx-1)+1:ll*(neg_idx-1),:) = omat_ll;
            end
            omat_jj(1:nrows(kk+1),:) = omat_kk;
        end
        for kk=1:numel(im_idcs)
            % Place the variable in all positions above Kvars(idx(cntr+1))
            im_idx = im_idcs(numel(im_idcs)-kk+1);
            Dmat_kk = [zeros(nvars_full-im_idx,im_idx), im_vals(numel(im_idcs)-kk+1)*eye(nvars_full-im_idx)];
            omat_kk = zeros(nrows(kk+numel(neg_idcs)+1),nvars_full);
            for ll=1:nrows(kk+numel(neg_idcs))
                omat_ll = repmat(omat_jj(ll,:),[nvars_full-im_idx,1]);
                Lmat_kk = [omat_ll(:,1:im_idx-1),tril(omat_ll(:,im_idx+1:end)),zeros(nvars_full-im_idx,1)];
                Umat_kk = [zeros(nvars_full-im_idx,im_idx),triu(omat_ll(:,im_idx+1:end),1)];
                omat_ll = Lmat_kk + Dmat_kk + Umat_kk;
                omat_kk((ll-1)*(nvars_full-im_idx)+1:ll*(nvars_full-im_idx),:) = omat_ll;
            end
            omat_jj(1:nrows(kk+numel(neg_idcs)+1),:) = omat_kk;
        end
        omat_full = [omat_full; omat_jj];
        param_idcs = [param_idcs; jj*ones(size(omat_jj,1),1)];
    end
    %omat_full = [repmat(idx_ii(1:cntr-1),[size(omat_ii,1),1]), omat_ii, repmat(idx_ii(cntr+1:end),[size(omat_ii,1),1])];

    % Augment the functional to include the full list of variables,
    % getting rid of duplicates
    [P,omat_unique] = uniquerows_integerTable(omat_full);    % P*omat_unique = omat
    P = sparse(param_idcs,1:numel(param_idcs),1,numel(Cii.params),numel(param_idcs))*P;
    Cii.params = Cii.params*P;
    Cii.omat = omat_unique;
    Cii.vars = vars_full;
    Cii.dom = Kdom;

    % % Remove duplicate terms from the functional operator
    state_idcs_m = midx+(0:numel(Fii.degmat)-1);
    state_idcs_m = repelem(state_idcs_m,Fii.degmat);
    state_idcs = [1:midx-1,state_idcs_m,state_idcs_m(end)+(1:numel(Kvars)-midx)];
    Cii = combine_terms(Cii,state_idcs);

    % Apply the augmented functional to the remaining factors
    if cntr>1
        Cii = mtimes_functional(Cii,Fx,cntr-1);
        %Cii = Cfun.C.ops{1};
    end

    % Add the kernels to the functional operator Cop
    if ii==1
        % Initialize an operator representing the functional
        Cop = Cii;
    else
        if ~isequal(Cop.vars,Cii.vars)
            error("Something is going wrong, the variables don't match...")
        end
        omat = [Cop.omat;Cii.omat];
        params_full = [Cop.params,Cii.params];
        [P,omat_unique] = uniquerows_integerTable(omat);    % P*omat_unique = omat
        Cparams = params_full*P;
        Cop.omat = omat_unique;
        Cop.params = Cparams;
    end
end



% Reorder the spatial variables to account for a possible
% reordering of the state variables,
%   [x1(s1)*x1(s2)*x2(s3)*x2(s4)]*[x1(s5)*x2(s6)]
%       --> [x1(s1)*x1(s2)*x1(s3)]*[x2(s4)*x2(s5)*x2(s6)]
if cntr==numel(Fx)
    varname_full = {};
    deg_full = [];
    varmat_full = [];
    pvarname_full = {};
    for ii=1:numel(Fx)
        varname_full = [varname_full; Fx{ii}.varname];
        deg_full = [deg_full, Fx{ii}.degmat]; 
        varmat_full = blkdiag(varmat_full,Fx{ii}.varmat);
        pvarname_full = [pvarname_full; Fx{ii}.pvarname];
    end
    % Establish a unique set of state variables
    [varname,old2new_idcs,new2old_idcs] = unique(varname_full);  % var_full = var_unique(new2old_idcs);
    degmat = deg_full*sparse((1:numel(new2old_idcs)),new2old_idcs,1);
    % Reorder variables in Cop to match new ordering of states
    var_idcs = mat2cell([[0,cumsum(deg_full(1:end-1))];cumsum(deg_full)],2,ones(1,numel(deg_full)));
    var_idcs = cellfun(@(a)(a(1)+1:a(2)),var_idcs,'UniformOutput',false);
    Cvars_new = Cop.vars;
    strt_idx = 0;
    Cvar_order = zeros(numel(Cop.vars),1);
    for jj=1:numel(varname)
        isvar_jj = new2old_idcs'==jj;
        var_idcs_jj = cell2mat(var_idcs(isvar_jj));
        nvars_jj = numel(var_idcs_jj);
        Cvars_new(strt_idx+(1:nvars_jj)) = Cop.vars(var_idcs_jj);
        Cvar_order(var_idcs_jj) = strt_idx+(1:nvars_jj);
        strt_idx = strt_idx + nvars_jj;
    end
    Cop.vars = Cvars_new;
    % Reorder columns of omat to match new ordering of variables
    omat_new = Cop.omat;
    for jj=1:numel(Cvar_order)
        omat_new(Cop.omat==jj) = Cvar_order(jj);
    end
    Cop.omat = omat_new;
    % Establish a unique set of independent variables
    [pvarname_unique,old2new_idcs_p] = unique(pvarname_full);
    varmat_unique = varmat_full(old2new_idcs,old2new_idcs_p);
end


% Also deal with possible multiplier term,
%       int_{a}^{b} Kop{3}(s)*(Rop1*u)(s)*(Rop2*u)(s) ds
% Note that quad2lin already accounts for possible reordering of the state
% variables
if numel(Kparams)>size(idx_mat,1)
    if numel(Fx)~=2
        error("Multiplier terms may appear only for quadratic monomial.")
    end
    Ktmp = Kparams(1,end);
    var1 = polynomial(Ktmp.varname);
    C_mult = quad2lin_term_old(Ktmp,Fx{1},Fx{2},Kop.dom,var1,Cop.vars);
    
    % % Match the variable names in C_mult with those used in Cop
    % varname_old = C_mult.params.varname;
    % varname_new = varname_old;
    % for jj=1:numel(C_mult.vars)
    %     old_var_idx = isequal(C_mult.vars(jj),polynomial(varname_old));
    %     if ~any(old_var_idx)
    %         continue
    %     end
    %     varname_new(old_var_idx) = Cop.vars(jj).varname;
    % end
    % C_mult.params.varname = varname_new;


    % Combine parameters involving the same integral
    omat = [Cop.omat;C_mult.omat];
    params_full = [Cop.params,C_mult.params];
    [P,omat_unique] = uniquerows_integerTable(omat);    % P*omat_unique = omat
    Cop.params = params_full*P;
    Cop.omat = omat_unique;
end

if cntr<numel(Fx)
    % If we're not done, return the operator rather than a polynomial
    Cfun = Cop;
else
    % Return a 'polyopvar' object representing the polynomial functional
    Cfun = polyopvar();
    Cfun.degmat = degmat;
    Cfun.varname = varname;
    Cfun.C = tensopvar();
    Cfun.C.ops{1} = Cop;
    Cfun.pvarname = pvarname_unique;
    Cfun.varmat = varmat_unique;
end

end
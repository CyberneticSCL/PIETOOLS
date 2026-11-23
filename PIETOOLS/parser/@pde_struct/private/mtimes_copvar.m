function PDE_out = mtimes_copvar(P,PDE_in)
% PDE_OUT = MTIMES_COPVAR(P,PDE_IN) declares the PDE terms corresponding to
% the PI operator P applied to the PDE variable or terms PDE_IN,
%   (P*x)_i(s^i) = sum_j sum_gamma int K_ij^gamma(s^i,s_dum) x_j(s_dum) ds_dum,
% where in each direction shared by the input and output of block (i,j)
% the parameter gamma takes a multiplier (gamma_k=1), an integral from the
% lower limit to s_k (gamma_k=2), or from s_k to the upper limit
% (gamma_k=3), and each input-only variable is integrated over its domain.
%
% INPUT
% - P:      'copvar' or 'sopvar' object, with total input dimension equal
%           to the size of PDE_in;
% - PDE_in: 'pde_struct' object representing a single PDE state or input
%           variable, or one or several free terms. The terms acted on by
%           column j of P can depend only on the input variables of that
%           column, each on the same domain as in P;
%
% OUTPUT
% - PDE_out: 'pde_struct' object representing the free terms P*PDE_in, one
%            equation per row of P, or a single equation if PDE_in is a
%            single equation;
%
% NOTES
% - Outputs and temporal derivatives of states cannot be multiplied with P.
% - A term partially integrated along a direction in which P integrates
%   is not supported.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes_copvar
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding DJ - 11/22/2026

% % % Process the inputs

% % Represent the operator as a container of 'sopvar' blocks.
if isa(P,'opvar')
    % The conversion requires the dummy variable to be named as var1 with
    % _dum appended; only the kernels R1 and R2 depend on it.
    dum_var = pvar([P.var1.varname{1},'_dum']);
    if ~isequal(P.var2,dum_var)
        P.R.R1 = subs(polynomial(P.R.R1),P.var2,dum_var);
        P.R.R2 = subs(polynomial(P.R.R2),P.var2,dum_var);
        P.var2 = dum_var;
    end
    P = opvar2copvar(P);
elseif isa(P,'sopvar')
    P = copvar({P});
end

% % Check that the terms can be multiplied with the operator.
% Outputs and temporal derivatives can only appear on the left-hand side.
for ii=1:numel(PDE_in.free)
    if ~isfield(PDE_in.free{ii},'term'),    continue;   end
    for jj=1:numel(PDE_in.free{ii}.term)
        if is_LHS_term(PDE_in.free{ii}.term{jj})
            error("Multiplication of PI operators with outputs or temporal derivatives of states is not supported.")
        end
    end
end
% The input dimension of the operator must match the size of the terms.
n_eqs = size(PDE_in,'free','vec_size_tot');
if sum(P.dim_in)~=n_eqs
    error("The input dimension of the operator should match the size of the PDE terms.")
elseif n_eqs==0
    PDE_out = PDE_in;
    return
elseif sum(P.dim_out)==0
    % Empty output.
    PDE_out = PDE_in;
    PDE_out.free = cell(0,1);
    return
end
% The spatial variables of the terms must exist on the same domains as in
% the operator.
objs = {'x','u','w'};
for kk=1:numel(objs)
    for ii=1:numel(PDE_in.(objs{kk}))
        comp_ii = PDE_in.(objs{kk}){ii};
        if ~isfield(comp_ii,'vars') || ~isfield(comp_ii,'dom') || isempty(comp_ii.vars)
            continue
        end
        vnames_ii = arrayfun(@(l) comp_ii.vars(l,1).varname{1},1:size(comp_ii.vars,1),'UniformOutput',false);
        [is_P_var,idcs_P] = ismember(vnames_ii,P.vars);
        if any(any(comp_ii.dom(is_P_var,:)~=P.dom(idcs_P(is_P_var),:)))
            error("The spatial variables of the PDE terms should exist on the same domains as in the operator.")
        end
    end
end


% % % Apply each block P.C{i,j} to the rows of PDE_in in column j, and
% % % add the results in each row.
M = numel(P.dim_out);   N = numel(P.dim_in);
ridcs_in = [0;cumsum(P.dim_in(:))];
PDE_rows = cell(M,1);
for j=1:N
    % Extract the terms in column j, using a selection matrix.
    if N==1
        PDE_j = PDE_in;
    else
        sel_j = zeros(P.dim_in(j),n_eqs);
        sel_j(:,ridcs_in(j)+1:ridcs_in(j+1)) = eye(P.dim_in(j));
        PDE_j = sel_j*PDE_in;
    end
    % Check that the terms depend only on input variables of column j.
    vnames_in = P.vars(P.space_in(j,:));
    for ii=1:numel(PDE_j.free)
        vars_ii = PDE_j.free{ii}.vars;
        if ~isempty(vars_ii) && ~all(ismember(vars_ii.varname,vnames_in))
            error("The PDE terms should depend only on spatial variables in the input space of the operator.")
        end
    end
    % Apply each nonzero block in column j.
    for i=1:M
        if isempty(P.C{i,j})
            continue
        end
        PDE_ij = mtimes_sopvar(P.C{i,j},PDE_j);
        if isempty(PDE_rows{i})
            PDE_rows{i} = PDE_ij;
        else
            PDE_rows{i} = PDE_rows{i}+PDE_ij;
        end
    end
end
% A row without nonzero blocks gives zero terms of the size of the row.
for i=1:M
    if isempty(PDE_rows{i})
        PDE_rows{i} = zeros(P.dim_out(i),n_eqs)*PDE_in;
    end
end
% A single input equation gives a single output equation, stacking the
% rows of P as for multiplication with a matrix; otherwise there is one
% equation per row of P.
if isscalar(PDE_in.free) && M>1
    ridcs_out = [0;cumsum(P.dim_out(:))];
    for i=1:M
        emb_i = zeros(ridcs_out(end),P.dim_out(i));
        emb_i(ridcs_out(i)+1:ridcs_out(i+1),:) = eye(P.dim_out(i));
        PDE_rows{i} = emb_i*PDE_rows{i};
    end
    PDE_out = PDE_rows{1};
    for i=2:M
        PDE_out = PDE_out+PDE_rows{i};
    end
else
    PDE_out = vertcat(PDE_rows{:});
end

end



function PDE_out = mtimes_sopvar(P,PDE_in)
% PDE_OUT = MTIMES_SOPVAR(P,PDE_IN) declares the PDE terms corresponding to
% the 'sopvar' object P applied to the free terms PDE_IN.
%
% INPUT
% - P:          'sopvar' object of dimension q x p;
% - PDE_in:     'pde_struct' object representing free terms of total size
%               p, depending only on variables in P.vars.in;
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the free terms P*PDE_in,
%               of total size q;

% % Partition the input variables of P: shared variables S3 (one direction
% % of the multi-index of P.params each, in sorted order) and input-only
% % variables S1 (integrated over their domain).
vars_S3 = intersect(P.vars.in,P.vars.out);
vars_S1 = setdiff(P.vars.in,P.vars.out);
[~,idcs] = ismember(vars_S3,P.vars.in);     dom_S3 = P.dom.in(idcs,:);
[~,idcs] = ismember(vars_S1,P.vars.in);     dom_S1 = P.dom.in(idcs,:);
n3 = numel(vars_S3);

% % Collect the variables on which the terms depend.
vnames_PDE = cell(0,1);
for ii=1:numel(PDE_in.free)
    vnames_PDE = [vnames_PDE; PDE_in.free{ii}.vars.varname];                %#ok<AGROW>
end

% % The kernels use the dummy variables of P.vars.in, named with _dum
% % appended; rename a dummy variable that the terms already use.
K = P.kernels;
vnames_dum = strcat(P.vars.in,'_dum');
vnames_used = [vnames_PDE; term_vnames(PDE_in)];
for l=1:numel(vnames_dum)
    vname_new = vnames_dum{l};
    while ismember(vname_new,vnames_used)
        vname_new = [vname_new,'_'];                                        %#ok<AGROW>
    end
    if ~strcmp(vname_new,vnames_dum{l})
        K = cellfun(@(Kk) subs(Kk,pvar(vnames_dum{l}),pvar(vname_new)),K,'UniformOutput',false);
        vnames_dum{l} = vname_new;
    end
end

% % Add the terms of each nonzero parameter, the kernel K_k times the
% % terms evaluated at the dummy variables, integrated over the dummy
% % variables.
PDE_out = [];
for k=1:numel(P.params)
    if isempty(P.params{k}) || ~any(P.params{k}(:))
        continue
    end
    % Multi-index of the parameter, direction 1 varying fastest.
    gamma = ones(1,n3);
    if n3==1
        gamma = k;
    elseif n3>1
        gamma_cell = cell(1,n3);
        [gamma_cell{:}] = ind2sub(3*ones(1,n3),k);
        gamma = [gamma_cell{:}];
    end
    % Variables integrated out: shared variables with gamma_k>1 from the
    % lower limit to s_k (gamma_k=2) or from s_k to the upper limit
    % (gamma_k=3), and input-only variables over their domain.
    is_int = gamma>1;
    vnames_int = [vars_S3(is_int), vars_S1];
    L = num2cell([dom_S3(is_int,1); dom_S1(:,1)]);
    U = num2cell([dom_S3(is_int,2); dom_S1(:,2)]);
    gamma_int = gamma(is_int);
    for l=1:numel(gamma_int)
        if gamma_int(l)==2
            U{l} = pvar(vnames_int{l});
        else
            L{l} = pvar(vnames_int{l});
        end
    end
    if isempty(vnames_int)
        % Multiplier only.
        PDE_k = K{k}*PDE_in;
    else
        % Evaluate the terms at the dummy variables, for the variables
        % they depend on, multiply with the kernel, and integrate.
        vars_int = cellfun(@pvar,vnames_int(:),'UniformOutput',false);
        vars_int = vertcat(vars_int{:});
        [~,idcs] = ismember(vnames_int,P.vars.in);
        vars_dum = cellfun(@pvar,vnames_dum(idcs)','UniformOutput',false);
        vars_dum = vertcat(vars_dum{:});
        is_subs = ismember(vnames_int,vnames_PDE);
        check_integrated(PDE_in,vnames_int(is_subs));
        PDE_k = PDE_in;
        if any(is_subs)
            PDE_k = subs(PDE_k,vars_int(is_subs),vars_dum(is_subs));
        end
        L = cellfun(@polynomial,L,'UniformOutput',false);
        U = cellfun(@polynomial,U,'UniformOutput',false);
        PDE_k = int(K{k}*PDE_k,vars_dum,vertcat(L{:}),vertcat(U{:}));
    end
    if isempty(PDE_out)
        PDE_out = PDE_k;
    else
        PDE_out = PDE_out+PDE_k;
    end
end
% A block with only zero parameters gives zero terms.
if isempty(PDE_out)
    PDE_out = zeros(P.dims(1),P.dims(2))*PDE_in;
end

end



function check_integrated(PDE_in,vnames)
% CHECK_INTEGRATED(PDE_IN,VNAMES) throws an error if a term in PDE_IN is
% partially integrated along one of the variables VNAMES, which cannot be
% integrated a second time with a separate kernel.
%
% INPUT
% - PDE_in:     'pde_struct' object representing free terms;
% - vnames:     cellstr of the names of the variables;

objs = {'x','u','w'};
for ii=1:numel(PDE_in.free)
    if ~isfield(PDE_in.free{ii},'term'),    continue;   end
    for jj=1:numel(PDE_in.free{ii}.term)
        trm = PDE_in.free{ii}.term{jj};
        for f=1:numel(trm)
            if ~isfield(trm(f),'I'),    continue;   end
            for kk=1:numel(objs)
                if ~isfield(trm(f),objs{kk}) || isempty(trm(f).(objs{kk}))
                    continue
                end
                obj_vars = PDE_in.(objs{kk}){trm(f).(objs{kk})}.vars;
                for l=1:numel(trm(f).I)
                    if ~isempty(trm(f).I{l}) && ~isdouble(polynomial(trm(f).I{l})) ...
                            && ismember(obj_vars(l,1).varname{1},vnames)
                        error("Multiplication of PI operators with terms partially integrated along the same variable is not supported.")
                    end
                end
            end
        end
    end
end

end



function vnames = term_vnames(PDE_in)
% VNAMES = TERM_VNAMES(PDE_IN) returns the names of the variables at which
% the terms in PDE_IN are evaluated, and of those in their coefficients.
%
% INPUT
% - PDE_in:     'pde_struct' object representing free terms;
%
% OUTPUT
% - vnames:     cellstr column of variable names, possibly repeated;

vnames = cell(0,1);
for ii=1:numel(PDE_in.free)
    if ~isfield(PDE_in.free{ii},'term'),    continue;   end
    for jj=1:numel(PDE_in.free{ii}.term)
        trm = PDE_in.free{ii}.term{jj};
        for f=1:numel(trm)
            if isfield(trm(f),'loc') && ~isempty(trm(f).loc)
                loc_f = polynomial(trm(f).loc);
                vnames = [vnames; loc_f.varname];                           %#ok<AGROW>
            end
            if isfield(trm(f),'C') && ~isempty(trm(f).C)
                C_f = polynomial(trm(f).C);
                vnames = [vnames; C_f.varname];                             %#ok<AGROW>
            end
        end
    end
end

end

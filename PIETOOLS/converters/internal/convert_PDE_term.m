function [Pblk,Px] = convert_PDE_term(term,Tops,PDE,opts)
% [Pblk,Px] = CONVERT_PDE_TERM(term,Tops,PDE,opts) converts a term of a
% PDE, a product of factors
%   term = fctr_1.*...*fctr_m,
%   fctr_k = int_{I{1}} ... int_{I{p}} C(s,theta)*D^{D}x_r(loc) dtheta
% (with x_r replaced by w_r or u_r for an input), to the corresponding
% term of the PIE, in the fundamental state xf and inputs wf, uf, using
% the operators D^{alpha}x_r = Tops.x*xf + Tops.w*wf + Tops.u*uf.
%
% INPUT
% - term:   m x 1 struct array of the m factors, with fields x, w or u,
%           D, loc, C, I (and delay), as the terms of a PDE;
% - Tops:   struct with fields
%           - x:    nx x 1 cell, with Tops.x{k}{alpha+1} the 1 x nV
%                   'copvar' (or 'tensopmat') xf --> D^{alpha}x_k, as
%                   Tops.x of convert_PIETOOLS_BVP with derivatives;
%           - w,u:  the same for wf and uf ([] if no such inputs);
%           - Iw,Iu:    (if the term involves inputs w_r or u_r) cells
%                   with Iw{r} the operator wf --> w_r, see info.Iw of
%                   convert_PIETOOLS_BVP ([] for inputs in the BCs);
% - PDE:    initialized 'pde_struct' object;
% - opts:   struct with fields nonlinear (true for a nonlinear PDE, the
%           only case in which a term can have several factors) and
%           outvars (cellstr, the variables of the equation of the term:
%           e.g. a constant term in an equation for x(s) depends on s);
%
% OUTPUT
% - Pblk:   struct with fields x, w, u, for a term with a single factor
%           a struct with fields B (cell with the 'sopvar' block acting
%           on each column of xf, wf or uf, [] for zero blocks) and cols
%           (the spaces of the columns), or [] if there is no
%           contribution;
% - Px:     for a nonlinear PDE, the contribution of xf to the term, as a
%           'polyopvar' object in the PIE state components x01_pie,...;
%           [] if there is none, and for a linear PDE;
%
% NOTES
% - Inputs can only appear linearly: a term with several factors cannot
%   involve inputs, or states that depend on inputs.
% - Delays are not supported. Integrals are supported over [a,s_k],
%   [s_k,b] or [a,b], not over a variable in which the state is
%   evaluated at a boundary, with C depending only on the s_j and on the
%   dummy variables of the integrated variables.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - convert_PDE_term
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
% DJ, 11/02/2026: Initial coding;
% DJ, 11/11/2026: Rename solve_PIETOOLS_BVP to convert_PIETOOLS_BVP;
% DJ (CC), 11/13/2026: Return the blocks first, and the polynomial only for
%                   a nonlinear PDE (of xf only); opts is a struct;

term = term(:);
nfactors = numel(term);
% The PIE is linear in wf and uf, so an input cannot be a factor in a
% product; inputs entering through a state (via the BCs) are caught below
if nfactors>1 && any(arrayfun(@(f) (isfield(f,'w') && ~isempty(f.w)) || (isfield(f,'u') && ~isempty(f.u)),term))
    error("convert_PDE_term:inputNonlinear","An input appears in a nonlinear or bilinear term, which is not supported.")
end
if nfactors>1 && ~opts.nonlinear
    error("A term of a linear PDE can only have a single factor; use the option 'nonlinear'.")
end

% % % Convert the factors to operators on the fundamental state and inputs
% Pfactors{k}.(x,w,u) holds the sopvar blocks of factor k acting on each
% column of xf, wf, uf. A state in a product that also depends on the
% inputs would make the inputs appear nonlinearly, which is not supported.
Pfactors = cell(nfactors,1);
for k=1:nfactors
    Pfactors{k} = convert_PDE_factor(term(k),Tops,PDE,reshape(opts.outvars,1,[]));
    if nfactors>1 && (~isempty(Pfactors{k}.w) || ~isempty(Pfactors{k}.u))
        error("convert_PDE_term:inputNonlinear","A state in a nonlinear term depends on the inputs, which then appear in a nonlinear or bilinear manner; this is not supported.")
    end
end
% Blocks are only returned for a single factor, as for a product they
% are not linear operators
Pblk = struct('x',[],'w',[],'u',[]);
if nfactors==1
    Pblk = Pfactors{1};
end

% % % For a nonlinear PDE, the product of the factors, as a polynomial in
% % % the components of xf (zero if a factor has no contribution)
Px = [];
if ~opts.nonlinear || any(cellfun(@(p) isempty(p.x),Pfactors))
    return
end
for k=1:nfactors
    % Declare the components X that the factor acts on (for nonzero blocks
    % only), named as the PIE state components (x01_pie,...), so that the
    % factor is the row B_cell acting on X, and the product a polynomial
    B_cell = Pfactors{k}.x.B;
    col_spaces = Pfactors{k}.x.cols;
    nz_idcs = find(~cellfun(@isempty,B_cell));
    state_names = arrayfun(@(j) sprintf('x%0*d_pie',numel(num2str(numel(B_cell))),j),nz_idcs,'UniformOutput',false)';
    [vars_X,uniq_idcs] = unique([col_spaces(nz_idcs).vars],'stable');
    vars_X = reshape(vars_X,1,[]);      dom_X = vertcat(col_spaces(nz_idcs).dom);      dom_X = dom_X(uniq_idcs,:);
    if isempty(vars_X)
        X = polyopvar(state_names,[col_spaces(nz_idcs).dim]');
    else
        X = polyopvar(state_names,[col_spaces(nz_idcs).dim]',vars_X,dom_X,cell2mat(arrayfun(@(j) ismember(vars_X,col_spaces(j).vars),nz_idcs','UniformOutput',false)));
    end
    F_factor = tensopmat(B_cell(nz_idcs))*X;
    if k==1
        Px = F_factor;
    else
        Px = Px.*F_factor;
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Pfactor = convert_PDE_factor(fctr,Tops,PDE,vars_eq)
% Pfactor = CONVERT_PDE_FACTOR(fctr,Tops,PDE,vars_eq) converts a single
% factor
%   fctr = int_{I} C(s,theta)*D^{D}x_r(loc) dtheta,
% or the same with x_r replaced by the input w_r or u_r, to operators on
% the fundamental state and inputs,
%   fctr = sum_j Pfactor.x.B{j}*xf_j + sum_j Pfactor.w.B{j}*wf_j
%               + sum_j Pfactor.u.B{j}*uf_j.
%
% INPUT
% - fctr:   'struct' defining a single factor, with fields x, w or u, D,
%           loc, C, I (and delay), see convert_PDE_term;
% - Tops:   struct with the operators mapping the fundamental state xf
%           and inputs wf, uf to the derivatives of the states, and to the
%           inputs of the PDE, see convert_PDE_term;
% - PDE:    initialized 'pde_struct' object;
% - vars_eq:    cellstr, the names of the variables of the equation in
%           which the factor appears;
%
% OUTPUT
% - Pfactor:  struct with fields x, w, u, each [] if the factor has no
%           contribution from xf, wf or uf, respectively, and otherwise
%           a struct with fields
%           - B:    1 x n cell of 'sopvar' blocks, with [] for the zero
%                   operators;
%           - cols: struct array with fields vars, dom, dim, the spaces
%                   of the n columns;
%
% NOTES
% - The factor is the composition of the operator K, defined by the
%   coefficient C and the integrals (see kernel_operator), with the
%   operators that give D^{D}x_r(loc), obtained from Tops.x, Tops.w and
%   Tops.u by evaluating the output variables with loc at their boundary
%   position (the 'subs' method of 'sopvar').
% Names of the primary variables s_k and the dummy variables theta_k
nvars = size(PDE.vars,1);
vars_primary = arrayfun(@(i) PDE.vars(i,1).varname{1},1:nvars,'UniformOutput',false);
vars_dummy = arrayfun(@(i) PDE.vars(i,2).varname{1},1:nvars,'UniformOutput',false);
% Fields of a factor are optional: has('D') etc. checks they are set
has = @(name) isfield(fctr,name) && ~isempty(fctr.(name));
if has('delay') && any(double(polynomial(fctr.delay(:))))
    error("Terms with delays are not supported.")
end

% % % Operators from the fundamental state and inputs to D^{D}x_r, or to
% % % the input w_r or u_r
% For a state, Tops_rows.x/.w/.u are the rows of Tops giving D^{D}x_r from
% xf, wf and uf. For an input, w_r or u_r is itself given by the row
% Tops.Iw{r} or Tops.Iu{r}, so only Tops_rows.w or Tops_rows.u is nonzero.
Tops_rows = struct('x',[],'w',[],'u',[]);
is_input = has('w') || has('u');
if is_input
    ty = 'w';       if has('u'),    ty = 'u';   end
    comp_idx = fctr.(ty);  obj_tab = PDE.([ty,'_tab']);
    fld_name = ['I',ty];
    if ~isfield(Tops,fld_name) || numel(Tops.(fld_name))<comp_idx || isempty(Tops.(fld_name){comp_idx})
        error("The operators Tops.%s{%d} mapping the fundamental inputs to the input of the PDE are not available; this input may be decomposed in a boundary condition, so that it cannot be used elsewhere.",fld_name,comp_idx)
    end
    Tops_rows.(ty) = Tops.(fld_name){comp_idx};
else
    comp_idx = 1;          if has('x'),    comp_idx = fctr.x;     end
    obj_tab = PDE.x_tab;
end
vars_idcs = find(logical(obj_tab(comp_idx,3:2+nvars)));             % variables of the component
vars_y = vars_primary(vars_idcs);
nvars_y = numel(vars_idcs);
if ~is_input
    % Derivative orders D in each variable of x_r; Tops.x{r}{D+1} maps
    % xf to D^{D}x_r (one cell index per variable), and similarly for w
    % and u
    D = zeros(1,nvars_y);    if has('D'),    D = fctr.D(:)';     end
    D_idcs = num2cell(D+1);
    if isempty(D_idcs)
        D_idcs = {1};
    end
    for obj_type={'x','w','u'}
        if isfield(Tops,obj_type{1}) && ~isempty(Tops.(obj_type{1}))
            R_cell = Tops.(obj_type{1}){comp_idx};
            if any(D+1>size(R_cell,1:max(nvars_y,1)))
                error("The order of the derivative of x_%d is higher than that to which the state is differentiable.",comp_idx)
            end
            Tops_rows.(obj_type{1}) = R_cell{D_idcs{:}};
            if isempty(Tops_rows.(obj_type{1})) && strcmp(obj_type{1},'x')
                error('convert_PDE_term:derivativeMissing',"The derivative D^{%s} of x_%d was not computed by convert_PIETOOLS_BVP: request it in the derivatives specified for convert_PIETOOLS_BVP.",mat2str(D),comp_idx) % DJ, 11/11/2026
            end
        end
    end
end

% % % The position at which the state is evaluated, and the integrals
% is_bdry(k): the state is evaluated at the boundary bdry_vals(k) in s_k,
% which is removed from the variables; Itype_cell{k} classifies the
% integral over s_k as 'as' ([a,s_k]), 'sb' ([s_k,b]) or 'ab' ([a,b])
is_bdry = false(1,nvars_y);       bdry_vals = zeros(1,nvars_y);
Itype_cell = repmat({''},1,nvars_y);
for k=1:nvars_y
    s_k = PDE.vars(vars_idcs(k),1);
    if has('loc') && ~is_input
        loc_k = polynomial(fctr.loc(k));
        is_bdry(k) = isdouble(loc_k);
        if is_bdry(k)
            bdry_vals(k) = double(loc_k);
        end
    end
    if has('I') && numel(fctr.I)>=k && ~isempty(fctr.I{k})
        if is_bdry(k)
            error("Integration over a variable at which the state is evaluated at a boundary is not supported.")
        end
        I_k = polynomial(fctr.I{k}(:));
        is_lower = isdouble(I_k(1)) && abs(double(I_k(1))-PDE.dom(vars_idcs(k),1))<1e-12;
        is_upper = isdouble(I_k(2)) && abs(double(I_k(2))-PDE.dom(vars_idcs(k),2))<1e-12;
        if is_lower && is_upper
            Itype_cell{k} = 'ab';
        elseif is_lower && isequal(I_k(2),s_k)
            Itype_cell{k} = 'as';
        elseif is_upper && isequal(I_k(1),s_k)
            Itype_cell{k} = 'sb';
        else
            error("Integrals over %s are only supported over [a,%s], [%s,b], or [a,b].",vars_primary{vars_idcs(k)},vars_primary{vars_idcs(k)},vars_primary{vars_idcs(k)})
        end
    end
end

% % % Evaluate at the boundary positions, and compose with the operator of
% % % the coefficient and the integrals
% K: operator for the coefficient C and the integrals. The factor is K*B_j,
% for the blocks B_j of Tops_rows.(ty_cur) evaluated at the boundary
% positions; blocks that vanish are dropped ([]), and col_spaces are the
% column spaces
C = polynomial(eye(obj_tab(comp_idx,2)));
if has('C')
    C = polynomial(fctr.C);
end
K = kernel_operator(C,vars_y(~is_bdry),Itype_cell(~is_bdry),vars_eq,vars_primary,vars_dummy,PDE.dom);
Pfactor = struct('x',[],'w',[],'u',[]);
for obj_type={'x','w','u'}
    obj = obj_type{1};
    P_row = Tops_rows.(obj);
    if isempty(P_row)
        continue
    end
    % Tops are 'copvar' (linear PDE) or 'tensopmat' (nonlinear PDE) rows
    if isa(P_row,'copvar')
        ncols = size(P_row.C,2);       B_cell = P_row.C(1,:);
        depmat_in = P_row.space_in;       vars_cols = P_row.vars;    dom_cols = P_row.dom;     dims_cols = P_row.dim_in;
    else
        ncols = size(P_row.ops,2);     B_cell = cell(1,ncols);
        for j=find(~cellfun(@isempty,P_row.ops(1,:)))
            B_cell{j} = P_row.ops{1,j}.ops{1};
        end
        depmat_in = P_row.depmat.in>0;    vars_cols = P_row.vars.in; dom_cols = P_row.dom.in;  dims_cols = P_row.dims.in;
    end
    col_spaces = arrayfun(@(j) struct('vars',{vars_cols(depmat_in(j,:))},'dom',dom_cols(depmat_in(j,:),:),'dim',dims_cols(j)),1:ncols);
    for j=find(~cellfun(@isempty,B_cell))
        if ~any(cellfun(@nnz,B_cell{j}.params(:)))
            B_cell{j} = [];  continue
        end
        if any(is_bdry)
            B_cell{j} = subs(B_cell{j},vars_y(is_bdry),bdry_vals(is_bdry));
        end
        B_cell{j} = K*B_cell{j};
        if ~any(cellfun(@nnz,B_cell{j}.params(:)))
            B_cell{j} = [];
        end
    end
    if ~all(cellfun(@isempty,B_cell))
        Pfactor.(obj) = struct('B',{B_cell},'cols',col_spaces);
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function K = kernel_operator(C,vars_y,Itype_cell,vars_eq,vars_primary,vars_dummy,dom)
% K = KERNEL_OPERATOR(C,VARS_Y,ITYPE_CELL,VARS_EQ,VARS_PRIMARY,VARS_DUMMY,
% DOM) returns the 'sopvar' object K representing
%   (K*y)(s) = int_{I_1} ... int_{I_p} C(s,theta)*y(theta_1,...) dtheta.
%
% INPUT
% - C:      'polynomial' of size m x n, the coefficient, which may depend
%           on s_j for all variables, and on the dummy variable theta_k of
%           the integrated variables vars_y{k} only;
% - vars_y: cellstr, the names of the variables of y, a function of
%           dimension n;
% - Itype_cell:    cell of the same size as vars_y, specifying the
%           integral over vars_y{k} as '' (no integration, so that y is
%           evaluated at s_k), 'as' (from a to s_k), 'sb' (from s_k to
%           b) or 'ab' (from a to b);
% - vars_eq:    cellstr, the names of the variables of the equation in
%           which the term appears, of which the output is a function;
% - vars_primary,vars_dummy:    cellstr, the names of the variables s_j of
%           the PDE, and of their dummy variables theta_j;
% - dom:    N x 2 array, the domain of the variables vars_primary;
%
% OUTPUT
% - K:      'sopvar' object mapping L2^n[S1,S3] --> L2^m[S2,S3], with S3
%           the variables of y that are in the output, S1 the variables
%           that are integrated out (over [a,b], with C independent of
%           s_k, and s_k not in vars_eq), and S2 the other variables that
%           C depends on, or that are in vars_eq (the output is a function
%           of these variables, constant if C is);
[~,idcs_primary] = ismember(vars_y,vars_primary);
vars_y_dummy = vars_dummy(idcs_primary);
is_int = ~cellfun(@isempty,Itype_cell);
% For an integral over [a,b] of a variable not in vars_eq, the integration
% variable may be written as the primary variable, e.g.
% int((1-s)*x,s,[0,1]); use the dummy instead
for k=find(strcmp(Itype_cell,'ab') & ~ismember(vars_y,vars_eq))
    C = polynomial(subs(C,pvar(vars_y{k}),pvar(vars_y_dummy{k})));
end
vars_C = C.varname(:)';
% C(s,theta) may only depend on s and on the dummy variables of integrals
if ~all(ismember(vars_C,[vars_primary,vars_y_dummy(is_int)]))
    error("The coefficient depends on variables that are not allowed: %s.",strjoin(vars_C(~ismember(vars_C,[vars_primary,vars_y_dummy(is_int)])),', '))
end
% Input and output variables, and their domains. A variable integrated
% over [a,b] disappears from the output unless it is a variable of the
% equation (vars_eq), in which case the result is constant in it.
vars_out_s =union(intersect(vars_C,vars_primary),vars_eq);          % variables of the output
is_S1 = strcmp(Itype_cell,'ab') & ~ismember(vars_y,vars_out_s);    % integrated out
S3 = sort(vars_y(~is_S1));
vars_out = sort([setdiff(vars_out_s,vars_y),S3]);
vars_in = sort(vars_y);
[~,idcs_out] = ismember(vars_out,vars_primary);     [~,idcs_in] = ismember(vars_in,vars_primary);
% Monomial bases and coefficients, with the dummy variable of each integral
% (a placeholder for variables without integral, in which the kernel is
% independent of the input variable, as the operator is a multiplier)
dummy_names = cell(1,numel(vars_in));
for k=1:numel(vars_in)
    is_yk = strcmp(vars_y,vars_in{k});
    if is_int(is_yk)
        dummy_names{k} = vars_y_dummy{is_yk};
    else
        dummy_names{k} = [vars_in{k},'_nodummy'];   % no integral
    end
end
quad_poly = quadPoly.polynomial2quadPoly(C,vars_out,dummy_names);
% Parameters, indexed by the slots of the variables in S3: multiplier (1),
% integral from a to s (2), from s to b (3), or both for a to b
slots_cell = cell(1,numel(S3));
slotmap_cell = {1,2,3,[2,3]};
for k=1:numel(S3)
    slots_cell{k} = slotmap_cell{strcmp(Itype_cell{strcmp(vars_y,S3{k})},{'','as','sb','ab'})};
end
% All slots selected by the integral types receive the same coefficients
% Q.C; the other slots are zero
size_params = [3*ones(1,numel(S3)),1];
params =repmat({sparse(size(quad_poly.C,1),size(quad_poly.C,2))},size_params);
if isempty(S3)
    params{1} = quad_poly.C;
else
    grid_cell = cell(1,numel(S3));
    [grid_cell{:}] = ndgrid(slots_cell{:});
    [params{sub2ind(size_params,grid_cell{:})}] = deal(quad_poly.C);
end
K = sopvar(params,struct('in',{vars_in},'out',{vars_out}),quad_poly.Zs,quad_poly.Zt,struct('in',dom(idcs_in,:),'out',dom(idcs_out,:)),size(C));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [Tops,Kops,Fop,info,Gops] = compose_Tmaps(maps,vars,dom,sizes,Delta,alphas,sources)
% [Tops,Kops,Fop,info,Gops]=COMPOSE_TMAPS(maps,vars,dom,sizes,Delta,
% alphas,sources) composes the solution operators of the 1D boundary value
% problems (BVPs) along the N directions of a group of states on the same
% variables, as a tensor product, after checking that they are compatible.
% This gives the operators of the ND fundamental state space Y, so that
% for all x_pde in the PDE domain D, with xf in Y and fundamental inputs
% wf and uf (the exogenous and actuator inputs in the BCs),
%       D^{alpha}x_pde_k = Tops.x{k}{alpha+1}*xf
%                          + Tops.w{k}{alpha+1}*wf
%                          + Tops.u{k}{alpha+1}*uf,
%       Kops.x*xf + Kops.w*wf + Kops.u*uf = 0.
%
% INPUT
% - maps:   1 x N struct array with fields Tops, Kops, Fop, info (the
%           outputs of solve_BVP_1D along direction i, involving all states
%           in the order of the group; empty for a direction in which no
%           state is differentiable), and optionally Qops (the fifth
%           output), forced (struct with logical arrays w, u, the
%           conditions along direction i with such inputs) and source
%           (struct with fields w, u, each with the fields pde, elem,
%           scale, dep of the field src of fund_inputs);
% - vars:   1 x N cellstr, the names of the variables;
% - dom:    N x 2 array, the intervals of the variables;
% - sizes:  1 x q array, the dimensions of the q states;
% - Delta:  q x N array, the orders of each state in each variable;
% - alphas: (optional) nal x N array of derivative orders alpha for which
%           to compute the operators in Tops, 'all' (default) for all
%           alpha<=Delta, or [] for none. The operator for alpha=0 (the
%           state itself) is always computed;
% - sources: (optional) true if all inputs in the BCs are determined by
%           other states (lifted states of lower dimension, see
%           solve_PIETOOLS_BVP), so that they satisfy the conditions along
%           rank-deficient directions, which are then passed through
%           (default false);
%
% OUTPUT
% All operators are 'copvar' objects, with the nV components xf_c of
% xf as columns.
% - Tops:   struct with fields x, w, u: q x 1 cells, with
%           x{k}{alpha_1+1,...,alpha_N+1} the 1 x nV operator
%           xf --> D^{alpha}x_pde_k for alpha<=Delta_k in alphas ([]
%           for the others), and w{k}, u{k} the same for wf, uf;
% - Kops:   struct with fields x, w, u, the operators of the constraints
%           Kops.x*xf + Kops.w*wf + Kops.u*uf = 0, with a row for each
%           constraint (no rows if there are none), and columns for the
%           components of xf, wf and uf. The inputs do not enter the
%           constraints, so Kops.w and Kops.u are zero;
% - Fop:    nV x q operator such that xf_c = sum_k Fop_{c,k}*
%           D^{Dorder_{c,k}}x_pde_k, with Dorder as in info.fund(c).Dorder;
% - info:   struct with fields
%           - vars:     the names of the variables;
%           - pass:     1 x N logical, true for the rank-deficient
%                       directions that are passed through (sources only);
%           - fund:     struct array describing each xf_c: vars, dom,
%                       dim, Dorder (q x N) and p (xf_c stacks the
%                       states for each auxiliary value if p=sum(sizes));
%           - inputs:   struct with fields w, u: struct arrays describing
%                       each component c of wf, uf: vars, dom, dim,
%                       pde, Dorder, loc, sel (component c is the
%                       derivative of order Dorder of the elements sel of
%                       the input pde of the group, evaluated at the lower
%                       boundary where loc is -1);
%           - rhs:      struct with fields w, u: for each direction with
%                       such inputs, the stacked right-hand sides g_i of
%                       its conditions (fields dir, vars, dom, nrow, rows,
%                       the rows with such inputs);
% - Gops:   struct with fields w, u: operators wf (uf) --> the
%           rows of g^w (g^u) with such inputs, one block row for each
%           element of info.rhs (fundamental inputs only, [] otherwise);
%
% NOTES
% - Let T_{1,i}, T_{0,i}, F_i, K_i, R_i be the inverse of the highest
%   derivative, the map from the auxiliary values, the auxiliary values,
%   the constraint and the derivative operators of the 1D BVP along s_i
%   (solve_BVP_1D). Let A be the directions with m_i=0 (admissible) and R
%   those with m_i>0 (rank deficient). Then xf consists of xf_I for
%   all subsets I of R, with D_i the highest derivative in s_i,
%       xf_I = (prod_{i in R\I} F_i)(prod_{i in I u A} D_i) x_pde,
%       x_pde = (prod_{j in A} T_{1,j})
%                       (sum_I prod_{i in R} T_i^{[I]} xf_I),
%   with T_i^{[I]}=T_{1,i} for i in I and T_{0,i} otherwise, and
%   K_i xf_I=0 for i in I. D^{alpha}x_pde follows by replacing each T by
%   its derivative.
% - The construction is valid if the BCs along each direction are order
%   compatible, the T_{1,j} of the admissible directions preserve the BCs
%   along all other directions, and all but one rank-deficient direction
%   are scalar (same conditions and order for every component); this is
%   tested at the start of the composition.
% - Inputs are supported in the conditions along admissible directions,
%   if all directions are admissible. Fundamental inputs are free functions
%   such that x_pde is in D, and every x_pde in D with compatible inputs is
%   of this form (see fund_inputs), so that the user need not enforce any
%   compatibility. The inputs w and u are decomposed separately.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - compose_Tmaps
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
% DJ (CC), 11/02/2026: Initial coding;

% % % Dimensions: N directions, q states of the group, ns scalar
% % % components. The 1D solutions are expressed in the states regrouped by
% % % order in direction i (maps(i).info.x_order); regroup_idcs{i} gives the
% % % position in x_pde of each regrouped component, to transform back to
% % % the ordering of the states.
vars = reshape(vars,1,[]);
N = numel(vars);        q = numel(sizes);       ns = sum(sizes);
cs_end = cumsum(sizes);     cs_start = cs_end-sizes+1;
d = max(Delta,[],1);                        % highest order in each direction
triv = arrayfun(@(mp) isempty(mp.info),maps);
m = zeros(1,N);                             % number of auxiliary values
regroup_idcs = cell(1,N);                            % position in u of the regrouped state
for i=find(~triv)
    m(i) = maps(i).info.m;
    regroup_idcs{i} = cell2mat(arrayfun(@(k) cs_start(k):cs_end(k),maps(i).info.x_order,'UniformOutput',false));
end


% % % Directions with inputs: the inputs need all directions to be admissible,
% % % except that, for sources, rank-deficient directions without inputs are
% % % passed through: the data of the inputs are functions of the variables,
% % % and need not be parametrized along these directions
if nargin<7
    sources = false;
end
pass = sources & (m>0);
itypes = {'w','u'};                         % types of inputs, w and u
forced = struct('w',false(1,N),'u',false(1,N));     % directions with inputs
if isfield(maps,'forced')
    for i=find(~arrayfun(@(mp) isempty(mp.forced),maps))
        for ty=itypes
            forced.(ty{1})(i) = any(maps(i).forced.(ty{1}));
        end
    end
end
any_forced = forced.w | forced.u;
if any(any_forced) && any(m>0) && (~sources || any(any_forced & m>0))
    error('compose_Tmaps:inputsRankDeficient',['The boundary conditions have inputs, but the directions ',strjoin(vars(m>0),', '),' are rank deficient; inputs in the boundary conditions are only supported if all directions are admissible.'])
end
% Rows without inputs of a type along a direction with inputs of that type:
% their right-hand side is zero, so they must be preserved too
rows_hom = repmat({{}},1,N);
for ty=itypes
    for i=find(forced.(ty{1}))
        rows_hom{i}{end+1} = find(~maps(i).forced.(ty{1}));
    end
end

% % % Check compatibility, and determine the roles of the directions:
% % % admissible (adm_idcs, m=0), rank deficient (rk_idcs, m>0), of which at
% % % most one non-scalar direction gen_idx (the general direction), and the
% % % scalar ones scalar_idcs
% Checks (see NOTES): order compatibility and scalar directions, and
% preservation of the BCs by the inverses T_{1,j}
scalar = false(1,N);

% % % Preservation of the BCs of direction i by T_{1,j} (tested below)
% % % needs the T_{1,j} in terms of matrix coefficients.
% Matrix coefficients (ns x ns, in the ordering of the states) of the
% kernel of the inverse T_1 of each admissible direction j, for each of
% the three partial integrals and each monomial in (s,theta)
T1_coefs = cell(1,N);
for j=find(~triv & m==0)
    T1_blk = maps(j).Tops.T1.C{1,1};
    nnZL_arr = prod(cellfun(@numel,T1_blk.ZL));    nnZR_arr = prod(cellfun(@numel,T1_blk.ZR));
    T1_coefs{j} = zeros(ns,ns,0);
    for iparam=1:numel(T1_blk.params)
        param = T1_blk.params{iparam};
        if nnz(param)==0
            continue
        end
        for izl=1:nnZL_arr
            for izr=1:nnZR_arr
                coef = full(param(izl:nnZL_arr:end,izr:nnZR_arr:end));
                if any(coef(:))
                    coef_state = zeros(ns);     coef_state(regroup_idcs{j},regroup_idcs{j}) = coef;
                    T1_coefs{j}(:,:,end+1) = coef_state;
                end
            end
        end
    end
end

% Order compatibility and scalar directions
Mmat_cell = cell(1,N);        Fcoef_cell = cell(1,N);
for i=find(~triv)
    info_i = maps(i).info;
    nrows = size(info_i.E,1);       d_i = numel(info_i.n)-1;
    F_flat = full(info_i.F.coefficient);
    if isempty(F_flat)
        F_arr = zeros(nrows,ns,0);
    else
        F_arr = reshape(F_flat.',nrows,ns,[]);
    end
    Mmat = [info_i.E,reshape(F_arr,nrows,[])];
    tol = 1e-10*max(1,norm(Mmat,'fro'));
    Mmat_cell{i} = Mmat;       Fcoef_cell{i} = F_arr;

    % Order compatibility: rank(M) = sum_C rank(M_C), rank(M_C) <= N_{i,C},
    % for M_C the columns of the states in a class C of states with the
    % same orders in all other directions
    if N==1
        class_id = ones(q,1);
    else
        [~,~,class_id] = unique(Delta(:,[1:i-1,i+1:N]),'rows');
    end
    col_class = class_id([info_i.lab_E,repmat(info_i.lab_F,1,size(F_arr,3))])';
    rank_sum = 0;
    for class_k=unique(class_id)'
        rank_class = rank(Mmat(:,col_class==class_k),tol);
        if rank_class>sum(sizes(class_id==class_k).*Delta(class_id==class_k,i)')
            error('compose_Tmaps:orderIncompatible',['The boundary conditions along ',vars{i},' are order incompatible.'])
        end
        rank_sum = rank_sum + rank_class;
    end
    if rank(Mmat,tol)~=rank_sum
        error('compose_Tmaps:orderIncompatible',['The boundary conditions along ',vars{i},' couple components that have different orders in another direction, which is not supported.'])
    end

    % Scalar direction: all components have the same order, and Mmat has
    % the row space of kron(e,I_ns), for e the rref of the first entries
    if all(Delta(:,i)==d_i)
        Mmat_rref = rref(Mmat,tol);
        Mmat_rref = Mmat_rref(any(abs(Mmat_rref)>tol,2),:);
        if mod(size(Mmat_rref,1),ns)==0
            e_scalar = Mmat_rref(1:ns:end,1:ns:end);
            scalar(i) = norm(Mmat_rref-kron(e_scalar,eye(ns)),'fro')<=10*tol;
        end
    end
end

% Preservation: T_{1,j}g satisfies the BCs of direction i for all g that
% do, iff, for each matrix coefficient M of T_{1,j}, the functionals Mmat
% with M acting on the columns lie in the row space of Mmat
for i=find(~triv)
    info_i = maps(i).info;
    nrows = size(info_i.E,1);       n = info_i.n;        d_i = numel(n)-1;
    Mmat = Mmat_cell{i};       F_arr = Fcoef_cell{i};
    tol = 1e-10*max(1,norm(Mmat,'fro'));
    offs = cumsum([0,n]);
    idx_E0 = cumsum([0,2*(1:d_i).*n(2:end)]);
    rank_Mmat = rank(Mmat,tol);
    for j=find(~triv & m==0 & (1:N)~=i)
        for k=1:size(T1_coefs{j},3)
            coef_reg = T1_coefs{j}(regroup_idcs{i},regroup_idcs{i},k);
            E_new = zeros(size(info_i.E));
            coef_diag = zeros(size(coef_reg));
            for jj=find(n>0)-1
                blk = offs(jj+1)+(1:n(jj+1));
                coef_diag(blk,blk) = coef_reg(blk,blk);
                if jj>0
                    cols = idx_E0(jj)+(1:2*jj*n(jj+1));
                    E_new(:,cols) = info_i.E(:,cols)*kron(eye(2*jj),coef_reg(blk,blk));
                end
            end
            F_new = zeros(size(F_arr));
            for r=1:size(F_arr,3)
                F_new(:,:,r) = F_arr(:,:,r)*coef_reg;
            end
            % M must not couple components of different orders in s_i
            if norm(coef_reg-coef_diag,'fro')>1e-10*max(1,norm(coef_reg,'fro')) || rank([Mmat;[E_new,reshape(F_new,nrows,[])]],tol)~=rank_Mmat
                if m(i)==0
                    error('compose_Tmaps:notConsistent',['The inverses of the (admissible) boundary conditions along ',vars{j},' and ',vars{i},' do not commute.'])
                else
                    error('compose_Tmaps:notPreserved',['The inverse of the boundary conditions along ',vars{j},' does not preserve the (rank-deficient) boundary conditions along ',vars{i},'.'])
                end
            end
            % The conditions without inputs must be preserved on their own,
            % so that the inputs along s_i are only constrained by their
            % own entries
            for rows_cell = rows_hom{i}
                rows_sel = rows_cell{1};
                if ~isempty(rows_sel) && rank([Mmat(rows_sel,:);[E_new(rows_sel,:),reshape(F_new(rows_sel,:,:),numel(rows_sel),[])]],tol)~=rank(Mmat(rows_sel,:),tol)
                    error('compose_Tmaps:inputRowsNotPreserved',['The inverse of the boundary conditions along ',vars{j},' does not preserve the conditions along ',vars{i},' that have no inputs, which is not supported.'])
                end
            end
        end
    end
end
adm_idcs = find(~triv & m==0);       rk_idcs = find(m>0);
nonscalar = rk_idcs(~scalar(rk_idcs));
if numel(nonscalar)>1
    error('compose_Tmaps:coupledRankDeficient',['The boundary conditions along ',strjoin(vars(nonscalar),', '),' are rank deficient, and couple components; at most one such direction is supported.'])
end
if isempty(nonscalar)
    gen_idx = rk_idcs(1:min(1,numel(rk_idcs)));   % general direction
else
    gen_idx = nonscalar;
end
scalar_idcs = rk_idcs(rk_idcs~=gen_idx);                             % scalar directions

% Inputs are only supported along admissible directions: the particular
% solution is built from T_{1,j} and Tp_i of these
if any(any_forced & (triv | m>0))
    error('compose_Tmaps:forcedRankDeficient',['The boundary conditions along ',strjoin(vars(any_forced & (triv | m>0)),', '),' involve inputs, but are rank deficient; inputs are only supported in the conditions along admissible directions.'])
end


% % % The derivative orders alpha (rows of alphas, with alpha=0 the first) for
% % % which to compute the operators, and the orders in each direction that
% % % this requires from the 1D operators
% alphas: the requested multi-indices, with alpha=0 first (all alpha<=d by
% default, with the first direction varying slowest)
if nargin<6 || strcmp(alphas,'all')
    alphas = zeros(1,0);
    for i=1:N
        alphas = [repmat(alphas,d(i)+1,1),kron((0:d(i))',ones(size(alphas,1),1))];
    end
else
    alphas_in = alphas;
    if ~isempty(alphas) && (size(alphas,2)~=N || any(alphas(:)<0) || any(alphas(:)~=round(alphas(:))) || any(any(alphas>d,2)))
        error('compose_Tmaps:alphas','The derivative orders must be nonnegative integers, N per row, not exceeding the orders of the states.')
    end
    alphas = zeros(1,N);
    if N>0 && ~isempty(alphas_in)
        [~,keep] = unique([alphas;alphas_in],'rows','stable');
        alphas = [alphas;alphas_in];       alphas = alphas(keep,:);
    end
end
nal = size(alphas,1);

% % % 1D operators in the ordering of the states (not regrouped by order),
% % % with zero rows for components of lower order than the derivative. Rp
% % % are the derivatives of the Taylor polynomial of the inputs
% ops_1D(i).R1{a+1}, R0{a+1}: a-th derivative of x_pde in terms of the
% data, with (R1) or without (R0) the auxiliary values; K, F: the
% constraint and the auxiliary-value operator of a rank-deficient
% direction.
ops_1D = struct('R1',cell(1,N),'R0',cell(1,N),'K',cell(1,N),'F',cell(1,N),'Rp',cell(1,N));
for i=find(~triv)
    n = maps(i).info.n;
    % Permutation back from the regrouped to the original ordering
    perm_back = zeros(1,ns);      perm_back(regroup_idcs{i}) = 1:ns;
    % Derivatives of order a of all components, in the order of the states,
    % for the orders a in the requested alphas
    for a=unique(alphas(:,i))'
        X = stack_rows(maps(i).Tops.R1,n,a,d(i),perm_back);
        ops_1D(i).R1{a+1} = X(:,perm_back);
        if m(i)>0
            ops_1D(i).R0{a+1} = stack_rows(maps(i).Tops.R0,n,a,d(i),perm_back);
        end
        if any_forced(i)
            ops_1D(i).Rp{a+1} = stack_rows(maps(i).Qops.Rp,n,a,d(i),perm_back);
        end
    end
    if m(i)>0
        ops_1D(i).K = maps(i).Kops.x.C{1,1}(:,perm_back);
        ops_1D(i).F = maps(i).Fop.C{1,1}(:,perm_back);
        % With one auxiliary value for each component of the general
        % direction, these are ordered as the regrouped components in the
        % 1D problem. xf is split in blocks per state, with the auxiliary
        % values in the order of the states (see solve_PIETOOLS_BVP), so
        % reorder them: the rows of F, and the columns of R0.
        if i==gen_idx && m(i)==ns && ~isequal(perm_back,1:ns)
            ops_1D(i).F = ops_1D(i).F(perm_back,:);
            for ia=find(~cellfun(@isempty,ops_1D(i).R0(:)))'
                ops_1D(i).R0{ia} = ops_1D(i).R0{ia}(:,perm_back);
            end
        end
    end
end

% Scalar directions: all components have the same BCs, so the 1D operators
% are X (x) I_ns, of which only the scalar X is used; the identity factor
% is reintroduced per component by kron when composing below
for i=scalar_idcs
    for R_name={'R1','R0'}
        for a=find(~cellfun(@isempty,ops_1D(i).(R_name{1})(:)))'
            ops_1D(i).(R_name{1}){a} = unkron(ops_1D(i).(R_name{1}){a},ns);
        end
    end
    % K only enters through K*xf=0, which an invertible left factor W
    % leaves unchanged. With mixed rows in the BCs (e.g. u1(0)+u2(0)=0,
    % u1(0)-u2(0)=0), K=P*(X (x) I_ns) for a mixing matrix P. Then row
    % reduce the coefficient blocks of K, and multiply K by the same W, to
    % get X (x) I_ns.
    try
        ops_1D(i).K = unkron(ops_1D(i).K,ns);
    catch ME
        K_op = ops_1D(i).K;
        nnZL_arr = prod(cellfun(@numel,K_op.ZL));       nnZR_arr = prod(cellfun(@numel,K_op.ZR));
        K_coefs = zeros(K_op.dims(1),0);
        for iparam=1:numel(K_op.params)
            param = K_op.params{iparam};
            if nnz(param)==0
                continue
            end
            for izl=1:nnZL_arr
                for izr=1:nnZR_arr
                    K_coefs = [K_coefs,full(param(izl:nnZL_arr:end,izr:nnZR_arr:end))];     %#ok<AGROW>
                end
            end
        end
        tol_K = 1e-10*max(1,norm(K_coefs,'fro'));
        if ~strcmp(ME.identifier,'sopvar:unkron:notKron') || isempty(K_coefs) || rank(K_coefs,tol_K)<size(K_coefs,1)
            rethrow(ME)
        end
        K_rref = rref([K_coefs,eye(size(K_coefs,1))],tol_K);
        ops_1D(i).K = unkron(K_rref(:,size(K_coefs,2)+1:end)*K_op,ns);
    end
    ops_1D(i).F = unkron(ops_1D(i).F,ns);
end


% % % Compose the operators for each component of the fundamental state,
% % % and each derivative alpha (row of alphas, with alpha=0 the first).
% % % There is a component xf_I for each subset I of the rank-deficient
% % % directions: in I, xf_I is defined by the highest derivative
% % % (T_{1,i}); outside I, by the auxiliary values (T_{0,i}), which are
% % % lost as variables of xf_I.
nrk = numel(rk_idcs);
ncomp = 2^nrk;
x_ops = cell(ncomp,nal);        F_comp = cell(ncomp,1);         K_blocks = cell(0,ncomp);
fund = struct('vars',cell(1,ncomp),'dom',[],'dim',[],'Dorder',[],'p',[]);
for c=1:ncomp
    % The subset I of rank-deficient directions (binary expansion of c-1),
    % the directions aux with auxiliary values, and the scalar ones among
    % these
    in_I = logical(bitget(c-1,1:nrk));
    I = rk_idcs(in_I);        aux = rk_idcs(~in_I);
    aux_scalar = aux(aux~=gen_idx);
    m_scalar = m(aux_scalar)/ns;                        % scalar auxiliary dimensions
    % Size p of the values along the general direction: all of the data if
    % gen_idx is in I, and its auxiliary values otherwise
    p = ns;
    if ~isempty(gen_idx) && ~any(I==gen_idx)
        p = m(gen_idx);
    end
    var_mask = true(1,N);      var_mask(aux) = false;     % variables of xf_I
    dim_I = prod(m_scalar)*p;                        % dimension of xf_I
    Dorder_I = Delta;         Dorder_I(:,aux) = 0;      % orders no longer a function of s_aux
    fund(c) = struct('vars',{vars(var_mask)},'dom',dom(var_mask,:),'dim',dim_I,'Dorder',Dorder_I,'p',p);

    % Maps xf_I --> D^{alpha}x_pde for the requested alpha, as stages
    % of 1D operators: the scalar directions in I, the scalar auxiliary
    % directions, the general direction, and the admissible directions,
    % chained by compose_sopvar
    for a=1:nal
        alpha_a = alphas(a,:);
        stages = struct('op',{},'var',{});
        curdim = dim_I;
        for s=I(ismember(I,scalar_idcs))
            stages(end+1) = struct('op',kron(ops_1D(s).R1{alpha_a(s)+1},speye(curdim)),'var',vars{s});  %#ok<AGROW>
        end
        for k=numel(aux_scalar):-1:1
            curdim = curdim/m_scalar(k);
            stages(end+1) = struct('op',kron(ops_1D(aux_scalar(k)).R0{alpha_a(aux_scalar(k))+1},speye(curdim)),'var',vars{aux_scalar(k)});  %#ok<AGROW>
        end
        if ~isempty(gen_idx)
            R_name = 'R0';
            if any(I==gen_idx)
                R_name = 'R1';
            end
            stages(end+1) = struct('op',ops_1D(gen_idx).(R_name){alpha_a(gen_idx)+1},'var',vars{gen_idx});                                 %#ok<AGROW>
        end
        for j=adm_idcs
            stages(end+1) = struct('op',ops_1D(j).R1{alpha_a(j)+1},'var',vars{j});                                      %#ok<AGROW>
        end
        x_ops{c,a} = compose_sopvar(stages,vars(var_mask),vars,dom,dim_I);
    end

    % Operator D^{Dorder}x_pde --> xf_I (the map Fop): the auxiliary
    % values of the general direction, then of the scalar directions, each
    % adding a factor to the left of the values
    stages = struct('op',{},'var',{});
    curdim = ns;
    if ~isempty(gen_idx) && ~any(I==gen_idx)
        stages(end+1) = struct('op',ops_1D(gen_idx).F,'var',vars{gen_idx});                                                  %#ok<AGROW>
        curdim = m(gen_idx);
    end
    for k=1:numel(aux_scalar)
        stages(end+1) = struct('op',kron(ops_1D(aux_scalar(k)).F,speye(curdim)),'var',vars{aux_scalar(k)});  %#ok<AGROW>
        curdim = m_scalar(k)*curdim;
    end
    F_comp{c} = compose_sopvar(stages,vars,vars,dom,ns);

    % Constraints K_i xf_I=0, for the general direction on the last
    % factor of the values, and for a scalar direction on all values.
    % K_blocks(r,c) is the r-th constraint on component c (empty if c has
    % none).
    for i=I
        if i==gen_idx
            X = kron(speye(prod(m_scalar)),ops_1D(i).K);
        else
            X = kron(ops_1D(i).K,speye(dim_I));
        end
        K_blocks(end+1,c) = {compose_sopvar(struct('op',X,'var',vars{i}),vars(var_mask),vars,dom,dim_I)};                       %#ok<AGROW>
    end
end


% % % Inputs of each type: operators of the fundamental inputs (fund_inputs),
% % % and rhs_info of the conditions they enter in; the conditions with such
% % % inputs are those with inputs of that type, the others homogeneous
% fund_inputs takes the structs group (the states and their layout in
% x_pde) and dirs (the 1D problem along each direction, with the entries
% of the boundary data), described in fund_inputs
group = struct('vars',{vars},'dom',{dom},'ns',{ns},'Delta',{Delta},'alphas',{alphas});
group.stpos = repelem((1:numel(sizes))',sizes(:));      % state of each component of u
group.cs_start = cs_start;
% dirs(i): the 1D BVP along direction i, with the entries of the boundary
% data ordered as the boundary values of solve_BVP_1D (for each order j,
% the coefficients 0,...,j-1 of its n(j+1) states)
dirs = repmat(struct('triv',false,'adm',false,'forced',false,'rows_in',[],'R1',[],'Rp',[],'Cw',[],'d',[],'pos',[],'ell',[],'src',[]),1,numel(maps));
for i=1:numel(maps)
    dirs(i).triv = triv(i) | pass(i);     dirs(i).adm = ~triv(i) & m(i)==0;
    dirs(i).R1 = ops_1D(i).R1;      dirs(i).Rp = ops_1D(i).Rp;
end
for i=find(any_forced)
    n = maps(i).info.n;     d_i = numel(n)-1;       n_offs = cumsum([0,n]);
    dirs(i).Cw = maps(i).Qops.Cw;
    dirs(i).d = d_i;
    dirs(i).pos = cell2mat(arrayfun(@(j) repmat(regroup_idcs{i}(n_offs(j+1)+(1:n(j+1)))',j,1),1:d_i,'UniformOutput',false)');
    dirs(i).ell = cell2mat(arrayfun(@(j) repelem((0:j-1)',n(j+1),1),1:d_i,'UniformOutput',false)');
end
rhs_info = struct('w',[],'u',[]);     inputs_info = rhs_info;        inputs_ops = rhs_info;        rhs_ops = rhs_info;
for ty=itypes
    obj_type = ty{1};
    rhs_info.(obj_type) = struct('dir',{},'vars',{},'dom',{},'nrow',{},'rows',{});
    for i=find(forced.(obj_type))
        rhs_info.(obj_type)(end+1) = struct('dir',i,'vars',{vars([1:i-1,i+1:N])},'dom',dom([1:i-1,i+1:N],:),'nrow',numel(maps(i).forced.(obj_type)),'rows',find(maps(i).forced.(obj_type)));
    end
    for i=find(any_forced)
        dirs(i).forced = forced.(obj_type)(i);     dirs(i).rows_in = maps(i).forced.(obj_type);
        dirs(i).src = maps(i).source.(obj_type);
    end
    rhs_ops.(obj_type) = cell(0,0);
    [inputs_ops.(obj_type),inputs_info.(obj_type),rhs_ops.(obj_type)] = fund_inputs(group,dirs);
end


% % % Declare the operators per state: select the rows of the stacked
% % % operators corresponding to each state k of the group, and collect the
% % % columns (components of xf, or of the inputs) in containers.
% % % R_state{k} then holds the derivatives of x_k, up to its own orders
% % % Delta(k,:).
% Containers for the operators of each state k, to be filled below
F_blocks = cell(ncomp,q);         R_state = cell(q,1);
R_inputs = struct('w',{cell(q,1)},'u',{cell(q,1)});
for k=1:q
    rows = cs_start(k):cs_end(k);
    select_rows = @(X) X;      select_cols = @(X) X;      % no indexing needed for one state
    if q>1
        select_rows = @(X) X(rows,':');            select_cols = @(X) X(':',rows);
    end
    F_blocks(:,k) = cellfun(select_cols,F_comp,'UniformOutput',false);
    R_state{k} = state_diffs(x_ops,Delta(k,:),alphas,select_rows);
    for ty=itypes
        obj_type = ty{1};
        if ~isempty(inputs_ops.(obj_type))
            R_inputs.(obj_type){k} = state_diffs(inputs_ops.(obj_type),Delta(k,:),alphas,select_rows);
        end
    end
end
% % % Assemble the outputs; the input operators are empty if there are no
% % % inputs of that type
Tops = struct('x',{R_state},'w',[],'u',[]);
Fop = copvar(F_blocks);
Gops = struct('w',[],'u',[]);
for ty=itypes
    obj_type = ty{1};
    if isempty(inputs_ops.(obj_type))
        continue
    end
    Tops.(obj_type) = R_inputs.(obj_type);
    % Map from the fundamental inputs to the right-hand sides of the
    % conditions (Gops), with zero blocks where a component does not enter
    for ii=1:numel(rhs_info.(obj_type))
        for c=1:size(inputs_ops.(obj_type),1)
            if isempty(rhs_ops.(obj_type){ii,c})        % direction in the set J of the component
                rhs_ops.(obj_type){ii,c} = mat2sopvar(zeros(numel(rhs_info.(obj_type)(ii).rows),inputs_info.(obj_type)(c).dim),struct('in',{inputs_info.(obj_type)(c).vars},'out',{rhs_info.(obj_type)(ii).vars}),struct('in',inputs_info.(obj_type)(c).dom,'out',rhs_info.(obj_type)(ii).dom));
            end
        end
    end
    Gops.(obj_type) = copvar(rhs_ops.(obj_type));
end
Kop = [];
if ~isempty(K_blocks)
    % Every column of the container needs a block to define its space
    blk_ref = K_blocks{1,find(~cellfun(@isempty,K_blocks(1,:)),1)};
    for c=find(all(cellfun(@isempty,K_blocks),1))
        K_blocks{1,c} = mat2sopvar(zeros(blk_ref.dims(1),fund(c).dim),struct('in',{fund(c).vars},'out',{blk_ref.vars.out}),struct('in',fund(c).dom,'out',blk_ref.dom.out));
    end
    Kop = copvar(K_blocks);
end
% The spaces of the components of xf, wf and uf
cols = struct('x',[],'w',[],'u',[]);
comps = struct('x',fund,'w',inputs_info.w,'u',inputs_info.u);
for ty={'x','w','u'}
    comps_ty = comps.(ty{1});
    if ~isstruct(comps_ty)
        comps_ty = struct('dim',{},'vars',{});
    end
    cols.(ty{1}) = struct('dims',reshape([comps_ty.dim],[],1),'depmat',false(numel(comps_ty),numel(vars)));
    for c=1:numel(comps_ty)
        cols.(ty{1}).depmat(c,:) = ismember(vars,comps_ty(c).vars);
    end
end
% Zero operators with the rows of Kop (none if there are no constraints)
empty_op = @(P,sp) zeros_copvar_sop(struct('out',P.dim_out,'in',sp.dims),struct('out',{arrayfun(@(i) P.vars(P.space_out(i,:)),1:numel(P.dim_out),'UniformOutput',false)},'in',{arrayfun(@(i) vars(sp.depmat(i,:)),1:numel(sp.dims),'UniformOutput',false)}),struct('vars',{vars},'dom',dom));
Kops = struct('x',Kop,'w',[],'u',[]);
if isempty(Kop)
    Kop = struct('dim_out',zeros(0,1),'vars',{cell(1,0)},'space_out',false(0,0));
    Kops.x = empty_op(Kop,cols.x);
end
Kops.w = empty_op(Kop,cols.w);
Kops.u = empty_op(Kop,cols.u);

% % % Information on the components of xf, and the inputs
info = struct('vars',{vars},'pass',pass,'fund',fund,'inputs',inputs_info,'rhs',rhs_info);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function R_stacked = stack_rows(Rops_1D,n,a,d_i,perm_back)
% R_stacked = STACK_ROWS(Rops_1D,n,a,d_i,perm_back) stacks the operators to
% the derivative of order a of the states of each order in direction i,
% into a single operator to the derivative of order a of all components.
%
% INPUT
% - Rops_1D:    (d_i+1) x (d_i+1) cell of 'copvar' objects, e.g. the
%           operators Tops.R1 of solve_BVP_1D, with Rops_1D{j+1,a+1}
%           mapping the 1D fundamental state to the derivative of order a
%           of the n(j+1) states of order j;
% - n:      1 x (d_i+1) array, with n(j+1) the number of states of order j;
% - a:      integer, the order of the derivative;
% - d_i:    integer, the highest order of the states in direction i;
% - perm_back:  1 x ns array, the position of each regrouped component in
%           the ordering of the states;
%
% OUTPUT
% - R_stacked:  'sopvar' object mapping the 1D fundamental state to the
%           derivative of order a of all ns components, in the ordering of
%           the states, with zeros for the states of order j<a, which have
%           no such derivative;
%
% NOTES
% - The 1D solution is organized by order j (n(j+1) states of order j), so
%   this gives the a-th derivative of all components in a single operator.
blk_ref = Rops_1D{d_i+1,1}.C{1,1};
blocks = cell(1,0);
for j=find(n>0)-1
    if j>=a
        blocks{end+1} = Rops_1D{j+1,a+1}.C{1,1};  %#ok<AGROW>
    else
        blocks{end+1} = mat2sopvar(zeros(n(j+1),blk_ref.dims(2)),blk_ref.vars,blk_ref.dom);  %#ok<AGROW>
    end
end
R_stacked = vertcat(blocks{:});
R_stacked = R_stacked(perm_back,:);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Rk = state_diffs(ops_stacked,Dk,alphas,select_rows)
% Rk = STATE_DIFFS(ops_stacked,Dk,alphas,select_rows) extracts the
% operators to the derivatives of a single state, up to its own orders,
% from those to the derivatives of the stacked states.
%
% INPUT
% - ops_stacked:    ncols x nal cell, with ops_stacked{c,a} the operator
%           from column c (of the fundamental state or inputs) to the
%           derivative of order alphas(a,:) of the stacked states;
% - Dk:    1 x N array, the orders up to which the state is differentiable
%           in each variable;
% - alphas: nal x N array, the orders of the derivatives in each variable
%           for the elements of ops_stacked;
% - select_rows:  function handle, selecting the rows of the state from an
%           operator on the stacked states;
%
% OUTPUT
% - Rk:    cell of size Dk+1, with Rk{alpha+1} a 1 x ncols 'copvar' object
%           mapping to D^{alpha}x for the state x, for all alpha<=Dk,
%           indexed linearly (first direction fastest);
ncols = size(ops_stacked,1);
Rk = cell([Dk+1,ones(1,max(0,2-numel(Dk)))]);
for a=find(all(alphas<=Dk,2))'
    blocks_row = cell(1,ncols);
    for c=1:ncols
        blocks_row{c} = select_rows(ops_stacked{c,a});
    end
    Rk{1+sum(alphas(a,:).*cumprod([1,Dk(1:end-1)+1]))} = copvar(blocks_row);
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

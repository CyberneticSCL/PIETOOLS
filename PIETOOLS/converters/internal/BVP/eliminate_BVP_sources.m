function [Rops,Gops,info] = eliminate_BVP_sources(PDE,src,Rops,Gops,info)
% [Rops,Gops,info] = ELIMINATE_BVP_SOURCES(PDE,src,Rops,Gops,info)
% eliminates the inputs introduced by lift_BVP_sources for the
% states of lower dimension (including ODE states) in the BCs, from the
% operators assembled in solve_PIETOOLS_BVP. A state x2 with the BC
% x2(s1,0)=x1(s1) is solved as x2 = T2*v2 + Tw2*wf_c, with wf_c the
% fundamental inputs derived from the BC data p=x1: derivatives of p, at
% the lower boundary in some variables. As p=x1 is itself a function of the
% fundamental state v1 (and of its inputs), so is wf_c,
%       wf_c = Qx*xf + Qw*wf + Qu*uf,
% in terms of the other fundamental inputs wf and uf, and
% x2 = (T2 + Tw2*Qx)*xf + ..., so that x2 depends on the fundamental
% state of x1, and wf_c is no longer an input of the PDE.
%
% INPUT
% - PDE:    original 'pde_struct' object, whose BCs involve the states of
%           lower dimension;
% - src:    struct array defining the new inputs of the BCs, from
%           lift_BVP_sources, the r-th input being w_{numel(PDE.w)+r};
% - Rops,Gops,info:    operators and info assembled in solve_PIETOOLS_BVP
%           for the PDE in which these states are replaced by the new
%           inputs;
%
% OUTPUT
% - Rops,Gops,info:    the same, for the states and inputs of PDE:
%           the columns of the operators of the exogenous inputs w that
%           represent the new inputs are replaced by operators acting on
%           the fundamental state and inputs of the groups of the states
%           of lower dimension, and are removed from Gops and info, with
%           the rows of Gops and info.rhs of the BCs that they enter in.
%
% NOTES
% - The fundamental input wf_c is the derivative D^{Dorder} of elements
%   sel of the new input p = sum C*D^{D}x_state(loc) (src), evaluated at
%   the lower boundary where loc is -1 (info.inputs). With the product
%   rule, it is a sum of terms C_k*D^{D+k}x_state(loc), which
%   convert_PDE_term evaluates using Rops of the state.
% - A group of BCs with new inputs may not also involve inputs of the PDE,
%   which could not be parameterized by free functions, error
%   solve_PIETOOLS_BVP:sourceAndInput (raised in solve_PIETOOLS_BVP).
% - The inputs are eliminated one by one, in an order in which each only
%   depends on the states of groups without new inputs left. A cyclic
%   dependence between groups gives the error
%   solve_PIETOOLS_BVP:cyclicSources.
% - The BCs involving the states of lower dimension are then checked: the
%   terms of their states must cancel for every xf, or the compatibility
%   (corner) conditions of the BCs would constrain xf, which is not
%   supported, error solve_PIETOOLS_BVP:sourceCompatibility. This is not
%   checked for groups with rank-deficient directions, whose BCs hold on
%   the constrained fundamental state only.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eliminate_BVP_sources
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

N = size(PDE.vars,1);       nw0 = numel(PDE.w);
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);
isvar = logical(PDE.x_tab(:,3:2+N));
nx = numel(PDE.x);

% The groups with sources: those with a BC that has a new input
groups_src = zeros(1,0);
for r=1:numel(src)
    groups_src(end+1) = find(arrayfun(@(g) ismember(src(r).BC,info.groups(g).BC_idcs),1:numel(info.groups)),1);  %#ok<AGROW>
end
groups_src = unique(groups_src);

% % % Eliminate the inputs of the sources one by one
while true
    % The columns of the exogenous inputs wf that are derived from the
    % sources
    inputs_w = info.inputs.w;        n_wf = numel(inputs_w);
    is_source = [inputs_w.pde]>nw0;
    if ~any(is_source)
        break
    end
    n_xf = size(Rops.x{1}{1}.C,2);      n_uf = 0;
    if ~isempty(Rops.u)
        n_uf = size(Rops.u{1}{1}.C,2);
    end

    % Find an input that only depends on states without sources left, and
    % eliminate it
    done = false;
    for in_idx=find(is_source)
        keep = setdiff(1:n_wf,in_idx);      % the inputs that remain
        % The maps from xf, wf and uf to the input wf_c (see
        % source_input_map), and, for the derivatives of the states in
        % directions that are passed through, to its derivatives (inc)
        in_c = inputs_w(in_idx);
        source = src(in_c.pde-nw0);
        ncols = struct('x',n_xf,'w',n_wf,'u',n_uf);
        input_map = source_input_map(PDE,source,in_c,Rops,zeros(1,N),ncols);
        maps_inc = containers.Map(mat2str(zeros(1,N)),input_map);
        if any(~cellfun(@isempty,input_map.w(is_source)))
            continue                                % depends on another source
        end
        % The maps to the states and their derivatives
        for k=1:nx
            % The directions of the group of the state that are passed
            % through (the data of the inputs are not parametrized there)
            vk = find(isvar(k,:));
            g_k = info.fund(find(arrayfun(@(f) ismember(k,f.x_idcs),info.fund),1)).group;
            pass_k = ismember(vars_pde,info.groups(g_k).pass);
            for lidx=find(~cellfun(@isempty,Rops.x{k}(:)))'
                ops = struct('x',Rops.x{k}{lidx},'w',Rops.w{k}{lidx},'u',[]);
                % The derivative of order alpha of the state in a direction
                % that is passed through is that of the data of the inputs
                inc = zeros(1,N);
                if any(pass_k(vk))
                    al = cell(1,numel(vk));     [al{:}] = ind2sub(size(Rops.x{k}),lidx);
                    inc(vk) = cell2mat(al)-1;       inc(~pass_k) = 0;
                end
                if ~isKey(maps_inc,mat2str(inc))
                    maps_inc(mat2str(inc)) = source_input_map(PDE,source,in_c,Rops,inc,ncols);
                end
                input_map = maps_inc(mat2str(inc));
                if ~isempty(Rops.u)
                    ops.u = Rops.u{k}{lidx};
                end
                % Substitute wf_c = Q*(xf,wf,uf) into each row
                % of ops, so that ops.x*xf + ops.w*wf + ops.u*uf
                % becomes (ops.x+ops.w(:,c)*Qx)*xf + (ops.w(:,others)
                % +ops.w(:,c)*Qw)*wf + (ops.u+ops.w(:,c)*Qu)*uf.
                % Each row needs a block to define its space, so rows left
                % with only zero blocks get an explicit one
                blocks = struct('x',[],'w',{ops.w.C},'u',[]);
                if ~isempty(ops.x)
                    blocks.x = ops.x.C;
                end
                if ~isempty(ops.u)
                    blocks.u = ops.u.C;
                end
                nrows = size(ops.w.C,1);
                ref_blks = cell(nrows,1);                           % spaces of the rows, for zero blocks
                for i=1:nrows
                    ref_blks{i} = ops.w.C{i,in_idx};
                    if isempty(ref_blks{i})
                        continue
                    end
                    for ty={'x','w','u'}
                        obj_type = ty{1};
                        for j=find(~cellfun(@isempty,input_map.(obj_type)))
                            if isempty(blocks.(obj_type))
                                error("solve_PIETOOLS_BVP:sourceInputs","The inputs of a state of lower dimension in a boundary condition have no counterpart in the operators.")
                            end
                            blk_new = ref_blks{i}*input_map.(obj_type){j};
                            if isempty(blocks.(obj_type){i,j})
                                blocks.(obj_type){i,j} = blk_new;
                            else
                                blocks.(obj_type){i,j} = blocks.(obj_type){i,j}+blk_new;
                            end
                        end
                    end
                end
                if ~isempty(ops.x)
                    ops.x = copvar(blocks.x);
                end
                if ~isempty(ops.u)
                    ops.u = copvar(blocks.u);
                end
                % Remove the column c from the operator of wf
                if isempty(keep)
                    ops.w = [];
                else
                    W_keep = blocks.w(:,keep);
                    for i=find(all(cellfun(@isempty,W_keep),2))'
                        W_keep{i,1} = mat2sopvar(zeros(ref_blks{i}.dims(1),inputs_w(keep(1)).dim),struct('in',{inputs_w(keep(1)).vars},'out',{ref_blks{i}.vars.out}),struct('in',inputs_w(keep(1)).dom,'out',ref_blks{i}.dom.out));
                    end
                    ops.w = copvar(W_keep);
                end
                Rops.x{k}{lidx} = ops.x;       Rops.w{k}{lidx} = ops.w;
                if ~isempty(Rops.u)
                    Rops.u{k}{lidx} = ops.u;
                end
            end
        end
        % Remove the input from the operators and the info
        if isempty(keep)
            Rops.w = [];
        end
        if ~isempty(Gops.w)
            if isempty(keep)
                Gops.w = [];
            else
                Gops.w = copvar(Gops.w.C(:,keep));
            end
        end
        info.inputs.w(in_idx) = [];
        done = true;
        break
    end
    if ~done
        error("solve_PIETOOLS_BVP:cyclicSources","The states of lower dimension in the boundary conditions depend on each other cyclically.")
    end
end

% The BCs of the groups with sources are not conditions with inputs of the
% PDE: remove their right-hand sides from the info and from Gops
is_src_rhs = ismember([info.rhs.w.group],groups_src);
if any(is_src_rhs)
    info.rhs.w(is_src_rhs) = [];
    if ~isempty(Gops.w)
        if all(is_src_rhs)
            Gops.w = [];
        else
            Gops.w = copvar(Gops.w.C(~is_src_rhs,:));
        end
    end
end
% % % Check that the BCs of the groups with sources hold for every xf:
% % % the terms of the states must cancel, as operators on xf.
% % % Otherwise, the compatibility conditions of the BCs, e.g. x1(0)=0 for
% % % x2(0,s2)=0 and x2(s1,0)=x1(s1), constrain xf, and constraints on
% % % xf are not generated.
amax = @(P) max([0;cellfun(@(p) max([0;abs(nonzeros(p))]),P.params(:))]);
% Per group with sources, BC terms of its states must cancel on xf
for g=unique(groups_src)
    if nnz([info.fund.group]==g)>1                       % rank-deficient directions
        continue
    end
    for k=info.groups(g).BC_idcs(:)'
        terms = PDE.BC{k}.term;
        if ~all(cellfun(@(t) isfield(t,'x'),terms))
            continue
        end   % has inputs
        vars_out = arrayfun(@(m) PDE.BC{k}.vars(m).varname{1},1:numel(PDE.BC{k}.vars),'UniformOutput',false);
        B_sum = cell(1,size(Rops.x{1}{1}.C,2));
        for j=1:numel(terms)
            [~,Pblk] = convert_PDE_term(terms{j},Rops,PDE,struct('nonlinear',false,'outvars',{vars_out}));
            if isempty(Pblk.x)
                continue
            end
            for col_idx=find(~cellfun(@isempty,Pblk.x.B))
                if isempty(B_sum{col_idx})
                    B_sum{col_idx} = Pblk.x.B{col_idx};
                else
                    B_sum{col_idx} = B_sum{col_idx}+Pblk.x.B{col_idx};
                end
            end
        end
        if any(cellfun(@(P) ~isempty(P) && amax(P)>1e-9,B_sum))
            error("solve_PIETOOLS_BVP:sourceCompatibility","The boundary conditions %s are not satisfied by the states coupled to states of lower dimension: these states must satisfy the compatibility (corner) conditions of the boundary conditions they enter, e.g. x1(0)=0 if x2(0,s2)=0 and x2(s1,0)=x1(s1). Constraints on the fundamental state are not generated.",mat2str(k))
        end
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function input_map = source_input_map(PDE,source,in_c,Rops,inc,ncols)
% input_map = SOURCE_INPUT_MAP(PDE,SOURCE,IN_C,ROPS,INC,NCOLS) gives the
% maps from xf, wf and uf to an input wf_c derived from a
% source, using the maps of the states in the source.
%
% INPUT
% - PDE:    'pde_struct' object of the original PDE;
% - source: element of src, see lift_BVP_sources;
% - in_c:   element of info.inputs.w, the input wf_c derived from the
%           source, see solve_PIETOOLS_BVP;
% - Rops:   struct with the maps to the states and their derivatives;
% - inc:    1 x N array, increase of the orders of the derivative of
%           wf_c in each variable, for the derivatives of the states in
%           directions that are passed through;
% - ncols:  struct with fields x, w, u, the numbers of columns of the
%           operators in Rops;
%
% OUTPUT
% - input_map:   struct with fields x, w, u, cells with the blocks of the
%           operator xf, wf or uf --> D^{Dorder+inc}p, with
%           p = sum C*D^{D}x_state(loc) the source, selecting the elements
%           sel, evaluated at the lower boundary in the variables where
%           loc is -1;
%
% NOTES
% - With the product rule, the derivative of order Dorder of C*D^{D}x is
%   sum_k binom(Dorder,k)*(D^{Dorder-k}C)*(D^{D+k}x), with k nonzero only
%   in the variables of the state in which it is not evaluated at a
%   boundary. The terms are evaluated by convert_PDE_term.
N = size(PDE.vars,1);
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);
isvar = logical(PDE.x_tab(:,3:2+N));
Dorder = in_c.Dorder+inc;
input_map = struct('x',{cell(1,ncols.x)},'w',{cell(1,ncols.w)},'u',{cell(1,ncols.u)});
for jt=1:numel(source.terms)
    term = source.terms(jt);
    vidx = find(isvar(term.state,:));           nvars_src = numel(vidx);
    free = arrayfun(@(m) ~isdouble(term.loc(m)),1:nvars_src);
    vfree = reshape(vidx(free),1,[]);
    if any(in_c.loc(vfree)==0 & ~ismember(vars_pde(vfree),in_c.vars))
        error("solve_PIETOOLS_BVP:sourceVars","A state of lower dimension in a boundary condition depends on a variable that the other states do not depend on.")
    end
    % The multi-indices k <= Dorder in the variables of the state
    nfree = numel(vfree);       k_arr = zeros(1,0);
    if nfree>0
        rng_k = arrayfun(@(q) 0:q,Dorder(vfree),'UniformOutput',false);
        k_grid = cell(1,nfree);     [k_grid{:}] = ndgrid(rng_k{:});
        k_arr = cell2mat(cellfun(@(g) g(:),k_grid,'UniformOutput',false));
    end
    for ik=1:size(k_arr,1)
        k_vec = zeros(1,nvars_src);       k_vec(free) = k_arr(ik,:);
        binom_coef = prod(arrayfun(@nchoosek,Dorder(vfree),k_arr(ik,:)));
        diff_orders = Dorder;          diff_orders(vfree) = diff_orders(vfree)-k_arr(ik,:);
        % Derivative of the coefficient, of the selected elements
        C_diff = term.C(in_c.sel,:);
        for q=find(diff_orders>0)
            C_diff = polynomial(diff(C_diff,pvar(vars_pde{q}),diff_orders(q)));
        end
        % Evaluation at the lower boundary of the coefficient and state in
        % the variables where wf_c is not a function
        loc_new = term.loc;
        for q=find(in_c.loc==-1)
            C_diff = polynomial(subs(C_diff,pvar(vars_pde{q}),PDE.dom(q,1)));
            m = find(vidx==q);
            if ~isempty(m) && free(m)
                loc_new(m) = polynomial(PDE.dom(q,1));
            end
        end
        if isempty(C_diff.coeff) || ~any(any(C_diff.coeff))
            continue
        end
        % The term of the state with the derivatives
        term_k = struct('x',term.state,'D',term.D+k_vec,'loc',loc_new,'C',binom_coef*C_diff,'I',{cell(nvars_src,1)},'delay',0);
        [~,Pblk] = convert_PDE_term(term_k,Rops,PDE,struct('nonlinear',false,'outvars',{reshape(in_c.vars,1,[])}));
        for ty={'x','w','u'}
            obj_type = ty{1};
            if isempty(Pblk.(obj_type))
                continue
            end
            for j=find(~cellfun(@isempty,Pblk.(obj_type).B))
                if isempty(input_map.(obj_type){j})
                    input_map.(obj_type){j} = Pblk.(obj_type).B{j};
                else
                    input_map.(obj_type){j} = input_map.(obj_type){j}+Pblk.(obj_type).B{j};
                end
            end
        end
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
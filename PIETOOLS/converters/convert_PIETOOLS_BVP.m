function [Tops,Kops,Dmap,info] = convert_PIETOOLS_BVP(PDE,diffs,opts)
% [Tops,Kops,Dmap,info] = CONVERT_PIETOOLS_BVP(PDE,diffs,opts) computes
% operators parameterizing all states of a PDE that satisfy its boundary
% conditions (BCs), in terms of the fundamental state xf (free of BCs,
% up to constraints Kops.x*xf=0) and fundamental inputs wf, uf:
%       x_pde = Tops.x*xf + Tops.w*wf + Tops.u*uf,
% and, if requested, the same for derivatives of the states. Conversely,
% Dmap defines xf, wf and uf in terms of the PDE. The dynamics are
% ignored.
%
% INPUT
% - PDE:    'pde_struct' object in the terms format;
% - diffs:  (optional) the derivatives of the states to compute: [] (the
%           default) for none, 'all' for all up to the order of each
%           state, 'terms' for those in the terms of the equations and
%           outputs (as convert_PDE_term needs), or an nx x 1 cell, with
%           diffs{k} 'all' or an array with a row of derivative orders
%           (in the variables of state k) for each derivative;
% - opts:   (optional) char or cellstr, with 'no_decouple' (solve all
%           BCs as one group), 'all_inputs' (also include the inputs that
%           are not in the BCs), 'blocks' (split xf per state, see
%           OUTPUT) and/or 'nonlinear' (return 'tensopmat' objects).
%           diffs and opts may be omitted or [], or given in either
%           order: a single input is diffs if it is 'all', 'terms' or
%           not text, and opts otherwise;
%
% OUTPUT
% The operators are 'copvar' objects, or 'tensopmat' objects for a
% nonlinear PDE or with 'nonlinear'. The components of xf (the columns of
% Tops and Kops, the rows of Dmap) are those of info.fund, which may stack
% the states of a group, or, with 'blocks' or 'tensopmat', a block for
% each state in each of these (info.blocks).
% - Tops:   struct with fields x, w, u: the maps from xf, wf and uf to
%           x_pde (zero if no such inputs are in the BCs). With diffs,
%           nx x 1 cells instead, with x{k}{alpha+1} (an index for each
%           variable of state k) the map to D^{alpha}x_pde_k for the
%           requested alpha (and alpha=0), and [] for the others;
% - Kops:   struct with fields x, w, u, the constraints
%           Kops.x*xf + Kops.w*wf + Kops.u*uf = 0, with a row for each
%           constraint (none if there are none). The inputs do not enter
%           the constraints, so Kops.w and Kops.u are zero;
% - Dmap:   struct with fields x, w, u, defining each component c of xf,
%           wf, uf as PIE.pde2pie does: the PDE state or input idx(c),
%           differentiated (orders diff(c,:)), evaluated at the lower
%           boundary where loc(c,:)=-1, with the operator C applied. For
%           x, idx(c)=0 if xf_c stacks several states, diff is
%           ncomp x nx x N and xf_c = sum_k C_{c,k}*D^{diff(c,k,:)}x_k;
%           for w and u, C selects the elements of the input;
% - info:   struct with fields
%           - space:    the spaces of xf, wf and uf: vars, dom, and x, w,
%                       u with fields dims and depmat (true where a
%                       component depends on a variable);
%           - groups:   for each group of coupled BCs, its BCs (BC_idcs),
%                       variables (vars), and the directions passed
%                       through (pass, see compose_Tmaps);
%           - fund:     for each xf_c: group, x_idcs, vars, dom, dim,
%                       Dorder and p, see compose_Tmaps;
%           - blocks:   for each block of xf: fund (the index c), x (the
%                       state), comp (its entries of xf_c), dim, vars,
%                       dom and Dorder;
%           - inputs:   struct with fields w, u, for each component of wf
%                       (uf): the PDE input pde, derivative Dorder,
%                       location loc and elements sel that it is (as in
%                       Dmap), and its group (0 if not in the BCs);
%           - input_maps:   struct with fields w, u, false if Tops.w
%                       (Tops.u) is zero as no such inputs are in the BCs;
%           - Iw, Iu:   with 'all_inputs', cells with Iw{r} the map
%                       wf --> w_r for an input that is its own
%                       fundamental input ([] if it is decomposed);
%           - surplus:  true if there are more BCs than the orders of the
%                       states require (see NOTES), so that the PDE may
%                       be overdetermined;
%
% NOTES
% - The BCs are split into groups sharing no state (decouple_BVPs), and
%   into 1D problems along each direction (split_BVPs), which are solved
%   (compute_Tmap_1D) and composed (compose_Tmaps, see there for the BCs
%   that are supported). More BCs along a direction than the order of
%   the states give constraints Kops.x*xf=0.
% - BCs that fix several variables in each term of a state with
%   variables (e.g. x(0,0)=0, or an integral over the full domain) are
%   imposed on the solution of the other BCs. A scalar BC that fixes an
%   auxiliary value without variables (as for periodic BCs) eliminates
%   it from xf; the others give constraints Kops.x*xf=0. Such BCs cannot
%   involve inputs, nor states with inputs in their BCs.
% - States of lower dimension (e.g. ODE states, or boundary values of
%   states on other variables, as in x4(a)=x3(d)) in the BCs of other
%   states act as inputs of these BCs: they are lifted to inputs
%   (lift_BVP_sources), and eliminated after solving
%   (eliminate_BVP_sources). They must satisfy the compatibility
%   conditions of the BCs they enter.
% - Inputs w, u in the BCs are decomposed into free fundamental inputs
%   (derivatives of the PDE inputs, see fund_inputs), so that the BC data
%   are compatible by construction. This needs all directions of the
%   group to be admissible (e.g. not periodic), and each BC row with an
%   input to have a single element of a single input as right-hand side.
%   In a group on several variables, BCs with inputs cannot involve
%   states of lower dimension.
% - A delayed input in a BC is treated as an input (its value is the BC
%   data); the delay is not represented in Dmap.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - convert_PIETOOLS_BVP
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
% DJ (CC), 11/08/2026: Lift the inputs in a BC with a state of lower dimension
%                   into its source, as fundamental inputs, and map each
%                   input that is its own fundamental input to the PDE input
%                   (Iin), also if in the BCs; add info.input_maps;
% DJ, 11/11/2026: Rename from solve_PIETOOLS_BVP;
% DJ (CC), 11/13/2026: Solve the 1D problems with compute_Tmap_1D, and
%                   simplify the code; impose BCs that fix several
%                   variables at once; add info.surplus;

% % % Check the inputs
msg = "Options should be specified as 'no_decouple', 'all_inputs', 'blocks' or 'nonlinear'.";
% With a single optional input, or in the order (opts,diffs), the inputs
% are told apart as options are text other than 'all' and 'terms'
if nargin<3
    opts = [];      
end
if nargin<2    
    diffs = [];     
end
is_text = @(a) ischar(a) || isstring(a) || iscellstr(a);
is_diffs = @(a) ~isempty(a) && (~is_text(a) || any(ismember(lower(cellstr(a)),{'all','terms'})));
if ~isempty(diffs) && ~is_diffs(diffs) && (isempty(opts) || is_diffs(opts))
    [diffs,opts] = deal(opts,diffs);
end
diffs = convertStringsToChars(diffs);
if isempty(opts)
    opts = '';
end
if ~is_text(opts)
    error(msg)
end
opts = lower(cellstr(opts));
opts = opts(~cellfun(@isempty,opts));
if ~all(ismember(opts,{'no_decouple','all_inputs','blocks','nonlinear'}))
    error(msg)
end
is_an_opt = @(o) any(strcmp(opts,o));
PDE = pde_struct(PDE);
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
nx = numel(PDE.x);
if nx==0
    error("The PDE does not involve any states.")
end
N = size(PDE.vars,1);
isvar = logical(PDE.x_tab(:,3:2+N));
Delta = PDE.x_tab(:,2+N+(1:N));

% % % The derivatives of the states to compute: none (only the maps to the
% % % states) by default, all, those in the terms of the PDE, or specified
% diffs_k{k}: 'all', or an array with a row of derivative orders for each
% derivative of state k (in its own variables); alpha=0 is always included
want_diffs = ~isempty(diffs);
diffs_k = cell(nx,1);
if want_diffs
    if ischar(diffs) && strcmpi(diffs,'all')
        diffs_k(:) = {'all'};
    elseif ischar(diffs) && strcmpi(diffs,'terms')
        eq_terms = cellfun(@(e) e.term(:),[PDE.x(:);PDE.z(:);PDE.y(:)],'UniformOutput',false);
        diffs_k = term_diffs(PDE,diffs_k,vertcat(eq_terms{:}));
    elseif iscell(diffs) && numel(diffs)==nx
        diffs_k = diffs(:);
    else
        error('convert_PIETOOLS_BVP:diffs',"The derivatives should be specified as 'all', 'terms', or an nx x 1 cell.")
    end
    for k=1:nx
        Dk = Delta(k,isvar(k,:));
        if ischar(diffs_k{k})
            bad = ~strcmpi(diffs_k{k},'all');
        else
            bad = ~isempty(diffs_k{k}) && (~isnumeric(diffs_k{k}) || size(diffs_k{k},2)~=numel(Dk) || any(diffs_k{k}(:)<0) || any(diffs_k{k}(:)~=round(diffs_k{k}(:))) || any(any(diffs_k{k}>Dk,2)));
        end
        if bad
            error('convert_PIETOOLS_BVP:diffs','The derivatives of state %d should be ''all'', or nonnegative integer orders, one for each of its variables, not exceeding the orders of the state.',k)
        end
    end
end


% % % Set aside the BCs that fix several variables at once in each term of
% % % a state with variables (by evaluating at a boundary, or integrating
% % % over the full domain), to impose them on the solution of the others
multi = false(numel(PDE.BC),1);
for k=1:numel(PDE.BC)
    nfix = zeros(1,0);
    for j=1:numel(PDE.BC{k}.term)
        tm = PDE.BC{k}.term{j};
        if ~isfield(tm,'x') || ~any(isvar(tm.x,:))
            continue
        end
        vidx = find(isvar(tm.x,:));     loc = polynomial(tm.loc);
        is_fixed = false(1,numel(vidx));
        for m=1:numel(vidx)
            I_m = [];
            if isfield(tm,'I') && numel(tm.I)>=m
                I_m = polynomial(tm.I{m});
            end
            is_fixed(m) = isdouble(loc(m)) || (numel(I_m)==2 && isdouble(I_m) && isequal(sort(double(I_m(:)))',PDE.dom(vidx(m),:)));
        end
        nfix(end+1) = nnz(is_fixed);                                        %#ok<AGROW>
    end
    multi(k) = ~isempty(nfix) && all(nfix>1);
end
msg_multi = "Boundary conditions that fix several variables at once are not supported for states with inputs in their boundary conditions.";
if any(cellfun(@(b) any(cellfun(@(tm) any(isfield(tm,{'w','u'})),b.term)),PDE.BC(multi)))
    error('convert_PIETOOLS_BVP:multiInput',msg_multi)
end
PDE_full = PDE;
PDE.BC(multi) = [];     PDE.BC_tab(multi,:) = [];

% % % Split the BCs into groups sharing no states
% Terms of states of lower dimension (e.g. ODE states) are first lifted to
% inputs, eliminated after solving, so that each group involves states on
% the same variables. Groups are solved independently, keeping the cost
% small; without decoupling, there is a single group.
[PDE_lift,src] = lift_BVP_sources(PDE);
nw0 = numel(PDE.w);
if is_an_opt('no_decouple')
    if size(unique(isvar,'rows'),1)>1
        error('convert_PIETOOLS_BVP:mixedDimension',"The option 'no_decouple' requires all states to depend on the same variables.")
    end
    PDEs = {PDE_lift};
    dec = struct('BC_idcs',{{1:numel(PDE.BC)}},'x_idcs',{{1:nx}},'w_idcs',{{1:numel(PDE_lift.w)}},'u_idcs',{{1:numel(PDE.u)}});
else
    [PDEs,dec] = decouple_BVPs(PDE_lift);
end
% Eliminating the sources needs all derivatives of the lifted states, and
% those in the BCs of the states they enter (see eliminate_BVP_sources);
% imposing the BCs set aside needs those in these BCs
diffs_comp = diffs_k;
multi_terms = cellfun(@(b) b.term(:),PDE_full.BC(multi),'UniformOutput',false);
diffs_comp = term_diffs(PDE,diffs_comp,vertcat(cell(0,1),multi_terms{:}));
for r=1:numel(src)
    for jt=1:numel(src(r).terms)
        diffs_comp{src(r).terms(jt).state} = 'all';
    end
end
for g=1:numel(PDEs)
    if any(dec.w_idcs{g}>nw0)
        BC_terms = cellfun(@(b) b.term(:),PDE.BC(dec.BC_idcs{g}(:)),'UniformOutput',false);
        BC_terms = vertcat(BC_terms{:});
        % The data of the inputs of lifted states are determined by other
        % states, so that no compatibility of the BCs with inputs of the
        % PDE can be ensured, which therefore cannot be mixed in. This
        % only holds for groups of several variables (corners): in a group
        % of one variable the inputs left in the BCs (not lifted into a
        % source) are free. The inputs of the sources are not inputs of
        % the PDE.
        BC_terms_l = cellfun(@(b) b.term(:),PDE_lift.BC(dec.BC_idcs{g}(:)),'UniformOutput',false);
        BC_terms_l = vertcat(BC_terms_l{:});
        has_input = cellfun(@(tm) any(arrayfun(@(f) (isfield(f,'w') && ~isempty(f.w) && f.w<=nw0) || (isfield(f,'u') && ~isempty(f.u)),tm)),BC_terms_l);
        if any(has_input) && size(PDEs{g}.vars,1)>1
            error('convert_PIETOOLS_BVP:sourceAndInput',"Boundary conditions of a state on several variables with both a state of lower dimension and an input are not supported.")
        end
        diffs_comp = term_diffs(PDE,diffs_comp,BC_terms);
    end
end

% % % Solve each group, and assemble the operators into those of the PDE
% Each group gives x_g = T_g*xf_g + ...; recoupling places these in the
% columns of xf = [xf_1;...;xf_G]. The inputs in the BCs are
% fundamental inputs, those derived from states of lower dimension being
% eliminated below.
sol = repmat(struct('Tops',[],'Kop',[],'Fop',[],'info',[]),1,numel(PDEs));
for g=1:numel(PDEs)
    % Derivative orders to compute for the states in the group
    diffs_g = diffs_comp(dec.x_idcs{g});
    if any(cellfun(@ischar,diffs_g))
        alphas_g = 'all';
    else
        alphas_g = vertcat(diffs_g{:});
    end
    try
        [sol(g).Tops,sol(g).Kop,sol(g).Fop,sol(g).info] = solve_BVP_group(PDEs{g},alphas_g,any(dec.w_idcs{g}>nw0));
    catch ME
        if strcmp(ME.identifier,'compose_Tmaps:inputsRankDeficient') && any(dec.w_idcs{g}>nw0)
            error('convert_PIETOOLS_BVP:sourceRankDeficient',"Boundary conditions with a state of lower dimension along a rank-deficient (e.g. periodic) direction are not supported.")
        end
        rethrow(ME)
    end
end

% % % Assemble the operators of the groups into those of the PDE
% The fundamental state of the PDE is xf = [xf_1;...;xf_G],
% stacking those of the groups, and the rows of the states of a group are
% zero in the columns of the other groups. The BCs, states and inputs of
% the groups are replaced by those of the PDE (using dec) in info. A
% 'copvar' needs a block in each row and column to define its spaces, so
% the blocks that no group defines are zero operators, which are only
% declared where needed for this.
nG = numel(PDEs);
% Zero operator with nr rows, mapping from the space s (fields vars, dom,
% dim) to the space with variables vars_out on domain dom_out
zero_op = @(nr,s,vars_out,dom_out) mat2sopvar(zeros(nr,s.dim),struct('in',{s.vars},'out',{vars_out}),struct('in',s.dom,'out',dom_out));

% % % Info on the groups, with the indices of the PDE instead of the group
info = struct('groups',struct('BC_idcs',cell(1,nG),'vars',[],'pass',[]));
fund = cell(1,nG);
for g=1:nG
    info.groups(g) = struct('BC_idcs',dec.BC_idcs{g},'vars',{sol(g).info.vars},'pass',{sol(g).info.vars(sol(g).info.pass)});
    fund{g} = sol(g).info.fund;
    [fund{g}.group] = deal(g);
    [fund{g}.x_idcs] = deal(dec.x_idcs{g});
end
info.fund = [fund{:}];
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);

% % % The maps from the fundamental state xf ('x') and the fundamental
% % % inputs wf and uf ('w', 'u') to the states and their
% % % derivatives, R. The rows of the states of a group are placed in the
% % % columns of that group, with zero blocks elsewhere.
% Whether a fundamental input f is the input r of the PDE itself (with
% table tab): all its elements, without derivative or boundary value, on
% its own variables; own_input(tab,r) is that fundamental input (group 0)
is_itself = @(f,tab,r) f.pde==r && ~any(f.Dorder) && ~any(f.loc) && isequal(sort(f.sel(:)'),1:tab(r,2)) && isequal(f.vars,vars_pde(logical(tab(r,3:2+N))));
own_input = @(tab,r) struct('vars',{vars_pde(logical(tab(r,3:2+N)))},'dom',PDE.dom(logical(tab(r,3:2+N)),:),'dim',tab(r,2),'pde',r,'Dorder',zeros(1,N),'loc',zeros(1,N),'sel',(1:tab(r,2))','group',0);
Rops = struct('x',[],'w',[],'u',[]);
Kops = struct('K',[]);
info.inputs = struct('w',[],'u',[]);
for ty={'x','w','u'}
    obj_type = ty{1};      is_state = strcmp(obj_type,'x');
    % Number of columns of each group, and the spaces of these columns
    ncols_g = zeros(1,nG);
    for g=1:nG
        if ~isempty(sol(g).Tops.(obj_type))
            ncols_g(g) = size(sol(g).Tops.(obj_type){1}{1}.C,2);
        end
    end
    off = [0,cumsum(ncols_g)];       ncols_tot = off(end);
    col_spaces = info.fund;
    % colmap{g}: the columns of the PDE that the columns of group g are
    colmap = arrayfun(@(g) off(g)+(1:ncols_g(g)),1:nG,'UniformOutput',false);
    if ~is_state
        col_spaces = cell(1,nG);
        for g=find(ncols_g>0)
            col_spaces{g} = sol(g).info.inputs.(obj_type);
            [col_spaces{g}.group] = deal(g);
            % The inputs of the PDE that the fundamental inputs derive
            % from, with the orders of the derivatives and the locations
            % in the variables of the PDE instead of those of the group
            [~,gpos] = ismember(info.groups(g).vars,vars_pde);
            for c=1:numel(col_spaces{g})
                col_spaces{g}(c).pde = dec.([obj_type,'_idcs']){g}(col_spaces{g}(c).pde);
                Dorder = zeros(1,N);    Dorder(gpos) = col_spaces{g}(c).Dorder;
                loc = zeros(1,N);       loc(gpos) = col_spaces{g}(c).loc;
                col_spaces{g}(c).Dorder = Dorder;       col_spaces{g}(c).loc = loc;
            end
        end
        col_spaces = [col_spaces{:}];
        if isempty(col_spaces)
            col_spaces = struct('vars',{},'dom',{},'dim',{},'pde',{},'Dorder',{},'loc',{},'sel',{},'group',{});
        else
            % An input shared by BCs of groups on different variables (see
            % decouple_BVPs) gives fundamental inputs in each of these
            % groups, which are the same input if they are the same
            % derivative at the same location of the same elements. Inputs
            % that are not identical, but depend on the same element of an
            % input of the PDE, are not supported
            col_group = [col_spaces.group];
            keys = cell(1,numel(col_spaces));
            for c=1:numel(col_spaces)
                f = col_spaces(c);
                keys{c} = sprintf('%d|%s|%s|%s',f.pde,mat2str(f.Dorder),mat2str(f.loc),mat2str(sort(f.sel(:))'));
            end
            [~,first_idx,col_id] = unique(keys,'stable');
            groups_c = arrayfun(@(u) unique(col_group(col_id==u)),1:numel(first_idx),'UniformOutput',false);
            for c1=1:numel(first_idx)
                for c2=c1+1:numel(first_idx)
                    f1 = col_spaces(first_idx(c1));     f2 = col_spaces(first_idx(c2));
                    if f1.pde==f2.pde && any(ismember(f1.sel,f2.sel)) && isempty(intersect(groups_c{c1},groups_c{c2}))
                        error('convert_PIETOOLS_BVP:sharedInputs','The input %s%d appears in the boundary conditions of states on different variables, in ways that cannot be written as the same fundamental input; this is not supported.',obj_type,f1.pde)
                    end
                end
            end
            colmap = arrayfun(@(g) col_id(col_group==g)',1:nG,'UniformOutput',false);
            col_spaces = col_spaces(first_idx);
            ncols_tot = numel(first_idx);
        end
        % The PDE inputs lifted into the sources (lift_BVP_sources) are
        % fundamental inputs too, with zero columns in the maps of the
        % groups: the existing fundamental input if that is the PDE input
        % itself (as for an input also in another BC), else one of its own
        % (group 0, like those not in the BCs below)
        tab_in = PDE.([obj_type,'_tab']);
        for lf = [struct('type',{},'idx',{},'C',{}),src.inputs]
            if ~strcmp(lf.type,obj_type) || any(arrayfun(@(f) is_itself(f,tab_in,lf.idx),col_spaces))
                continue
            elseif any([col_spaces.pde]==lf.idx)
                error('convert_PIETOOLS_BVP:sourceAndInput','An input in the boundary conditions of a state of lower dimension is also decomposed in other boundary conditions, in a way that is not the input itself; this is not supported.')
            end
            col_spaces(end+1) = own_input(tab_in,lf.idx);   %#ok<AGROW>
            ncols_tot = ncols_tot+1;
        end
        info.inputs.(obj_type) = col_spaces;
    end
    if ncols_tot==0
        continue
    end
    % R{k}{alpha+1}: maps to D^alpha x_k, which only depends on the columns
    % of the group of x_k (the other columns are zero)
    R_cell = cell(nx,1);
    for g=1:nG
        for k_grp=1:numel(dec.x_idcs{g})
            R_group = sol(g).Tops.x{k_grp};
            Rk = cell(size(R_group));
            for lidx=find(~cellfun(@isempty,R_group(:)))'
                blk_ref = R_group{lidx}.C{1,1};
                row = cell(1,ncols_tot);
                if ~isempty(colmap{g})
                    row(colmap{g}) = sol(g).Tops.(obj_type){k_grp}{lidx}.C(1,:);
                end
                for c=find(cellfun(@isempty,row))
                    row{c} = zero_op(blk_ref.dims(1),col_spaces(c),blk_ref.vars.out,blk_ref.dom.out);
                end
                Rk{lidx} = copvar(row);
            end
            R_cell{dec.x_idcs{g}(k_grp)} = Rk;
        end
    end
    Rops.(obj_type) = R_cell;

    if is_state
        % F: the rows of the fundamental state of each group, in its states
        F_blocks = cell(ncols_tot,nx);
        for g=1:nG
            F_blocks(off(g)+(1:ncols_g(g)),dec.x_idcs{g}) = sol(g).Fop.C;
        end
        Fop = copvar(F_blocks);
        % K: the constraints of each group, in the columns of its state;
        % the columns without a block (of other groups) are zero
        K_blocks = cell(0,ncols_tot);
        for g=1:nG
            if ~isempty(sol(g).Kop)
                K_blocks(end+(1:size(sol(g).Kop.C,1)),off(g)+(1:ncols_g(g))) = sol(g).Kop.C;
            end
        end
        if ~isempty(K_blocks)
            blk_ref = K_blocks{1,find(~cellfun(@isempty,K_blocks(1,:)),1)};
            for c=find(all(cellfun(@isempty,K_blocks),1))
                K_blocks{1,c} = zero_op(blk_ref.dims(1),col_spaces(c),blk_ref.vars.out,blk_ref.dom.out);
            end
            Kops.K = copvar(K_blocks);
        end
    end
end
% Eliminate the inputs that represent states of lower dimension, by
% composing with the maps of these states
if ~isempty(src)
    [Rops,info] = eliminate_BVP_sources(PDE,src,Rops,info);
end
info.surplus = any(arrayfun(@(s) s.info.surplus,sol));
bc_keep = find(~multi);
for g=1:nG
    info.groups(g).BC_idcs = reshape(bc_keep(info.groups(g).BC_idcs),size(info.groups(g).BC_idcs));
end

% % % Impose the BCs set aside: on the maps Rops.x, each is a block row R
% % % of operators with R*xf=0 (convert_PDE_term). The scalar BCs (without
% % % variables) are solved for the auxiliary values without variables of
% % % states with variables (as for periodic BCs) where possible, giving
% % % xf_c=Q_c*xf for these (pinned) components c, which are eliminated by
% % % xf=P*xf_free. The other rows (nonzero) are constraints, added to K.
if any(multi)                                                               % DJ (CC), 11/13/2026
    fund = info.fund;       ncol = numel(fund);
    idx_m = find(multi);    nm = numel(idx_m);
    Tops_m = struct('x',{Rops.x},'w',{Rops.w},'u',{Rops.u});
    R = cell(nm,ncol);      vars_m = cell(nm,1);        sz_m = zeros(nm,1);
    for r=1:nm
        BC = PDE_full.BC{idx_m(r)};
        vars_m{r} = vars_pde(logical(PDE_full.BC_tab(idx_m(r),3:2+N)));
        sz_m(r) = BC.size;
        for j=1:numel(BC.term)
            Pblk = convert_PDE_term(BC.term{j},Tops_m,PDE_full,struct('nonlinear',false,'outvars',{vars_m{r}}));
            if ~isempty(Pblk.w) || ~isempty(Pblk.u)
                error('convert_PIETOOLS_BVP:multiInput',msg_multi)
            elseif isempty(Pblk.x)
                continue
            end
            for c=find(~cellfun(@isempty,Pblk.x.B))
                if isempty(R{r,c})
                    R{r,c} = Pblk.x.B{c};
                else
                    R{r,c} = R{r,c}+Pblk.x.B{c};
                end
            end
        end
        % Zero blocks, so that each block row defines its space
        dom_r = PDE.dom(ismember(vars_pde,vars_m{r}),:);
        for c=find(cellfun(@isempty,R(r,:)))
            R{r,c} = zero_op(sz_m(r),fund(c),vars_m{r},dom_r);
        end
    end
    is_sc = cellfun(@isempty,vars_m);
    % S{c}: the scalar BCs stacked, acting on xf_c; on the candidates for
    % pinning (cand) these are matrices, the columns of A. Row reduction
    % W*A=rref(A) gives a pivot row for each pinned entry, which has the
    % identity on this entry and no other pinned entry; components of which
    % only some entries are pivots cannot be pinned, and are not
    % candidates.
    cand = find(arrayfun(@(f) isempty(f.vars) && ~isempty(info.groups(f.group).vars),fund));
    m_sc = sum(sz_m(is_sc));
    S = cell(1,ncol);
    for c=1:ncol
        S{c} = vertcat(R{is_sc,c});
    end
    W = eye(m_sc);      jb = zeros(1,0);        cols = zeros(1,0);
    dims_c = arrayfun(@(c) fund(c).dim,cand);
    comp_of = zeros(1,0);
    if ~isempty(cand)
        comp_of = repelem(1:numel(cand),dims_c);
    end
    keep = true(1,numel(cand));     npiv = zeros(1,numel(cand));
    if m_sc>0 && ~isempty(cand)
        A = cell2mat(cellfun(@(c) full(S{c}.params{1}),num2cell(cand),'UniformOutput',false));
        tol = 1e-10*max(1,norm(A,'fro'));
        while true
            cols = find(keep(comp_of));
            [RW,jb] = rref([A(:,cols),eye(m_sc)],tol);
            jb = jb(jb<=numel(cols));
            npiv = accumarray([comp_of(cols(jb))';numel(cand)],[ones(numel(jb),1);0])';
            part = keep & npiv>0 & npiv<dims_c;
            if ~any(part)
                break
            end
            keep(part) = false;
        end
        W = RW(:,numel(cols)+1:end);
    end
    pinned = cand(npiv>0);
    WS = cellfun(@(X) W*X,S,'UniformOutput',false);
    % Q{i,c}: xf_{pinned(i)} = sum_c Q{i,c}*xf_c, from the pivot rows of
    % its entries (in order)
    Q = cell(numel(pinned),ncol);
    for i=1:numel(pinned)
        E = eye(m_sc);
        E = E(cand(comp_of(cols(jb)))==pinned(i),:);
        Q(i,:) = cellfun(@(X) -E*X,WS,'UniformOutput',false);
    end
    % The constraints: the scalar rows other than the pivot rows, and the
    % BCs with variables; zero rows (BCs implied by the others) are dropped
    E = eye(m_sc);      E = E(numel(jb)+1:end,:);
    K_rows = R(~is_sc,:);
    if ~isempty(E)
        K_rows = [cellfun(@(X) E*X,WS,'UniformOutput',false);K_rows];
    end
    coefs = cellfun(@(X) max([0;cellfun(@(p) max([0;abs(full(p(:)))]),X.params(:))]),K_rows);
    tol_K = 1e-10*max([1;coefs(:)]);
    is_zero = all(coefs<=tol_K,2);
    if ~isempty(E) && ~is_zero(1)
        % Zero rows of the scalar block, by entry: the params of an operator
        % without output variables have a row for each entry
        nz = false(size(E,1),1);
        for c=1:ncol
            for p=K_rows{1,c}.params(:)'
                nz = nz | any(abs(full(p{1}))>tol_K,2);
            end
        end
        E = eye(numel(nz));     E = E(nz,:);
        K_rows(1,:) = cellfun(@(X) E*X,K_rows(1,:),'UniformOutput',false);
    end
    K_rows = K_rows(~is_zero,:);
    info.surplus = info.surplus || ~isempty(pinned) || ~isempty(K_rows);

    % Eliminate the pinned components, xf=P*xf_free, from the maps to the
    % states, the constraints, and the map F from the states to xf
    free = setdiff(1:ncol,pinned);
    if ~isempty(pinned)
        Pc = cell(ncol,numel(free));
        for j=1:numel(free)
            f = fund(free(j));
            Pc{free(j),j} = mat2sopvar(eye(f.dim),struct('in',{f.vars},'out',{f.vars}),struct('in',f.dom,'out',f.dom));
        end
        Pc(pinned,:) = Q(:,free);
        P = copvar(Pc);
        for k=1:nx
            for l=find(~cellfun(@isempty,Rops.x{k}(:)))'
                Rops.x{k}{l} = Rops.x{k}{l}*P;
            end
        end
        if ~isempty(Kops.K)
            Kops.K = Kops.K*P;
        end
        Fop = copvar(Fop.C(free,:));
        info.fund = fund(free);
    end
    if ~isempty(K_rows)
        K_new = copvar(K_rows);
        if ~isempty(pinned)
            K_new = K_new*P;
        end
        if isempty(Kops.K)
            Kops.K = K_new;
        else
            Kops.K = [Kops.K;K_new];
        end
    end
end
% Keep only the requested derivatives (the others were needed internally)
if want_diffs
    for k=1:nx
        keep = false(size(Rops.x{k}));
        keep(1) = true;
        if ischar(diffs_k{k})
            keep(:) = true;
        elseif any(isvar(k,:)) && ~isempty(diffs_k{k})
            Dk = Delta(k,isvar(k,:));
            keep(1+diffs_k{k}*cumprod([1,Dk(1:end-1)+1])') = true;
        end
        for types={'x','w','u'}
            if ~isempty(Rops.(types{1}))
                Rops.(types{1}){k}(~keep) = {[]};
            end
        end
    end
end

% % % The blocks of the fundamental state: for each xf_c that stacks the
% % % states in its group (p equals the total size), a block for each state
% % % with components ordered as (auxiliary values, states)
% A xf_c with auxiliary values stacks the states as [xf_c(aux 1);...;
% xf_c(aux m)], so state k has entries (j-1)*ns+off(k)+(1:size_k). If
% p~=ns, xf_c cannot be split per state and stays one block, with the
% Dorder of its state if it has one, and zero otherwise.
sizes = PDE.x_tab(:,2)';
info.blocks = struct('fund',{},'x',{},'comp',{},'dim',{},'vars',{},'dom',{},'Dorder',{});
for c=1:numel(info.fund)
    f = info.fund(c);
    [~,gpos] = ismember(info.groups(f.group).vars,vars_pde);
    ns = sum(sizes(f.x_idcs));      off = [0,cumsum(sizes(f.x_idcs))];
    if f.p==ns
        for k=1:numel(f.x_idcs)
            % Indices of state k in xf_c, and the orders of its
            % derivative (in the variables of the group) that xf_c
            % represents
            comp = reshape((0:f.dim/ns-1)'*ns+off(k)+(1:sizes(f.x_idcs(k))),1,[]);
            D = zeros(1,N);     D(gpos) = f.Dorder(k,:);
            info.blocks(end+1) = struct('fund',c,'x',f.x_idcs(k),'comp',comp,'dim',numel(comp),'vars',{f.vars},'dom',f.dom,'Dorder',D);
        end
    else
        D = zeros(1,N);
        if isscalar(f.x_idcs)
            D(gpos) = f.Dorder;
        end
        info.blocks(end+1) = struct('fund',c,'x',f.x_idcs,'comp',1:f.dim,'dim',f.dim,'vars',{f.vars},'dom',f.dom,'Dorder',D);
    end
end

% % % Add the inputs that are not in the BCs, if requested
% These inputs enter the PDE only through the equations, so they become
% fundamental inputs of their own, with zero columns in the maps to the
% state.
if is_an_opt('all_inputs')
    % Spaces of the states (the rows of the new zero columns)
    dom_pde = struct('vars',{vars_pde},'dom',PDE.dom);
    spaces_x = arrayfun(@(k) vars_pde(logical(PDE.x_tab(k,3:2+N))),1:nx,'UniformOutput',false);
    % The exogenous (w) and actuator (u) inputs are treated identically
    for ty={'w','u'}
        obj_type = ty{1};
        tab = PDE.([obj_type,'_tab']);
        nin = size(tab,1);
        inputs = info.inputs.(obj_type);
        n_fund = numel(inputs);

        % % % The inputs without a fundamental input (those not in the
        % % % BCs, nor lifted into a source), as additional fundamental
        % % % inputs (group 0, after the others) with zero columns in the
        % % % maps to the state, as they only appear in the equations
        inter = setdiff(1:nin,[inputs.pde]);
        for r=inter
            inputs(end+1) = own_input(tab,r);                               %#ok<AGROW>
        end
        % Append zero columns to the maps from the fundamental inputs to
        % the states, so that their columns match the full list of inputs.
        % With no decomposed inputs, Rops.(t) is empty and stays so.
        % zero_cols(k): zero operator from the new inputs to the states k
        if n_fund>0 && ~isempty(inter)
            new_in = inputs(n_fund+1:end);
            zero_cols = @(k) zeros_copvar_sop(struct('out',PDE.x_tab(k,2),'in',[new_in.dim]'),struct('out',{spaces_x(k)},'in',{{new_in.vars}}),dom_pde);
            for k=1:numel(Rops.(obj_type))
                Zcols_k = zero_cols(k);
                for l=find(~cellfun(@isempty,Rops.(obj_type){k}(:)))'      % computed derivatives
                    Rops.(obj_type){k}{l} = [Rops.(obj_type){k}{l},Zcols_k];
                end
            end
        end

        % % % Maps from the fundamental inputs to the inputs of the PDE
        % Iin{r}: the row of operators with w_r = Iin{r}*wf, the identity
        % on the fundamental input of w_r if it is the only one and is
        % w_r itself (for the inputs of their own, also those lifted into
        % sources, and an input in the BCs that is itself fundamental,
        % which may then also be used in the equations), and zero
        % elsewhere; empty for the inputs decomposed in the BCs, which
        % cannot be used elsewhere (see convert_PDE_term)
        Iin = cell(nin,1);
        for r=1:nin
            c_r = find([inputs.pde]==r);
            if numel(c_r)~=1 || ~is_itself(inputs(c_r),tab,r)
                continue
            end
            dep_mask = logical(tab(r,3:2+N));
            Mat = zeros(tab(r,2),sum([inputs.dim]));
            Mat(:,sum([inputs(1:c_r-1).dim])+(1:tab(r,2))) = eye(tab(r,2));
            Iin{r} = mat2copvar_sop(Mat,struct('out',tab(r,2),'in',[inputs.dim]'),struct('out',{{vars_pde(dep_mask)}},'in',{{inputs.vars}}),dom_pde);
        end
        info.(['I',obj_type]) = Iin;
        info.inputs.(obj_type) = inputs;
    end
end

% % % The maps to the states: single operators (the first of the
% % % derivatives of each state, stacked) if no derivatives were
% % % requested, or the cells Rops of the maps to the requested derivatives
if want_diffs
    Tops = struct('x',{Rops.x},'w',{Rops.w},'u',{Rops.u});
else
    Tops = struct('x',[],'w',[],'u',[]);
    for ty={'x','w','u'}
        if ~isempty(Rops.(ty{1}))
            rows = cellfun(@(Rk) Rk{1},Rops.(ty{1}),'UniformOutput',false);
            Tops.(ty{1}) = vertcat(rows{:});
        end
    end
end

% % % Return the operators as 'tensopmat' objects for a nonlinear PDE, with
% % % a block for each state, and as 'copvar' objects for a linear PDE
% A nonlinear PDE needs a block for each state: the columns of the maps
% from xf and the rows of F are split according to info.blocks (inputs
% are already single columns). A linear PDE keeps the merged objects.
use_nonlin = PDE.is_nonlinear || is_an_opt('nonlinear');
is_split = use_nonlin || is_an_opt('blocks');
if is_split
    split_op = @(P,rsplit,csplit) split_copvar(P,rsplit,csplit);
    if use_nonlin
        split_op = @(P,rsplit,csplit) tensopmat(split_copvar(P,rsplit,csplit,'tensopmat'));
    end
    % csplit{c}: the index sets of the components of xf_c forming the
    % blocks
    csplit = arrayfun(@(c) {{info.blocks([info.blocks.fund]==c).comp}},1:numel(info.fund));
    % K*xf=0
    if ~isempty(Kops.K)
        Kops.K = split_op(Kops.K,[],csplit);
    end
    % x_pde --> part of xf that is not free (rows split, columns are PDE
    % states)
    Fop = split_op(Fop,csplit(:),[]);
    % Maps from xf (split in columns), wf and uf to x_pde (or its
    % derivatives, where computed: the other entries of Rops.x{k}{alpha+1}
    % are [] and stay so)
    col_split = struct('x',{csplit},'w',[],'u',[]);
    for types={'x','w','u'}
        obj_type = types{1};
        if isempty(Rops.(obj_type))
            continue
        elseif want_diffs
            for k=1:numel(Rops.(obj_type))
                for l=find(~cellfun(@isempty,Rops.(obj_type){k}(:)))'
                    Tops.(obj_type){k}{l} = split_op(Tops.(obj_type){k}{l},[],col_split.(obj_type));
                end
            end
        else
            Tops.(obj_type) = split_op(Tops.(obj_type),[],col_split.(obj_type));
        end
    end
    % Maps of the fundamental inputs to the PDE inputs (Iw, Iu), [] for
    % the inputs in the BCs
    if is_an_opt('all_inputs')
        for fn={'Iw','Iu'}
            for r=find(~cellfun(@isempty,info.(fn{1}))')
                info.(fn{1}){r} = split_op(info.(fn{1}){r},[],[]);
            end
        end
    end
end

% % % The spaces of the PDE and of the fundamental state and inputs; the
% % % components of xf are the blocks if the operators are split, and
% % % those of info.fund otherwise
space = struct('vars',{vars_pde},'dom',reshape(double(PDE.dom),N,2));
comps_x = info.fund;
if is_split
    comps_x = info.blocks;
end
for types={'x','w','u'}
    comps = comps_x;
    if ~strcmp(types{1},'x')
        comps = info.inputs.(types{1});
    end
    % depmat(c,j) is true if the component c depends on the variable j
    space.(types{1}) = struct('dims',reshape([comps.dim],[],1),'depmat',false(numel(comps),N));
    for c=1:numel(comps)
        space.(types{1}).depmat(c,:) = ismember(vars_pde,comps(c).vars);
    end
end

% % % The map Dmap from the PDE to the fundamental state and inputs: xf_c
% % % = sum_k C_{c,k}*D^{diff(c,k,:)}x_pde_k, with C the operator Fop, and
% % % the inputs as the PDE inputs (see the header for those in the BCs)
ncomp = numel(comps_x);
diff_arr = zeros(ncomp,nx,N);      idx_x = zeros(ncomp,1);
for c=1:ncomp
    if is_split
        f = info.fund(comps_x(c).fund);     states = comps_x(c).x;
    else
        f = info.fund(c);                   states = f.x_idcs;
    end
    [~,gpos] = ismember(info.groups(f.group).vars,vars_pde);
    for j=find(ismember(f.x_idcs,states))
        diff_arr(c,f.x_idcs(j),gpos) = f.Dorder(j,:);
    end
    if isscalar(states)
        idx_x(c) = states;
    end
end
Dmap = struct('x',struct('idx',idx_x,'diff',diff_arr,'loc',zeros(ncomp,N),'C',Fop),'w',[],'u',[]);
for types={'w','u'}
    obj_type = types{1};
    inputs = info.inputs.(obj_type);
    Dmap.(obj_type) = struct('idx',[],'diff',zeros(0,N),'loc',zeros(0,N),'C',[]);
    if ~isempty(inputs)
        % The PDE input of each fundamental input, differentiated, and
        % evaluated at the lower boundary of the variables that the
        % fundamental input does not depend on; C selects the elements
        % of this PDE input that the fundamental input is
        Dmap.(obj_type).idx = [inputs.pde]';
        Dmap.(obj_type).diff = reshape([inputs.Dorder],N,[])';
        Dmap.(obj_type).loc = reshape([inputs.loc],N,[])';
        C_blocks = cell(numel(inputs));
        for c=1:numel(inputs)
            dim_pde = PDE.([obj_type,'_tab'])(inputs(c).pde,2);
            C_blocks{c,c} = mat2sopvar(full(sparse(1:inputs(c).dim,inputs(c).sel,1,inputs(c).dim,dim_pde)),struct('in',{inputs(c).vars},'out',{inputs(c).vars}),struct('in',inputs(c).dom,'out',inputs(c).dom));
        end
        Dmap.(obj_type).C = copvar(C_blocks);
        if is_split
            Dmap.(obj_type).C = split_op(Dmap.(obj_type).C,[],[]);
        end
    end
end

% Whether the inputs enter the states, or Tops.w, Tops.u are zero (below)
info.input_maps = struct('w',~isempty(Tops.w),'u',~isempty(Tops.u));
% % % The constraints Kx*xf + Kw*wf + Ku*uf = 0, with a row for each
% % % constraint (none if there are none); the inputs do not enter the
% % % constraints (yet), so Kw and Ku are zero. Tops.w and Tops.u are zero
% % % if there are no such inputs in the BCs. zero_map(P,sp) is the zero
% % % operator with the rows of P and the columns of the space sp
Kx = Kops.K;
Kops = struct('x',Kx,'w',[],'u',[]);
if use_nonlin
    [vars_srt,ord] = sort(vars_pde);
    zero_map = @(P,sp) mat2tensopmat(zeros(sum(P.dims.out),sum(sp.dims)),struct('in',sp.dims,'out',P.dims.out),...
        struct('in',{vars_srt},'out',{P.vars.out}),struct('in',space.dom(ord,:),'out',P.dom.out),struct('in',double(sp.depmat(:,ord)),'out',P.depmat.out));
    no_rows = struct('dims',struct('out',zeros(0,1)),'vars',struct('out',{cell(1,0)}),'dom',struct('out',zeros(0,2)),'depmat',struct('out',zeros(0,0)));
else
    zero_map = @(P,sp) zeros_copvar_sop(struct('out',P.dim_out,'in',sp.dims),...
        struct('out',{arrayfun(@(i) P.vars(P.space_out(i,:)),1:numel(P.dim_out),'UniformOutput',false)},'in',{arrayfun(@(i) vars_pde(sp.depmat(i,:)),1:numel(sp.dims),'UniformOutput',false)}),...
        struct('vars',{vars_pde},'dom',space.dom));
    no_rows = struct('dim_out',zeros(0,1),'vars',{cell(1,0)},'space_out',false(0,0));
end
if isempty(Kx)
    Kx = no_rows;
    Kops.x = zero_map(Kx,space.x);
end
for ty={'w','u'}
    sp = space.(ty{1});
    Kops.(ty{1}) = zero_map(Kx,sp);
    if ~isempty(Tops.(ty{1}))
        continue
    end
    if ~want_diffs
        Tops.(ty{1}) = zero_map(Tops.x,sp);
        continue
    end
    Tops.(ty{1}) = cell(size(Tops.x));
    for k=1:numel(Tops.x)
        Tops.(ty{1}){k} = cell(size(Tops.x{k}));
        for l=find(~cellfun(@isempty,Tops.x{k}(:)))'
            Tops.(ty{1}){k}{l} = zero_map(Tops.x{k}{l},sp);
        end
    end
end
info.space = space;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %




%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function diffs_k = term_diffs(PDE,diffs_k,terms)
% diffs_k = TERM_DIFFS(PDE,DIFFS_K,TERMS) adds the derivative orders of
% the states in the terms to the lists of derivatives of each state.
%
% INPUT
% - PDE:        initialized 'pde_struct' object;
% - diffs_k:   nx x 1 cell, with diffs_k{k} either 'all' or an array with
%               a row of derivative orders for each derivative of state k;
% - terms:      cell of terms of the PDE, each a struct array of factors
%               with fields x, w or u, and D;
%
% OUTPUT
% - diffs_k:   the same cell, with the derivatives D of each factor of
%               state k in the terms added (unless diffs_k{k} is 'all');
N = size(PDE.vars,1);
for j=1:numel(terms)
    for fctr = terms{j}(:)'
        if ~isfield(fctr,'x') || isempty(fctr.x) || ischar(diffs_k{fctr.x})
            continue
        end
        D = zeros(1,nnz(PDE.x_tab(fctr.x,3:2+N)));
        if isfield(fctr,'D') && ~isempty(fctr.D)
            D = fctr.D(:)';
        end
        if isempty(diffs_k{fctr.x}) || ~ismember(D,diffs_k{fctr.x},'rows')
            diffs_k{fctr.x} = [diffs_k{fctr.x};D];
        end
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [Tops,Kop,Fop,info] = solve_BVP_group(PDE,alphas,sources)
% [Tops,Kop,Fop,info] = SOLVE_BVP_GROUP(PDE,alphas,sources) computes
% the operators parameterizing the states that satisfy the BCs of one
% group of coupled BCs (see decouple_BVPs): the BCs are split into 1D
% problems along each direction (split_BVPs), these are solved by
% compute_Tmap_1D, and the solutions are composed by compose_Tmaps.
%
% INPUT
% - PDE:    initialized 'pde_struct' whose states all depend on the same
%           variables. Only the BCs are used;
% - alphas: the derivative orders to compute, as in compose_Tmaps;
% - sources: true if the inputs in the BCs represent states of lower
%           dimension, see compose_Tmaps;
%
% OUTPUT
% - Tops,Kop,Fop,info:   as returned by compose_Tmaps, for the
%           states, BCs and inputs of PDE;
N = size(PDE.vars,1);
sizes = PDE.x_tab(:,2)';
isvar = logical(PDE.x_tab(:,3:2+N));
Delta = PDE.x_tab(:,2+N+(1:N));         % order of x_k in the i-th variable
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);
dir_idcs = 1:N;
if isempty(PDE.BC)
    % States without BCs must not be differentiable, and have xf=x
    if any(Delta(:)>0)
        error('convert_PIETOOLS_BVP:noBCs','A state that is differentiable in some direction has no boundary conditions along it.')
    end
    dir_idcs = find(any(isvar,1));
else
    % The BCs along each direction, with the other variables as parameters
    dirs = split_BVPs(PDE);
    % The states must all depend on the same variables, so that the
    % right-hand sides of the BCs (forcing) only involve inputs
    if ~all(isvar(:))
        error('convert_PIETOOLS_BVP:mixedDimension','Boundary conditions that couple states on different variables are not supported.')
    end
end

% % % Solve the 1D problem along each direction in which a state is
% % % differentiable (split_BVPs checks that there are BCs along it)
maps = repmat(struct('Tops',[],'Rops',[],'Kops',[],'Fop',[],'Qops',[],'info',[],'forced',[],'source',[]),1,numel(dir_idcs));
for i=find(max(Delta(:,dir_idcs),[],1)>0)
    % forced.w(r) (forced.u) is true if row r of the BCs along direction i
    % has inputs w (u) on its right-hand side; the maps of these (Qops of
    % compute_Tmap_1D) are only needed then
    nrow = dirs(i).nrow;
    for ty={'w','u'}
        obj_type = ty{1};
        has_input = cellfun(@(f) any(cellfun(@(x) isfield(x,obj_type),f)),dirs(i).forcing);
        maps(i).forced.(obj_type) = repelem(has_input,nrow);
        % The element of an input that is the right-hand side of each
        % condition (src.pde=0 if it is not one element of one input)
        tab = PDE.([obj_type,'_tab']);
        face = true(1,N);           face(i) = false;
        src = struct('pde',zeros(sum(nrow),1),'elem',zeros(sum(nrow),1),'scale',ones(sum(nrow),1),'dep',false(sum(nrow),N));
        off = 0;
        for k=1:numel(nrow)
            terms = dirs(i).forcing{k};
            terms = terms(cellfun(@(tm) isfield(tm,obj_type),terms));
            % Coefficients as numbers: constant polynomials are
            % converted, and a coefficient depending on the variables
            % (as for lifted states) is not a single element
            C_num = cell(size(terms));
            for jj=1:numel(terms)
                C_num{jj} = terms{jj}.C;
                if ~isnumeric(C_num{jj})
                    try
                        C_num{jj} = double(C_num{jj});
                    catch
                        C_num{jj} = NaN(size(C_num{jj}));
                    end
                end
            end
            for r=1:nrow(k)
                % Elements with a nonzero coefficient in this row
                nz = cellfun(@(C) find(C(r,:)~=0),C_num,'UniformOutput',false);
                hit = find(~cellfun(@isempty,nz));
                if isscalar(hit) && isscalar(nz{hit}) && ~isnan(C_num{hit}(r,nz{hit})) && all(cellfun(@isempty,terms{hit}.I)) && all(face(logical(tab(terms{hit}.(obj_type),2+(1:N)))))
                    src.pde(off+r) = terms{hit}.(obj_type);
                    src.dep(off+r,:) = logical(tab(terms{hit}.(obj_type),2+(1:N)));
                    src.elem(off+r) = nz{hit};
                    src.scale(off+r) = C_num{hit}(r,nz{hit});
                end
            end
            off = off+nrow(k);
        end
        maps(i).source.(obj_type) = src;
    end
    % The orders of the derivatives in this direction to compute
    orders = 0:max(Delta(:,i));
    if ~ischar(alphas)
        orders = 0;
        if ~isempty(alphas)
            orders = unique([0;alphas(:,i)])';
        end
    end
    % Errors are prefixed with the direction, to locate the BCs
    try
        if any(maps(i).forced.w) || any(maps(i).forced.u)
            [maps(i).Tops,maps(i).Rops,maps(i).Kops,maps(i).Fop,maps(i).Qops] = compute_Tmap_1D(dirs(i).n,dirs(i).E,dirs(i).F,PDE.dom(i,:),PDE.vars(i,:),orders);
        else
            [maps(i).Tops,maps(i).Rops,maps(i).Kops,maps(i).Fop] = compute_Tmap_1D(dirs(i).n,dirs(i).E,dirs(i).F,PDE.dom(i,:),PDE.vars(i,:),orders);
            maps(i).forced = [];
        end
    catch ME
        msg = ['Boundary conditions along ',vars_pde{i},': ',ME.message];
        if isempty(ME.identifier)
            error('convert_PIETOOLS_BVP:direction','%s',msg)
        else
            error(ME.identifier,'%s',msg)
        end
    end
    maps(i).info = dirs(i);
end

% % % Compose the solutions along the different directions
[Tops,Kop,Fop,info] = compose_Tmaps(maps,vars_pde(dir_idcs),PDE.dom(dir_idcs,:),sizes,Delta(:,dir_idcs),alphas,sources);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

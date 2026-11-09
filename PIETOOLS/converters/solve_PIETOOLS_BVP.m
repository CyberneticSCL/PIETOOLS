function [Tops,Kops,Dmap,info] = solve_PIETOOLS_BVP(PDE,diffs,opts)
% [Tops,Kops,Dmap,info] = SOLVE_PIETOOLS_BVP(PDE,diffs,opts) computes
% operators parameterizing all states of an ND PDE that satisfy its
% boundary conditions (BCs): with the fundamental state xf (free of BCs,
% except for the constraints Kops.x*xf=0) and the fundamental inputs wf,
% uf (the right-hand sides of the BCs), for all x_pde in the domain D
%       x_pde = Tops.x*xf + Tops.w*wf + Tops.u*uf,
% and, if requested, for the derivatives of the states
%       D^{alpha}x_pde_k = Tops.x{k}{alpha+1}*xf
%                          + Tops.w{k}{alpha+1}*wf
%                          + Tops.u{k}{alpha+1}*uf.
% Conversely, Dmap defines xf and wf, uf in terms of the PDE
% states and inputs. The BCs are split into groups that share no state
% (decouple_BVPs), then into 1D problems along each direction (split_BVPs),
% solved by solve_BVP_1D, and composed (compose_Tmaps). The dynamics are
% ignored.
%
% INPUT
% - PDE:    'pde_struct' object in the terms format. The BCs may involve
%           states of fewer variables, and inputs w and u, see NOTES;
% - diffs:  (optional) the derivatives of the states to compute. [] (the
%           default) for none, 'all' for all derivatives up to the order
%           of each state, 'terms' for those in the terms of the equations
%           and outputs of PDE (as needed by convert_PDE_term), or an
%           nx x 1 cell with diffs{k} for state k either 'all', or an
%           array with a row alpha=(alpha_1,...,alpha_Nk) of derivative
%           orders for each derivative, in the variables of x_pde_k;
% - opts:   (optional) string or cell of strings: 'no_decouple' (treat all
%           BCs as one group), 'all_inputs' (also include the inputs that
%           are not in the BCs), 'blocks' (see OUTPUT),
%           and 'nonlinear' (return 'tensopmat' objects).
%           Either of diffs and opts may be omitted or [], so that e.g.
%           solve_PIETOOLS_BVP(PDE,'terms') and
%           solve_PIETOOLS_BVP(PDE,'blocks') both work: a single input is
%           diffs if it is 'all', 'terms' or not text, and opts otherwise.
%           The order (opts,diffs) is also recognized;
%
% OUTPUT
% All operators are 'copvar' objects for a linear PDE, and 'tensopmat'
% objects for a nonlinear PDE or the option 'nonlinear'. The components
% xf_c of xf are the columns of Tops and Kops, and the rows of Dmap:
% a 'copvar' merges the states of a group in one component (info.fund),
% but with the option 'blocks' (and for 'tensopmat') splits it into a block
% for each state (info.blocks).
% - Tops:   struct with fields x, w, u: the maps xf --> x_pde,
%           wf --> x_pde and uf --> x_pde (w, u are zero operators if
%           the BCs have no such inputs, without columns if there are
%           no inputs at all). Without diffs, these are single
%           operators. With diffs, nx x 1 cells, with
%           x{k}{alpha_1+1,...,alpha_N+1} the
%           1 x nV operator to D^{alpha}x_pde_k for the requested alpha
%           (and alpha=0), and [] for the others, and w{k}, u{k} the same
%           for wf, uf;
% - Kops:   struct with fields x, w, u, the operators of the constraints
%           Kops.x*xf + Kops.w*wf + Kops.u*uf = 0, with a row for each
%           constraint, as operators of the same class as Tops. Without
%           constraints, they have no rows. The inputs do not enter the
%           constraints (yet), so Kops.w and Kops.u are zero operators
%           (without rows if there are no constraints);
% - Dmap:   struct with fields x, w, u, defining the components of xf,
%           wf and uf in terms of the PDE, as PIE.pde2pie does for
%           a PIE: component c is obtained from the PDE component idx(c)
%           by taking the derivative of order diff(c,:), evaluating at the
%           lower boundary where loc(c,:)=-1, and applying the operator C.
%           Each field is a struct with fields
%           - idx:  ncomp x 1 array, the index of that PDE state or input
%                   (0 for a component of xf stacking several states);
%           - diff: ncomp x N array, the orders of the derivative in the
%                   variables of the PDE, for x ncomp x nx x N, with
%                   diff(c,k,:) for the state k in xf_c;
%           - loc:  ncomp x N array, -1 for the variables of the PDE input
%                   that component c does not depend on, in which it
%                   is evaluated at the lower boundary, 0 otherwise (0
%                   for x);
%           - C:    for x, the operator with a row for each component and
%                   a column for each state, so that
%                   xf_c = sum_k C_{c,k}*D^{diff(c,k,:)}x_pde_k, the
%                   bounded part of the map (the auxiliary values); for w
%                   and u, a block diagonal operator selecting the
%                   elements of the PDE input that component c of wf
%                   (uf) consists of, so that, for wf,
%                   wf_c = C_{c,c}*D^{diff(c,:)}w_{idx(c)}(loc);
%           Each fundamental input is thus a derivative of an input of the
%           PDE, evaluated at a boundary in the variables that it does not
%           depend on (a PDE input can give several fundamental inputs,
%           e.g. its value and derivatives at the lower boundary, and its
%           highest derivative). For the inputs in the BCs, these are only
%           the derivatives needed to parametrize the compatible right-hand
%           sides, see fund_inputs. The other BC data
%           (e.g. the value of a second input at a corner) follow from
%           these, and are not independent of them;
% - info:   struct with fields
%           - space:    space of the fundamental
%                       state and inputs, a struct with fields vars (names
%                       of the variables), dom (N x 2 array) and x, w, u,
%                       each with fields dims (ncomp x 1 array) and depmat
%                       (ncomp x N logical, true if the component depends
%                       on the variable). The components of space are those
%                       of Tops and Dmap;
%           - groups:   struct array, with for each group its BCs
%                       (BC_idcs) and variables (vars);
%           - fund:     struct array, with for each xf_c its group,
%                       x_idcs, vars, dom, dim, Dorder and p;
%           - blocks:   struct array, with for each block of xf: fund
%                       (index c), x (the PDE state), comp (entries of
%                       xf_c), dim, vars, dom, Dorder (derivative of the
%                       state in each variable);
%           - inputs, rhs:  the components of wf and uf, with for
%                       each the PDE input pde, derivative Dorder,
%                       location loc and elements sel that it is (see
%                       Dmap), and the BC right-hand sides they enter,
%                       see compose_Tmaps;
%           - input_maps: struct with fields w, u, true if the inputs enter
%                       the states (Tops.w, Tops.u), false if these are
%                       zero operators, as no inputs are in the BCs;
%           - Gops:     struct with fields w, u: operators mapping wf
%                       (uf) to the right-hand sides of the BCs with
%                       such inputs (fundamental inputs only, [] else),
%                       and with 'all_inputs' fields Iw, Iu: cells with
%                       Iw{r} the operator wf --> w_r for the input w_r
%                       that is its own fundamental input ([] for those
%                       decomposed in the BCs);
% NOTES
% - See compose_Tmaps for the construction, and conditions for support.
% - The number of BCs along a direction should be at least the order of the
%   states along it. Surplus BCs (beyond this order) are constraints on the
%   fundamental state, which are included in Kops.x*xf=0, see compose_Tmaps
%   for the conditions (e.g. no inputs in the BCs along such a direction).
% - Computing only the derivatives needed (diffs) avoids an operator for
%   each of the up to prod(Delta_k+1) derivatives of each state.
% - States of lower dimension in the BCs of other states (e.g. ODE states)
%   act as inputs to those BCs: they are lifted (lift_BVP_sources), solved
%   as fundamental inputs, and eliminated (eliminate_BVP_sources). They
%   must satisfy the compatibility conditions of the BCs they enter. All
%   their derivatives are computed, and those of the states in these BCs.
%   Rank-deficient (e.g. periodic) directions of the states they enter,
%   without inputs, are passed through: the data of the inputs are the
%   source functions in these directions.
% - Inputs w and u in the BCs are decomposed into free fundamental inputs
%   wf, uf, so that the BC data are compatible by construction and
%   the user need not enforce any compatibility. This requires all
%   directions of the group of BCs to be admissible: inputs in the BCs of a
%   group with a rank-deficient (e.g. periodic) direction, in a direction
%   with rank-deficient or no conditions, or in the BCs of a group of
%   several variables that also involves a state of lower dimension (in a
%   group of one variable, with a scalar source, the inputs in the BC of the
%   source are lifted into it), raise an error (the states of
%   lower dimension themselves may enter along admissible directions of a
%   group with rank-deficient directions). 
% - The fundamental inputs are derivatives of the PDE inputs (see Dmap),
%   chosen in fund_inputs. This needs the right-hand side of each BC row
%   with an input to be a single element of a single input, that depends
%   on some of the variables of the BCs other than the one it is along (a
%   function of fewer variables is constant in the others, which
%   restricts the other inputs).
% - A boundary value of a state in the BC of a state on other variables,
%   e.g. x4(a)=x3(d), acts as a state of lower dimension (see
%   lift_BVP_sources).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - solve_PIETOOLS_BVP
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
        error('solve_PIETOOLS_BVP:diffs',"The derivatives should be specified as 'all', 'terms', or an nx x 1 cell.")
    end
    for k=1:nx
        Dk = Delta(k,isvar(k,:));
        if ischar(diffs_k{k})
            bad = ~strcmpi(diffs_k{k},'all');
        else
            bad = ~isempty(diffs_k{k}) && (~isnumeric(diffs_k{k}) || size(diffs_k{k},2)~=numel(Dk) || any(diffs_k{k}(:)<0) || any(diffs_k{k}(:)~=round(diffs_k{k}(:))) || any(any(diffs_k{k}>Dk,2)));
        end
        if bad
            error('solve_PIETOOLS_BVP:diffs','The derivatives of state %d should be ''all'', or nonnegative integer orders, one for each of its variables, not exceeding the orders of the state.',k)
        end
    end
end


% % % Split the BCs into groups sharing no states
% Terms of states of lower dimension (e.g. ODE states) are first lifted to
% inputs, eliminated after solving, so that each group involves states on
% the same variables. Groups are solved independently, keeping the cost
% small; without decoupling, there is a single group.
[PDE_lift,src] = lift_BVP_sources(PDE);
if isempty(src)
    PDE_lift = PDE;
end
nw0 = numel(PDE.w);
if is_an_opt('no_decouple')
    PDEs = {PDE_lift};
    dec = struct('BC_idcs',{{1:numel(PDE.BC)}},'x_idcs',{{1:nx}},'w_idcs',{{1:numel(PDE_lift.w)}},'u_idcs',{{1:numel(PDE.u)}});
else
    [PDEs,dec] = decouple_BVPs(PDE_lift);
end
% Eliminating the sources needs all derivatives of the lifted states, and
% those in the BCs of the states they enter (see eliminate_BVP_sources)
diffs_comp = diffs_k;
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
        % PDE can be ensured, which therefore cannot be mixed in. This only
        % holds for groups of several variables (corners): in a group of one
        % variable the inputs left in the BCs (not lifted into a source) are
        % free. The inputs of the sources are not inputs of the PDE.
        % % % Start of changes DJ (CC), 11/08/2026
        BC_terms_l = cellfun(@(b) b.term(:),PDE_lift.BC(dec.BC_idcs{g}(:)),'UniformOutput',false);
        BC_terms_l = vertcat(BC_terms_l{:});
        has_input = cellfun(@(tm) any(arrayfun(@(f) (isfield(f,'w') && ~isempty(f.w) && f.w<=nw0) || (isfield(f,'u') && ~isempty(f.u)),tm)),BC_terms_l);
        % % % End of changes DJ (CC), 11/08/2026
%         has_input = cellfun(@(tm) any(arrayfun(@(f) (isfield(f,'w') && ~isempty(f.w)) || (isfield(f,'u') && ~isempty(f.u)),tm)),BC_terms); % DJ (CC), 11/08/2026 (was)
%         if any(has_input)                                                 % DJ (CC), 11/08/2026 (was)
        if any(has_input) && size(PDEs{g}.vars,1)>1                         % DJ (CC), 11/08/2026
            error('solve_PIETOOLS_BVP:sourceAndInput','The boundary conditions %s involve both a state of lower dimension and an input; this is not supported, as the input could not be parameterized by free functions.',mat2str(dec.BC_idcs{g}(:)'))
        end
        diffs_comp = term_diffs(PDE,diffs_comp,BC_terms);
    end
end

% % % Solve each group, and assemble the operators into those of the PDE
% Each group gives x_g = T_g*xf_g + ...; recoupling places these in the
% columns of xf = [xf_1;...;xf_G]. The inputs in the BCs are
% fundamental inputs, those derived from states of lower dimension being
% eliminated below.
sol = repmat(struct('Tops',[],'Kops',[],'Fop',[],'info',[],'Gops',[]),1,numel(PDEs));
for g=1:numel(PDEs)
    % Derivative orders to compute for the states in the group
    diffs_g = diffs_comp(dec.x_idcs{g});
    if any(cellfun(@ischar,diffs_g))
        alphas_g = 'all';
    else
        alphas_g = vertcat(diffs_g{:});
    end
    try
        [sol(g).Tops,sol(g).Kops,sol(g).Fop,sol(g).info,sol(g).Gops] = solve_BVP_group(PDEs{g},alphas_g,any(dec.w_idcs{g}>nw0));
    catch ME
        if strcmp(ME.identifier,'compose_Tmaps:inputsRankDeficient') && any(dec.w_idcs{g}>nw0)
            error('solve_PIETOOLS_BVP:sourceRankDeficient','The boundary conditions %s involve states of lower dimension, in a direction that is rank deficient (e.g. periodic); this is not supported.',mat2str(dec.BC_idcs{g}(:)'))
        end
        rethrow(ME)
    end
end

% % % Assemble the operators of the groups into those of the PDE
% The fundamental state of the PDE is xf = [xf_1;...;xf_G],
% stacking those of the groups, and the rows of the states of a group are
% zero in the columns of the other groups. The BCs, states and inputs of
% the groups are replaced by those of the PDE (using dec) in info.groups
% and info.rhs, so that the right-hand sides of the BCs in info.rhs are
% expressed in the inputs w,u of the PDE. A 'copvar' needs a block in each
% row and column to define its spaces, so the blocks that no group defines
% are zero operators, which are only declared where needed for this.
nG = numel(PDEs);
% Zero operator with nr rows, mapping from the space s (fields vars, dom,
% dim) to the space with variables vars_out on domain dom_out
zero_op = @(nr,s,vars_out,dom_out) mat2sopvar(zeros(nr,s.dim),struct('in',{s.vars},'out',{vars_out}),struct('in',s.dom,'out',dom_out));

% % % Info on the groups, with the indices of the PDE instead of the group
info = struct('groups',struct('BC_idcs',cell(1,nG),'vars',[],'pass',[]));
fund = cell(1,nG);
for g=1:nG
    info.groups(g) = struct('BC_idcs',dec.BC_idcs{g},'vars',{sol(g).info.groups.vars},'pass',{sol(g).info.groups.pass});
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
% Whether a fundamental input f is the input r of the PDE itself, with table tab: all its
% elements, without derivative or boundary value, on its own variables
is_itself = @(f,tab,r) f.pde==r && ~any(f.Dorder) && ~any(f.loc) && isequal(sort(f.sel(:)'),1:tab(r,2)) && isequal(f.vars,vars_pde(logical(tab(r,3:2+N)))); % DJ (CC), 11/08/2026
% The fundamental input that is the input r of the PDE itself (group 0), see is_itself
own_input = @(tab,r) struct('vars',{vars_pde(logical(tab(r,3:2+N)))},'dom',PDE.dom(logical(tab(r,3:2+N)),:),'dim',tab(r,2),'pde',r,'Dorder',zeros(1,N),'loc',zeros(1,N),'sel',(1:tab(r,2))','group',0); % DJ (CC), 11/08/2026
Rops = struct('x',[],'w',[],'u',[]);
Kops = struct('K',[]);                      Gops = struct('w',[],'u',[]);
info.inputs = struct('w',[],'u',[]);        info.rhs = struct('w',[],'u',[]);
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
            info.rhs.(obj_type) = struct('dir',{},'vars',{},'dom',{},'nrow',{},'rows',{},'group',{},'BC_idcs',{},'forcing',{});
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
                        error('solve_PIETOOLS_BVP:sharedInputs','The input %s%d appears in the boundary conditions of states on different variables, in ways that cannot be written as the same fundamental input; this is not supported.',obj_type,f1.pde)
                    end
                end
            end
            colmap = arrayfun(@(g) col_id(col_group==g)',1:nG,'UniformOutput',false);
            col_spaces = col_spaces(first_idx);
            ncols_tot = numel(first_idx);
        end
        % % % Start of changes DJ (CC), 11/08/2026
        % The PDE inputs lifted into the sources (lift_BVP_sources) are fundamental
        % inputs too, with zero columns in the maps of the groups: the existing
        % fundamental input if that is the PDE input itself (as for an input also in
        % another BC), else one of its own (group 0, like those not in the BCs below)
        tab_in = PDE.([obj_type,'_tab']);
        for lf = [struct('type',{},'idx',{},'C',{}),src.inputs]
            if ~strcmp(lf.type,obj_type) || any(arrayfun(@(f) is_itself(f,tab_in,lf.idx),col_spaces))
                continue
            elseif any([col_spaces.pde]==lf.idx)
                error('solve_PIETOOLS_BVP:sourceAndInput','An input in the boundary conditions of a state of lower dimension is also decomposed in other boundary conditions, in a way that is not the input itself; this is not supported.')
            end
            col_spaces(end+1) = own_input(tab_in,lf.idx);   %#ok<AGROW>
            ncols_tot = ncols_tot+1;
        end
        % % % End of changes DJ (CC), 11/08/2026
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
            if size(sol(g).Kops.x,1)>0
                K_blocks(end+(1:size(sol(g).Kops.x.C,1)),off(g)+(1:ncols_g(g))) = sol(g).Kops.x.C;
            end
        end
        if ~isempty(K_blocks)
            blk_ref = K_blocks{1,find(~cellfun(@isempty,K_blocks(1,:)),1)};
            for c=find(all(cellfun(@isempty,K_blocks),1))
                K_blocks{1,c} = zero_op(blk_ref.dims(1),col_spaces(c),blk_ref.vars.out,blk_ref.dom.out);
            end
            Kops.K = copvar(K_blocks);
        end
    else
        % The BCs, and the terms of the inputs on their right-hand sides,
        % with the indices of the BCs and inputs of the PDE
        rhs = cell(1,nG);
        for g=find(ncols_g>0)
            rhs{g} = sol(g).info.rhs.(obj_type);
            for c=1:numel(rhs{g})
                rhs{g}(c).group = g;
                rhs{g}(c).BC_idcs = dec.BC_idcs{g}(rhs{g}(c).BC_idcs);
                for kk=1:numel(rhs{g}(c).forcing)
                    for jj=1:numel(rhs{g}(c).forcing{kk})
                        forcing_term = rhs{g}(c).forcing{kk}{jj};
                        forcing_term.(obj_type) = dec.([obj_type,'_idcs']){g}(forcing_term.(obj_type));
                        rhs{g}(c).forcing{kk}{jj} = forcing_term;
                    end
                end
            end
        end
        info.rhs.(obj_type) = [rhs{:}];
        % G: the maps to the right-hand sides of the BCs, with a row for
        % each element of rhs, zero for the inputs of the other groups
        nrows_G = zeros(1,nG);
        for g=1:nG
            if ~isempty(sol(g).Gops.(obj_type))
                nrows_G(g) = size(sol(g).Gops.(obj_type).C,1);
            end
        end
        if any(nrows_G)
            G_off = [0,cumsum(nrows_G)];
            G_blocks = cell(G_off(end),ncols_tot);
            for g=find(nrows_G>0)
                G_blocks(G_off(g)+(1:nrows_G(g)),colmap{g}) = sol(g).Gops.(obj_type).C;
            end
            rhs_all = [rhs{nrows_G>0}];
            for r=1:size(G_blocks,1)
                for c=find(cellfun(@isempty,G_blocks(r,:)))
                    G_blocks{r,c} = zero_op(numel(rhs_all(r).rows),col_spaces(c),rhs_all(r).vars,rhs_all(r).dom);
                end
            end
            Gops.(obj_type) = copvar(G_blocks);
        end
    end
end
% Eliminate the inputs that represent states of lower dimension, by
% composing with the maps of these states
if ~isempty(src)
    [Rops,Gops,info] = eliminate_BVP_sources(PDE,src,Rops,Gops,info);
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

        % % % The inputs that are not in the BCs (those in the right-hand
        % % % sides of the BCs were decomposed into fundamental inputs
        % % % above), as additional fundamental inputs with zero columns in
        % % % the maps to the state
        % These inputs do not enter the BCs, so they do not affect x
        % through R; they are their own fundamental input (group 0),
        % listed after the decomposed ones, and only appear in the
        % equations.
        % % % Start of changes DJ (CC), 11/08/2026
        % The inputs without a fundamental input: not in the BCs, nor lifted into
        % a source (was: not in the right-hand sides info.rhs of the BCs, which gave
        % the same inputs, except for those lifted, that have a fundamental input)
        inter = setdiff(1:nin,[inputs.pde]);
        % % % End of changes DJ (CC), 11/08/2026
        for r=inter
%             dep_mask = logical(tab(r,3:2+N));                             % DJ (CC), 11/08/2026 (was)
%             inputs(end+1) = struct('vars',{vars_pde(dep_mask)},'dom',PDE.dom(dep_mask,:),'dim',tab(r,2),'pde',r,'Dorder',zeros(1,N),'loc',zeros(1,N),'sel',(1:tab(r,2))','group',0);   %#ok<AGROW> % DJ (CC), 11/08/2026 (was)
            inputs(end+1) = own_input(tab,r);                               % DJ (CC), 11/08/2026
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
        % on its own fundamental input, and zero elsewhere (empty for
        % inputs decomposed in the BCs, other than as the input itself, which
        % cannot be used elsewhere, see convert_PDE_term)
        Iin = cell(nin,1);
        % % % Start of changes DJ (CC), 11/08/2026
        % The identity on the fundamental input if that is the only one of w_r and is w_r
        % itself: for the inputs of their own, also those lifted into sources, and for an
        % input in the BCs that is itself fundamental, which may then also be used in the
        % equations (was: only the inputs of their own, Iin{r} empty for those in the BCs)
%         for p=1:numel(inter)                                              % DJ (CC), 11/08/2026 (was)
%             r = inter(p);                                                 % DJ (CC), 11/08/2026 (was)
        for r=1:nin
            c_r = find([inputs.pde]==r);
            if numel(c_r)~=1 || ~is_itself(inputs(c_r),tab,r)
                continue
            end
            dep_mask = logical(tab(r,3:2+N));
            Mat = zeros(tab(r,2),sum([inputs.dim]));
%             Mat(:,sum([inputs(1:n_fund+p-1).dim])+(1:tab(r,2))) = eye(tab(r,2)); % DJ (CC), 11/08/2026 (was)
            Mat(:,sum([inputs(1:c_r-1).dim])+(1:tab(r,2))) = eye(tab(r,2));     % DJ (CC), 11/08/2026
            Iin{r} = mat2copvar_sop(Mat,struct('out',tab(r,2),'in',[inputs.dim]'),struct('out',{{vars_pde(dep_mask)}},'in',{{inputs.vars}}),dom_pde);
        end
        % % % End of changes DJ (CC), 11/08/2026
        Rops.(['I',obj_type]) = Iin;
        info.inputs.(obj_type) = inputs;
    end
end

% % % The maps to the states: single operators (the first of the derivatives
% % % of each state, stacked) if no derivatives were requested, or the cells
% % % Rops of the maps to the requested derivatives
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
    % Maps of the fundamental inputs to the BC right-hand sides (Gops), and
    % to the PDE inputs (Iw, Iu)
    for types={'w','u'}
        obj_type = types{1};
        if ~isempty(Gops.(obj_type))
            Gops.(obj_type) = split_op(Gops.(obj_type),[],[]);
        end
        % The maps to the PDE inputs are [] for the inputs in the BCs
        if isfield(Rops,['I',obj_type])
            for r=find(~cellfun(@isempty,Rops.(['I',obj_type]))')
                Rops.(['I',obj_type]){r} = split_op(Rops.(['I',obj_type]){r},[],[]);
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
        if ~isstruct(comps)
            comps = struct('dim',{},'vars',{});
        end
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
    if ~isstruct(inputs)
        inputs = struct('dim',{},'vars',{});
    end
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

% % % The maps from the fundamental inputs to the inputs of the PDE
if isfield(Rops,'Iw')
    Gops.Iw = Rops.Iw;      Gops.Iu = Rops.Iu;
end
% Whether the inputs enter the states, or Tops.w, Tops.u are the zero operators below
info.input_maps = struct('w',~isempty(Tops.w),'u',~isempty(Tops.u));        % DJ (CC), 11/08/2026
% % % The constraints Kx*xf + Kw*wf + Ku*uf = 0, as operators of the class
% % % of Tops with a row for each constraint: Kops.x is K, or has no rows if
% % % there are no constraints, and the inputs do not enter the constraints
% % % (yet)
% % % Kops.w and Kops.u are zero operators with the rows of Kops.x. Tops.w
% % % and Tops.u are zero operators with the rows of Tops.x, if there are no
% % % inputs in the BCs. The columns are the components of space.w, space.u
if use_nonlin
    [vars_srt,ord] = sort(vars_pde);
    Kx = Kops.K;
    Kops = struct('x',Kx,'w',[],'u',[]);
    if isempty(Kx)
        % No constraints: no rows
        Kx = struct('dims',struct('out',zeros(0,1)),'vars',struct('out',{cell(1,0)}),'dom',struct('out',zeros(0,2)),'depmat',struct('out',zeros(0,0)));
        Kops.x = mat2tensopmat(zeros(0,sum(space.x.dims)),struct('in',space.x.dims,'out',zeros(0,1)),struct('in',{vars_srt},'out',{cell(1,0)}),...
            struct('in',space.dom(ord,:),'out',zeros(0,2)),struct('in',double(space.x.depmat(:,ord)),'out',zeros(0,0)));
    end
    for ty={'w','u'}
        sp = space.(ty{1});
        in_dom = space.dom(ord,:);      in_dep = double(sp.depmat(:,ord));
        Kops.(ty{1}) = mat2tensopmat(zeros(sum(Kx.dims.out),sum(sp.dims)),struct('in',sp.dims,'out',Kx.dims.out),...
            struct('in',{vars_srt},'out',{Kx.vars.out}),struct('in',in_dom,'out',Kx.dom.out),struct('in',in_dep,'out',Kx.depmat.out));
        if ~isempty(Tops.(ty{1}))
            continue
        end
        if ~want_diffs
            P = Tops.x;
            Tops.(ty{1}) = mat2tensopmat(zeros(sum(P.dims.out),sum(sp.dims)),struct('in',sp.dims,'out',P.dims.out),...
                struct('in',{vars_srt},'out',{P.vars.out}),struct('in',in_dom,'out',P.dom.out),struct('in',in_dep,'out',P.depmat.out));
            continue
        end
        Tops.(ty{1}) = cell(size(Tops.x));
        for k=1:numel(Tops.x)
            Tops.(ty{1}){k} = cell(size(Tops.x{k}));
            for l=find(~cellfun(@isempty,Tops.x{k}(:)))'
                P = Tops.x{k}{l};
                Tops.(ty{1}){k}{l} = mat2tensopmat(zeros(sum(P.dims.out),sum(sp.dims)),struct('in',sp.dims,'out',P.dims.out),...
                    struct('in',{vars_srt},'out',{P.vars.out}),struct('in',in_dom,'out',P.dom.out),struct('in',in_dep,'out',P.depmat.out));
            end
        end
    end
else
    Kx = Kops.K;
    Kops = struct('x',Kx,'w',[],'u',[]);
    if isempty(Kx)
        % No constraints: no rows
        Kx = struct('dim_out',zeros(0,1),'vars',{cell(1,0)},'space_out',false(0,0));
        Kops.x = zeros_copvar_sop(struct('out',zeros(0,1),'in',space.x.dims),struct('out',{cell(1,0)},'in',{arrayfun(@(i) vars_pde(space.x.depmat(i,:)),1:numel(space.x.dims),'UniformOutput',false)}),struct('vars',{vars_pde},'dom',space.dom));
    end
    for ty={'w','u'}
        sp = space.(ty{1});
        in_spaces = arrayfun(@(i) vars_pde(sp.depmat(i,:)),1:numel(sp.dims),'UniformOutput',false);
        out_spaces = arrayfun(@(i) Kx.vars(Kx.space_out(i,:)),1:numel(Kx.dim_out),'UniformOutput',false);
        Kops.(ty{1}) = zeros_copvar_sop(struct('out',Kx.dim_out,'in',sp.dims),struct('out',{out_spaces},'in',{in_spaces}),struct('vars',{vars_pde},'dom',space.dom));
        if ~isempty(Tops.(ty{1}))
            continue
        end
        if ~want_diffs
            P = Tops.x;
            out_spaces = arrayfun(@(i) P.vars(P.space_out(i,:)),1:numel(P.dim_out),'UniformOutput',false);
            Tops.(ty{1}) = zeros_copvar_sop(struct('out',P.dim_out,'in',sp.dims),struct('out',{out_spaces},'in',{in_spaces}),struct('vars',{vars_pde},'dom',space.dom));
            continue
        end
        Tops.(ty{1}) = cell(size(Tops.x));
        for k=1:numel(Tops.x)
            Tops.(ty{1}){k} = cell(size(Tops.x{k}));
            for l=find(~cellfun(@isempty,Tops.x{k}(:)))'
                P = Tops.x{k}{l};
                out_spaces = arrayfun(@(i) P.vars(P.space_out(i,:)),1:numel(P.dim_out),'UniformOutput',false);
                Tops.(ty{1}){k}{l} = zeros_copvar_sop(struct('out',P.dim_out,'in',sp.dims),struct('out',{out_spaces},'in',{in_spaces}),struct('vars',{vars_pde},'dom',space.dom));
            end
        end
    end
end
info.space = space;
info.Gops = Gops;
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
function [Tops,Kops,Fop,info,Gops] = solve_BVP_group(PDE,alphas,sources)
% [Tops,Kops,Fop,info,Gops] = SOLVE_BVP_GROUP(PDE,alphas,sources)
% computes the operators parameterizing all states satisfying the BCs of a
% single group of coupled BCs, such as the groups returned by
% decouple_BVPs: the BCs are split into 1D problems along each direction
% (split_BVPs), which are solved using solve_BVP_1D, and their solution
% operators are composed into ND operators by compose_Tmaps.
%
% INPUT
% - PDE:    initialized 'pde_struct' object defining a PDE in the terms
%           format, whose states all depend on the same variables. Only the
%           BCs are used;
% - alphas: the derivative orders to compute for Tops, as in compose_Tmaps
%           ('all' for all);
% - sources: true if the inputs in the BCs are derived from states of lower
%           dimension, see compose_Tmaps;
%
% OUTPUT
% - Tops,Kops,Fop,Gops: operators as in compose_Tmaps;
% - info:   struct with fields groups, fund, inputs, rhs as in
%           solve_PIETOOLS_BVP, for a single group (group index 1, states
%           x_idcs=1:nx). The fields inputs and rhs have subfields w and u
%           for the exogenous and actuator inputs, and rhs.w and rhs.u
%           have the BCs and the terms of w, resp. u, of the right-hand
%           sides of the BCs along each direction with such inputs,
%           BC_idcs and forcing.
% % % Extract the sizes, variables and orders of differentiability of the
% % % states of the group (Delta(k,i): order of x_k in the i-th variable)
nx = numel(PDE.x);
N = size(PDE.vars,1);
sizes = PDE.x_tab(:,2)';
isvar = logical(PDE.x_tab(:,3:2+N));
Delta = PDE.x_tab(:,2+N+(1:N));
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);
dir_idcs = 1:N;
if isempty(PDE.BC)
    % States without BCs must not be differentiable, and have xf=u
    if any(Delta(:)>0)
        error('solve_PIETOOLS_BVP:noBCs','A state that is differentiable in some direction has no boundary conditions along it.')
    end
    dir_idcs = find(any(isvar,1));
    maps = repmat(struct('Tops',[],'Kops',[],'Fop',[],'info',[],'Qops',[],'forced',[],'source',[]),1,numel(dir_idcs));
else
    % Split into 1D problems along each direction, with the BCs on the
    % states as functions of s_i, and the other variables as parameters,
    % and solve each of them
    [PDEs_1D,info_1D] = split_BVPs(PDE);
    % The states must all depend on the same variables, so that the
    % right-hand sides of the BCs (forcing) only involve inputs
    if ~all(isvar(:))
        error('solve_PIETOOLS_BVP:mixedDimension','The states in a group of coupled boundary conditions should all depend on the same variables.')
    end
    maps = repmat(struct('Tops',[],'Kops',[],'Fop',[],'info',[],'Qops',[],'forced',[],'source',[]),1,N);
    for i=1:N
        if max(Delta(:,i))==0 && isempty(PDEs_1D{i}.BC)
            continue                                % no derivatives, so no BCs
        end
        % forced.w(r) (forced.u) is true if row r of the BCs along
        % direction i has inputs w (u) on its right-hand side; the maps of
        % these (sixth output of solve_BVP_1D) are only needed then
        nrow = cellfun(@(c) c.size,PDEs_1D{i}.BC);
        for ty={'w','u'}
            has_input = cellfun(@(f) iscell(f) && any(cellfun(@(x) isfield(x,ty{1}),f)),info_1D.forcing{i});
            maps(i).forced.(ty{1}) = repelem(has_input,nrow(:)');
            % The element of an input that is the right-hand side of each
            % condition (src.pde=0 if it is not one element of one input)
            obj_type = ty{1};
            tab = PDE.([obj_type,'_tab']);
            face = true(1,N);           face(i) = false;
            src = struct('pde',zeros(sum(nrow),1),'elem',zeros(sum(nrow),1),'scale',ones(sum(nrow),1),'dep',false(sum(nrow),N));
            off = 0;
            for k=1:numel(nrow)
                terms = info_1D.forcing{i}{k};
                if iscell(terms)
                    terms = terms(cellfun(@(tm) isfield(tm,obj_type),terms));
                else
                    terms = {};
                end
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
        % Errors are prefixed with the direction, to locate the BCs. The 1D
        % problem is solved for the requested orders of the derivatives in
        % this direction, with input maps if some BC has inputs
        diffs_i = 'all';
        if ~ischar(alphas)
            diffs_i = 0;
            if ~isempty(alphas)
                diffs_i = unique([0;alphas(:,i)])';
            end
        end
        try
            if any(maps(i).forced.w) || any(maps(i).forced.u)
                [maps(i).Tops,maps(i).Kops,maps(i).Fop,maps(i).info,maps(i).Qops] = solve_BVP_1D(PDEs_1D{i},diffs_i);
            else
                [maps(i).Tops,maps(i).Kops,maps(i).Fop,maps(i).info] = solve_BVP_1D(PDEs_1D{i},diffs_i);
                maps(i).forced = [];
            end
        catch ME
            msg = ['Boundary conditions along ',vars_pde{i},': ',ME.message];
            if isempty(ME.identifier)
                error('solve_PIETOOLS_BVP:direction','%s',msg)
            else
                error(ME.identifier,'%s',msg)
            end
        end
    end
end

% % % Compose the solutions along the different directions: each maps the
% % % 1D fundamental state to the states and their derivatives in that
% % % direction
[Tops,Kops,Fop,info_comp,Gops] = compose_Tmaps(maps,vars_pde(dir_idcs),PDE.dom(dir_idcs,:),sizes,Delta(:,dir_idcs),alphas,sources);
% Info on the group, in terms of the states and BCs of the group; the
% indices are mapped to those of the PDE in solve_PIETOOLS_BVP
info = struct('groups',struct('BC_idcs',1:numel(PDE.BC),'vars',{info_comp.vars},'pass',{info_comp.vars(info_comp.pass)}));
info.fund = info_comp.fund;
[info.fund.group] = deal(1);
[info.fund.x_idcs] = deal(1:nx);

% % % The inputs, and the right-hand sides of the BCs they enter in, with
% % % the BCs of each direction and the terms of this type on their
% % % right-hand sides (the BCs may have terms of both w and u)
% Label the inputs and right-hand sides with their group and forcing terms
info.inputs = info_comp.inputs;     info.rhs = info_comp.rhs;
for ty={'w','u'}
    obj_type = ty{1};
    if ~isempty(info.inputs.(obj_type))
        [info.inputs.(obj_type).group] = deal(1);
    end
    for c=1:numel(info.rhs.(obj_type))
        info.rhs.(obj_type)(c).group = 1;
        info.rhs.(obj_type)(c).BC_idcs = info_1D.BC_idcs{info.rhs.(obj_type)(c).dir};
        forcing_terms = info_1D.forcing{info.rhs.(obj_type)(c).dir};
        for kk=1:numel(forcing_terms)
            forcing_terms{kk} = forcing_terms{kk}(cellfun(@(x) isfield(x,obj_type),forcing_terms{kk}));
        end
        info.rhs.(obj_type)(c).forcing = forcing_terms;
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

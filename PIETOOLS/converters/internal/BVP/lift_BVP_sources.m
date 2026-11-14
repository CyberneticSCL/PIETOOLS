function [PDE_out,src] = lift_BVP_sources(PDE)
% [PDE_out,src] = LIFT_BVP_SOURCES(PDE) replaces the terms in the BCs of
% PDE that involve a state of lower dimension than another state in the
% same BC (e.g. ODE states, or boundary values of states on other
% variables) by a term of a new exogenous input, so that each group of
% coupled BCs only involves states on the same variables. E.g. the BC
%       x2(s1,0) = x1(s1)
% of x2(s1,s2) becomes x2(s1,0) = p(s1), for a new input p (the "source"
% x1). The BC is then solved as a BC of x2 with an input, and p is
% eliminated afterwards, see eliminate_BVP_sources.
%
% INPUT
% - PDE:    initialized 'pde_struct' object in the terms format;
%
% OUTPUT
% - PDE_out:    initialized 'pde_struct' object with the states and BCs
%           of PDE (no equations), in which the lifted terms are replaced
%           by terms of the new inputs w_{numel(PDE.w)+r}, r=1,...,
%           numel(src) (or dropped, if their coefficient is zero). PDE
%           itself if no term is lifted;
% - src:    struct array, defining each new input p_r as the sum of the
%           terms of one BC, with fields
%           - BC:     the index of the BC in PDE;
%           - terms:  struct array with for each term fields state (the
%                     index of the state), D (orders of the derivative in
%                     each of its p variables), loc (1 x p 'polynomial',
%                     where it is evaluated), C (dim x size of the state
%                     'polynomial' coefficient) and I (p x 1 cell, the
%                     limits of the integral over each variable, [] if
%                     none), so that p_r = sum_{terms} C*D^{D}x_state(loc);
%           - vars, dom, dim:   the variables (pvar array, those of the
%                     terms), domain and size of p_r;
%           - inputs: struct array with the input terms of the BC lifted
%                     into the source, with fields type ('w' or 'u'), idx
%                     (the input) and C (constant), adding sum C*w_idx to
%                     p_r; empty if none, see NOTES;
%
% NOTES
% - A term of a state is of lower dimension if the variables of another
%   state in the BC strictly contain its variables. Integrals or delays
%   of such a term are not supported.
% - A term that fixes several variables of its state (a corner value, or
%   an integral over the full domain in several variables) cannot be a
%   condition along a single direction (see split_BVPs), so it is lifted
%   whatever the dimension of its state, if the BC involves other states:
%   e.g. x1(0)=int_{0}^{1}int_{0}^{1}x2 is a BC of x1(s1) with the
%   source int x2 of x2(s1,s2), constant in the integrated variables.
% - Of the other terms, the first is that of the constrained state. A
%   term of a state on other variables is also lifted if it only depends
%   on variables of the constrained state, e.g. x4(a)=x3(d) for states
%   x4(s1), x3(s2).
% - The terms lifted from a BC form one input, for the rows of the BC
%   that they enter, so that each row has a single element of a single
%   input as right-hand side. If the source has no variables and the
%   constrained state has one, the input terms of the BC (of inputs
%   without variables, with constant coefficients) are lifted into the
%   source too; other inputs stay in the BCs, which is not supported in
%   a group with a source on several variables (convert_PIETOOLS_BVP).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lift_BVP_sources
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
% DJ (CC), 11/08/2026: Lift the input terms of a BC into a scalar source of a
%                   1-D state, so that a row with a state of lower dimension
%                   and an input has a single input as right-hand side;
% DJ, 11/11/2026: Rename solve_PIETOOLS_BVP to convert_PIETOOLS_BVP;
% DJ (CC), 11/13/2026: Lift terms fixing several variables of a state;
%                   drop lifted terms with a zero coefficient; return PDE
%                   if nothing is lifted; shorter error messages;

N = size(PDE.vars,1);
nw = numel(PDE.w);
isvar = logical(PDE.x_tab(:,3:2+N));
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);
has_int = @(term) isfield(term,'I') && any(~cellfun(@isempty,term.I));
has_delay = @(term) isfield(term,'delay') && any(double(polynomial(term.delay(:))));
% Rows of the matrix polynomial C with a nonzero entry
rows_of = @(C) reshape(any(reshape(any(full(C.coeff)~=0,1),size(C)),2),[],1);

% % % Find the terms of lower dimension in each BC, and replace them
src = struct('BC',{},'terms',{},'vars',{},'dom',{},'dim',{},'inputs',{});
BC = PDE.BC;
changed = false;
for k=1:numel(BC)
    terms = BC{k}.term;
    state_terms = find(cellfun(@(term) isfield(term,'x'),terms));      state_terms = state_terms(:)';
    nst = numel(state_terms);
    isvar_states = isvar(cellfun(@(term) term.x,terms(state_terms)),:);          % variables of the states

    % For each state term, the variables in which it is a function (not
    % evaluated at a boundary or integrated), and the limits of its
    % integrals, of which those over the full domain, with the variables
    % that it fixes (see split_BVPs). A term fixing several variables is a
    % scalar functional of its state, which cannot be a condition along a
    % single direction. If the BC has another state, the BC constrains
    % that state, and the term is lifted, whatever the dimension of its
    % state.
    free_states = false(nst,N);     is_multi = false(nst,1);
    I_term = cell(nst,1);           full_int = cell(nst,1);
    for ii=1:nst
        term = terms{state_terms(ii)};
        vidx = find(isvar_states(ii,:));
        loc = polynomial(term.loc);
        is_bdry = arrayfun(@(m) isdouble(loc(m)),1:numel(vidx));
        I_term{ii} = cell(numel(vidx),1);       full_int{ii} = false(numel(vidx),1);
        for m=1:numel(vidx)
            if isfield(term,'I') && numel(term.I)>=m && ~isempty(term.I{m})
                I_term{ii}{m} = term.I{m};
                I_lim = polynomial(term.I{m});
                full_int{ii}(m) = isdouble(I_lim) && numel(I_lim)==2 && isequal(sort(double(I_lim(:))'),sort(PDE.dom(vidx(m),:)));
            end
        end
        free_states(ii,vidx) = ~is_bdry;
        if has_int(term)
            free_states(ii,:) = false;
        end
        is_multi(ii) = nnz(is_bdry(:) | full_int{ii})>1;
    end
    if all(is_multi)
        is_multi(:) = false;        % no state to constrain, left to split_BVPs
    end
    % Lower dimension: the variables are strictly contained in those of
    % another state in the BC that can be constrained by it
    is_lower = false(nst,1);
    for ii=1:nst
        is_lower(ii) = is_multi(ii) || any(all(isvar_states(~is_multi,:)>=isvar_states(ii,:),2) & any(isvar_states(~is_multi,:)>isvar_states(ii,:),2));
    end
    % The first of the other terms is that of the constrained state
    ii_main = find(~is_lower,1);
    lift_ii = zeros(1,0);
    for ii=1:nst
        term = terms{state_terms(ii)};
        % A state on other variables than the constrained state, whose term
        % is a function of variables of the term of the constrained state
        % only, is also of lower dimension (e.g. a boundary value, as in
        % x4(a)=x3(d))
        is_other = ~is_lower(ii) && ii~=ii_main && ~isequal(isvar_states(ii,:),isvar_states(ii_main,:)) && ~any(free_states(ii,:) & ~free_states(ii_main,:)) && ~has_int(term) && ~has_delay(term);
        if ~is_lower(ii) && ~is_other
            continue
        end
        % Integrals are only supported over the full domain, in a term that
        % fixes several variables
        if has_delay(term)
            error("convert_PIETOOLS_BVP:sourceTerm","Delays of a state of lower dimension in a boundary condition are not supported.") % DJ, 11/11/2026
        elseif has_int(term) && ~(is_multi(ii) && all(full_int{ii} | cellfun(@isempty,I_term{ii})))
            error("convert_PIETOOLS_BVP:sourceTerm","Integrals of a state of lower dimension in a boundary condition are only supported over its full domain, in terms that fix two or more variables.")
        end
        lift_ii(end+1) = ii;                                                    %#ok<AGROW>
    end
    if isempty(lift_ii)
        continue
    end
    changed = true;

    % The new input is the sum of the lifted terms, for the rows of the BC
    % that they enter in, a function of the variables in which the states
    % are not evaluated at a boundary or integrated, and of those of the
    % coefficients. Terms with a zero coefficient are dropped.
    nrow = BC{k}.size;
    rows_src = false(nrow,1);
    vars_src = false(1,N);
    st = struct('state',{},'D',{},'loc',{},'C',{},'I',{});
    for ii=lift_ii
        term = terms{state_terms(ii)};
        C = polynomial(term.C);
        rows_C = rows_of(C);
        if ~any(rows_C)                                                     % DJ (CC), 11/13/2026
            continue
        end
        vidx = find(isvar_states(ii,:));
        loc = polynomial(term.loc);
        D = zeros(1,numel(vidx));
        if isfield(term,'D') && ~isempty(term.D)
            D = term.D(:)';
        end
        % The source is constant in the variables over which the state is
        % integrated, also if the coefficient depends on them
        free = arrayfun(@(m) ~isdouble(loc(m)),1:numel(vidx)) & ~full_int{ii}(:)';
        vars_src(vidx(free)) = true;
        vars_src(ismember(vars_pde,setdiff(C.varname,vars_pde(vidx(full_int{ii}))))) = true;
        rows_src = rows_src | rows_C;
        st(end+1) = struct('state',term.x,'D',D,'loc',loc,'C',C,'I',{I_term{ii}});      %#ok<AGROW>
    end
    if isempty(st)                  % only zero terms: no source
        terms(state_terms(lift_ii)) = [];
        BC{k}.term = terms;
        continue
    end
    % The input terms are lifted into a scalar source of a 1-D state
    % (inputs without variables, constant coefficients), see NOTES. Then
    % each row of the BC has a single input as right-hand side, the
    % source; the lifted inputs enter the source, and are eliminated with
    % it (eliminate_BVP_sources). Their rows are rows of the source.
    inp = struct('type',{},'idx',{},'C',{});    in_idx = zeros(1,0);
    if ~any(vars_src) && nnz(isvar_states(ii_main,:))==1
        for jj=1:numel(terms)
            tm = terms{jj};
            if numel(tm)~=1 || ~(isfield(tm,'w') || isfield(tm,'u')) || has_int(tm) || has_delay(tm)
                continue
            end
            ty = 'w';   if isfield(tm,'u'), ty = 'u'; end
            Cj = polynomial(tm.C);
            if any(PDE.([ty,'_tab'])(tm.(ty),3:2+N)) || any(full(Cj.degmat(:)))
                continue                    % input with variables, or coefficient of a variable
            end
            inp(end+1) = struct('type',ty,'idx',tm.(ty),'C',Cj);                        %#ok<AGROW>
            in_idx(end+1) = jj;                                                         %#ok<AGROW>
            rows_src = rows_src | rows_of(Cj);
        end
    end
    rows = find(rows_src);
    for jt=1:numel(st)
        st(jt).C = st(jt).C(rows,:);
    end
    for ji=1:numel(inp)
        inp(ji).C = inp(ji).C(rows,:);
    end
    vidx_src = find(vars_src);
    src(end+1) = struct('BC',k,'terms',st,'vars',PDE.vars(vidx_src,1),'dom',reshape(PDE.dom(vidx_src,:),[],2),'dim',numel(rows),'inputs',inp); %#ok<AGROW>
    E = zeros(nrow,numel(rows));
    E(sub2ind(size(E),rows(:),(1:numel(rows))')) = 1;
    terms([state_terms(lift_ii),in_idx]) = [];
    terms{end+1} = struct('w',nw+numel(src),'C',E,'I',{cell(numel(vidx_src),1)},'delay',0);        %#ok<AGROW>
    BC{k}.term = terms;
end
PDE_out = PDE;
if ~changed
    return
end

% % % Declare the BVP: the states, inputs (old and new), and modified BCs
% copy_fields copies the fields fields of each component in the cell c
copy_fields = @(c,fields) cellfun(@(o) cell2struct(cellfun(@(f) o.(f),fields,'UniformOutput',false),fields,2),c,'UniformOutput',false);
PDE_s = struct('vars',PDE.vars,'dom',PDE.dom);
if PDE.has_delay
    PDE_s.tau = PDE.tau;
end
PDE_s.x = copy_fields(PDE.x,{'size','vars','dom','diff'});
PDE_s.w = [copy_fields(PDE.w,{'size','vars','dom'});arrayfun(@(s) struct('size',s.dim,'vars',s.vars,'dom',s.dom),src(:),'UniformOutput',false)];
PDE_s.u = copy_fields(PDE.u,{'size','vars','dom'});
PDE_s.BC = cellfun(@(b) struct('size',b.size,'term',{b.term}),BC,'UniformOutput',false);
PDE_out = initialize(pde_struct(PDE_s),true);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

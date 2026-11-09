function [PDE_out,src] = lift_BVP_sources(PDE)
% [PDE_out,src] = LIFT_BVP_SOURCES(PDE) replaces, in the BCs of PDE, the
% terms that involve a state of lower dimension than another state in the
% same BC (including ODE states, and boundary values of states on other
% variables) by a term of a new exogenous input, so that each group of
% coupled BCs only involves states that depend on the same
% variables, as required by solve_BVP_group. For example, the BC
%       x2(s1,0) = x1(s1)
% of a state x2(s1,s2), with a state x1(s1), becomes x2(s1,0) = p(s1) for a
% new input p, representing the "source" x1. The BC is then the boundary
% condition of x2 with an input, the solution of which is composed with the
% operators of x1 to eliminate p, see eliminate_BVP_sources.
%
% INPUT
% - PDE:    'pde_struct' object defining a PDE in the terms format, see
%           also "@pde_struct/initialize";
%
% OUTPUT
% - PDE_out:    'pde_struct' object with the same spatial variables, states
%           and BCs as PDE, and without equations, in which the terms of
%           states of lower dimension in the BCs are replaced by terms of
%           the new inputs w, numbered numel(PDE.w)+(1:numel(src)). Empty
%           ([]) if no BC involves a state of lower dimension;
% - src:    struct array with an element for each new input
%           r=1,...,numel(src), defining it as the sum of the terms of the
%           states of lower dimension of one BC, with fields
%           - BC:     index of the BC in PDE;
%           - terms:  struct array with for each term the fields
%                     - state:  index of the state x_state;
%                     - D:      1 x p array of the orders of the
%                               derivatives of the state in each of its p
%                               variables;
%                     - loc:    1 x p 'polynomial' array, at which the
%                               state is evaluated;
%                     - C:      the coefficient, a dim x size(x_state)
%                               'polynomial' matrix,
%                     so that p_r = sum_{terms} C*D^{D}x_state(loc);
%           - vars,dom,dim:   variables (pvar array), domain, and size of
%                     the input, the variables being those on which the
%                     terms depend;
%           - inputs: struct array with the input terms of the BC lifted
%                     into the source, with fields type ('w' or 'u'), idx
%                     (the input), and C (constant, dim x size of the
%                     input), so that p_r = sum C*D^{D}x_state(loc) +
%                     sum C*w_idx; empty if none, see NOTES;
%
% NOTES
% - A term of a state is considered to be of lower dimension if the
%   variables of another state in the same BC strictly contain its
%   variables. Terms with integrals or delays of such a state are not
%   supported.
% - Of the other terms, the first is that of the constrained state. The
%   term of a state on other variables than the constrained state is also
%   lifted, if it is a function of variables of the term of the constrained
%   state only, e.g. a boundary value, as in x4(a)=x3(d) for states x4(s1)
%   and x3(s2).
% - All terms lifted from a BC are combined in one input, for the rows of
%   the BC that they enter in, so that each row has a single element of a
%   single input as right-hand side, whatever the coefficients C of the
%   terms. The input depends on the variables of the terms and
%   coefficients.
% - If the source has no variables and the constrained state depends on one
%   variable, the input terms of the BC, of inputs without variables and
%   with constant coefficients, are also lifted into the source, so that a
%   row with a state of lower dimension and an input still has a single
%   input as right-hand side (see eliminate_BVP_sources). Without this,
%   or with more variables, such inputs stay in the BCs, which is not
%   supported in the same group as a source (solve_PIETOOLS_BVP).
%

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

PDE = pde_struct(PDE);
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
N = size(PDE.vars,1);
nw = numel(PDE.w);
isvar = logical(PDE.x_tab(:,3:2+N));
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);
has_int_delay = @(term) (isfield(term,'I') && any(~cellfun(@isempty,term.I))) || (isfield(term,'delay') && any(double(polynomial(term.delay(:)))));
% Rows of the matrix polynomial C with a nonzero entry
rows_of = @(C) reshape(any(reshape(any(full(C.coeff)~=0,1),size(C)),2),[],1); % DJ (CC), 11/08/2026

% % % Find the terms of lower dimension in each BC, and replace them
% src = struct('BC',{},'terms',{},'vars',{},'dom',{},'dim',{});             % DJ (CC), 11/08/2026 (was)
src = struct('BC',{},'terms',{},'vars',{},'dom',{},'dim',{},'inputs',{});   % DJ (CC), 11/08/2026
BC = PDE.BC;
for k=1:numel(BC)
    terms = BC{k}.term;
    state_terms = find(cellfun(@(term) isfield(term,'x'),terms));      state_terms = state_terms(:)';
    isvar_states = isvar(cellfun(@(term) term.x,terms(state_terms)),:);          % variables of the states
    % Variables on which each term is a function: those of the state, in
    % which it is not evaluated at a boundary or integrated
    free_states = false(numel(state_terms),N);
    for ii=1:numel(state_terms)
        term = terms{state_terms(ii)};
        vidx = find(isvar(term.x,:));
        loc = polynomial(term.loc);
        free_states(ii,vidx) = arrayfun(@(m) ~isdouble(loc(m)),1:numel(vidx));
        if isfield(term,'I') && any(~cellfun(@isempty,term.I))
            free_states(ii,:) = false;
        end
    end
    % Lower dimension: the variables are strictly contained in those of
    % another state in the BC
    is_lower = false(numel(state_terms),1);
    for ii=1:numel(state_terms)
        is_lower(ii) = any(all(isvar_states>=isvar_states(ii,:),2) & any(isvar_states>isvar_states(ii,:),2));
    end
    % The first of the other terms is that of the constrained state
    ii_main = find(~is_lower,1);
    lift_idx = zeros(1,0);
    for ii=1:numel(state_terms)
        term = terms{state_terms(ii)};
        % A state on other variables than the constrained state, whose term
        % is a function of variables of the term of the constrained state
        % only, is also of lower dimension (e.g. a boundary value, as in
        % x4(a)=x3(d))
        is_other = ~is_lower(ii) && ii~=ii_main && ~isequal(isvar_states(ii,:),isvar_states(ii_main,:)) && ~any(free_states(ii,:) & ~free_states(ii_main,:)) && ~has_int_delay(term);
        if ~is_lower(ii) && ~is_other
            continue
        end
        if has_int_delay(term)
            error("solve_PIETOOLS_BVP:sourceTerm","Terms with an integral or a delay of a state of lower dimension in a boundary condition are not supported.")
        end
        lift_idx(end+1) = state_terms(ii);                                         %#ok<AGROW>
    end
    if isempty(lift_idx)
        continue
    end

    % The new input is the sum of the lifted terms, for the rows of the BC
    % that they enter in, a function of the variables in which the states
    % are not evaluated at a boundary, and of those of the coefficients
    nrow = BC{k}.size;
    rows_src = false(nrow,1);
    vars_src = false(1,N);
    st = struct('state',{},'D',{},'loc',{},'C',{});
    for jj=lift_idx
        term = terms{jj};
        vidx = find(isvar(term.x,:));
        loc = polynomial(term.loc);
        D = zeros(1,numel(vidx));
        if isfield(term,'D') && ~isempty(term.D)
            D = term.D(:)';
        end
        C = polynomial(term.C);
        free = arrayfun(@(m) ~isdouble(loc(m)),1:numel(vidx));
        vars_src(vidx(free)) = true;
        vars_src(ismember(vars_pde,C.varname)) = true;
%         rows_src = rows_src | reshape(any(reshape(any(full(C.coeff)~=0,1),size(C)),2),[],1); % DJ (CC), 11/08/2026 (was)
        rows_src = rows_src | rows_of(C);                                   % DJ (CC), 11/08/2026
        st(end+1) = struct('state',term.x,'D',D,'loc',loc,'C',C);                  %#ok<AGROW>
    end
    % % % Start of changes DJ (CC), 11/08/2026
    % The input terms are lifted into a scalar source of a 1-D state (inputs
    % without variables, constant coefficients), see NOTES. Then each row of
    % the BC has a single input as right-hand side, the source; the lifted
    % inputs enter the source, and are eliminated with it, with the states
    % (eliminate_BVP_sources). Their rows are rows of the source.
    inp = struct('type',{},'idx',{},'C',{});    in_idx = zeros(1,0);
    if ~any(vars_src) && nnz(isvar_states(ii_main,:))==1
        for jj=1:numel(terms)
            tm = terms{jj};
            if numel(tm)~=1 || ~(isfield(tm,'w') || isfield(tm,'u')) || has_int_delay(tm)
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
    % % % End of changes DJ (CC), 11/08/2026
    rows = find(rows_src);
    for jt=1:numel(st)
        st(jt).C = st(jt).C(rows,:);
    end
    for ji=1:numel(inp)                                                     % DJ (CC), 11/08/2026
        inp(ji).C = inp(ji).C(rows,:);                                      % DJ (CC), 11/08/2026
    end                                                                     % DJ (CC), 11/08/2026
    vidx_src = find(vars_src);
%     src(end+1) = struct('BC',k,'terms',st,'vars',PDE.vars(vidx_src,1),'dom',reshape(PDE.dom(vidx_src,:),[],2),'dim',numel(rows)); %#ok<AGROW> % DJ (CC), 11/08/2026 (was)
    src(end+1) = struct('BC',k,'terms',st,'vars',PDE.vars(vidx_src,1),'dom',reshape(PDE.dom(vidx_src,:),[],2),'dim',numel(rows),'inputs',inp); %#ok<AGROW> % DJ (CC), 11/08/2026
    E = zeros(nrow,numel(rows));
    E(sub2ind(size(E),rows(:),(1:numel(rows))')) = 1;
%     terms(lift_idx) = [];                                                 % DJ (CC), 11/08/2026 (was)
    terms([lift_idx,in_idx]) = [];                                          % DJ (CC), 11/08/2026
    terms{end+1} = struct('w',nw+numel(src),'C',E,'I',{cell(numel(vidx_src),1)},'delay',0);        %#ok<AGROW>
    BC{k}.term = terms;
end
if isempty(src)
    PDE_out = [];
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

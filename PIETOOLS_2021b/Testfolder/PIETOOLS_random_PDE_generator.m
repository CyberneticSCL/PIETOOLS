function PDE = PIETOOLS_random_PDE_generator(dim,max_diff,n_comps,max_size,max_terms,max_deg,add_BC_terms)
% PDE = PIETOOLS_random_PDE_generator(dim,max_diff,n_comps,max_size,max_terms,max_deg,add_BC_terms)
% generates a ``random'' pde_struct objecting representing a PDE in the
% PIETOOLS PDE terms format.
%
% INPUT
% - dim:        1x1 integer the spatial dimension of the system. At this
%               time, only dim=0, dim=1 or dim=2 systems can be analysed in
%               PIETOOLS. Note that the dependence of each state component
%               on each of the dim variables will be randomly set, so even
%               using dim=2, the PDE may describe only a 1D system, or even
%               an ODE.
% - max_diff:   1 x dim array of integers specifying the maximal order of 
%               differentiability of the state components in each of the
%               spatial variables that appear in the system. Note that this
%               is a maximal order of differentiability; the actual order
%               of differentiability of each state component will be
%               set randomly.
% - n_comps:    1x5 array of integers specifying how many state (1),
%               exogenous input (2), actuator input (3), regulated output
%               (4), and observed output (5) components must appear in the
%               system.
% - max_size:   1x5 array of integers specifying the maximal number of
%               variables in each state (1), exogenous input (2), actuator
%               input (3), regulated output (4) and observed output (5)
%               component that appears. Note that the actual size of the
%               components will be set randomly, with different components
%               potentially containing different numbers of variables.
% - max_terms:  1x3 array of integers specifying how many terms may appear
%               in the PDE equation for the state (1), the regulated
%               outputs (2), and the actuator outputs (3). Note that the
%               actual number of terms in each equation will be randomly
%               set.
% - max_deg:    1x1 integer specifying the maximal degree of the monomials 
%               appearing in the coefficients describing each of the terms
%               in the PDE (and output equations). Here we define the
%               degree of a monomial as the sum of the degrees of each
%               variable appearing in it.
% - add_BC_terms:   1x1 integer specifying how many terms may be added to
%                   the BCs in order to increase complexity. Simple BCs of
%                   the form 0 = partial_{s}^{k} x(s=a) or 
%                   0 = partial_{s}^{k} x(s=b) for some k below the order
%                   of differentiability of x wrt s are set initially,
%                   aiming to ensure that the PDE is well-posed, and
%                   conversion to a PIE is possible. Setting
%                   add_BC_terms=n>0, at most n additional terms MAY be 
%                   added to each BC, increasing complexity, but also
%                   increasing the chance that the system is no longer
%                   well-posed.
%
% OUTPUT:
% - PDE:    A pde_struct object, specifying a PDE in the PIETOOLS
%           terms-based format.
%
% NOTES
% Although setting add_BC_terms=0 should ensure the PDE to be well-posed,
% we do not guarantee that this indeed the case (we have not proved this).
% Certainly, adding more (exotic) terms to the BCs by setting
% add_BC_terms>0 increases the probability that the system cannot be
% converted to a PIE, especially in the 2D case.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_random_PDE_generator
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 08/23/2022
%

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% % % Extract and check the inputs.
%rng('default');
%rng(12345);

arguments
    dim {mustBeMember(dim,[0;1;2])} = 2;
    max_diff {mustBeNonnegative} = [3,3];
    n_comps {mustBeNonnegative} = [randi(3),randi([0,2]),randi([0,2]),randi([0,2]),randi([0,2])];
    max_size {mustBeNonnegative} = [randi(3),randi(2,[1,4])];
    max_terms {mustBeNonnegative} = 4;
    max_deg {mustBeNonnegative} = 2;
    add_BC_terms {mustBeNonnegative} = 0;
end

% Check that the maximal orders of differentiability have been properly
% specified.
if numel(max_diff)==1
    max_diff = max_diff*ones(1,dim);
elseif numel(max_diff)~=dim
    error('The maximal order of differentiability should be specified as a 1 x nvars array, where nvars is the number of spatial variables in the system.')
end

% Check that the maximal number of components has been properly
% specified.
if numel(n_comps)==1
    n_comps = n_comps*ones(1,5);
elseif numel(n_comps)==3
    n_comps = [n_comps(1), n_comps(2)*ones(1,2), n_comps(2)*ones(1,2)];
elseif numel(n_comps)~=5
    error('The maximal size of the components should be specified as a 1x5 array, providing sizes for [x, w, u, z, y].')
else
    n_comps = n_comps(:)';
end

% Check that the maximal size of the components has been properly
% specified.
if numel(max_size)==1
    max_size = max_size*ones(1,5);
elseif numel(max_size)==3
    max_size = [max_size(1), max_size(2)*ones(1,2), max_size(2)*ones(1,2)];
elseif numel(max_size)~=5
    error('The maximal size of the components should be specified as a 1x5 array, providing sizes for [x, w, u, z, y].')
else
    max_size = max_size(:)';
end    

% Check that the maximal number of terms in each equation has been properly
% specified.
if numel(max_terms)==1
    max_terms = max_terms*ones(1,3);
elseif numel(max_terms)==2
    max_terms = [max_terms(1), max_terms(2)*ones(1,2)];
elseif numel(max_terms)~=3
    error('The maximal size of the components should be specified as a 1x5 array, providing sizes for [x, w, u, z, y].')
else
    max_terms = max_terms(:)';
end



%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% % % Initialize a PDE, and define the PDE components
% Set variables and their domain.
nvars = dim;
pvar s1 s2 theta1 theta2
vars = [s1,theta1; s2,theta2];
vars = vars(1:nvars,:);
dom = randi(3,[2,1])-2;
dom = [dom, dom+randi(3,[2,1])];
dom = dom(1:nvars,:);

% Initialize an empty PDE in the proposed variables.
PDE = pde_struct();
PDE.vars = vars;
PDE.dom = dom;

% Establish the number of state, input, and output components.
nx = n_comps(1);
nw = n_comps(2);     nz = n_comps(4);
nu = n_comps(3);     ny = n_comps(5);

% Define the state components (size, variables, order of
% differentiability).
x_tab = zeros(nx,2+2*nvars);
for ii=1:nx
    % Set the size of the component.
    sz_ii = randi(max_size(1));
    
    % Set the number of variables on which the component depends.
    dep_ii = zeros(1,nvars);    % Binary indices, indicating whether the component depends on each variable
    diff_ii = zeros(1,nvars);   % Order of differentiability of the component wrt each variable
    nvars_ii = randi(nvars+1) - 1;
    if nvars_ii~=0
        var_indcs = randperm(nvars);
        var_indcs = var_indcs(1:nvars_ii);
        max_diff_ii = max_diff(var_indcs);
        dep_ii(var_indcs) = 1;
        for kk=1:nvars_ii
            diff_ii(var_indcs(kk)) = randi([0,max_diff_ii(kk)]);
        end
    end
    
    x_tab(ii,1) = ii;
    x_tab(ii,2) = sz_ii;
    x_tab(ii,3:2+nvars) = dep_ii;
    x_tab(ii,3+nvars:2+2*nvars) = diff_ii;
end
    
% Define the inputs and outputs (really just their sizes)
u_tab = zeros(nu,2+nvars);  y_tab = zeros(ny,2+nvars);
w_tab = zeros(nw,2+nvars);  z_tab = zeros(nz,2+nvars);

w_tab(:,1) = (1:nw);    w_tab(:,2) = randi(max_size(2),[nw,1]);
u_tab(:,1) = (1:nu);    u_tab(:,2) = randi(max_size(3),[nu,1]);
z_tab(:,1) = (1:nz);    z_tab(:,2) = randi(max_size(4),[nz,1]);
y_tab(:,1) = (1:ny);    y_tab(:,2) = randi(max_size(5),[ny,1]);

% Add the component information to the PDE.
PDE.x_tab = x_tab;  
PDE.w_tab = w_tab;          PDE.z_tab = z_tab;
PDE.u_tab = u_tab;          PDE.y_tab = y_tab;

% Set the state components, inputs, and outputs in the PDE structure.
objs = {'x','u','w','y','z'};
for objc = objs
    obj = objc{1};
    obj_tab = PDE.([obj,'_tab']);
    ncomps = size(obj_tab,1);
    for ii=1:ncomps
        PDE.(obj){ii}.size = obj_tab(ii,2);
        has_vars_Rcomp = logical(obj_tab(ii,3:2+nvars));
        PDE.(obj){ii}.vars = vars(has_vars_Rcomp,:);
        PDE.(obj){ii}.dom = dom(has_vars_Rcomp,:);
        if strcmp(obj,'x')
            diff_full = x_tab(ii,3+nvars:2+2*nvars);
            PDE.(obj){ii}.diff = diff_full(has_vars_Rcomp);
        end
    end
end



%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% % % Define terms in the PDE
objs = {'x','z','y'};
for kk = 1:length(objs)
obj = objs{kk};
obj_tab = PDE.([obj,'_tab']);
ncomps = size(obj_tab,1);
for ii=1:ncomps
    % Extract information regarding the component for which we are defining
    % an equation.
    nr = obj_tab(ii,2);     % Number of elements of the component
    has_vars_Rcomp = logical(obj_tab(ii,3:2+nvars));
    
    % Create a random number of new terms to describe the equation of
    % component ii.
    nterms_ii = randi(max_terms(kk));
    PDE.(obj){ii}.term = cell(1,nterms_ii);

    for jj=1:nterms_ii
        % Build a new random term.
        term_jj = build_random_term(PDE,has_vars_Rcomp,nr,max_deg);        
        % Add the new term to the PDE.
        PDE.(obj){ii}.term{jj} = term_jj;
    end 
end
end



%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% % % Construct BCs for the state components.

% First, make a list of all maximal orders of differentiability.
diff_vals = PDE.x_tab(:,3+nvars:2+2*nvars)';
diff_vals = diff_vals(:);
% For each degree, keep track of which state component is differentiated,
% and wrt which variable.
comp_nums = kron(PDE.x_tab(:,1),ones(nvars,1));
var_indcs = repmat((1:nvars)',[nx,1]);

% Make a full list of all allowed degrees of the derivative wrt each
% variable.
nBCs = sum(diff_vals);
diff_list = zeros(nBCs,1);
var_idx_list = zeros(nBCs,1);
comp_list = zeros(nBCs,1);

rnum = 0;
for kk=1:length(diff_vals)
    % If a state component is differentiable up to degree N wrt variable s,
    % then we need N-1 BCs for this state component in the variable s,
    % where the order of the derivative cannot exceed N-1.
    diff_list_kk = (0:diff_vals(kk)-1);
    diff_list(diff_list_kk+rnum+1) = diff_list_kk;
    var_idx_list(diff_list_kk+rnum+1) = var_indcs(kk);
    comp_list(diff_list_kk+rnum+1) = comp_nums(kk);
    rnum = rnum + length(diff_list_kk);
end

% Make a table describing all standard boundary conditions, indicating
% 1. Which state component is set equal to zero.
% 2. Which variable is evaluated at an upper or lower boundary
% 3. To what order is the state differentiated wrt this variable
% 4. Whether an upper boundary (1) or lower boundary (0) is considered.
BC_tab_0 = [comp_list, var_idx_list, diff_list, zeros(nBCs,1)];
BC_tab_1 = [comp_list, var_idx_list, diff_list, ones(nBCs,1)];
BC_tab = [BC_tab_0; BC_tab_1];

% Implement nBCs boundary conditions, choosing from those described by
% BC_tab.
for ii=1:nBCs
    
    % Establish for which state component we are imposing a boundary
    % condition, wrt which variable we are imposing a BC, and what is the
    % maximal order of the derivative of the state component wrt this
    % variable in this BC.
    Lcomp = comp_list(ii);
    nr = PDE.x_tab(Lcomp,2);
    var_idx_ii = var_idx_list(ii);
    diff_max_ii = diff_list(ii);
    has_vars_Rcomp = logical(PDE.x_tab(Lcomp,3:2+nvars));
    
    % Determine which of the BCs in BC_tab involve the appropriate state
    % component, spatial variable, and derivative order.
    BC_opts_indcs = find(BC_tab(:,1)==Lcomp & BC_tab(:,2)==var_idx_ii & BC_tab(:,3)<=diff_max_ii);
    % Choose one of the allowed BCs to impose.
    BC_choice_indx = BC_opts_indcs(randi(length(BC_opts_indcs)));
    
    % Establish the order of the derivative associated to the BC.
    Dval = zeros(1,nvars);
    Dval(var_idx_ii) = BC_tab(BC_choice_indx,3);
    Dval = Dval(has_vars_Rcomp);
    
    % Establish which spatial boundary is associated to the BC.
    use_upper_bndry = BC_tab(BC_choice_indx,4);
    locval = vars(:,1)';
    locval(var_idx_ii) = dom(var_idx_list(ii),use_upper_bndry+1);
    locval = locval(has_vars_Rcomp);
    
    % Set the first term to describe the desired boundary condition.
    PDE.BC{ii}.term{1}.x = Lcomp;
    PDE.BC{ii}.term{1}.D = Dval;
    PDE.BC{ii}.term{1}.loc = locval;
    PDE.BC{ii}.term{1}.I = cell(sum(has_vars_Rcomp),1);
    PDE.BC{ii}.term{1}.C = eye(nr);

    % Remove the implemented BC from the list of possible BCs;
    BC_tab(BC_choice_indx,:) = [];
    
    % Add some more random terms if desired.
    if add_BC_terms
        % Add some empty terms to the BC.
        nterms_ii = randi(add_BC_terms)+1;
        PDE.BC{ii}.term = [PDE.BC{ii}.term, cell(1,nterms_ii-1)];
        % Check on which variables the BC function 0=F(s) depends.
        has_vars_Lcomp = has_vars_Rcomp;
        has_vars_Lcomp(var_idx_ii) = false;
        for jj=2:nterms_ii
            % Build a new random term.
            term_jj = build_random_term(PDE,has_vars_Lcomp,nr,max_deg,true);        
            % Add the new term to the PDE.
            PDE.BC{ii}.term{jj} = term_jj;
        end
    end
    
end

% Finally, initialize the PDE structure.
PDE = initialize_PIETOOLS_PDE(PDE,true);


end



%% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % 
function term_jj = build_random_term(PDE,has_vars_Lcomp,nr,max_deg,is_BC)
% Build a random term for the PDE.
%
% INPUT
% - PDE:        A pde_struct object, defining the PDE to which the term
%               will be added (though we don't know which equation)
% - has_vars_Lcomp: 1 x nvars logical array indicating for each variable in
%                   PDE.vars whether the term is allowed to depend on this
%                   variable.
% - nr:         1x1 integer indicating the number of elements (rows) the
%               term should consist of.
% - max_deg:    1x1 integer specifying the maximally allowed degree of the
%               monomials appearing in the coefficient associated to this
%               term. The coefficient will be a randomly generate
%               polynomial of the appropriate size.
% - is_BC:      logical value indicating whether the term will be added to
%               a BC.
%
% OUTPUT
% - term_jj     A struct specifying a term according to the PIETOOLS PDE
%               terms-based format.
%
% NOTES
% In defining the term, the function will randomly assign a state component
% or input component to this term, randomly establish a spatial derivative
% for the term, randomly establish whether the term is evaluated at any
% boundary (and if so, at which one), randomly choose to perform no
% integration, partial integration, or full integration, and randomly
% assign coefficients (or a kernel in case of integration).
% In doing so, the function makes sure the term does not depend on any
% variables that are not allowed by "has_vars_Lcomp", and makes sure the
% order of differentiation, boundary evaluation, etc. match those allowed
% for the considered (state) component.

if nargin<5
    is_BC = false;
end

% Extract variables and their domain.
vars = PDE.vars;
dom = PDE.dom;
var1_Lcomp = vars(has_vars_Lcomp,1);
nvars = size(vars,1);

% Establish the number of state and input components.
nx = size(PDE.x_tab,1);
nu = size(PDE.u_tab,1);
nw = size(PDE.w_tab,1);

term_jj = struct();
% Enhance the probability of using a state component.
if randi([0,1])
    Robj_indx = randi(nx);
else
    Robj_indx = randi(nx+nw+nu);
end
if Robj_indx>nx
    % % The term involves an input.
    if Robj_indx>(nx+nw)
        % The input is an actuator input.
        Robj_indx = Robj_indx - nx - nw;
        term_jj.u = Robj_indx;
        nc = PDE.u_tab(Robj_indx,2);
    else
        % The input is an exogenous input.
        Robj_indx = Robj_indx - nx;
        term_jj.w = Robj_indx;
        nc = PDE.w_tab(Robj_indx,2);
    end
    if isempty(var1_Lcomp)
        % Create a random matrix.
        Cval = randi(max_deg+1,[nr,nc])-1;
    else
        % Create a random polynomial.
        Cval = PIETOOLS_random_poly_generator(nr,nc,var1_Lcomp,max_deg);
        if max(abs(Cval.coeff))==0
            % Avoid adding a zero term to the PDE.
            Cval = eye(nr,nc);
        end
    end
    term_jj.C = Cval;
else
    % % % The term involves a state component.
    % Establish the size of the state component.
    nc = PDE.x_tab(Robj_indx,2);
    % Determine on which variables the component depends.
    has_vars_Rcomp = logical(PDE.x_tab(Robj_indx,3:2+nvars));
    nvars_Rcomp = sum(has_vars_Rcomp);
    var1_Rcomp = vars(has_vars_Rcomp,1);
    dom_Rcomp = dom(has_vars_Rcomp,:);
    var_dom_Rcomp = [var1_Rcomp,dom_Rcomp];
    % Determine to what degree the component is differentiable wrt
    % each variable.
    diff_Rcomp = PDE.x_tab(Robj_indx,3+nvars:2+2*nvars);

    % % Randomly generate a derivative, integral, and position at
    % % which to evaluate the state.
    Dval = zeros(1,nvars_Rcomp);
    Idom = cell(nvars_Rcomp,1);
    locval = var1_Rcomp';
    % For the integrals and position, we have a fixed number of
    % options for each variable.
    int_type = randi(4,[nvars_Rcomp,1])-1;
    loc_type = randi(3,[nvars_Rcomp,1])-1;
    for kk=1:nvars_Rcomp
        % Set the order of the derivative wrt var kk.
        Dval(kk) = randi(diff_Rcomp(kk)+1) - 1;

        % Set the domain of integration.
        if int_type(kk)==0
            % No integration is performed wrt var kk.
            Idom{kk} = [];
        elseif int_type(kk)==1
            % Partial integral int_{a}^{s}.
            Idom{kk} = [dom_Rcomp(kk,1),var1_Rcomp(kk)];
        elseif int_type(kk)==2
            % Partial integral int_{s}^{b}.
            Idom{kk} = [var1_Rcomp(kk),dom_Rcomp(kk,2)];
        else
            % Full integral int_{a}^{b}.
            Idom{kk} = dom_Rcomp(kk,:);
        end
        % Set the position at which to evaluate the state.
        if int_type(kk)<3 && ...
                loc_type(kk)==0 && ~ismember(var1_Rcomp(kk).varname,var1_Lcomp.varname)
            % If the term is not evaluated at a boundary for variable kk,
            % but the LHS component does not depend on variable kk, we have
            % to perform full integration to get rid of the variable
            % dependence.
            int_type(kk) = 3;
            Idom{kk} = dom_Rcomp(kk,:);
        elseif int_type(kk)<3 && ~ismember(var1_Rcomp(kk).varname,var1_Lcomp.varname)
            % If the state is alreeady evaluated at a boundary, partial
            % integrals don't make sense.
            int_type(kk) = 0;
            Idom{kk} = [];
        end
        if int_type(kk)>0
            % If we're integrating wrt var kk, don't also evaluate
            % the state at a boundary wrt var kk.
            locval(kk) = var1_Rcomp(kk);
        elseif ((loc_type(kk)>0)+Dval(kk))>diff_Rcomp(kk)
            % For a state differentiable up to degree N wrt s, we
            % can only evaluate derivatives d^n/ds^n x(s) at s=a or
            % s=b if n is strictly smaller than N.
            locval(kk) = var1_Rcomp(kk);
            if int_type(kk)<3 && ~ismember(var1_Rcomp(kk).varname,var1_Lcomp.varname)
                % If the term is not evaluated at a boundary for variable kk,
                % but the LHS component does not depend on variable kk, we have
                % to perform full integration to get rid of the variable
                % dependence.
                int_type(kk) = 3;
                Idom{kk} = dom_Rcomp(kk,:);
            end
        else
            locval(kk) = var_dom_Rcomp(kk,loc_type(kk)+1);
        end
    end

    % % % Finally, generate random coefficients.
    % First determine which primary and dummy variables may appear
    % in the coefficient matrix.

    % Primary variables can only be used if either the LHS
    % component depends on them, or if they are integrated out from
    % the RHS component (using a type 3 integral).
    % Note that in the BCs, kernels for full integrals must also be
    % specified using dummy variables.
    use_primary_vars = has_vars_Lcomp;
    use_primary_vars(has_vars_Rcomp) = use_primary_vars(has_vars_Rcomp) | ((int_type==3)' & ~is_BC);
    % Dummy variables can be used only if we're using partial
    % integrals, i.e. if both the LHS and RHS depend on this
    % variable, and we are integrating the RHS.
    use_dummy_vars = has_vars_Lcomp & has_vars_Rcomp;
    use_dummy_vars(has_vars_Rcomp) = use_dummy_vars(has_vars_Rcomp) & (int_type>0)';

    % Set the full list of variables we can use.
    varlist = [vars(use_primary_vars,1); vars(use_dummy_vars,2)];
    % To reduce complexity, don't retain all variables per se.
    if ~isempty(varlist)
        retain_indcs = logical(randi([0,1],[length(varlist),1]));
        varlist = varlist(retain_indcs);
    end
    % Generate a random polynomial in these variables.
    if isempty(varlist)
        % Create a random matrix.
        Cval = randi(max_deg+1,[nr,nc])-1;
    else
        % Create a random polynomial.
        Cval = PIETOOLS_random_poly_generator(nr,nc,varlist,max_deg);
        if max(abs(Cval.coeff))==0
            % Avoid adding a zero term to the PDE.
            Cval = eye(nr,nc);
        end
    end

    % Set all properties of the term.
    term_jj.x = Robj_indx;
    term_jj.D = Dval;
    term_jj.loc = locval;
    term_jj.I = Idom;
    term_jj.C = Cval;
end
        
end
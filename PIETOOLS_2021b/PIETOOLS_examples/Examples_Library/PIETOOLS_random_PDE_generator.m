function PDE = PIETOOLS_random_PDE_generator(dim,max_diff,n_comps,max_size,max_terms,max_deg)
%PIETOOLS_RANDOM_PDE_GENERATOR Summary of this function goes here
%   Detailed explanation goes here


%rng('default');
% rng(12345);

arguments
    dim {mustBeMember(dim,[0;1;2])} = 2;
    max_diff {mustBeNonnegative} = 3;
    n_comps {mustBeNonnegative} = [randi(3),randi([0,2]),randi([0,2]),randi([0,2]),randi([0,2])];
    max_size {mustBeNonnegative} = [randi(3),randi(2,[1,4])];
    max_terms {mustBeNonnegative} = 4;
    max_deg {mustBeNonnegative} = 2;
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
    



nvars = dim;
pvar s1 s2 theta1 theta2
vars = [s1,theta1; s2,theta2];
dom = randi(3,[2,1])-2;
dom = [dom, dom+randi(3,[2,1])];

PDE = pde_struct();
PDE.vars = vars;
PDE.dom = dom;

% Establish the actual dimension of the domain
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
        dep_ii(var_indcs) = 1;
        diff_ii(var_indcs) = randi(max_diff+1,[1,nvars_ii]) - 1;
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

PDE.x_tab = x_tab;  
PDE.u_tab = u_tab;
PDE.w_tab = w_tab;
PDE.y_tab = y_tab;
PDE.z_tab = z_tab;


% Set the state components, inputs, and outputs in the PDE structure.
objs = {'x','u','w','y','z'};
for objc = objs
    obj = objc{1};
    obj_tab = PDE.([obj,'_tab']);
    ncomps = size(obj_tab,1);
    for ii=1:ncomps
        PDE.(obj){ii}.size = obj_tab(ii,2);
        has_vars_Lcomp = logical(obj_tab(ii,3:2+nvars));
        PDE.(obj){ii}.vars = vars(has_vars_Lcomp,:);
        PDE.(obj){ii}.dom = dom(has_vars_Lcomp,:);
        if strcmp(obj,'x')
            diff_full = x_tab(ii,3+nvars:2+2*nvars);
            PDE.(obj){ii}.diff = diff_full(has_vars_Lcomp);
        end
    end
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% % % Define terms in the PDE
for objc = {'x','y','z'}
obj = objc{1};
obj_tab = PDE.([obj,'_tab']);
ncomps = size(obj_tab,1);
for ii=1:ncomps
    % Extract information regarding the component for which we are defining
    % an equation.
    nr = obj_tab(ii,2);     % Number of elements of the component
    has_vars_Lcomp = logical(obj_tab(ii,3:2+nvars));
    
    % Create a random number of new terms to describe the equation of
    % component ii.
    nterms_ii = randi(max_terms);
    PDE.(obj){ii}.term = cell(1,nterms_ii);

    for jj=1:nterms_ii
        % Build a new random term.
        term_jj = build_random_term(PDE,has_vars_Lcomp,nr,max_deg);        
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
    has_vars_Lcomp = logical(PDE.x_tab(Lcomp,3:2+nvars));
    
    % Determine which of the BCs in BC_tab involve the appropriate state
    % component, spatial variable, and derivative order.
    BC_opts_indcs = find(BC_tab(:,1)==Lcomp & BC_tab(:,2)==var_idx_ii & BC_tab(:,3)<=diff_max_ii);
    % Choose one of the allowed BCs to impose.
    BC_choice_indx = BC_opts_indcs(randi(length(BC_opts_indcs)));
    
    % Establish the order of the derivative associated to the BC.
    Dval = zeros(1,nvars);
    Dval(var_idx_ii) = BC_tab(BC_choice_indx,3);
    Dval = Dval(has_vars_Lcomp);
    
    % Establish which spatial boundary is associated to the BC.
    use_upper_bndry = BC_tab(BC_choice_indx,4);
    locval = vars(:,1)';
    locval(var_idx_ii) = dom(var_idx_list(ii),use_upper_bndry+1);
    locval = locval(has_vars_Lcomp);
    
    % Set the first term to describe the desired boundary condition.
    PDE.BC{ii}.term{1}.x = Lcomp;
    PDE.BC{ii}.term{1}.D = Dval;
    PDE.BC{ii}.term{1}.loc = locval;
    PDE.BC{ii}.term{1}.I = cell(sum(has_vars_Lcomp),1);
    PDE.BC{ii}.term{1}.C = eye(nr);

    % Remove the implemented BC from the list of possible BCs;
    BC_tab(BC_choice_indx,:) = [];
    
    % Add some more random terms if desired.
    toggle = false;
    if toggle
        nterms_ii = randi(max_terms);
        PDE.BC{ii}.term = [PDE.BC{ii}.term, cell(1,nterms_ii-1)];
        for jj=2:nterms_ii
            % Build a new random term.
            term_jj = build_random_term(PDE,has_vars_Lcomp,nr,max_deg);        
            % Add the new term to the PDE.
            PDE.BC{ii}.term{jj} = term_jj;
        end
    end
    
end

PDE = initialize_PIETOOLS_PDE(PDE,true);


end



%%
function term_jj = build_random_term(PDE,has_vars_Lcomp,nr,max_deg)

vars = PDE.vars;
dom = PDE.dom;
var1_Lcomp = vars(has_vars_Lcomp,1);
nvars = size(vars,1);

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
    use_primary_vars = has_vars_Lcomp;
    use_primary_vars(has_vars_Rcomp) = use_primary_vars(has_vars_Rcomp) | (int_type==3)';
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
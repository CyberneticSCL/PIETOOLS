function [PDE,Gvar_order] = initialize(PDE,suppress_summary)
% PDE = initialize(PDE) takes a PDE data structure and
% checks that all the necessary fields are appropriately specified, and
% assigns a default value to all the optional fields that have not been
% specified.
% 
% INPUT
%   PDE:              "pde_struct" class object.
%   suppress_summary: Logical index. If 'false' or not specified, a summary
%                     of the observed state, input, output and BC
%                     information will be printed in the command window at
%                     the end of the initialization process.
%
% OUTPUT
%   PDE: "pde_struct" class object describing the same system as the input.
%   Gvar_order:     New order of the variables in PDE.vars, assuming
%                   a field PDE.vars was already specified in the input
%                   PDE: PDE_out.vars = PDE_in.vars(Gvar_order,:);
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% The PDE structure has the following fields:
%
% - PDE.dim:    Integer "p" specifying the dimension of the
%   (dependent) spatial domain of the system. Currently, only p<=2 is
%               supported.
% - PDE.vars:   A px2 pvar (polynomial class) object specifying
%   (dependent) the p spatial variables that appear in the system in the
%               first column, and the associated dummy variables in the
%               second column.
% - PDE.dom:    A px2 "double" array specifying the interval on
%   (dependent) which each of p spatial variables exists.
%
% - PDE.tau:    A qx2 "polynomial" array specifying temporal variables in
%   (mandatory) the first column, and the maximal value these delays can
%               assume in the second column.
%
% - PDE.x:      A cell with each element i defining a state 
%   (mandatory) component x_i and associated differential equation
%               \dot{x}_i = ....
% - PDE.y:      A cell with each element i defining an observed
%   (optional)  output y_i and associated equation y_i = .... If not
%               specified, it is assumed that the system has no observed
%               outputs.
% - PDE.z:      A cell with each element i defining a regulated
%   (optional)  output z_i and associated equation z_i = .... If not
%               specified, it is assumed that the system has no regulated
%               outputs.
% - PDE.u:      A cell with each element i defining an actuator
%   (optional)  input u_i that appears in the system. If not specified, it
%               is assumed that the system has no actuator inputs.
% - PDE.w:      A cell with each element i defining an exogenous
%   (optional)  input w_i that appears in the system. If not specified, it
%               is assumed that the system has no exogenous inputs.
% - PDE.BC:     A cell with each element i defining a boundary
%   (mandatory) condition 0 = ... of the system.
%
% For each state component, input, and output, the size of this component,
% as well as the number of variables it depends on and their associated 
% domain are specified in the same manner, through fields "size", "vars" 
% and "dom". In addition, the order of differentiability of state
% components is indicated through "diff", and the order of the temporal
% derivative in the LHS of the PDE is indicated through "tdiff". That is,
% PDE.x{i} has fields:
%
% - PDE.x{i}.size:  An integer n_i indicating the number of
%   (optional)      variables in the ith state component x_i. If
%                   not specified, a size is determined based on the size
%                   of coefficient matrices in the different equations
%                   defining or relying on this state. If specified, an
%                   error will be thrown if the size of these coefficient
%                   matrices does not match the specified size.
% - PDE.x{i}.vars:  A p_ix2 "pvar" (polynomial class) object
%   (mandatory)     specifying the p_i spatial variables that the
%                   component depends on in the first column, and the
%                   associated dummy variables in the second column. If not
%                   specified, the component is assumed to be
%                   finite-dimensional, unless global variables PDE.vars
%                   are specified. In that case, any state component is
%                   assumed to vary in these global variables, but inputs
%                   and outputs are still assumed to be finite-dimensional.
% - PDE.x{i}.dom:   A p_ix2 "double" array specifying the
%   (mandatory)     spatial interval on which each of the p_i variables
%                   that the component depends on exists.
% - PDE.x{i}.diff:  A 1xp_i array of integers specifying for
%   (optional)      each of the p_i variables on which the state component
%                   depends to what order the component is differentiable
%                   with respect to this variable. The field should not be
%                   specified for inputs, outputs, or BCs. If no "diff" is
%                   specified, the function will look through all the terms
%                   in each equation that depends on the component x_i, and
%                   use the maximal order of any derivative taken of this
%                   component as the maximal order of differentiability.
%
% - PDE.x{i}.tdiff: A scalar integer r_i indicating the order
%   (optional)      of the temporal derivative in the PDE d_t^{r_i} x_i=...
%                   associated to the ith state component. Should not be
%                   specified for the inputs, outputs, or BCs. Defaults to
%                   1 if not specified.
%
% - PDE.x{i}.term:  A cell with each element defining a term in
%   (mandatory)     the (differential) equation associated to the
%                   component. Should not be specified for the inputs
%                   PDE.w{i} or PDE.u{i}.
%
% To specify the PDEs, output equations, and boundary conditions, a cell
% "term" is used, each element describing a term in the equation. Each of
% these terms is structured in the same way as
%
% term_j = int_{I{1}(1)}^{I{1}(2) ... int_{I{p_r}(1)}^{I{p_r}(2)}
%             [ C(s_1,...,s_{p_i},theta_1,...,theta_{p_r}) * 
%                   d_{s_{1}}^D(1) ... d_{s_{p_r}}^{D(p_r)}
%                       x_r(loc(1),...,loc(p_r)) ]dtheta_{p_r} ... dtheta_1
%
% where x_r could be replaced by w_r or u_r. Such a term can then be
% specified through the fields
%
% -    term{j}.x;   An integer r specifying which state component x_r
%   OR term{j}.w;   or input w_r or u_r the term is defined by. Only one of
%   OR term{j}.u:   these fields can be specified, and one of these fields
%   (mandatory)     must be specified. If none of these fields is
%                   specified, and the systems pertains only a single state
%                   component, the value defaults to the index of this
%                   state component (1).
% - term{j}.D:      A 1xp_r integer array specifying the order of the
%   (optional)      derivative of the state component x_r in each of the
%                   p_r variables that this state component depends on.
%                   Should not be specified if the term involves an input
%                   w_r or u_r (we cannot take derivatives of inputs). If
%                   not specified, it defaults to D=zeros(1,p_r), taking no
%                   derivative.
% - term{j}.loc:    A 1xp_r "pvar" (polynomial class) or "double" array
%   (optional)      specifying the spatial position at which to evaluate
%                   the state. For any variable s_k in [a,b] on which the
%                   state depends, only loc(k)=s_k, loc(k)=a or loc(k)=b is
%                   allowed. Should not be specified if the term involves
%                   an input w_r or u_r. For BCs, this field is mandatory,
%                   as otherwise it is unclear at which boundary the state
%                   is considered. In the other equations, if no "loc" is
%                   specified, it defaults to loc=[s_{1},...,s_{p_r}],
%                   where (s_{1},...,s_{p_r}) are the spatial variables on
%                   which the state depends.
% - term{j}.C:      A n_ixn_r "polynomial" or "double" array, specifying
%   (optional)      the coefficient matrix with which the component x_r,
%                   w_r or u_w is multiplied, as it contributes to the
%                   considered equation. If it is polynomial, it can only
%                   vary in the state variables on which the LHS component
%                   x_i, y_i, or z_i depends. If integration is performed,
%                   it can also depend on the dummy variables associated
%                   to the spatial variables on which the RHS component
%                   x_r, w_r or u_r depends. If not specified, it defaults
%                   to an identity matrix.
% - term{j}.I:      A cell of p_r elements, the kth of which must be a 1x2
%   (optional)      "pvar" or "double" array specifying on which domain to
%                   integrate the RHS component with respect to the kth
%                   variable on which it depends. If the kth element is
%                   empty, no integration is performed along this
%                   direction. Otherwise, for the kth variables s_k in
%                   [a,b], only I{k}=[a,s_k], I{k}=[s_k,b] or I{k}=[a,b]
%                   can be specified. Any other domain of integration will
%                   not be allowed. If no I is specified, it defaults to a
%                   p_rx1 cell of which each element is empty, performing
%                   no integration.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% The initialization function verifies that all mandatory fields are
% proprly specified, and assigns a value to the optional fields which the
% user has not specified, as well as to the dependent fields. In addition,
% it builds the tables x_tab, y_tab, z_tab, u_tab, w_tab and BC_tab. The
% number of rows in each of these tables matches the number of associated
% components (e.g. x_tab has one row for each element of PDE.x), and for a 
% system with p global variables, the number of columns is given by 2+p.
% The first column specifies the index associated to each component.
% The second column specifies the size of this (vector-valued) component.
% The remaining columns j+2 for j=1:p are logical indices indicating
% whether the considered component varies in the jth global variable.
% The x_tab has p additional columns, indicating for each variable j=1:p to
% what order the state component is differentiable in this spatial
% variable.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - initialize_PIETOOLS_PDE
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 07/07/2022


% % % --------------------------------------------------------------- % % %
% % % Check the PDE structure, and initialize any fields that the user has
% % % not set.

% If no input is provided, initialize an empty system
if nargin==0
    PDE = pde_struct;
    return
elseif nargin==1
    suppress_summary = false;
end
% Convert the input to a "pde_struct" object.
if ~isa(PDE,'struct') && ~isa(PDE,'pde_struct')
    error('Input must be a ''struct'' or ''pde_struct'' class object.')
elseif isa(PDE,'struct')
    PDE = pde_struct(PDE);
end
ncomps_x = length(PDE.x);

% % Re-order the cell elements.
PDE.x = [PDE.x(:)];
PDE.y = [PDE.y(:)];
PDE.z = [PDE.z(:)];
PDE.u = [PDE.u(:)];
PDE.w = [PDE.w(:)];
PDE.BC = [PDE.BC(:)];

% Sort the rows of the tables.
[~,new_order_x] = sort(PDE.x_tab(:,1));
PDE.x_tab = PDE.x_tab(new_order_x,:);
[~,new_order_y] = sort(PDE.y_tab(:,1));
PDE.y_tab = PDE.y_tab(new_order_y,:);
[~,new_order_z] = sort(PDE.z_tab(:,1));
PDE.z_tab = PDE.z_tab(new_order_z,:);
[~,new_order_u] = sort(PDE.u_tab(:,1));
PDE.u_tab = PDE.u_tab(new_order_u,:);
[~,new_order_w] = sort(PDE.w_tab(:,1));
PDE.w_tab = PDE.w_tab(new_order_w,:);
[~,new_order_BC] = sort(PDE.BC_tab(:,1));
PDE.BC_tab = PDE.BC_tab(new_order_BC,:);

% % % --------------------------------------------------------------- % % %
% % % Build a full list of the spatial variables that appear, and
% % % their domain.

% % First check if a set of global variables and associated domain has been
% % specified.
if ~isempty(PDE.dom) && ~isempty(PDE.vars)
    % If global variables have been specified, make sure they have been
    % appropriately specified.
    global_dom = PDE.dom;
    global_vars = PDE.vars;
    if size(global_dom,2)~=2
        error(['The global domain has not been appropriately specified: '...
                '"dom" must be a px2 array specifying the interval on which each of the p spatial variables exists.'])
    elseif any(global_dom(:,1) >= global_dom(:,2))
        error(['The global domain has not been appropriately specified: '...
                'lower boundaries of spatial domain cannot exceed the upper boundaries.'])
    end
    if size(global_vars,2)>2
        error(['The global variables have not been appropriately specified: '...
                '"vars" should be a px2 array, specifying the p primary and dummy spatial variables.'])
    elseif size(global_vars,1)~=size(global_dom,1)
        if size(global_vars,1)==1 && size(global_vars,2)==size(global_dom,1)
            % Only primary spatial variables have been specified
            global_vars = global_vars';
        else
            error(['The global variables have not been appropriately specified :'...
                    'the number of variables does not match the dimension of the specified domain.'])
        end
    end
    % If variable names are specified, convert to pvar objects.
    if ~ispvar(global_vars) && ~iscellstr(global_vars)
        error(['The global variables have not been appropriately specified: '...
                'spatial variables must be specified as pvar (polynomial class) objects.'])
    elseif iscellstr(global_vars)
        global_vars = polynomial(global_vars);
    end
    if length(global_vars.varname)~=prod(size(global_vars))
        error(['Global variables have not been appropriately specified: '...
                'spatial variables must be distinct.'])
    end
    % If no dummy variables have been specified, do not set any dummy
    % variables just yet.
    if size(global_vars,2)==1
        global_vars = [global_vars,zeros(size(global_vars,1),1)];
    end
elseif isempty(PDE.dom) && ~isempty(PDE.vars)
    % We do not allow variables to be specified without domain.
    error(['A domain must be specified for the spatial variables.'])
elseif isempty(PDE.vars) && ~isempty(PDE.dom)
    % We do not allow a domain to be specified without variables.    
    error(['No spatial variables have been specified.'])
else
    % If no global variables are available, initialize an empty set.
    global_dom = zeros(0,2);
    global_vars = polynomial(zeros(0,2));
end
PDE.dom = global_dom;
PDE.vars = global_vars;
    
% % Next, determine for each state, input, and output what spatial variables
% % it depends on, along with their domain.
[PDE,dom_list_x,var_list_x,comp_list_x] = extract_vars(PDE,'x');
[PDE,dom_list_w,var_list_w,comp_list_w] = extract_vars(PDE,'w');
[PDE,dom_list_u,var_list_u,comp_list_u] = extract_vars(PDE,'u');
[PDE,dom_list_z,var_list_z,comp_list_z] = extract_vars(PDE,'z');
[PDE,dom_list_y,var_list_y,comp_list_y] = extract_vars(PDE,'y');

% Combine the lists to get a full list of variable names and their
% associated domains.
dom_list = [dom_list_x; dom_list_w; dom_list_u; dom_list_z; dom_list_y];
var_list = [var_list_x; var_list_w; var_list_u; var_list_z; var_list_y];
comp_list = [comp_list_x; comp_list_w; comp_list_u; comp_list_z; comp_list_y];

% Keep track of which rows in the full list correspond to states, inputs,
% and outputs.
nvars_arr = [size(comp_list_x,1); size(comp_list_w,1); size(comp_list_u,1); size(comp_list_z,1); size(comp_list_y,1)];
nnvars_arr = cumsum([0; nvars_arr]);    


% % Establish a unique set of primary variables from the full set.
[var1_list_unique,id1,varnum_list] = unique(var_list(:,1),'stable'); % list_unique = list(id1);  list = list_unique(varnum_list)
nvars = length(var1_list_unique);
if nvars>2 && ~suppress_summary
    warning(['Currently, PIETOOLS supports only problems with at most two distinct spatial variables.'...
            ' Analysis of the returned PDE structure will not be possible.'...
            ' Try running "combine_vars" to reduce the dimensionality of your problem.']);
end

% % Loop over the unique variables, and make sure each is coupled to just
% % one domain and (at most) one dummy variable.
for varnum=1:nvars
    % Find rows in var_list associated to unique spatial variable varnum.
    var_locs = find(varnum_list==varnum);   
    % Check that each spatial variable is coupled with just one interval on
    % which it exists.
    if any(any(dom_list(var_locs(2:end,:),:) - dom_list(var_locs(1:end-1,:),:)))
        error(['The same spatial variable seems to exist on different domains for different states/inputs/outputs;'...
                ' please make sure that the domain of any one variable across different components x{i}, w{i}, u{i}, z{i} and y{i} is identical.'])
    end
    % For the given variable varnum, check the third column of var_list to
    % establish whether any dummy variables have been coupled to this 
    % variable
    has_var2_ii = cell2mat(var_list(var_locs,3));
    if ~any(has_var2_ii)
        % If there is no dummy variable associated to the spatial
        % variable varnum, define a new dummy variable.
        var2_name_ii = [var1_list_unique{varnum},'_dum'];
        var_list(var_locs,2) = {var2_name_ii};
    else
        % Otherwise, make sure only one dummy variable is coupled to the
        % spatial variable.
        var2_list_ii = var_list(var_locs,2);
        var2_list_ii = unique(var2_list_ii(has_var2_ii));
        if length(var2_list_ii)>1
            error(['The same spatial variable seems to be coupled with different dummy variables across different states/inputs/outputs;'...
                    ' please make sure that the dummy variable associated to any one variable across different components x{i}, w{i}, u{i}, z{i} and y{i} is the same.'])
        else
            % Make sure that the dummy variable is linked to the spatial
            % variable varnum throughout the var_list
            var_list(var_locs,2) = var2_list_ii;
        end
    end
end

% Make sure that each dummy variable is coupled to at most one spatial
% variable.
var2_names = var_list(id1,2);
if length(unique(var2_names))~=length(var2_names)
    error(['Different spatial variables seem to be coupled to the same dummy variable;'...
            ' please use a distinct dummy variable for each distinct primary spatial variable.'])
end

% % Set a (unique) set of global variables, and their associated domain,
% % and add these to the PDE structure.
global_vars = polynomial(var_list(id1,1:2));
global_dom = dom_list(id1,:);

Gvar_order = (1:size(global_vars,1));
if size(PDE.vars,1)==size(global_vars,1)
    for kk=1:size(global_vars,1)
        % var_order is such that vars_new = vars_old(var_order)
        Gvar_order(kk) = find(isequal(global_vars(kk,1),PDE.vars),1,'first');
    end
end

PDE.dom = global_dom;
PDE.vars = global_vars;


% % Check whether the temporal delay variables have been properly
% % specified.
if size(PDE.tau,2)~=2
    error('Delay variables "PDE.tau" should be specified as a qx2 array of pvars in the first array, and the associated maximal delay in the second column.')
elseif size(PDE.tau,1)>1
    if ~ispvar(PDE.tau(:,1)) || length(PDE.tau(:,1))~=length(PDE.tau(:,1).varname)
        error('The first column of "PDE.tau" should be given by single pvar objects representing temporal delay (dummy) variables.')
    elseif length(PDE.tau(:,1).varname)~=length(unique(PDE.tau(:,1).varname))
        error('The variables in the first column of "PDE.tau" should be unique.');
    end
    if ~isa(PDE.tau(:,2),'double') && ~isdouble(PDE.tau(:,2))
        error('The second column of "PDE.tau" should be real values providing the maximal value which the variables in the first column exist.')
    end
    PDE.tau(:,2) = polynomial(abs(double(PDE.tau(:,2))));
end


% % % --------------------------------------------------------------- % % %
% % % Now that we have established a unique set of variables, we build
% % % tables for each state, input, and output, keeping track of which
% % % of the global variables they depend on, as well as their size.

% % First, establish for each state, input, and output which rows in the
% % var_list and dom_list correspond to this component.
n_comps = [numel(PDE.x); numel(PDE.w); numel(PDE.u); numel(PDE.z); numel(PDE.y)];
comp_indcs = cell(5,1);
for j=1:numel(comp_indcs)
    % For j=1:5, indices correspond to the rows in var_list associated with
    % x (j=1), w(j=2), u (j=3), z (j=4), and y (j=5).
    comp_indcs{j} = nnvars_arr(j)+1:nnvars_arr(j+1);
end

% % For the states, inputs, and outputs, we build a table with the number of
% % rows equal to the number of components of this state, input and output.
tab_cell = cell(5,1);   % {x_tab (j=1); w_tab (j=2); u_tab (j=3); z_tab (j=4); y_tab (j=5)}.
objs = {'x','w','u','z','y'};
for j=1:numel(tab_cell)
    % The first two columns in the table list the index associated to this
    % component, and respectively the size of the component.
    tab_j = zeros(n_comps(j),2 + nvars);
    tab_j(:,1) = 1:n_comps(j);
    if size(PDE.([objs{j},'_tab']),1)==size(tab_j,1)
        tab_j(:,2) = PDE.([objs{j},'_tab'])(:,2);
    end
    %tab_j(:,2) = 1;     % Assume each component is scalar-valued for now
    
    % The remaining nvars columns are binary indices, indicating for each
    % varnum=1:nvars whether the component varies in the variable varnum.
    if isempty(comp_indcs{j})
        % There are no rows in var_list associated to object j.
        tab_cell{j} = tab_j;
        continue
    end
    % Extract variable names and associated indices assoicated to this
    % object (the state, input, or output).
    comp_list_j = comp_list(comp_indcs{j});
    varnum_list_j = varnum_list(comp_indcs{j});    
    for varnum=1:nvars
        var_locs = varnum_list_j==varnum;   % Rows in var_list associated to spatial variable varnum
        row_nums = comp_list_j(var_locs);   % Components that vary in variable varnum
        tab_j(row_nums,2+varnum) = 1;       % Use binary index to indicate dependence on variable varnum of component j
    end
    % Store the new table.
    tab_cell{j} = tab_j;
end

% Separate and store the tables in the PDE structure
PDE.x_tab = tab_cell{1};
PDE.w_tab = tab_cell{2};
PDE.u_tab = tab_cell{3};
PDE.z_tab = tab_cell{4};
PDE.y_tab = tab_cell{5};

% Also build a table for the boundary conditions, where once again, each
% row corresponds to a particular condition, with the first column
% providing an index associated to this condition, and the second column
% a number of equations (size) associated to this condition. The remaining
% columns are filled with binary indices, indicating for each of the global
% variables whether the function f(s) in the BC 0=f(s) depends on this
% variable.
BC_tab = zeros(numel(PDE.BC),2+nvars);
BC_tab(:,1) = 1:numel(PDE.BC);
%BC_tab(:,2) = 1;    % Assume each BC concerns just a single equation for now
PDE.BC_tab = BC_tab;


% Using the tables, set the variables on which each state, input, and
% output depends.
PDE = set_vars(PDE,'x');
PDE = set_vars(PDE,'y');
PDE = set_vars(PDE,'z');
PDE = set_vars(PDE,'u');
PDE = set_vars(PDE,'w');


% % For the state components, we also store the order of differentiability
% % in each spatial variable, in diff_tab.
diff_tab = zeros(numel(PDE.x),nvars);
BC_diff_tab = zeros(numel(PDE.BC),nvars);

% For the different equations for x, y, z, and the BCs, we use 
% "get_diff_tab" to loops over the terms in each of these equations, 
% checking the order of the derivatives of the state components to
% establish an order of differentiability of each component.
% The functions also determine a size for each state, input, output, and
% BC, and stores it in the second column of the appropriate table 
% PDE.(*)_tab. 
[PDE,diff_tab] = get_diff_tab(PDE,'x',diff_tab,Gvar_order);
[PDE,diff_tab] = get_diff_tab(PDE,'z',diff_tab,Gvar_order);
[PDE,diff_tab] = get_diff_tab(PDE,'y',diff_tab,Gvar_order);
[PDE,diff_tab] = get_diff_tab(PDE,'BC',diff_tab,Gvar_order);


% % Finally, set sizes of the objects.
% For any state component of which the size is not clear, assume the
% component to be scalar.
PDE.x_tab(PDE.x_tab(:,2)==0,2) = 1;
% For the inputs, outputs, and BCs of which the size is unclear, if all
% state components are the same size, assume the other objects to be of
% this same size as well
if all(PDE.x_tab(2:end,2)==PDE.x_tab(1:end-1,2))
    shared_size = PDE.x_tab(1,2);
    PDE.y_tab(PDE.y_tab(:,2)==0,2) = shared_size;
    PDE.z_tab(PDE.z_tab(:,2)==0,2) = shared_size;
    PDE.u_tab(PDE.u_tab(:,2)==0,2) = shared_size;
    PDE.w_tab(PDE.w_tab(:,2)==0,2) = shared_size;
    %PDE.BC_tab(PDE.BC_tab(:,2)==0,2) = shared_size;
else
    % Otherwise, assume they are all scalar.
    PDE.y_tab(PDE.y_tab(:,2)==0,2) = 1;
    PDE.z_tab(PDE.z_tab(:,2)==0,2) = 1;
    PDE.u_tab(PDE.u_tab(:,2)==0,2) = 1;
    PDE.w_tab(PDE.w_tab(:,2)==0,2) = 1;
    %PDE.BC_tab(PDE.BC_tab(:,2)==0,2) = 1;
end


% % % --------------------------------------------------------------- % % %
% % % Now that we have all the necessary information on our states, inputs
% % % and outputs, we loop over the different terms in the different
% % % equations (except for the BCs) and make sure all the mandatory fields
% % % are appropriately specified. Any optional fields that are no
% % % specified are assigned the associated default value.

[PDE,diff_tab] = check_terms(PDE,'x',diff_tab);
[PDE,diff_tab] = check_terms(PDE,'z',diff_tab);
[PDE,diff_tab] = check_terms(PDE,'y',diff_tab);

% Append the order of differentiability to the x_tab.
PDE.x_tab = [PDE.x_tab, diff_tab];


% % % --------------------------------------------------------------- % % %
% % % Finally: The boundary conditions!
% % % For the boundary conditions, we perform very similar checks as for
% % % the previous equations. However, for the BCs, a field ".loc" MUST be
% % % specified in each term involving a state component, indicating at
% % % what (boundary) position this state component should be evaluated.
% % % Also, we determine for each condition whether it is finite or
% % % infinite-dimensional, and in what variables it might vary, and keep 
% % % track of this in the table PDE.BC_tab.

BC_state_tab = zeros(numel(PDE.BC),numel(PDE.x));

for ii=1:numel(PDE.BC)
    
    % Set the number of BCs that this field describes.
    if isfield(PDE.BC{ii},'size') && PDE.BC{ii}.size ~= PDE.BC_tab(ii,2)
        error(['The specified size "BC{',num2str(ii),'}.size" does not match the observed number of BC equations of "BC{',num2str(ii),'};',...
                ' please adjust the specified size or the dimensions of the matrices "C" in the appropriate terms.'])
    end
    PDE.BC{ii}.size = PDE.BC_tab(ii,2);
    
    [PDE,BC_diff_tab_ii,BC_state_indcs] = check_BC_terms(PDE,ii);
    BC_state_tab(ii,BC_state_indcs) = true;
    BC_diff_tab(ii,:) = BC_diff_tab_ii;
    
    % Make sure that the boundary condition is indeed a boundary condition.
    if all(PDE.BC_tab(ii,3:2+nvars))
        error(['The boundary condition "BC{',num2str(ii),'}", of the form 0=F(s), is defined by a function F(s) that varies on the interior of the spatial domain.',...
                ' Please use boundary positions in "BC{',num2str(ii),'}.term{j}.loc" and integrals "BC{',num2str(ii),'}.term{j}.I" to make sure the boundary condition function F(s) does not vary along all spatial directions.'])
    else
        % Keep track of which variables s the function f(s) in the BC
        % 0=f(s) depends on, and their associated domain.
        PDE.BC{ii}.vars = global_vars(logical(PDE.BC_tab(ii,3:2+nvars)),:);
        PDE.BC{ii}.dom = global_dom(logical(PDE.BC_tab(ii,3:2+nvars)),:);
    end
    
    if isfield(PDE.BC{ii},'eq_num')
        % The field eq_num was only for internal use...
        PDE.BC{ii} = rmfield(PDE.BC{ii},'eq_num');
    end
    % Also order the fields specifying each BC.
    PDE.BC{ii} = orderfields(PDE.BC{ii},{'size','vars','dom','term'});    
end

% Append the orders of differentiability to the BC_tab.
PDE.BC_tab = [PDE.BC_tab, BC_diff_tab];

% For each order of differentiability of the state components, a boundary
% condition must be specified for the system to be well-posed.
required_num_BCs_arr = sum(diff_tab,2);         % Required number of BCs for each state component
observed_num_BCs_arr = sum(BC_state_tab,1)';    % Observed number of BCs which contain each state component
if any(observed_num_BCs_arr < required_num_BCs_arr) || sum(PDE.BC_tab(:,2)) < PDE.x_tab(:,2)'*required_num_BCs_arr
    warning(['The number of boundary conditions specified for certain state components seems to be insufficient to ensure a well-posed system;',...
                ' conversion to a PIE will likely not be possible.',...
                ' Please ensure that for each state component, the number of boundary conditions along any spatial dimension matches the order of differentiability of this component along this dimension.'])
end
% Should we also check that sufficient BCs are specified along each spatial
% dimension?


% % % --------------------------------------------------------------- % % %
% % % With that, we have finished the initialization. We return the
% % % structure as a "PDE_struct" class object, and display some results.

% % First, clean up the different objects.
objs = {'x','y','z','u','w','BC'};
for kk = 1:numel(objs)
    obj = objs{kk};
    for ii = 1:numel(PDE.(obj))
        if isfield(PDE.(obj){ii},'eq_num')
            % The field eq_num was only for internal use...
            PDE.(obj){ii} = rmfield(PDE.(obj){ii},'eq_num');
        end
        if isfield(PDE.(obj){ii},'var_order')
            % The field var_order was only for internal use...
            PDE.(obj){ii} = rmfield(PDE.(obj){ii},'var_order');
        end
        % Also order the fields of the LHS component.
        if strcmp(obj,'x')
            PDE.(obj){ii} = orderfields(PDE.(obj){ii},{'size','vars','dom','diff','tdiff','term'});
        elseif strcmp(obj,'y') || strcmp(obj,'z') || strcmp(obj,'BC')
            PDE.(obj){ii} = orderfields(PDE.(obj){ii},{'size','vars','dom','term'});
        else
            if ~isfield(PDE.(obj){ii},'size')
                % Set number of input variables.
                PDE.(obj){ii}.size = PDE.([obj,'_tab'])(ii,2);
            end
            PDE.(obj){ii} = orderfields(PDE.(obj){ii},{'size','vars','dom'});
        end
    end
end

% % Finally, if desired, display an overview of how many state components,
% % inputs, outputs, and BCs were encountered, what their sizes are, and on
% % what variables they depend.
if ~suppress_summary
    print_initialization_summary(PDE,'x',ncomps_x);
    print_initialization_summary(PDE,'u');
    print_initialization_summary(PDE,'w');
    print_initialization_summary(PDE,'y');
    print_initialization_summary(PDE,'z');
    print_initialization_summary(PDE,'BC');
end

% Indicate that the PDE has been initialized.
PDE.is_initialized = true;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [PDE,dom_list,var_list,comp_list] = extract_vars(PDE,obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PDE,dom_list,var_list,comp_list] = extract_vars(PDE,obj)
% checks component "comp" of a PDE structure "PDE", and verifies that the
% variables that this component depends on, as well as their associated
% domain, have been properly specified. The observed variables and their
% associated domain are returned along with the updated PDE structure.
%
% INPUTS:
% - PDE:    A "struct" or "pde_struct" class object defining a PDE.
% - obj:    A "char" object 'x', 'y', 'z', 'u' or 'w', specifying for which
%           object of the PDE we are checking the variables.
%
% OUTPUTS:
% - PDE:        The same structure as the input "PDE", now with fields
%               "vars" and "dom" specified for all elements PDE.(comp){ii}.
% - dom_list:   An nx2 "double" array collecting the spatial domains 
%               "PDE.(obj){ii}.dom" for all elements ii.
% - var_list:   An nx3 cell listing all the observed variable names
%               "PDE.(comp){ii}.vars(:,1)" for all elements ii in the first
%               column, the observed associated dummy variable names
%               "PDE.(obj){ii}.vars(:,2)" in the second column, and a
%               logical index indicating whether a dummy variable has in
%               fact been specified or not.
% - comp_list:  An nx1 integer array specifying for each of the observed
%               variable names from which element ii of "PDE.(comp)" it was
%               extracted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize an empty list of variable names and associated domains.
var_list = cell(0,3);
dom_list = zeros(0,2);
comp_list = zeros(0,1);     % Keep track of which component is associated to each variable

% % Loop over the equations for each component of "comp".
for ii=1:numel(PDE.(obj))
    PDE_comp = PDE.(obj){ii};
    if ~isfield(PDE_comp,'dom') && ~isfield(PDE_comp,'vars')
        % Set default values.
        if strcmp(obj,'x') && ~isempty(PDE.dom)
            % For state components default variables/domain to the global
            % values.
            PDE_comp.vars = PDE.vars;
            PDE_comp.dom = PDE.dom;
        else
            % For inputs and outputs, assume finite-dimensionality if no
            % variables have been specified.
            PDE.(obj){ii}.vars = polynomial(zeros(0,2));
            PDE.(obj){ii}.dom = zeros(0,2);
            continue
        end
    elseif isfield(PDE_comp,'dom') && ~isfield(PDE_comp,'vars')
        % A lack of variables is only allowed if the domain indicates the
        % component to be finite-dimensional.
        if isempty(PDE_comp.dom)
            PDE.(obj){ii}.vars = polynomial(zeros(0,2));
            continue
        else
            error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                    ' no spatial variables have been specified.'])
        end
    elseif isfield(PDE_comp,'vars') && ~isfield(PDE_comp,'dom')
        % A lack of domain is only allowed if the vars indicate the
        % component to be finite-dimensional.
        if isempty(PDE_comp.vars)
            PDE.(obj){ii}.dom = zeros(0,2);
            continue
        else
            error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                    ' no spatial domain has been specified.'])
        end
    else
        % Check that the domain has been properly specified
        if isempty(PDE_comp.dom)
            PDE_comp.dom = zeros(0,2);
        end
        if size(PDE_comp.dom,2)~=2
            error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                    ' "dom" must be a px2 array specifying the interval on which each of the p spatial variables exists.'])
        elseif any(PDE_comp.dom(:,1) >= PDE_comp.dom(:,2))
            error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                    ' lower boundaries of spatial domain cannot exceed the upper boundaries.'])
        end
        % Check that the variables have been properly specified
        if size(PDE_comp.vars,2)>2
            error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                    ' "vars" should be a px2 array, specifying the p primary and dummy spatial variables.'])
        elseif size(PDE_comp.vars,1)~=size(PDE_comp.dom,1)
            if size(PDE_comp.vars,1)==1 && size(PDE_comp.vars,2)==size(PDE_comp.dom,1)
                % Only primary spatial variables are specified.
                PDE_comp.vars = PDE_comp.vars';
            else
                error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                        ' the number of specified variables in "vars" does not match the dimension of the spatial domain in "dom".'])
            end
        end
        % If variable names (rather than variables) have been specified,
        % convert to polynomial objects.
        if ~isempty(PDE_comp.vars) && ~ispvar(PDE_comp.vars) && ~iscellstr(PDE_comp.vars)
            error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                    ' spatial variables must be specified as pvar (polynomial class) objects.'])
        elseif iscellstr(PDE_comp.vars)
            PDE_comp.vars = polynomial(PDE_comp.vars);
        end
        % Make sure the specified variables are distinct.
        if ~isempty(PDE_comp.vars) && length(PDE_comp.vars.varname)~=prod(size(PDE_comp.vars)) % would use numel, but not available for polynomials...
            error(['Component ',obj,'{',num2str(ii),'} is not appropriately specified:'...
                    ' spatial variables must be distinct.'])
        end
    end
    % Add the proposed domain to the list
    dom_list = [dom_list; PDE_comp.dom];
    % Add the proposed variable names to  the list.
    % > Unfortunately, we have to extract the variable names
    % individually, as the order of vars.varname may not match the
    % order of the actual variables. <
    if ~isempty(PDE_comp.vars) && size(PDE_comp.vars,2)==2 && ~isdouble(PDE_comp.vars(2))
        % If dummy variables are available add these to the list as well.
        for pp=1:size(PDE_comp.vars,1)
            var_list = [var_list; [PDE_comp.vars(pp,1).varname, PDE_comp.vars(pp,2).varname, {true}]];
        end
    elseif ~isempty(PDE_comp.vars)
        % Otherwise, set variables as "0", and use "false" to indicate that
        % no dummy variables have been specified for these primary
        % variables.
        for pp=1:size(PDE_comp.vars,1)
            var_list = [var_list; [PDE_comp.vars(pp,1).varname, 0, {false}]];
        end
    end
    % Keep track of which component these variables correspond to.
    comp_list = [comp_list; ii*ones(size(PDE_comp.vars,1),1)];
    PDE.(obj){ii} = PDE_comp;
end
    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function PDE = set_vars(PDE,obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE = set_vars(PDE,obj)
% sets the properties of each component of the state, input, or output
% "obj", such as the variables it depends on.
%
% INPUTS:
% - PDE:        A "struct" or "pde_struct" class object defining a PDE.
% - obj:        A "char" object 'x', 'u', 'w', 'y', or 'z', specifying for
%               which objects we are assigning the vars.
%
% OUTPUTS:
% - PDE:        The same structure as the input "PDE", now with
%               "PDE.obj{ii}.vars", "PDE.obj{ii}.dom" 
%               specified for each appropriate element ii.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the global variables and their associated domain.
global_vars = PDE.vars;
global_dom = PDE.dom;

% Extract the table associated to the considered object.
Lobj_tab = PDE.([obj,'_tab']);

% Loop over the components, assigning the values of each field.
n_comps = numel(PDE.(obj));
for ii=1:n_comps
%     % First, set the size of the component
%     if isfield(PDE.(obj){ii},'size') && PDE.(obj){ii}.size ~= PDE.([obj,'_tab'])(ii,2)
%         error(['The specified size "',obj,'{',ii,'}.size" does not match the observed size of component "',obj,'{',num2str(ii),'};',...
%                 ' please adjust the specified size or the dimensions of the matrices "C" in the appropriate terms.'])
%     end
    
    % Next, check which of the global variables our LHS component (state
    % or output) actually depends on.
    has_vars_Lcomp = logical(Lobj_tab(ii,3:end));
    nvars_Lcomp = sum(has_vars_Lcomp);
    Lvars_new = global_vars(has_vars_Lcomp,:);
    
    % Set the variables and domain of the component
    if ~isfield(PDE.(obj){ii},'vars') || isempty(PDE.(obj){ii}.vars)
        var_order = (1:nvars_Lcomp);
    else
        % If variables have already been specified, we have to keep track
        % of the order in which they were originally specified.
        Lvars_old = PDE.(obj){ii}.vars(:,1);
        var_order = zeros(1,nvars_Lcomp);
        for kk=1:nvars_Lcomp
            % var_order is such that Lvars_new = Lvars_old(var_order)
            var_order(kk) = find(isequal(Lvars_new(kk,1),Lvars_old),1,'first');
        end
    end
    PDE.(obj){ii}.vars = Lvars_new;
    PDE.(obj){ii}.var_order = var_order;  
    PDE.(obj){ii}.dom = global_dom(has_vars_Lcomp,:); 
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [PDE,diff_tab] = get_diff_tab(PDE,obj,diff_tab,Gvar_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PDE,diff_tab] = get_diff_tab(PDE,obj,diff_tab)
% checks the different equations for each of the components of state, 
% output, or BC "obj", to establish the required order of differentiability
% for each state component x{ii}, returned as diff_tab.
%
% INPUTS:
% - PDE:        A "struct" or "pde_struct" class object defining a PDE.
% - obj:        A "char" object 'x', 'y', 'z', or 'BC', specifying for
%               which equations of the PDE we are checking the terms.
% - diff_tab:   A nx x p integer array specifying for each of the nx
%               elements of PDE.x (i.e. each state component) to what order
%               it must be differentiable with respect to each of the p
%               global variables PDE.vars.
%
% OUTPUTS:
% - PDE:        The same structure as the input "PDE", now with the terms
%               in the equations for "PDE.(obj){ii}" updated to explicitly
%               include orders of derivative if not specified (all orders
%               default to zero), as well as spatial positions (default to
%               interior of domain) if not specified. The tables x_tab
%               through BC_tab are also updated to include new sizes for
%               the different states, inputs, outputs, and BCs.
% - diff_tab:   A nx x p integer array specifying for each of the nx
%               elements of PDE.x (i.e. each state component) to what order
%               it must be differentiable with respect to each of the p
%               global variables PDE.vars. Orders of differentiability will
%               not decrease with respect to the input, but may increase if
%               in some term of some equation, a higher order spatial
%               derivative is taken of some state component than the
%               current order of differentiability for this state component
%               in the input "diff_tab".
%
% For example, if diff_tab(i,k) = 3, but a term involving 
% \partial_{s_k}^{5} x_i is encountered, diff_tab will be updated such that
% diff_tab(i,k) = 5. If then another term involving \partial_{s_k}^{2} x_i
% is encountered, diff_tab is not adjusted, as the maximal required order
% of differentiability does not change.
%
% NOTE: For any variable s\in[a,b], if x_i is differentiable wrt s up to
% order d, at s=a or s=b, it is taken to be differentiable only up to order
% d-1. For example, if a term involving partial_{s_k}^2 x_i(s_{k}) is 
% encountered, the order of differentiability will be set such that
% diff_tab(i,k)>=2. However, if a term involving 
% partial_{s_k}^2 x_i(s_{k}=a_{k}) is encountered, the order of 
% differentiability will be set such that diff_tab(i,k)>=3. This also means
% that, if a term involving simply x_i(s_k=a_k) is encountered, state
% component x_i is taken to be at least first order differentiable wrt s_k.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the number of global variables, and the table associated to the
% considered object x, y, z, or BC.
nvars = size(diff_tab,2);
Lobj_tab = PDE.([obj,'_tab']);

% % Loop over the equations for each component of "comp".
for ii=1:numel(PDE.(obj))
    % First assign an equation number to the component.
    PDE.(obj){ii}.eq_num = ii;
    % This is only for internal use, and will be removed at the end.
    
    
    % If an explicit order of differentiability has already been specified,
    % set this order in our diff_tab.
    if strcmp(obj,'x') && isfield(PDE.x{ii},'diff')
        % Establish the global variables on which the component depends.
        has_var_Lstate = logical(Lobj_tab(ii,3:2+nvars));
        nvars_Lstate = sum(has_var_Lstate);
        % Make sure the order of differentiability is appropriately
        % specified.
        if size(PDE.x{ii}.diff,1)~=1 && size(PDE.x{ii}.diff,2)==1
            PDE.x{ii}.diff = PDE.x{ii}.diff';
        end
        if isempty(PDE.x{ii}.diff) && nvars_Lstate==0
            PDE.x{ii}.diff = zeros(1,0);
        elseif size(PDE.x{ii}.diff,1)~=1
            error(['The order of differentiability "x{',num2str(ii),'}.diff" is not appropriately specified;',...
                    ' the field should be specified as a 1xp array indicating the order of the derivative in each of the p variables on which the considered state "x{',num2str(ii),'}" depends.'])
        elseif size(PDE.x{ii}.diff,2)==nvars
            % Orders of differentiation have been specified with respect to
            % all of the global variables.
            if any(PDE.x{ii}.diff(1,~has_var_Lstate))
                error(['State componen "x{',num2str(ii),'}" appears to be differentiable with respect to variables it does not depend on;',...
                        ' please specify orders of differentiability "x{',num2str(ii),'}.diff" only for variables on which the component depends.']);
            else
                % Retain derivative orders only for variables on which
                % the state actually depends.
                PDE.x{ii}.diff = PDE.x{ii}.diff(Gvar_order);    % account for re-ordering of vars
                PDE.x{ii}.diff = PDE.x{ii}.diff(1,has_var_Lstate);
            end
        elseif size(PDE.x{ii}.diff,2)~=nvars_Lstate
            % The number of derivative orders does not match the number
            % of variables.
            error(['The order of differentiability "x{',num2str(ii),'}.diff" is not appropriately specified;',...
                    ' please make sure the number of derivative orders matches the number of variables on which the considered state component depends.']);
        else
            PDE.x{ii}.diff = PDE.x{ii}.diff(PDE.x{ii}.var_order);    % account for re-ordering of vars
        end
        % Update the order of differentiability.
        diff_tab(ii,has_var_Lstate) = max(diff_tab(ii,has_var_Lstate),PDE.x{ii}.diff);
    end
    
    % Also check if a size has been specified for the component.
    if isfield(PDE.(obj){ii},'size')
        Lobj_tab(ii,2) = PDE.(obj){ii}.size;
    end
    
    % % Next, loop over the terms in the PDE, and check the order of the
    % % derivative of the state component considered in each term.
    if ~isfield(PDE.(obj){ii},'term')
        warning(['No equation has been specified for ',obj,'{',num2str(ii),'};',...
                    ' to specify terms in the differential equation for this component, use the field ',obj,'{',num2str(ii),'}.term']);
        PDE.(obj){ii}.term = cell(0);
        continue
    else
        % Make sure the terms are specified as a 1xn cell.
        PDE.(obj){ii}.term = PDE.(obj){ii}.term(:)';
    end
    for jj=1:numel(PDE.(obj){ii}.term)
        % Assign a term number to the term (only for internal use in the
        % "check_terms" function), and keep track of the term name for
        % error messages and warnings.
        PDE.(obj){ii}.term{jj}.term_num = jj;
        term_jj = PDE.(obj){ii}.term{jj};
        term_name = [obj,'{',num2str(ii),'}.term{',num2str(jj),'}'];
        
        % Establish whether the term involves a state component or input.
        if ~isfield(term_jj,'x') && ~isfield(term_jj,'w') && ~isfield(term_jj,'u')
            if numel(PDE.x)==1
                % If there is only one state component, assume the term
                % involves this state.
                term_jj.x = 1;
            else
                error(['It is unclear which state component or input term ',obj,'{',num2str(ii),'}.term{',num2str(jj),'} involves;',...
                        ' please explicitly specify a state component through "term{',num2str(jj),'}.x", or an input through "term{',num2str(jj),'}.w" or "term{',num2str(jj),'}.u".'])
            end
        elseif (isfield(term_jj,'x') + isfield(term_jj,'w') + isfield(term_jj,'u'))~=1
            error(['Term ',obj,'{',num2str(ii),'}.term{',num2str(jj),'} appears to involve both a state component and an input;',...
                    ' please specify at most one of the fields "term{',num2str(jj),'}.x", "term{',num2str(jj),'}.w", or "term{',num2str(jj),'}.u".'])
        end
        if isfield(term_jj,'x')
            Robj = 'x';
            is_x_Rcomp = true;
        elseif isfield(term_jj,'w')
            Robj = 'w';
            is_x_Rcomp = false;
        else
            Robj = 'u';
            is_x_Rcomp = false;
        end
        
        % Determine which state component or input the term depends on.
        Rindx = term_jj.(Robj)(:);
        nR = length(Rindx);
        if ~isa(Rindx,'double')
            error(['Field "',term_name,'.',Robj,'" is not appropriately specified;',...
                    ' please provide a single positive integer to indicate which component the considered term involves.'])
        end
        if isempty(Rindx) && numel(PDE.(Robj))==1
            % If no compoennt is specified, but only one component of the
            % considered kind (state or input) is available, assume this to
            % be the desired component.
            Rindx = 1;
            nR = 1;
        elseif isempty(Rindx) || any(Rindx<=0)
            error(['Field "',term_name,'.',Robj,'" is not appropriately specified;',...
                    ' please provide a single positive integer to indicate which component the considered term involves.'])
        elseif any(Rindx>numel(PDE.(Robj)))
            error(['The component "',term_name,'.',Robj,'" does not exist;',...
                    ' please initialize the component by setting "PDE.',Robj,'{',num2str(jj),'}.vars" and "PDE.',Robj,'{',num2str(jj),'}.dom".'])
        end
        
        % Establish which of the global variables each component depends on.
        has_vars_Rcomp = PDE.([Robj,'_tab'])(Rindx,3:end);
        if nR>1
            % If multiple components are specified, make sure that they 
            % depend on the same spatial variables.
            if any(any(has_vars_Rcomp(2:end,:) - has_vars_Rcomp(1:end-1,:)))
                error(['Term "',term_name,'" is not appropriate;',...
                        ' inclusion of multiple components ".',Robj,'" is supported only if all components vary in the same spatial variables.'])
            else
                has_vars_Rcomp = has_vars_Rcomp(1,:);
            end
        end
        has_vars_Rcomp = logical(has_vars_Rcomp);
        nvars_Rcomp = sum(has_vars_Rcomp);
        Rvar_order = zeros(nR,nvars_Rcomp);
        for kk=1:size(Rindx,1)
            Rvar_order(kk,:) = PDE.(Robj){Rindx(kk)}.var_order;  % new_vars = old_vars(var_order)
        end
        
        % % Finally, if a derivative is specified, update the maximal
        % % order of differentiability of the involved state component.
        % First check if derivatives are appropriately specified.
        if is_x_Rcomp && (~isfield(term_jj,'D') || isempty(term_jj.D))
            % If no derivative is specified, assume no derivative is
            % desired (default to 0).
            Dval = zeros(nR,nvars_Rcomp);
        elseif is_x_Rcomp
            % Otherwise, check that the number of orders matches the number
            % of variables.
            Dval = term_jj.D;
            if ~isa(Dval,'double') || any(Dval(:)<0)
                error(['The order of differentiation "',obj,'{',num2str(ii),'}.term{',num2str(jj),'}.D" is not appropriately specified;',...
                        ' derivative orders must be specified as a type "double" array of nonnegative integers.']);
            end
            if size(term_jj.D,2)==nvars
                % Derivatives have been specified with respect to all of
                % the global variables.
                if any(any(term_jj.D(:,~has_vars_Rcomp)))
                    error(['Term "',obj,'{',num2str(ii),'}.term{',num2str(jj),'}" seems to differentiate state component x{',num2str(Rindx),'} with respect to a variable it does not depend on;',...
                            ' please specify orders of the derivative only for variables on which the component depends.']);
                else
                    % Retain derivative orders only for variables on which
                    % the state actually depends.
                    Dval = Dval(:,Gvar_order);
                    Dval = Dval(:,has_vars_Rcomp);
                end
            elseif size(term_jj.D,2)~=nvars_Rcomp
                % The number of derivative orders does not match the number
                % of variables
                error(['The order of differentiation "',obj,'{',num2str(ii),'}.term{',num2str(jj),'}.D" is not appropriately specified;',...
                        ' please make sure the number of derivative orders matches the number of variables on which the considered state component depends.']);
            else
                if size(Rvar_order,1)==1
                    Dval = Dval(:,Rvar_order);
                elseif size(Dval,1)==1
                    Dval = repmat(Dval,[size(Rvar_order,1),1]);
                    for kk=1:size(Rvar_order,1)
                        Dval(kk,:) = Dval(kk,Rvar_order(kk,:));
                    end
                elseif size(Dval,1)==size(Rvar_order,1)
                    for kk=1:size(Rvar_order,1)
                        Dval(kk,:) = Dval(kk,Rvar_order(kk,:));
                    end
                end                    
            end            
            % Also make sure the number of specified derivatives matches
            % the number of proposed states.
            if size(Dval,1)==1 && nR>1
                Dval = repmat(Dval,[nR,1]);
            elseif size(Dval,1)>1 && nR==1
                nR = size(term_jj.D,1);
                Rindx = repmat(Rindx,[nR,1]);
            elseif size(Dval,1)~=nR
                error(['The number of state components in "',term_name,'.x" should match the number of derivatives "',term_name,'.D".'])
            end
        else
            % The term involves an input --> no derivative is allowed.
            if isfield(term_jj,'D') && any(term_jj.D)
                error(['Term "',term_name,'" is not appropriate;',...
                        ' derivatives of inputs are not supported.'])
            elseif isfield(term_jj,'D')
                PDE.(obj){ii}.term{jj} = rmfield(PDE.(obj){ii}.term{jj},'D');
            end
        end
        
        % Next, check if a spatial position is appropriately specified.
        if is_x_Rcomp && (~isfield(term_jj,'loc') || isempty(term_jj.loc))
            % % For BCs, a spatial position MUST be specified.
            %if strcmp(obj,'BC') && any(has_vars_Rcomp) && ~isfield(term_jj,'I') && ~isfield(term_jj,'int')
            %    error(['BC term "',term_name,'" is not appropriately specified;',...
            %            ' a spatial position "',term_name,'.loc" at which to evaluate the state is required when specifying BCs.'])
%             if strcmp(obj,'BC') && isfield(term_jj,'I') || isfield(term_jj,'int')
%                 % If a full integral is used, it's okay if no position is
%                 % specified.
%                 if isfield(term_jj,'int') && ~isfield(term_jj,'I')
%                     term_jj.I = term_jj.int;
%                     PDE.(obj){ii}.term{jj}.I = term_jj.I;
%                     PDE.(obj){ii}.term{jj} = rmfield(PDE.(obj){ii}.term{jj},'int');
%                 end
%                 if ~isa(term_jj.I,'cell')
%                     error(['BC term "',term_name,'" is not appropriately specified;',...
%                             ' the field "',term_name,'.I" should be specified as a cell with each element describing the domain of integration along the associated variable.'])
%                 end
%                 % Check whether full integral is taken along any spatial
%                 % direction.
%                 if numel(term_jj.I)~=nvars_Rcomp
%                     error(['The integral "',term_name,'.I" is not properly specified;',...
%                             ' the number of elements should match the number of variables on which the component "',term_name,'.',Robj,'{',num2str(Rindx),'}" depends.']);
%                 end
% %                 use_full_int = false;
% %                 for kk=1:numel(term_jj.I(Rvar_order))
% %                     if ~isempty(term_jj.I{kk}) && (isa(term_jj.I{kk},'double') || (isa(term_jj.I{kk},'polynomial') && isdouble(term_jj.I{kk})))
% %                         use_full_int = true;
% %                         break
% %                     end
% %                 end
% %                 if ~use_full_int
% %                     error(['BC term "',term_name,'" is not appropriately specified;',...
% %                         ' a spatial position "',term_name,'.loc" at which to evaluate the state is required when specifying BCs.'])
% %                 end
% %             elseif strcmp(obj,'BC')
% %                 PDE.(obj){ii}.term{jj}.loc = zeros(1,0);
%             end
            % Otherwise, if no position is specified, evaluate the state on 
            % the interior of the domain.
            Rloc = repmat(PDE.vars(has_vars_Rcomp,1)',[nR,1]);
        elseif is_x_Rcomp
            % If a position is specified, make sure that it is appropriate.
            Rloc = term_jj.loc;
            if ~isa(Rloc,'double') && ~isa(Rloc,'polynomial')
                error(['The spatial position "',term_name,'.loc" is not appropriately specified;',...
                        ' the position should be specified as an array of type "double" or "polynomial".'])
            end
            % Check that the dimension matches the number of spatial
            % variables on which the considered state component depends.
            if size(Rloc,2)==nvars
                % Location specified for global vars --> retain only Rvars
                Rloc = Rloc(:,Gvar_order);
                Rloc = Rloc(:,has_vars_Rcomp);
            elseif size(Rloc,2)~=nvars_Rcomp
                error(['The spatial position "',term_name,'.loc" is not appropriately specified;',...
                        ' the number of elements should match the number of variables on which "',term_name,'.',Robj,'" depends.'])
            else
                if size(Rvar_order,1)==1
                    Rloc = Rloc(:,Rvar_order);
                elseif size(Rloc,1)==1
                    Rloc = repmat(Rloc,[size(Rvar_order,1),1]);
                    for kk=1:size(Rvar_order,1)
                        Rloc(kk,:) = Rloc(kk,Rvar_order(kk,:));
                    end
                elseif size(Rloc,1)==size(Rvar_order,1)
                    for kk=1:size(Rvar_order,1)
                        Rloc(kk,:) = Rloc(kk,Rvar_order(kk,:));
                    end
                end                
            end
            % Also make sure the number of positions matches the number of
            % considered state components.
            if size(Rloc,1)==1 && nR>1
                Rloc = repmat(Rloc,[nR,1]);
            elseif size(Rloc,1)>1 && nR==1
                % Repeat the Rindx and derivative orders if the user asks
                % the state to be evaluated at multiple positions.
                Rloc = term_jj.loc;
                nR = size(Rloc,1);
                Rindx = repmat(Rindx,[nR,1]);
                Dval = repmat(Dval,[nR,1]);
            elseif size(Rloc,1)>1 && nR>1
                error(['The number of state components in "',term_name,'.x" should match the number of derivatives "',term_name,'.D".'])
            end
            Rloc = polynomial(Rloc);
        elseif isfield(term_jj,'loc')
            % The term involves an input --> no spatial position may be 
            % specified.
            if  (length(term_jj.loc)~=nvars && length(term_jj.loc)~=nvars_Rcomp) || ...
                (length(term_jj.loc)==nvars && ~all(isequal(polynomial(term_jj.loc(Gvar_order)),PDE.vars(:,1)') | isequal(polynomial(term_jj.loc(Gvar_order)),PDE.vars(:,2)'))) || ...
                (length(term_jj.loc)==nvars_Rcomp && ~all(isequal(polynomial(term_jj.loc(Rvar_order)),PDE.vars(has_vars_Rcomp,1)') | isequal(polynomial(term_jj.loc(Rvar_order)),PDE.vars(has_vars_Rcomp,2)')))
                % A spatial position has been provided which is not on the
                % interior of the domain.
                error(['Term "',term_name,'" is not appropriate;',...
                        ' no spatial position ".loc" cannot be specified for terms involving an input.'])
            elseif isfield(term_jj,'loc')
                PDE.(obj){ii}.term{jj} = rmfield(PDE.(obj){ii}.term{jj},'loc');
                continue
            end
        end
        
        % Also, if a coefficient matrix is provided, use this matrix to
        % establish a size of the LHS component, and RHS component.
        if isfield(term_jj,'C')
            if ~isa(term_jj.C,'double') && ~isa(term_jj.C,'polynomial')
                error(['The coefficients "',term_name,'.C" are not appropriately specified;',...
                        ' coefficients should be specified as an array of type "double" or "polynomial".'])
            elseif isfield(PDE.(obj){ii},'size') && PDE.(obj){ii}.size~=size(term_jj.C,1)
                error(['The number of rows of the matrix "',term_name,'.C" does not match the specified size "',obj,'{',num2str(ii),'}.size".'])
            else
                Lobj_tab(ii,2) = size(term_jj.C,1);
            end
            % If only a single RHS component is provided, the number of
            % columns of C should match the size of this component.
            if nR==1 && isfield(PDE.(Robj){Rindx},'size') && PDE.(Robj){Rindx}.size~=size(term_jj.C,2)
                error(['The number of columns of the matrix "',term_name,'.C" does not match the specified size "',Robj,'{',num2str(Rindx),'}.size".'])
            elseif nR==1
                PDE.([Robj,'_tab'])(Rindx,2) = size(term_jj.C,2);
            elseif nR==size(term_jj.C,2)
                PDE.([Robj,'_tab'])(Rindx,2) = 1;
            end
        end
        
        if ~is_x_Rcomp
            % At this point, if the term does not involve a state
            % component, we can move on to the next equation.
            continue
        end
        
        % Finally, update the order of differentiability of the considered
        % state components.
        for kk=1:nR
            % If the state is evaluated at a boundary, the order of
            % the derivative can be 1 less than the maximal order of 
            % differentiability.
            isboundary_var = false(1,nvars_Rcomp);
            for ll=1:nvars_Rcomp
                isboundary_var(ll) = isdouble(Rloc(kk,ll));
            end
            diff_order_kk = Dval(kk,:) + isboundary_var;
            diff_tab(Rindx(kk),has_vars_Rcomp) = max(diff_tab(Rindx(kk),has_vars_Rcomp),diff_order_kk);
        end
        % Update the specified term.
        PDE.(obj){ii}.term{jj}.(Robj) = Rindx;
        PDE.(obj){ii}.term{jj}.D = Dval;
        PDE.(obj){ii}.term{jj}.loc = Rloc;
    end
end

% Store the updated table with potentially new size information for the
% state, output, or BCs.
PDE.([obj,'_tab']) = Lobj_tab;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [PDE,diff_tab] = check_terms(PDE,obj,diff_tab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PDE,diff_tab] = check_terms(PDE,obj,diff_tab)
% checks the different terms in the different equations for each of the 
% components of the state, output, or BC "obj", to check that all mandatory
% fields have been properly specified, and to assign a default value to all
% optional fields that the user has not specified.
%
% INPUTS:
% - PDE:        A "struct" or "pde_struct" class object defining a PDE.
% - obj:        A "char" object 'x', 'y', or 'z', specifying for
%               which equations of the PDE we are checking the terms.
% - diff_tab:   A nx x p integer array specifying for each of the nx
%               elements of PDE.x (i.e. each state component) to what order
%               it must be differentiable with respect to each of the p
%               global variables PDE.vars.
%
% OUTPUTS:
% - PDE:        The same structure as the input "PDE", now with all fields
%               in the structures "PDE.(obj){ii}.term{jj}" for all ii and
%               jj specified. In addition, any term corresponding to
%               multiple terms in the actual equation will be split up into
%               separate terms, for the PDEs, if higher order temporal
%               derivatives are requested, additional state components and
%               PDEs will be added to represent the system using only
%               first-order temporal derivatives.
% - diff_tab:   Similar array as in the input, only with potential
%               additional rows. These rows correspond to the state
%               components added to account for higher order temporal
%               derivatives in the PDE, but the spatial differentiability
%               of all of these added state components will be zero wrt all
%               spatial variables.
%
% For example, if \partial_t^3 x_i = eq_i(x_1,...,x_{nx}) is specified, 
% then this will be split into
%   \partial_t x_i      = x_{nx+1}
%   \partial_t x_{nx+1} = x_{nx+2}
%   \partial_t x_{nx+2} = eq_i(x_1,...,x_{nx})
% where the added state components x_{nx+1} and x_{nx+2} are assumed not to
% be differentiable with respect to any spatial variable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the global variables and their associated domain.
global_vars = PDE.vars;
global_dom = PDE.dom;
nvars = size(global_vars,1);

% Extract the table associated to the considered object x, y, or z.
Lobj_tab = PDE.([obj,'_tab']);

% Loop over the equations, noting that equations may be added if higher 
% order temporal derivatives are requested.
n_eqs = numel(PDE.(obj));
ii = 1;
while ii<=n_eqs
    
    if isfield(PDE.(obj){ii},'eq_num')
        eq_num_str = num2str(PDE.(obj){ii}.eq_num);
    else
        eq_num_str = num2str(ii);
    end
    
    % Set the size of the object.
    PDE.(obj){ii}.size = PDE.([obj,'_tab'])(ii,2);
    
    % First, check which of the global variables our LHS component (state
    % or output) actually depends on.
    has_vars_Lcomp = logical(Lobj_tab(ii,3:end));

    % For the state components, also check the order of differentiability,
    % and requested temporal derivative.
    if strcmp(obj,'x')
        % Check if an order of differentiability in each of the spatial
        % variables has been specified, and set it otherwise.
        if isfield(PDE.x{ii},'diff') && any(PDE.x{ii}.diff~=diff_tab(ii,has_vars_Lcomp))
            warning(['The specified order of differentiability "PDE.x{',eq_num_str,'}.diff" is smaller than the observed order of this state component in the PDE, and will therefore be increased.',...
                        ' If you wish to retain the original order of differentiability, please make sure that the order of any derivative of the state does not exceed this specified order.']);
        end
        PDE.x{ii}.diff = diff_tab(ii,has_vars_Lcomp);
        
        % Check if a temporal derivative is specified.
        if ~isfield(PDE.x{ii},'tdiff')
            PDE.x{ii}.tdiff = 1;
        elseif PDE.x{ii}.tdiff>1
            % Higher order temporal derivatives will need to be processed
            % before conversion to a PIE is possible.
            % This can be done using the function "expand_tderivatives".
            PDE.has_hotd = true;
        end
    end
    
    % % Now we dive into the equation associated to the considered 
    % % component.    
    % Loop over the terms, and make sure they are properly specified.
    n_terms = numel(PDE.(obj){ii}.term);
    jj = 1;
    while jj<=n_terms
        term_jj = PDE.(obj){ii}.term{jj};
        % For the error messages, keep track of the name of this term.
        % We use "eq_num" and "term_num" since the number of equations and
        % terms may change as we initialize.
        if isfield(PDE.(obj){ii}.term{jj},'term_num')
            term_num = PDE.(obj){ii}.term{jj}.term_num;
        else
            term_num = jj;
        end
        term_name = [obj,'{',eq_num_str,'}.term{',num2str(term_num),'}'];
        
        
        % % % Assing default values to the empty fields.
        % % First check which state component or input is involved.
        if isfield(term_jj,'x')
            Robj = 'x';
            is_x_Robj = true;
        elseif isfield(term_jj,'w')
            Robj = 'w';
            is_x_Robj = false;
        else
            Robj = 'u';
            is_x_Robj = false;
        end
        Rindx = term_jj.(Robj)(:);
        nR = length(Rindx);
        % Establish which spatial variables each component depends on.
        has_vars_Rcomp = logical(PDE.([Robj,'_tab'])(Rindx(1),3:end));
        nvars_Rcomp = sum(has_vars_Rcomp);
        % Note that we already checked all the different components listed
        % in Rindx to depend on the same variables in "get_diff".
        
        % % Extract a derivative (was already initialized in "get_diff").
        if is_x_Robj && (~isfield(term_jj,'D') || isempty(term_jj.D))
            % If no derivative is specified, assume no derivative is
            % desired.
            Dval = zeros(nR,nvars_Rcomp);
        elseif is_x_Robj
            Dval = term_jj.D;
        end
        
        % % Extract a spatial position (was already initialized in 
        % % "get_diff").
        if is_x_Robj && (~isfield(term_jj,'loc') || isempty(term_jj.loc))
            % If no position is specified, evaluate state on the interior
            % of the domain.
            Rloc = repmat(global_vars(has_vars_Rcomp,1)',[nR,1]);
        elseif is_x_Robj
            Rloc = term_jj.loc;
        end
        
        % % Initialize delay.
        if ~isfield(term_jj,'delay')
            delay = 0;
            delay_varname = cell(0,1);
        else
            delay = term_jj.delay;
            if any(size(delay)~=[1,1]) || ~(isa(delay,'double') || isdouble(delay) || ispvar(delay))
                error(['The delay "',term_name,'.delay" is not properly specified;',...
                        ' please specify a scalar double or pvar value.'])
            end
            if isa(delay,'double') || isdouble(delay)
                delay = abs(double(delay));
                delay_varname = cell(0,1);     % The coefficients are allowed to depend on the delay
                if delay~=0
                    PDE.has_delay = true;
                end
            elseif ispvar(delay) && length(delay.varname)>1
                error(['The delay "',term_name,'.delay" is not properly specified;',...
                        ' distributed delays should be specified using a single scalar variable.'])
            elseif ~ismember(delay.varname{1},PDE.tau.varname)
                error(['The delay "',term_name,'.delay" is not properly specified;',...
                        ' any distributed delay variable should be specified in the field "PDE.tau", along with the domain on which it exists.'])
            else
                delay_varname = delay.varname;
                PDE.has_delay = true;
            end
        end
        % We won't be doing anything with the delay here, we process delays
        % in a different function "expand_delays".
        PDE.(obj){ii}.term{jj}.delay = delay;
        
        % % Initialize coefficients.
        if ~isfield(term_jj,'C') || isempty(term_jj.C)
            % If no coefficient are specified, assume each term is scaled
            % with a factor 1.
            if any(PDE.([Robj,'_tab'])(Rindx,2)~=Lobj_tab(ii,2))
                error(['The size of the component "',obj,'{',eq_num_str,'}" does not match that of the components specified in "',term_name,'.',Robj,'";',...
                        ' please specify a matrix "',term_name,'.C" describing how each component on the right-hand side contributes to the equation for "',obj,'{',eq_num_str,'}".'])
            else
                Cval = repmat(eye(Lobj_tab(ii,2)),[1,nR]);
            end
        else
            % Otherwise, verify that the size of the matrix makes sense.
            Cval = term_jj.C;
            if ~isa(Cval,'double') && ~isa(Cval,'polynomial')
                error(['The coefficients "',term_name,'.C" are not appropriately specified;',...
                        ' coefficients should be specified as an array of type "double" or "polynomial".'])
            elseif size(Cval,1)~=Lobj_tab(ii,2)
                error(['The number of rows of the coefficient matrix "',term_name,'.C" does not match the number of equations "',obj,'{',eq_num_str,'}.size" it contributes to.'])
            elseif size(Cval,2)~=sum(PDE.([Robj,'_tab'])(Rindx,2))
                error(['The number of columns of the coefficient matrix "',term_name,'.C" does not match the number of state variables it maps.'])
            end
        end
        
        % % Initialize an empty integral (for now) if no integral is 
        % % specified.
        if isfield(term_jj,'int') && ~isfield(term_jj,'I')
            term_jj.I = term_jj.int;
            PDE.(obj){ii}.term{jj} = rmfield(PDE.(obj){ii}.term{jj},'int');
        end
        if isfield(term_jj,'I')
            Ival = term_jj.I(:);
            if ~isa(Ival,'cell')
                error(['The field "',term_name,'.I" is not appropriately specified;',...
                        ' the field should be a cell with each element k providing the desired domain of integration in variable k of component "',Robj,'{',num2str(Rindx),'}" depends.'])
            end
        else
            Ival = cell(0);
        end
        
        % % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % % % If multiple terms are included, split into separate elements.
        C_colnums = cumsum([0;PDE.([Robj,'_tab'])(Rindx,2)]);
        if nR>1
            % % We squeeze the new terms between the current one and the
            % % following ones.
            PDE.(obj){ii}.term = [PDE.(obj){ii}.term(1:jj), repmat(PDE.(obj){ii}.term(jj),[1,nR-1]), PDE.(obj){ii}.term(jj+1:end)];
            if is_x_Robj
                % Make sure each term indeed describes just one term.
                for kk=1:nR-1
                    C_cols = C_colnums(kk+1)+1 : C_colnums(kk+2);
                    PDE.(obj){ii}.term{jj+kk}.(Robj) = Rindx(kk+1);
                    PDE.(obj){ii}.term{jj+kk}.D = Dval(kk+1,:);
                    PDE.(obj){ii}.term{jj+kk}.loc = Rloc(kk+1,:);
                    PDE.(obj){ii}.term{jj+kk}.C = Cval(:,C_cols);
                    PDE.(obj){ii}.term{jj+kk}.I = Ival;
                    PDE.(obj){ii}.term{jj+kk}.delay = delay;
                    PDE.(obj){ii}.term{jj+kk}.term_num = term_num;
                end
            else
                % Make sure each term indeed describes just one term.
                for kk=1:nR-1
                    C_cols = C_colnums(kk+1)+1 : C_colnums(kk+2);
                    PDE.(obj){ii}.term{jj+kk}.(Robj) = Rindx(kk+1);
                    PDE.(obj){ii}.term{jj+kk}.C = Cval(:,C_cols);
                    PDE.(obj){ii}.term{jj+kk}.I = Ival;
                    PDE.(obj){ii}.term{jj+kk}.delay = delay;
                    PDE.(obj){ii}.term{jj+kk}.term_num = term_num;
                end
            end
            % Update the number of terms.
            n_terms = n_terms + nR-1;
        end
        % % Adjust the current element "term{jj}" to describe only the 
        % % first term.
        Rindx = Rindx(1);
        sz_R = PDE.([Robj,'_tab'])(Rindx,2);
        Cval = Cval(:,1:sz_R);
        if is_x_Robj
            Dval = Dval(1,:);
            Rloc = Rloc(1,:);
            PDE.(obj){ii}.term{jj}.(Robj) = Rindx;
            PDE.(obj){ii}.term{jj}.D = Dval;
            PDE.(obj){ii}.term{jj}.loc = Rloc;
            PDE.(obj){ii}.term{jj}.C = Cval;
            PDE.(obj){ii}.term{jj}.I = Ival;
            PDE.(obj){ii}.term{jj}.delay = delay;
        else
            PDE.(obj){ii}.term{jj}.(Robj) = Rindx;
            PDE.(obj){ii}.term{jj}.C = Cval;
            PDE.(obj){ii}.term{jj}.I = Ival;
            PDE.(obj){ii}.term{jj}.delay = delay;
        end
        
        
        % % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % % % Having initialized all the fields, and ensured that the
        % % % current term truly describes just one term, we now make sure
        % % % everything is properly specified, and fill in any remaining
        % % % gaps (most notably, the integrals).
        
        % % First, make sure that the specified location is appropriate.
        % Extract the variables associated to the RHS component, and their
        % associated domains
        Rvars = global_vars(has_vars_Rcomp,:);
        Rdom = global_dom(has_vars_Rcomp,:);
        Rvar_order = PDE.(Robj){Rindx}.var_order;   % new_vars = old_vars(var_order)
        % Keep track of which variables the term actually varies in
        % (without integration).
        isvariable_term_jj = true(1,nvars_Rcomp);   % Element kk will be false if variable kk is evaluated at a boundary.
        if is_x_Robj
            for kk=1:nvars_Rcomp
                % For each variable s\in[a,b], check that it is evaluated
                % at either s=s, s=a or s=b. 
                Rloc_kk = Rloc(kk);
                if isa(Rloc_kk,'double') || isa(Rloc_kk,'polynomial') && isdouble(Rloc_kk)
                    % s=a or s=b.
                    Rloc_kk = double(Rloc_kk);
                    if ~ismember(Rloc_kk,Rdom(kk,:))
                        error(['The spatial position "',term_name,'.loc" is not appropriate;',...
                        ' spatial position must match a boundary of the domain, or the appropriate spatial variable.'])
                    else
                        % Indicate that the term does not vary along this
                        % dimension.
                        isvariable_term_jj(kk) = false;
                    end
                elseif ispvar(Rloc_kk)
                    % s=s.
                    if ~isequal(Rloc_kk,Rvars(kk,1)) && ~isequal(Rloc_kk,Rvars(kk,2))
                        error(['The spatial position "',term_name,'.loc" is not appropriate;',...
                                ' for any variable s in [a,b], PIETOOLS can only evaluate s=a, s=b, or s=s.'])
                    else
                        % Replace any dummy variables with primary
                        % variables.
                        Rloc_kk = Rvars(kk,1);
                        PDE.(obj){ii}.term{jj}.loc(kk) = Rloc_kk;
                    end
                else
                    error(['The spatial position "',term_name,'.loc" is not appropriate;',...
                                ' spatial position must be specified as "double" or "pvar" (polynomial class) object.'])                
                end
            end
        end
        
        % % Establish which variables the term depends on but the LHS does
        % % not, and must therefore be integrated out.
        must_int_var = ~(has_vars_Lcomp(has_vars_Rcomp)) & isvariable_term_jj;
        % Keep track of along which variables integration is actually
        % performed.
        is_int_var = must_int_var;
        if ~isfield(term_jj,'I') || isempty(term_jj.I)
            % If no integral is specified, but an integral is required,
            % throw an error.
            if any(must_int_var)
                error(['The term "',term_name,'" appears to depend on a variable that the component "',obj,'{',num2str(ii),'}" does not depend on. This is not suppoprted.',...
                        ' Please use the field "',term_name,'.I" to introduce an integral, or the field "',term_name,'.loc" to evaluate the term at a particular, so as to remove the dependence on any "illegal" variables.'])
            end
            % Otherwise, if no integral is specified, assume that no
            % integral is  desired.
            Ival = cell(nvars_Rcomp,1);
            %Ival(must_int_var,1) = mat2cell(Rdom(must_int_var,:),ones(sum(must_int_var),1));    % Integrate over full domain
            PDE.(obj){ii}.term{jj}.I = Ival;
        else
            % If integrals are specified, make sure that they are
            % appropriate.
            if numel(Ival)>nvars_Rcomp
                error(['The field "',term_name,'.I" is not appropriately specified;',...
                        ' the number of elements should match the number of spatial variables on which component "',term_name,'.',Robj,'{',num2str(Rindx),'}" depends.'])
            end
            % If insufficient domains are specified, assume remaining
            % domains are empty.
            if numel(Ival)<nvars_Rcomp
                Ival = [Ival;cell(nvars_Rcomp-numel(Ival),1)];
            end
            Ival = Ival(Rvar_order);    % Reorder based on new ordering of variables.
            for kk = 1:numel(Ival)
                % Check for all spatial dimensions that the integral is
                % properly specified and appropriate.
                if isempty(Ival{kk}) && ~must_int_var(kk)
                    continue
                elseif isempty(Ival{kk}) 
                    % If integration must be performed, but is not
                    % specified, throw an error
                    error(['The term "',term_name,'" appears to depend on a variable that the component "',obj,'{',eq_num_str,'}" does not depend on. This is not suppoprted.',...
                            ' Please use the field "',term_name,'.I" to introduce an integral, or the field "',term_name,'.loc" to evaluate the state at a particular, so as to remove the dependence on any "illegal" variables.'])
                    %Ival{kk} = Rdom(kk,:);
                elseif ~isempty(Ival{kk}) && ~isvariable_term_jj(kk)
                    % Avoid integration in variables along which the
                    % considered component does not depend.
                    error(['The proposed integral "',term_name,'.I{',num2str(kk),'}" is not allowed;',...
                            ' the considered component "',term_name,'.',Robj,'{',num2str(Rindx),'}" at the specified position "',term_name,'.loc" does not vary along the proposed dimension of integration.'])
                end
                if (~isa(Ival{kk},'double') && ~isa(Ival{kk},'polynomial')) || ~all(size(Ival{kk})==[1,2])
                    error(['The field "',term_name,'.I" is not appropriately specified;',...
                            ' each element k must be a 1x2 array of type "double" or "polynomial", indicating the domain of integration in variable k of component "',Robj,'{',num2str(Rindx),'}" depends.'])
                end
                Ival_kk = polynomial(Ival{kk});
                % For variables we have to get rid of, we must integrate
                % along the entire domain.
                if must_int_var(kk) && ~all(isequal(Ival_kk,Rdom(kk,:)))
                    error(['The proposed integral "',term_name,'.I{',num2str(kk),'}" is not appropriate;',...
                            ' integration over the full domain of variable "',Rvars(kk,1).varname{1},'" is required as component "',obj,'{',eq_num_str,'}" does not depend on this variable.'])
                end
                % Check that the integral is taken over one of the allowed
                % domains.
                if all(isequal(Ival_kk,polynomial(Rdom(kk,2:-1:1)))) || all(isequal(Ival_kk,[Rvars(kk,1),Rdom(kk,1)])) || all(isequal(Ival_kk,[Rdom(kk,2),Rvars(kk,1)]))
                    % Integration is performed over mirror image of allowed
                    % domain (for whatever reason);
                    Cval = -Cval;
                    PDE.(obj){ii}.term{jj}.C = Cval;
                    Ival{kk} = Ival{kk}(1,end:-1:1);
                elseif ~(all(isequal(Ival_kk,polynomial(Rdom(kk,:)))) || all(isequal(Ival_kk,[Rdom(kk,1),Rvars(kk,1)])) || all(isequal(Ival_kk,[Rvars(kk,1),Rdom(kk,2)])))
                    error(['The proposed integral "',term_name,'.I{',num2str(kk),'}" is not allowed;',...
                            ' for a variable s in [a,b], integration is only allowed over [a,s], [s,b], or [a,b].'])
                end   
                is_int_var(kk) = true;  % Indicate that we are integrating along this spatial direction.
            end
            PDE.(obj){ii}.term{jj}.I = Ival;       
        end
        
        % % Finally, we check that the coefficients are appropriately
        % % specified.
        if isa(Cval,'polynomial') && ~isdouble(Cval)
            % Establish variable names that the term is allowed to have:
            
            % Cval can have variables that appear in the LHS component.
            Lcomp_varname = global_vars(has_vars_Lcomp,1).varname;            
            % Currently, opvar and opvar2d objects only use dummy variables
            % in partial integrals, not integrals over the full domain. As
            % such, dummy variables in Cval corresponding to spatial
            % dimensions that are integrated out should be
            % replaced with primary variables as well.
            if any(must_int_var)
                Cval = subs(Cval,Rvars(must_int_var,2),Rvars(must_int_var,1));
                Rcomp_varname_1 = Rvars(must_int_var,1).varname;
            else
                Rcomp_varname_1 = {};
            end
            % For variables that both the LHS component and RHS component
            % share, but for which a (partial) integral is used, a dummy 
            % variable may appear in Cval.
            Rcomp_varname_2 = Rvars(is_int_var & ~must_int_var,2).varname;
            
            % Check that the polynomial does not depend on any "illegal"
            % variables.
            if any(~ismember(Cval.varname,[Lcomp_varname; Rcomp_varname_1; Rcomp_varname_2; delay_varname]))
                error(['The proposed polynomial "',term_name,'.C" is not appropriate;',...
                        ' the (spatial) variables it depends on do not make sense with the variables of the component "',obj,'{',num2str(ii),'}"  it contributes to, and/or the component "',Robj,'{',num2str(Rindx),'}" it maps.'])
            end
            
            % Also check if a particular state variable in the component
            %  may not need to be differentiated.
            if is_x_Robj && any(Dval==diff_tab(Rindx,has_vars_Rcomp)) && size(Cval,2)>1
                for kk=1:size(Cval,2)
                    if all(isequal(Cval(:,kk),zeros(size(Cval,1),1)))
                        if kk==1
                            kkth = [num2str(kk),'st'];
                        elseif kk==2
                            kkth = [num2str(kk),'nd'];
                        elseif kk==3
                            kkth = [num2str(kk),'rd'];
                        else
                            kkth = [num2str(kk),'th'];
                        end
                        warning(['The ',kkth,' state variable in component "x{',num2str(Rindx),'}" appears not to contribute to the term "',term_name,'", but will still be differentiated up to degrees "',term_name,'.D".',...
                                    ' If this state variable does not need to be differentiable up to the same order as the other variables in state component component "x{',num2str(Rindx),'}", then please add it as a separate component to the PDE structure.']) 
                    end
                end
            end
            iszero_Cval = false;
        else
            % If not polynomial, only check if a particular state variable 
            % in the component may not need to be differentiated.
            Cval = double(Cval);
            if is_x_Robj && any(Dval==diff_tab(Rindx,has_vars_Rcomp)) && size(Cval,2)>1
                for kk=1:size(Cval,2)
                    if ~any(Cval(:,kk))
                        if kk==1
                            kkth = [num2str(kk),'st'];
                        elseif kk==2
                            kkth = [num2str(kk),'nd'];
                        elseif kk==3
                            kkth = [num2str(kk),'rd'];
                        else
                            kkth = [num2str(kk),'th'];
                        end
                        warning(['The ',kkth,' state variable in component "x{',num2str(Rindx),'}" appears not to contribute to the term "',term_name,'", but will still be differentiated up to degrees "',term_name,'.D".',...
                                    ' If this state variable does not need to be differentiable up to the same order as the other variables in state component component "x{',num2str(Rindx),'}", then please add it as a separate component to the PDE structure.']) 
                    end
                end
            end
            % Check if coefficients are all zero.
            iszero_Cval = all(all(Cval==0));
        end
        % If the coefficients are all zero, we can remove this term.
        if iszero_Cval
            PDE.(obj){ii}.term = [PDE.(obj){ii}.term(1:jj-1),PDE.(obj){ii}.term(jj+1:end)];
            n_terms = n_terms - 1;
            continue
        else
            PDE.(obj){ii}.term{jj}.C = Cval;
        end
        
        
        % % With that, the term should be good. We order the fields, and 
        % % move on to the next term.        
        if isfield(PDE.(obj){ii}.term{jj},'term_num')
            % The field term_num was only for internal use...
            PDE.(obj){ii}.term{jj} = rmfield(PDE.(obj){ii}.term{jj},'term_num');
        end
        if is_x_Robj
            PDE.(obj){ii}.term{jj} = orderfields(PDE.(obj){ii}.term{jj},{Robj,'D','loc','C','I','delay'});
        else
            PDE.(obj){ii}.term{jj} = orderfields(PDE.(obj){ii}.term{jj},{Robj,'C','I','delay'});
        end
        jj = jj+1; 
    end    
    ii = ii+1;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [PDE,BC_diff_tab,BC_state_indcs] = check_BC_terms(PDE,eq_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PDE,BC_state_tab_ii] = check_BC_terms(PDE,eq_num)
% checks the different terms in the BC specified by "eq_num", to check that 
% all mandatory fields have been properly specified, and to assign a 
% default value to all optional fields that the user has not specified. In
% addition, "PDE.BC_tab" is updated to indicate on which variables s the
% function f(s) defining the boundary condition 0=f(s) depends.
%
% INPUTS:
% - PDE:        A "struct" or "pde_struct" class object defining a PDE.
% - eq_num:     Integer specifying which equation "PDE.BC{eq_num}" to
%               check.
%
% OUTPUTS:
% - PDE:        The same structure as the input "PDE", now with all fields
%               in the structures "PDE.BC{eq_num}.term{jj}" for all jj
%               specified. In addition, the ith row of the array 
%               "PDE.BC_tab" is updated to indicate for each spatial
%               variable whether the function F_i(s) defining the BC
%               0=F_i(s) depends on this spatial variable.
% - BC_diff_tab:    A 1xp array indicating for the specified
%                   boundary condition to what order the RHS function
%                   F_i(s) defining the BC 0=F_i(s) is differentiable in
%                   each of the p global variables s1,...,sp.
% - BC_state_indcs: An integer array specifying which state components
%                   appear in the considered boundary condition.
%
% The main difference between this function and the other "check_terms"
% function is that, where we know the states and outputs to depend on
% particular spatial variables, and can therefore check that each term also
% depends on these same variables (and does not introduce other variables),
% fot the BCs, no variables are specified. Instead, we check for each of
% the specified terms what variables it depends on, and then say that the
% BC depends on the union of all of these variables.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the global variables and their associated domains.
global_vars = PDE.vars;
global_dom = PDE.dom;
nvars = size(global_vars,1);

% Keep track of which state components appear in the BC.
diff_tab = PDE.x_tab(:,3+nvars:2+2*nvars);
BC_diff_tab = zeros(1,nvars);
BC_state_indcs = [];

% Loop over the terms, noting that terms may be added if some element 
% PDE.BC{eq_num}.term{jj} actually describes multiple terms.
jj = 1;
n_terms = numel(PDE.BC{eq_num}.term);
while jj<=n_terms
    term_jj = PDE.BC{eq_num}.term{jj};
    
    % For the error messages, keep track of the name of the term.
    if isfield(PDE.BC{eq_num},'eq_num')
        eq_num_str = num2str(PDE.BC{eq_num}.eq_num);
    else
        eq_num_str = num2str(eq_num);
    end
    if isfield(PDE.BC{eq_num}.term{jj},'term_num')
        term_num = PDE.BC{eq_num}.term{jj}.term_num;
    else
        term_num = jj;
    end
    term_name = ['BC{',eq_num_str,'}.term{',num2str(term_num),'}'];
    
    
    % % % Assign a default value to all the non-specified fields.    
    % % First check which state component or input is involved.
    if isfield(term_jj,'x')
        Robj = 'x';
        is_x_Robj = true;
    elseif isfield(term_jj,'w')
        Robj = 'w';
        is_x_Robj = false;
    else
        Robj = 'u';
        is_x_Robj = false;
    end
    Rindx = term_jj.(Robj)(:);
    nR = length(Rindx);
    % Establish which spatial variables each component depends on.
    has_vars_Rcomp = logical(PDE.([Robj,'_tab'])(Rindx(1),3:2+nvars));
    nvars_Rcomp = sum(has_vars_Rcomp);
    Rvar_order = PDE.(Robj){Rindx}.var_order;
    % Note that we already checked all the different components listed
    % in Rindx to depend on the same variables in "get_diff".
        
    % % Extract a derivative (was already initialized in "get_diff").
    if is_x_Robj && (~isfield(term_jj,'D') || isempty(term_jj.D))
        % If no derivative is specified, assume no derivative is
        % desired.
        Dval = zeros(nR,nvars_Rcomp);
    elseif is_x_Robj
        Dval = term_jj.D;
    end

    % % Extract a spatial position (was already initialized in 
    % % "get_diff").
    if is_x_Robj
        Rloc = term_jj.loc;
    end
    
    % % Initialize delay.
    if ~isfield(term_jj,'delay')
        delay = 0;
        delay_varname = cell(0,1);
    else
        delay = term_jj.delay;
        if any(size(delay)~=[1,1]) || ~(isa(delay,'double') || isdouble(delay) || ispvar(delay))
            error(['The delay "',term_name,'.delay" is not properly specified;',...
                    ' please specify a scalar double or pvar value.'])
        end
        if isa(delay,'double') || isdouble(delay)
            delay = abs(double(delay));
            delay_varname = cell(0,1);     % The coefficients are allowed to depend on the delay
            if delay~=0
                PDE.has_delay = true;
            end
        elseif ispvar(delay) && length(delay.varname)>1
            error(['The delay "',term_name,'.delay" is not properly specified;',...
                    ' distributed delays should be specified using a single scalar variable.'])
        elseif ~ismember(delay.varname{1},PDE.tau.varname)
            error(['The delay "',term_name,'.delay" is not properly specified;',...
                    ' any distributed delay variable should be specified in the field "PDE.tau", along with the domain on which it exists.'])
        else
            delay_varname = delay.varname;
            PDE.has_delay = true;
        end
    end
    % We won't be doing anything with the delay here, we process delays
    % in a different function "expand_delays".
    PDE.BC{eq_num}.term{jj}.delay = delay;
    
    % % Initialize coefficients.
    if ~isfield(term_jj,'C') || isempty(term_jj.C)
        % If no coefficients are specified, assume each term is scaled
        % with a factor 1.
        if PDE.BC_tab(eq_num,2)==0
            % If we don't know the size of the BC, assume it is equal to
            % the size of the considered state component/input
            PDE.BC_tab(eq_num,2) = PDE.([Robj,'_tab'])(Rindx,2);
            PDE.BC{eq_num}.size = PDE.BC_tab(eq_num,2);
            Cval = repmat(eye(PDE.BC_tab(eq_num,2)),[1,nR]);
        elseif any(PDE.([Robj,'_tab'])(Rindx,2)~=PDE.BC_tab(eq_num,2))
            error(['The size of the component "BC{',eq_num_str,'}" does not match that of the components specified in "',term_name,'.',Robj,'";',...
                ' please specify a matrix "',term_name,'.C" describing how each component on the right-hand side contributes to the equation for "BC{',eq_num_str,'}".'])
        else
            Cval = repmat(eye(max(PDE.BC_tab(eq_num,2),1)),[1,nR]);
        end
    else
        % Otherwise, verify that the dimensions of the matrix makes sense.
        Cval = term_jj.C;        
        if PDE.BC_tab(eq_num,2)==0
            PDE.BC_tab(eq_num,2) = size(Cval,1);
            PDE.BC{eq_num}.size = PDE.BC_tab(eq_num,2);
        end
        if ~isa(Cval,'double') && ~isa(Cval,'polynomial')
            error(['The coefficients "',term_name,'.C" are not appropriately specified;',...
                ' coefficients should be specified as an array of type "double" or "polynomial".'])
        elseif size(Cval,1)~=PDE.BC_tab(eq_num,2)
            error(['The number of rows of the coefficient matrix "',term_name,'.C" does not match the number of BC equations "BC{',eq_num_str,'}.size" it contributes to.'])
        elseif size(Cval,2)~=sum(PDE.([Robj,'_tab'])(Rindx,2))
            error(['The number of columns of the coefficient matrix "',term_name,'.C" does not match the number of state variables it maps.'])
        end
        if size(Cval,2)~=sum(PDE.([Robj,'_tab'])(Rindx,2))
            error(['The number of columns of the coefficient matrix "',term_name,'.C" does not match the number of state variables it maps.'])
        end
    end
    
    % % Initialize an empty integral (for now) if no integral is
    % % specified.
    if isfield(term_jj,'int') && ~isfield(term_jj,'I')
        term_jj.I = term_jj.int;
        PDE.BC{eq_num}.term{jj} = rmfield(PDE.BC{eq_num}.term{jj},'int');
    end
    if isfield(term_jj,'I')
        Ival = term_jj.I(:);
        if ~isa(Ival,'cell')
            error(['The field "',term_name,'.I" is not appropriately specified;',...
                ' the field should be a cell with each element k providing the desired domain of integration in variable k of component "',Robj,'{',num2str(Rindx),'}" depends.'])
        end
    else
        Ival = cell(nvars_Rcomp,1);
    end
    
    % % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % % % If multiple terms are included, split into separate elements.
    C_colnums = cumsum([0;PDE.([Robj,'_tab'])(Rindx,2)]);
    if nR>1
        % % We squeeze the new terms between the current one and the
        % % following ones.
        PDE.BC{eq_num}.term = [PDE.BC{eq_num}.term(1:jj), repmat(PDE.BC{eq_num}.term(jj),[1,nR-1]), PDE.BC{eq_num}.term(jj+1:end)];
        if is_x_Robj
            % Make sure each term indeed describes just one term.
            for kk=1:nR-1
                C_cols = C_colnums(kk+1)+1 : C_colnums(kk+2);
                PDE.BC{eq_num}.term{jj+kk}.(Robj) = Rindx(kk+1);
                PDE.BC{eq_num}.term{jj+kk}.D = Dval(kk+1,:);
                PDE.BC{eq_num}.term{jj+kk}.loc = Rloc(kk+1,:);
                PDE.BC{eq_num}.term{jj+kk}.C = Cval(:,C_cols);
                PDE.BC{eq_num}.term{jj+kk}.I = Ival;
                PDE.BC{eq_num}.term{jj+kk}.delay = delay;
                PDE.BC{eq_num}.term{jj+kk}.term_num = term_num;
            end
        else
            % Make sure each term indeed describes just one term.
            for kk=1:nR-1
                C_cols = C_colnums(kk+1)+1 : C_colnums(kk+2);
                PDE.BC{eq_num}.term{jj+kk}.(Robj) = Rindx(kk+1);
                PDE.BC{eq_num}.term{jj+kk}.C = Cval(:,C_cols);
                PDE.BC{eq_num}.term{jj+kk}.I = Ival;
                PDE.BC{eq_num}.term{jj+kk}.delay = delay;
                PDE.BC{eq_num}.term{jj+kk}.term_num = term_num;
            end
        end
        n_terms = n_terms + nR-1;
    end
    % % Adjust the current element to describe only the first term.
    Rindx = Rindx(1);
    sz_R = PDE.([Robj,'_tab'])(Rindx,2);
    Cval = Cval(:,1:sz_R);
    if is_x_Robj
        Dval = Dval(1,:);
        Rloc = Rloc(1,:);
        PDE.BC{eq_num}.term{jj}.(Robj) = Rindx;
        PDE.BC{eq_num}.term{jj}.D = Dval;
        PDE.BC{eq_num}.term{jj}.loc = Rloc;
        PDE.BC{eq_num}.term{jj}.C = Cval;
        PDE.BC{eq_num}.term{jj}.I = Ival;
        PDE.BC{eq_num}.term{jj}.delay = delay;
    else
        PDE.BC{eq_num}.term{jj}.(Robj) = Rindx;
        PDE.BC{eq_num}.term{jj}.C = Cval;
        PDE.BC{eq_num}.term{jj}.I = Ival;
        PDE.BC{eq_num}.term{jj}.delay = delay;
    end
    
    
    % % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % % % Having initialized all the fields, and ensured that the
    % % % current term truly describes just one term, we now make sure
    % % % everything is properly specified, and fill in any remaining
    % % % gaps.
    
    % % First, make sure that the specified location is appropriate.
    % Extract the variables associated to the RHS component, and their
    % associated domains
    Rvars = global_vars(has_vars_Rcomp,:);
    Rdom = global_dom(has_vars_Rcomp,:);
    % Keep track of which variables the term actually varies in
    % (without integration).
    isvariable_term_jj = true(1,nvars_Rcomp);
    if is_x_Robj
        for kk=1:nvars_Rcomp
            % For each of the spatial variables, check if an appropriate
            % spatial position has been specified.
            Rloc_kk = Rloc(kk);
            if isa(Rloc_kk,'double') || isa(Rloc_kk,'polynomial') && isdouble(Rloc_kk)
                Rloc_kk = double(Rloc_kk);
                if ~ismember(Rloc_kk,Rdom(kk,:))
                    error(['The spatial position "',term_name,'.loc" is not appropriate;',...
                        ' for any variable s in [a,b], PIETOOLS can only evaluate s=a, s=b, or s=s.'])
                else
                    % Indicate that the term does not vary along this
                    % dimension.
                    isvariable_term_jj(kk) = false;
                end
            elseif ispvar(Rloc_kk)
                if ~isequal(Rloc_kk,Rvars(kk,1)) && ~isequal(Rloc_kk,Rvars(kk,2))
                    error(['The spatial position "',term_name,'.loc" is not appropriate;',...
                        ' for any variable s in [a,b], PIETOOLS can only evaluate s=a, s=b, or s=s.'])
                else
                    % Replace any dummy variables with primary
                    % variables.
                    Rloc_kk = Rvars(kk,1);
                    PDE.BC{eq_num}.term{jj}.loc(kk) = Rloc_kk;
                end
            else
                error(['The spatial position "',term_name,'.loc" is not appropriate;',...
                    ' spatial position must be specified as "double" or "pvar" (polynomial class) object.'])
            end
        end
    end
    
    % % Check the integrals, and establish which variables the boundary
    % % condition depends on.
    % Keep track of wrt which variables some integration is performed.
    is_int_var_Rcomp = false(1,nvars_Rcomp);
    if ~isfield(term_jj,'I') || isempty(term_jj.I)
        % If no integral is specified, assume no integration is desired.
        Ival = cell(nvars_Rcomp,1);
        PDE.BC{eq_num}.term{jj}.I = Ival;
    else
        % If integrals are specified, make sure that they are
        % appropriate.
        if numel(Ival)>nvars_Rcomp
            error(['The field "',term_name,'.I" is not appropriately specified;',...
                ' the number of elements should match the number of spatial variables on which component "',term_name,'.',Robj,'{',num2str(Rindx),'}" depends.'])
        end
        % If insufficient domains are specified, assume remaining
        % domains are empty.
        if numel(Ival)<nvars_Rcomp
            Ival = [Ival;cell(nvars_Rcomp-numel(Ival),1)];
        end
        Ival = Ival(Rvar_order);    % Reorder based on new ordering of variables.
        % Keep track of which along which directions a full integral is
        % performed.
        %is_full_int = false(1,numel(Ival));
        for kk = 1:numel(Ival)
            if isempty(Ival{kk})
                continue
            elseif ~isvariable_term_jj(kk)
                % Avoid integration in variables along which the
                % considered component does not depend.
                error(['The proposed integral "',term_name,'.I{',num2str(kk),'}" is not allowed;',...
                    ' the considered component "',term_name,'.',Robj,'{',num2str(Rindx),'}" at the specified position "',term_name,'.loc" does not vary along the proposed dimension of integration.'])
            end
            if (~isa(Ival{kk},'double') && ~isa(Ival{kk},'polynomial')) || ~all(size(Ival{kk})==[1,2])
                error(['The field "',term_name,'.I" is not appropriately specified;',...
                    ' each element k must be a 1x2 array of type "double" or "polynomial", indicating the domain of integration in variable k of component "',Robj,'{',num2str(Rindx),'}" depends.'])
            end
            Ival_kk = polynomial(Ival{kk});
            % Check that the integral is taken over one of the allowed
            % domains.
            if all(isequal(Ival_kk,Rdom(kk,:)))
                % The integral will discard any spatial dependence on the kkth
                % variable.
                isvariable_term_jj(kk) = false;
                %is_full_int(kk) = true;
            elseif all(isequal(Ival_kk,polynomial(Rdom(kk,2:-1:1))))
                % Integration is performed over mirror image of full domain
                % (for whatever reason).
                Cval = -Cval;
                PDE.BC{eq_num}.term{jj}.C = Cval;
                Ival{kk} = Ival{kk}(1,end:-1:1);
                isvariable_term_jj(kk) = false;
            elseif all(isequal(Ival_kk,[Rvars(kk,1),Rdom(kk,1)])) || all(isequal(Ival_kk,[Rdom(kk,2),Rvars(kk,1)]))
                % Integration is performed over mirror image of partial domain
                % (for whatever reason).
                Cval = -Cval;
                PDE.BC{eq_num}.term{jj}.C = Cval;
                Ival{kk} = Ival{kk}(1,end:-1:1);
            elseif ~(all(isequal(Ival_kk,[Rdom(kk,1),Rvars(kk,1)])) || all(isequal(Ival_kk,[Rvars(kk,1),Rdom(kk,2)])))
                error(['The proposed integral "',term_name,'.I{',num2str(kk),'}" is not allowed;',...
                    ' for a variable s in [a,b], integration is only allowed over [a,s], [s,b], or [a,b].'])
            end
            is_int_var_Rcomp(kk) = true;  % Indicate that we are integrating this variable.
        end
        PDE.BC{eq_num}.term{jj}.I = Ival;
    end
    % Logical indices indicating which of the global variables are
    % integrated.
    is_int_var_full = false(1,nvars);
    is_int_var_full(has_vars_Rcomp) = is_int_var_Rcomp;
    %is_full_int_full = false(1,nvars);
    %is_full_int_full(has_vars_Rcomp) = is_full_int;
    
    % Logical indices indicating which of the global variables the term
    % actually varies in (are not boundary positions or integrated out).
    isvariable_term_jj_full = false(1,nvars);
    isvariable_term_jj_full(has_vars_Rcomp) = isvariable_term_jj;
    
    % % Finally, we check that the coefficients are appropriately
    % % specified.
    if isa(Cval,'polynomial') && ~isdouble(Cval)
        % Make sure the polynomial does not introduce any new variables.
        if any(~ismember(Cval.varname,[global_vars.varname;delay_varname]))
            error(['The proposed polynomial "',term_name,'.C" is not appropriate;',...
                ' some of the (spatial) variables it depends on do not appear in any state component, input, or output.'])
        end
        
        % Check for each global variable and associated dummy variable if
        % it appears in Cval.
        has_var1_Cval = false(1,nvars);
        has_var2_Cval = false(1,nvars);
        for kk=1:nvars
            has_var1_Cval(kk) = ismember(global_vars(kk,1).varname,Cval.varname);
            has_var2_Cval(kk) = ismember(global_vars(kk,2).varname,Cval.varname);
        end
        
        % Check that no dummy variables are used along directions of no
        % integration.
        if any(has_var2_Cval & ~is_int_var_full)
            error(['The proposed polynomial "',term_name,'.C" is not appropriate;',...
                ' dummy variables are allowed only in directions along which the term is integrated.'])
        end
        % Check if any variable is (accidentally) introduced which the
        % integrals would otherwise remove.
        if any(has_var1_Cval & ~has_var2_Cval & is_int_var_full)
            warning(['The proposed polynomial "',term_name,'.C" introduces a dependence of the term on a variable which the integral removes;',...
                ' if you wish the polynomial to act as kernel, use a dummy variable, which can be specified in the second column of "',Robj,'{',num2str(Rindx),'}.vars".'])
        end
        
        % Any primary variable that the polynomial Cval depends on will be
        % one the term depends on.
        isvariable_term_jj_full = isvariable_term_jj_full | has_var1_Cval;
    end
    PDE.BC{eq_num}.term{jj}.C = Cval;
    
    % Add any variable dependence that the term introduces to the table.
    PDE.BC_tab(eq_num,3:2+nvars) = PDE.BC_tab(eq_num,3:2+nvars) | isvariable_term_jj_full;
    
    % Also keep track of which state variables appear in the BC, and what
    % the maximal order of differentiability is of the states.
    if is_x_Robj
        BC_state_indcs = [BC_state_indcs;Rindx];
        BC_diff_tab = max(BC_diff_tab,diff_tab(Rindx,:));
    end
    
    % % With that, the term should be good. We order the fields, and
    % % move on to the next term.
    if isfield(PDE.BC{eq_num}.term{jj},'term_num')
        % The field term_num was only for internal use...
        PDE.BC{eq_num}.term{jj} = rmfield(PDE.BC{eq_num}.term{jj},'term_num');
    end
    if is_x_Robj
        PDE.BC{eq_num}.term{jj} = orderfields(PDE.BC{eq_num}.term{jj},{Robj,'D','loc','C','I','delay'});
    else
        PDE.BC{eq_num}.term{jj} = orderfields(PDE.BC{eq_num}.term{jj},{Robj,'C','I','delay'});
    end
    jj = jj+1;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_initialization_summary(PDE,obj,ncomps_x)
% print_initialization_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning how many
% components "PDE.(obj)" were encountered, on what variables each of these
% components depends, and (iff obj=='x') to what order the component is
% differentiable in each of the spatial variables on which it depends.
%
% INPUTS:
% - PDE:        A "struct" or "pde_struct" class object defining a PDE.
% - obj:        char 'x', 'y', 'z', 'u', 'w', or 'BC', indicating for which
%               field the information is desired.
% - ncomps_x:   Integer scalar specifying how many state components were
%               present in the PDE before initialization. Needs only be
%               specified if obj=='x'.
%
% OUTPUTS:
% Displays information in the command window concerning the number of
% observed components of type obj, along with the size of each of these
% components, what variables they depend on, and (if obj=='x') to what
% order they are differentiable in each of these variables.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set a name associated to each object.
if strcmp(obj,'x')
    object_name = 'state component';
elseif strcmp(obj,'y')
    object_name = 'observed output';
elseif strcmp(obj,'z')
    object_name = 'regulated output';
elseif strcmp(obj,'u')
    object_name = 'actuator input';
elseif strcmp(obj,'w')
    object_name = 'exogenous input';
elseif strcmp(obj,'BC')
    object_name = 'boundary condition';
end
% If there are multiple components of the object, we add an 's' at the end.
num_comps = numel(PDE.(obj));
if num_comps>1
    add_s = 's';
else
    add_s = '';
end

% Use UNICODE to add subscript indices to different components.
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
% Determine how many subscripts will be needed.
n_digits = 1+floor(log(num_comps)/log(10));


% Determine which variables appear in the different components of this
% object.
nvars = size(PDE.vars,1);
global_vars_obj = PDE.vars(any(PDE.([obj,'_tab'])(:,3:2+nvars),1),:);
if isempty(global_vars_obj)
    global_varnames = '(t)';
else
    global_varnames =  global_vars_obj(1,1).varname{1};
    for kk=2:size(PDE.vars,1)
        global_varnames = [global_varnames,',',PDE.vars(kk,1).varname{1}];
    end
    global_varnames = ['(t,',global_varnames,')'];
end

% Estimate the size of the string of characters denoting the components
% "PDE.(obj){ii}".
LHS_length_max = 1+length(obj)+n_digits+length(global_varnames) + 3; % e.g. ' x13(t,s1,s2,s3)', 
if num_comps==1
    LHS_length_max = LHS_length_max - 1;    % No subscript index.
end

if num_comps==0
    % If there are no components for this object, display nothing.
    %fprintf(['Encountered no ',object_name,add_s,'.\n']);  % Should we indicate that no object of this type was encountered?
    return
else
    % Indicate how many components of the considered object were
    % encountered.
    if strcmp(obj,'x')
        fprintf(['\n','Encountered ',num2str(ncomps_x),' ',object_name,add_s,': \n']);
    else
        fprintf(['\n','Encountered ',num2str(num_comps),' ',object_name,add_s,': \n']);
    end
end
% For each of the components, display its size, and which variables it
% depends on.
for ii=1:num_comps
    comp_ii = PDE.(obj){ii};
    
    % Establish the names of the variables on which the component depends.
    if isempty(comp_ii.vars)
        varnames_ii = '';
        varnames_ii_t = 't';    % All components may vary in time
    else
        % Set a comma separated list.
        varnames_ii = comp_ii.vars(1,1).varname{1};
        for kk=2:size(comp_ii.vars,1)
            varnames_ii = [varnames_ii,',',comp_ii.vars(kk,1).varname{1}];
        end
        varnames_ii_t = ['t,',varnames_ii]; % All components may vary in time
    end
    
    % Establish the (subscript) index for the component.
    LHS_length = length(varnames_ii_t);
    if numel(PDE.(obj))==1
        % There is only one state component --> no need to give index
        Lcomp_idx = '';
    elseif numel(PDE.(obj))<=9
        % The state number consists of a single decimal.
        Lcomp_idx = sub_num{ii+1};
        LHS_length = LHS_length + 1;
    else
        % The state number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(ii)')+1)');
        LHS_length = LHS_length + length(num2str(ii));
    end
    % Set the name of the component, including its depdence on spatial
    % variables.
    if strcmp(obj,'BC')
        LHS_name = [' F',Lcomp_idx,'(',varnames_ii_t,') = 0'];
        LHS_length = length(obj) + LHS_length + 7;
    else
        LHS_name = [' ',obj,Lcomp_idx,'(',varnames_ii_t,')'];
        LHS_length = length(obj) + LHS_length + 3;
    end
    
    % If the state component is one that has been added to impose the
    % higher order temporal derivative, indicate to which component the
    % added state corresponds.
    if strcmp(obj,'x') && ii>ncomps_x
        if (numel(PDE.(obj))-ncomps_x)~=1
            add_s2 = 's';
        else
            add_s2 = '';
        end
        % Indicate that the remaining state components were all added
        % during initialization.
        if ii==ncomps_x+1
            fprintf(['\n','Added ',num2str(numel(PDE.(obj))-ncomps_x),' ',object_name,add_s2,': \n']);
        end
        % For the added state components, indicate of which state component
        % they are the temporal derivative.
        Rcomp_idx = cell2mat(sub_num(str2num(num2str(PDE.x{ii}.term{1}.x)')+1)');
        xdot = '\x1E8B';    % UNICODE \dot{x}
        RHS = [' := ',xdot,Rcomp_idx,'(',varnames_ii_t,')'];
    else
        RHS = [',',repmat(' ',[1,LHS_length_max-LHS_length])];
    end
    
    % Determine the order of differentiability in each spatial variable for
    % the state component
    if strcmp(obj,'x')
        D = comp_ii.diff;
        %has_vars_comp = logical(PDE.([obj,'_tab'])(ii,2+1:2+nvars));
        %D = D(has_vars_comp);
        if isempty(D)
            % If the component does not vary in space, indicate that it is
            % finite dimensional.
            fprintf([LHS_name,RHS,' of size ',num2str(comp_ii.size),', finite-dimensional;\n'])
        else
            % Otherwise, indicate the order of differentiability for the state
            diff_order_ii = num2str(D(1));
            for kk=2:length(D)
                diff_order_ii = [diff_order_ii,',',num2str(D(kk))];
            end

            fprintf([LHS_name,RHS,' of size ',num2str(comp_ii.size),', differentiable up to order (',diff_order_ii,') in variables (',varnames_ii,');\n'])
        end
    else
        % If the object is not a state component, there is no order of
        % differentiability to specify.
        fprintf([LHS_name,RHS,' of size ',num2str(comp_ii.size),';\n'])
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
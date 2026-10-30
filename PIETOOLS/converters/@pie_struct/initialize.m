function PIE = initialize(PIE)
% PIE = INITIALIZE(PIE) takes a "pie_struct" object and checks that all the
% necessary fields are appropriately specified, and assigns a default value
% to all the optional fields that have not been specified, so that the PIE
% is well-defined and consistent.
%
% INPUT
%   PIE:    "pie_struct" class object, with operators specified as 'opvar',
%           'opvar2d', or 'copvar' objects.
%
% OUTPUT
%   PIE:    "pie_struct" class object describing the same system as the
%           input, with the operators T, Tw, Tu, A, Bw, Bu, Cz, Dzw, Dzu,
%           Cy, Dyw, Dyu, K, Kw, and Ku of type 'copvar', with the function
%           spaces PIE.space specified for each of the components, and with
%           the field PIE.is_initialized set to true.
%
% NOTES
% - The operators T and A are mandatory. Any other of the operators that
%   has not been specified is initialized as a zero operator between the
%   spaces specified by the other operators.
% - Operators specified as 'opvar' or 'opvar2d' objects are converted to
%   'copvar' objects, which do not involve dummy variables. Dummy variables
%   of the operators are first renamed to follow the convention of the PIE,
%   adding '_dum' to the name of the spatial variable, see also
%   PIE.space.vars.
% - The variables and domain of each operator must match those of the PIE,
%   and the output (input) spaces of all operators mapping to (from) the
%   same type of component (x, w, u, z, or y) must match.
% - If PIE.space.(type) has not been specified for the state (x), inputs
%   (w and u), or outputs (z and y), the components are determined from the
%   spaces of the operators, with one component for each space. If it has
%   been specified, the components are checked against the operators.
% - The operators K, Kw, and Ku defining the algebraic constraint must all
%   map to the same space. If none of them has been specified, they are
%   initialized as operators without any rows (no output space), from the
%   spaces of the state x, the inputs w, and the actuators u, respectively.
%   If one of them has rows, any of the others without rows is treated as
%   unspecified, and initialized as a zero operator.
% - The fields of PIE.pde2pie must have a row for each of the components in
%   the associated function space, if they are specified.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - initialize
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding SS - 08/01/2022
% Adjust and add as class file function, DJ - 11/01/2022
% DJ, 12/07/2024: Use new default vars (s1,s2) and (s1_dum,s2_dum);
% DJ, 10/29/2026: Update to new pie_struct structure: operators converted
%                   to copvar, spaces and maps between PDE and PIE checked;

% Names of the operators, with the output (row) space of the operator
% specified by the row type, and the input (column) space by the column.
opnames = {'T', 'Tw', 'Tu';                                                 % DJ, 10/29/2026
           'A', 'Bw', 'Bu';
           'Cz','Dzw','Dzu';
           'Cy','Dyw','Dyu';
           'K', 'Kw', 'Ku'};
row_type = {'x';'x';'z';'y';'k'};
col_type = {'x','w','u'};
is_unset = @(P) isnumeric(P) && isempty(P);
% If one of the operators for the constraint has rows, any of the others
% without rows is treated as unspecified.
has_rows = cellfun(@(n) size(PIE.(n),1)>0,opnames(5,:));
for n = opnames(5,~has_rows)
    if any(has_rows)
        PIE.(n{1}) = [];
    end
end

% % % Establish the spatial variables and domain of the PIE.
% The state operator T is mandatory, and determines the variables.
T = PIE.T;
if is_unset(T)
    error("PIE.T is a mandatory element of PIE structure and cannot be defaulted to zero");
elseif isa(T,'opvar') || isa(T,'opvar2d')
    % Adapt variables and domain of the PIE to those of the T operator.
    PIE.space.vars = arrayfun(@(ii) T.var1(ii).varname{1},1:size(T.var1,1),'UniformOutput',false);
    PIE.space.dom = T.I;
elseif isa(T,'copvar')
    if isempty(PIE.space.vars)
        PIE.space.vars = T.vars;
        PIE.space.dom = T.dom;
    elseif isempty(PIE.space.dom)
        % Variables specified, but no domain
        % --> use domain of T, or [0,1] for variables not in T.
        dom = repmat([0,1],numel(PIE.space.vars),1);
        [is_T_var,loc] = ismember(PIE.space.vars,T.vars);
        dom(is_T_var,:) = T.dom(loc(is_T_var),:);
        PIE.space.dom = dom;
    end
else
    error("PIE.T should be specified as an 'opvar', 'opvar2d', or 'copvar' object.")
end
vnames = PIE.space.vars;
dom = PIE.space.dom;
nvars = numel(vnames);
if size(dom,1)~=nvars
    error("The number of domains specified for the PIE does not match the number of spatial variables.")
end


% % % Check that the operators are appropriately specified, and convert
% % % them to copvar objects.
for k = 1:numel(opnames)
    fname = opnames{k};
    Pop = PIE.(fname);
    if is_unset(Pop)
        % The operator has not been specified.
        if any(strcmp(fname,{'T','A'}))
            error(['PIE.',fname,' is a mandatory element of PIE structure and cannot be defaulted to zero']);
        end
        continue
    elseif isa(Pop,'opvar') || isa(Pop,'opvar2d')
        % Convert to copvar, using the dummy variables of the PIE. An
        % operator without any L2 space does not involve any variables.
        if any(any(Pop.dim(2:end,:)))
            op_vars = reshape(Pop.var1.varname,1,[]);
            op_dom = double(Pop.I);
        else
            op_vars = cell(1,0);
            op_dom = zeros(0,2);
        end
        Pop = op2copvar_sop(pie_rename_dummy(Pop));
    elseif isa(Pop,'copvar')
        op_vars = Pop.vars;
        op_dom = Pop.dom;
    else
        error(['PIE.',fname," should be specified as an 'opvar', 'opvar2d', or 'copvar' object."])
    end
    % The variables and domain of the operator must match those of the PIE.
    [is_pie_var,loc] = ismember(op_vars,vnames);
    if ~all(is_pie_var)
        error(['The spatial variables of the operator ',fname,' do not match those of the PIE.'])
    elseif any(any(op_dom~=dom(loc,:)))
        error(['The spatial domain of the operator ',fname,' does not match the spatial domain of the PIE.']);
    end
    PIE.(fname) = Pop;
end


% % % Establish the function spaces of the state, inputs and outputs from
% % % the operators, checking that the operators are consistent.
% For each type of component, collect the output (row) and input (column)
% spaces of all the operators involving this component, in the form of a
% depmat indicating on which variables each space depends, and the number of
% components in the space.
spaces = struct('x',[],'w',[],'u',[],'z',[],'y',[],'k',[]);
for ii=1:numel(opnames)
    fname = opnames{ii};
    if is_unset(PIE.(fname))
        continue
    end
    [r,c] = ind2sub(size(opnames),ii);
    for side = {'out','in'}
        if strcmp(side{1},'out')
            obj = row_type{r};
            depmat_op = PIE.(fname).space_out;
            dims = PIE.(fname).dim_out(:);
        else
            obj = col_type{c};
            depmat_op = PIE.(fname).space_in;
            dims = PIE.(fname).dim_in(:);
        end
        % The registry of the copvar is sorted, so permute to the order of
        % the variables in the PIE.
        [~,loc] = ismember(PIE.(fname).vars,vnames);
        depmat = false(size(depmat_op,1),nvars);
        depmat(:,loc) = logical(depmat_op);
        if isempty(spaces.(obj))
            spaces.(obj) = struct('depmat',depmat,'dims',dims,'name',[fname,' (',side{1},'put)']);
        elseif ~isequal(pie_space_canon(depmat,dims),pie_space_canon(spaces.(obj).depmat,spaces.(obj).dims))
            error(['The spaces of the operator ',fname,' do not match those of the operator ',spaces.(obj).name,', as suggested for the ',obj,' components.'])
        end
    end
end
% The operator T should be a map between symmetric spaces. We know its
% output spaces match those of x, so check its input spaces.
[~,loc] = ismember(PIE.T.vars,vnames);
depmat = false(size(PIE.T.space_in,1),nvars);
depmat(:,loc) = logical(PIE.T.space_in);
if ~isequal(pie_space_canon(depmat,PIE.T.dim_in),pie_space_canon(spaces.x.depmat,spaces.x.dims))
    error("PIE.T operator must be a map between symmetric spaces");
end
% Components for which no operators were specified have no spaces.
for obj = fieldnames(spaces)'
    if isempty(spaces.(obj{1}))
        spaces.(obj{1}) = struct('depmat',false(0,nvars),'dims',zeros(0,1),'name','');
    end
end


% % % Initialize all remaining operators as zero operators between the
% % % appropriate spaces. Without any constraint, the operators K, Kw, Ku
% % % have no rows.
for ii=1:numel(opnames)
    fname = opnames{ii};
    if ~is_unset(PIE.(fname))
        continue
    end
    [r,c] = ind2sub(size(opnames),ii);
    depmat_out = spaces.(row_type{r}).depmat;
    depmat_in = spaces.(col_type{c}).depmat;
    % The registry of the copvar consists of the variables on which at
    % least one of the spaces depends, sorted alphabetically.
    used = any(depmat_out,1) | any(depmat_in,1);
    [reg,idx] = sort(reshape(vnames(used),1,[]));
    dom_used = dom(used,:);
    depmat_out = depmat_out(:,used);    depmat_out = depmat_out(:,idx);
    depmat_in = depmat_in(:,used);      depmat_in = depmat_in(:,idx);
    meta = struct('vars',{reg},'dom',dom_used(idx,:),...
                  'space_out',depmat_out,'space_in',depmat_in,...
                  'dim_out',spaces.(row_type{r}).dims,'dim_in',spaces.(col_type{c}).dims);
    PIE.(fname) = copvar(cell(size(depmat_out,1),size(depmat_in,1)),meta);
end


% % % Set the function spaces of the PIE, and check that those that have
% % % already been specified match the spaces of the operators.
for obj_type = {'x','w','u','z','y'}
    obj = obj_type{1};
    depmat = spaces.(obj).depmat;
    dims = spaces.(obj).dims;
    sp = PIE.space.(obj);
    if isempty(sp.dims)
        % Use one component for each space of the operators.
        PIE.space.(obj).dims = dims(dims>0);
        PIE.space.(obj).depmat = depmat(dims>0,:);
    elseif size(sp.depmat,2)~=nvars || size(sp.depmat,1)~=numel(sp.dims)
        error(['The size of PIE.space.',obj,'.depmat does not match the number of components and the number of spatial variables.'])
    elseif any(sp.dims<=0) || any(sp.dims~=round(sp.dims))
        error(['The sizes PIE.space.',obj,'.dims of the components should be positive integers.'])
    elseif ~isequal(pie_space_canon(sp.depmat,sp.dims),pie_space_canon(depmat,dims))
        error(['The components in PIE.space.',obj,' do not match the spaces of the operators.'])
    end
    PIE.space.(obj).depmat = logical(PIE.space.(obj).depmat);
end


% % % Check that the maps between the PDE and PIE match the PIE components,
% % % and refer to components in the PDE.
for obj_type = {'x','w','u','z','y'}
    obj = obj_type{1};
    ncomps = numel(PIE.space.(obj).dims);
    for field = fieldnames(PIE.pde2pie.(obj))'
        if ~isempty(PIE.pde2pie.(obj).(field{1})) && size(PIE.pde2pie.(obj).(field{1}),1)~=ncomps
            error(['The number of rows of PIE.pde2pie.',obj,'.',field{1},' does not match the number of components in PIE.space.',obj,'.'])
        end
    end
    % For the state, the columns of idx refer to PDE states, inputs, and
    % actuators, respectively.
    ref_types = {obj};
    if strcmp(obj,'x')
        ref_types = {'x','w','u'};
    end
    for jj=1:numel(ref_types)
        sp_pde = PIE.space_pde.(ref_types{jj});
        if size(sp_pde.depmat,1)~=numel(sp_pde.dims) || (~isempty(sp_pde.ID) && numel(sp_pde.ID)~=numel(sp_pde.dims))
            error(['The fields of PIE.space_pde.',ref_types{jj},' should have a row for each PDE component.'])
        elseif ~isempty(sp_pde.dims) && ~isempty(PIE.pde2pie.(obj).idx) && any(PIE.pde2pie.(obj).idx(:,jj)>numel(sp_pde.dims))
            error(['PIE.pde2pie.',obj,'.idx refers to components not in PIE.space_pde.',ref_types{jj},'.'])
        end
    end
end

PIE.is_initialized = true;

end

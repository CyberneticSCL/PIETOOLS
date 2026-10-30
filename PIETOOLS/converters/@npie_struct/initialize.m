function PIE = initialize(PIE)
% PIE = initialize(PIE) takes an "npie_struct" object and checks that all
% the necessary fields are appropriately specified, and assigns a default
% value to all the optional fields that have not been specified, so that
% the PIE is well-defined and consistent.
%
% INPUT
%   PIE:    "npie_struct" class object, with operators T, Tw, Tu, Bw, Bu,
%           Dzw, Dzu, Dyw, and Dyu specified as 'tensopmat' objects, and
%           the polynomials f, gz, and gy specified as 'polyopvar' objects.
%
% OUTPUT
%   PIE:    "npie_struct" class object describing the same system as the
%           input, with the function spaces PIE.space specified for each of
%           the components, any operators that were not specified
%           initialized as zero operators, and the field
%           PIE.is_initialized set to true.
%
% NOTES
% - The operators T and f are mandatory. Any other of the operators that
%   has not been specified is initialized as a zero operator between the
%   spaces specified by the other operators.
% - Operators other than the polynomials that are not of type 'tensopmat',
%   e.g. 'opvar' or 'sopvar' objects, are converted to 'tensopmat' objects.
%   The polynomials f, gz, and gy must be specified as 'polyopvar' objects,
%   and are never converted.
% - The output function spaces of all operators mapping to the same type of
%   component (x, z, y, or the constraint) must match, so that e.g. T, Tw,
%   Tu, f, Bw, and Bu all map to the same space. Similarly, the input
%   spaces of all operators mapping from the same type of component (w or
%   u) must match, as must the input spaces of T, Bw, Bu, K, .... The state
%   variables of the polynomials f, gz, and gy must be among the components
%   named in PIE.space.x.varname, with the same sizes and dependency on the
%   spatial variables, but a polynomial need not depend on all components.
% - If gz or gy is not specified, it is initialized as a zero polynomial.
% - The variables and domain of each operator must match those of the PIE.
% - If PIE.space.(type) has not been specified for the state (x), inputs
%   (w and u), or outputs (z and y), the components are determined from the
%   spaces of the operators, with one component for each space. If it has
%   been specified, the components are checked against the operators.
% - The components in PIE.space and PIE.space_pde that have not been named
%   (varname empty) are assigned default names. The names must be unique,
%   and sorted for the state in PIE.space (see 'npie_struct').
% - The operators K, Kw, and Ku defining the algebraic constraint must all
%   map to the same space. If none of them has been specified, they are
%   initialized as operators without any rows (no output space), from the
%   spaces of the state x, the inputs w, and the actuators u, respectively.
%   If one of them has rows, any of the others without rows is treated as
%   unspecified, and initialized as a zero operator.
% - If specified, the 'tensopmat' objects in PIE.pie2pde must map from the
%   components of the PIE, and PIE.pde2pie.(type).C to the components of
%   the PIE, and the fields of PIE.pde2pie must have a row for each of the
%   components in the associated function space.

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
% DJ, 10/29/2026: Initial coding, based on the initialize function of
%                   pie_struct;

% Names of the operators, with the output (row) space of the operator
% specified by the row type, and the input (column) space by the column.
% The last row involves the operators for the algebraic constraint on the
% state, mapping to a space "k" that is not a component of the PIE.
opnames = {'T', 'Tw', 'Tu';
           'f', 'Bw', 'Bu';
           'gz','Dzw','Dzu';
           'gy','Dyw','Dyu';
           'K', 'Kw', 'Ku'};
row_type = {'x';'x';'z';'y';'k'};
col_type = {'x','w','u'};
poly_names = {'f','gz','gy'};
is_unset = @(P) isnumeric(P) && isempty(P);
% If one of the operators for the constraint has rows, any of the others
% without rows is treated as unspecified.
has_rows = cellfun(@(n) size(PIE.(n),1)>0,opnames(5,:));
for n = opnames(5,~has_rows)
    if any(has_rows)
        PIE.(n{1}) = [];
    end
end

% % % Check that the operators are specified, and convert them to tensopmat
% % % objects. The polynomials must be polyopvar objects, as the monomials
% % % in the state must be known to e.g. compute derivatives along the PIE.
for k=1:numel(opnames)
    fname = opnames{k};
    Pop = PIE.(fname);
    if is_unset(Pop)
        if any(strcmp(fname,{'T','f'}))
            error(['PIE.',fname,' is a mandatory element of PIE structure and cannot be defaulted to zero']);
        end
    elseif any(strcmp(fname,poly_names))
        if ~isa(Pop,'polyopvar')
            error(['PIE.',fname,' should be specified as a ''polyopvar'' object (was: ',class(Pop),').'])
        end
    elseif ~isa(Pop,'tensopmat')
        % Convert e.g. opvar and sopvar objects to tensopmat
        PIE.(fname) = tensopmat(Pop);
    end
end


% % % Establish the spatial variables and domain of the PIE.
% The operator T is mandatory, and determines the variables.
T = PIE.T;
if isempty(PIE.space.vars)
    [vnames,idx] = unique([T.vars.in,T.vars.out]);
    dom = [T.dom.in; T.dom.out];
    PIE.space.vars = reshape(vnames,1,[]);
    PIE.space.dom = dom(idx,:);
elseif isempty(PIE.space.dom)
    % Variables specified, but no domain
    % --> use domain of T, or [0,1] for variables not in T.
    [is_T_var,loc] = ismember(PIE.space.vars,[T.vars.in,T.vars.out]);
    dom = [T.dom.in; T.dom.out];
    PIE.space.dom = repmat([0,1],numel(PIE.space.vars),1);
    PIE.space.dom(is_T_var,:) = dom(loc(is_T_var),:);
end
vnames = PIE.space.vars;
dom = PIE.space.dom;
nvars = numel(vnames);
if size(dom,1)~=nvars
    error("The number of domains specified for the PIE does not match the number of spatial variables.")
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
    Pop = PIE.(fname);
    if is_unset(Pop)
        continue
    end
    [r,c] = ind2sub(size(opnames),ii);
    for side = {'in','out'}
        side = side{1};
        % The variables and domain must match those of the PIE.
        [is_pie_var,loc] = ismember(Pop.vars.(side),vnames);
        if ~all(is_pie_var)
            error(['The spatial variables of the operator ',fname,' do not match those of the PIE.'])
        elseif any(any(Pop.dom.(side)~=dom(loc,:)))
            error(['The spatial domain of the operator ',fname,' does not match the spatial domain of the PIE.']);
        end
        if strcmp(side,'out')
            type = row_type{r};
        else
            if any(strcmp(fname,poly_names))
                continue    % states of the polynomials are checked by name below
            end
            type = col_type{c};
        end
        % Depmat indicating on which variables each space depends, in the
        % order of the variables in the PIE.
        depmat = false(numel(Pop.dims.(side)),nvars);
        depmat(:,loc) = logical(Pop.depmat.(side));
        dims = Pop.dims.(side)(:);
        if isempty(spaces.(type))
            spaces.(type) = struct('depmat',depmat,'dims',dims,'name',[fname,' (',side,'put)']);
        elseif ~isequal(pie_space_canon(depmat,dims),pie_space_canon(spaces.(type).depmat,spaces.(type).dims))
            error(['The spaces of the operator ',fname,' do not match those of the operator ',spaces.(type).name,', as suggested for the ',type,' components.'])
        end
    end
end
% Components for which no operators were specified have no spaces.
for type = fieldnames(spaces)'
    if isempty(spaces.(type{1}))
        spaces.(type{1}) = struct('depmat',false(0,nvars),'dims',zeros(0,1),'name','');
    end
end


% % % Initialize all remaining operators as zero operators between the
% % % appropriate spaces. Without any constraint, the operators K, Kw, and Ku
% % % have no rows. Zero polynomials are initialized below, as they require
% % % the names of the state components.
for ii=1:numel(opnames)
    fname = opnames{ii};
    if ~is_unset(PIE.(fname)) || any(strcmp(fname,poly_names))
        continue
    end
    [r,c] = ind2sub(size(opnames),ii);
    sp.out = spaces.(row_type{r});
    sp.in = spaces.(col_type{c});
    dims = struct();    vars = struct();    dom_op = struct();    depmat_op = struct();
    for side = {'in','out'}
        side = side{1};
        % The variables are those on which at least one of the spaces
        % depends, sorted alphabetically.
        used = any(sp.(side).depmat,1);
        [vars.(side),idx] = sort(reshape(vnames(used),1,[]));
        dom_used = dom(used,:);
        dom_op.(side) = dom_used(idx,:);
        depmat_used = sp.(side).depmat(:,used);
        depmat_op.(side) = double(depmat_used(:,idx));
        dims.(side) = sp.(side).dims;
    end
    PIE.(fname) = mat2tensopmat(sparse(sum(dims.out),sum(dims.in)),dims,vars,dom_op,depmat_op);
end


% % % Set the function spaces of the PIE, and check that those that have
% % % already been specified match the spaces of the operators.
for obj_type = {'x','w','u','z','y'}
    obj = obj_type{1};
    depmat = spaces.(obj).depmat;
    dims = spaces.(obj).dims;
    sp_pie = PIE.space.(obj);
    if isempty(sp_pie.dims)
        % Use one component for each space of the operators.
        PIE.space.(obj).dims = dims(dims>0);
        PIE.space.(obj).depmat = depmat(dims>0,:);
    elseif size(sp_pie.depmat,2)~=nvars || size(sp_pie.depmat,1)~=numel(sp_pie.dims)
        error(['The size of PIE.space.',obj,'.depmat does not match the number of components and the number of spatial variables.'])
    elseif any(sp_pie.dims<=0) || any(sp_pie.dims~=round(sp_pie.dims))
        error(['The sizes PIE.space.',obj,'.dims of the components should be positive integers.'])
    elseif ~isequal(pie_space_canon(sp_pie.depmat,sp_pie.dims),pie_space_canon(depmat,dims))
        error(['The components in PIE.space.',obj,' do not match the spaces of the operators.'])
    end
    PIE.space.(obj).depmat = logical(PIE.space.(obj).depmat);
end


% % % Set the names of the components of the PIE and PDE, assigning default
% % % names if they have not been specified.
for sp_name = {'space','space_pde'}
    for obj_type = {'x','w','u','z','y'}
        obj = obj_type{1};
        sp = PIE.(sp_name{1}).(obj);
        n = numel(sp.dims);
        if ~isfield(sp,'varname'),  sp.varname = cell(0,1);     end
        names = sp.varname(:);
        if isempty(names) && n>0
            % Zero-padded numbers, so that the names are sorted in the
            % order of the components.
            fmt = ['%0',num2str(numel(num2str(n))),'d'];
            names = cell(n,1);
            for ii=1:n
                if strcmp(sp_name{1},'space')
                    names{ii} = [obj,num2str(ii,fmt),'_pie'];
                else
                    names{ii} = [obj,num2str(ii)];
                end
            end
        elseif ~iscellstr(names) || numel(names)~=n
            error(['PIE.',sp_name{1},'.',obj,'.varname should be a cellstr with a name for each of the components.'])
        elseif numel(unique(names))~=n
            error(['The names in PIE.',sp_name{1},'.',obj,'.varname should be unique.'])
        elseif strcmp(sp_name{1},'space') && strcmp(obj,'x') && ~issorted(names)
            error('The names in PIE.space.x.varname should be sorted alphabetically, as the components of PIE.x are sorted by name.')
        end
        PIE.(sp_name{1}).(obj).varname = names;
    end
end

% % % Check the polynomials f, gz, and gy in terms of the state variables, and
% % % initialize any that were not specified as zero polynomials.
xsp = PIE.space.x;
[vsort,vidx] = sort(vnames);
for k=1:numel(poly_names)
    fname = poly_names{k};
    Pop = PIE.(fname);
    ii = find(strcmp(opnames(:,1),fname));
    sp_out = PIE.space.(row_type{ii});
    if is_unset(Pop)
        % Linear zero polynomial in all state variables, with variables and
        % domain sorted as for any polyopvar.
        dims = struct('in',xsp.dims,'out',sp_out.dims);
        vars_op = struct('in',{vsort},'out',{vsort});
        dom_op = struct('in',dom(vidx,:),'out',dom(vidx,:));
        depmat_op = struct('in',xsp.depmat(:,vidx),'out',double(sp_out.depmat(:,vidx)));
        Pop = polyopvar();
        Pop.C = mat2tensopmat(sparse(sum(dims.out),sum(dims.in)),dims,vars_op,dom_op,depmat_op);
        Pop.Z = (1:numel(xsp.dims))';
        Pop.varname = xsp.varname;      Pop.dims = dims;
        Pop.vars = vars_op;             Pop.dom = dom_op;
        Pop.depmat = depmat_op;         Pop.chkval = true;
        PIE.(fname) = Pop;
        continue
    end
    if ~Pop.chkval
        Pop = initialize(Pop);
        PIE.(fname) = Pop;
    end
    [is_x,loc] = ismember(Pop.varname,xsp.varname);
    [is_var,vloc] = ismember(Pop.vars.in,vnames);
    if ~all(is_x)
        error(['The state variables of PIE.',fname,' must be components named in PIE.space.x.varname.'])
    elseif ~all(is_var)
        error(['The spatial variables of PIE.',fname,' do not match those of the PIE.'])
    end
    depmat_in = false(numel(loc),nvars);
    depmat_in(:,vloc) = logical(Pop.depmat.in);
    if ~isequal(Pop.dims.in(:),xsp.dims(loc)) || ~isequal(depmat_in,xsp.depmat(loc,:))
        error(['The state variables of PIE.',fname,' do not match the sizes and dependency on the spatial variables of the components in PIE.space.x.'])
    end
end


% % % Check that the maps between the PDE and PIE match the PIE components,
% % % and refer to components in the PDE.
for obj_type = {'x','w','u','z','y'}
    obj = obj_type{1};
    dims_pie = PIE.space.(obj).dims(:);
    dims_pde = PIE.space_pde.(obj).dims(:);
    % The map from PIE to PDE has the PIE components as input, and the PDE
    % components as output, if specified.
    X_pde = PIE.pie2pde.(obj);
    if ~isempty(X_pde.ops) && ~isequal(X_pde.dims.in(:),dims_pie)
        error(['The input space of PIE.pie2pde.',obj,' does not match the components in PIE.space.',obj,'.'])
    elseif ~isempty(X_pde.ops) && ~isempty(dims_pde) && ~isequal(X_pde.dims.out(:),dims_pde)
        error(['The output space of PIE.pie2pde.',obj,' does not match the components in PIE.space_pde.',obj,'.'])
    end
    % The map from PDE to PIE has the PIE components as output.
    for field = fieldnames(PIE.pde2pie.(obj))'
        val = PIE.pde2pie.(obj).(field{1});
        if strcmp(field{1},'C')
            if ~isempty(val.ops) && ~isequal(val.dims.out(:),dims_pie)
                error(['The output space of PIE.pde2pie.',obj,'.C does not match the components in PIE.space.',obj,'.'])
            end
        elseif ~isempty(val) && size(val,1)~=numel(dims_pie)
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

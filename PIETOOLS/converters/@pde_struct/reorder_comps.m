function [PDE,comp_order] = reorder_comps(PDE,obj,suppress_summary)
% This function reorders the rows of tables PDE.x_tab, PDE.w_tab,
% PDE.u_tab, PDE.y_tab and PDE.z_tab, in preparation of converting the PDE
% to a PIE.
%
% INPUTS:
% - PDE:    A struct or pde_struct type object, defining a PDE in the terms
%           format (see also the "@pde_struct/initialize" function).
% - obj:    'x', 'y', 'z', 'u', 'w', or 'BC', indicating for what object to
%           reorder the components. Defaults to 'all', in which case all
%           components of all objects are reordered.
% - suppress_summary:   Logical value indicating whether to suppress the
%                       summary of the reorder process, defaults to false.
%
%
% OUTPUTS:
% - PDE:    A structure of the same type as the input, with the table
%           PDE.([obj'_tab']) and the fields PDE.(obj) adjusted to match
%           the new order of the components of obj. If obj='x', obj='u' or
%           obj='w', the terms in each of the PDE equations are also
%           adjusted to make sure each term still refers to the appropriate
%           component of PDE.x, PDE.u, or PDE.w.
% - comp_order: An nx1 array indicating the new order of the components of
%               obj, so that PDE_new.(obj){j} = PDE_old.(obj){comp_order(j)}.
%               If multiple objects are specified. comp_order will be a
%               struct with field "comp_order.obj" specifying how the
%               components of "obj" have been reordered.
%
% NOTES:
% The components are seperated first based on which variables they depend
% on, combining e.g. the state components x_{j} into a full state
% x = [x0; x1(s1); x2(s2); x3(s1,s2)];
% For each function space, the components x_{j} within x1, x2 and x3 are
% then ordered based on their order of differentiability wrt each of the
% variables they depend on. For example, letting x_{i,j}(s1,s2) within x2
% be differentiable up to degree i wrt s1 and degree j wrt s2, they will be
% ordered as
% x3 = [x_{0,0}; x_{1,0}; ...; x_{N1,0}; x_{0,1}; x_{1,1}; ...; x_{N1,N2}];
% Finally, if different state components are differentiable up to the same
% degree, they will be ordered based on their original order in the PDE
% structure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Initial coding DJ - 08/09/2022
%

% % % Check the input arguments
if nargin==1
    % If no particular object is specified, reorder all components
    obj = 'all';
    suppress_summary = false;
elseif nargin==2
    suppress_summary = false;    
end

% If multiple objects are specified, run the function for each object
% separately.
if isa(obj,'char') && strcmp(obj,'all')
    obj = {'x','w','u','z','y','BC'};
end
if isa(obj,'cell')
    comp_order = struct();
    for k=1:numel(obj)
        obj_k = obj{k};
        [PDE,comp_order_obj] = reorder_comps(PDE,obj_k,suppress_summary);
        comp_order.(obj_k) = comp_order_obj;
    end
    return
elseif ~ismember(obj,{'x','w','u','z','y','BC'})
    error('The second argument must be one of {''x'',''w'',''u'',''z'',''y'',''BC'',''all''}.')
end


% % % Run the actual function.
% Extract PDE information.
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
% If there is no reordering to be done, just return.
if numel(PDE.(obj))==0 || numel(PDE.(obj))==1
    comp_order = ones(numel(PDE.(obj)),1);
    return
end
nvars = PDE.dim;

% Re-order the rows of obj_tab such that:
% - the index of the components increases fastest;
% - the order of differentiability wrt each of the nvars variables
%   increase next, with the order of differentiability wrt the first
%   variable increasing fastest and the order of differentiability wrt the
%   last variable increases slowest;
% - the dependence of the component on each of the nvars variables
%   increase last, with the dependence on the first variable increasing
%   fastest and the dependence on the last variable increasing slowest.
% For inputs, there is no order of differentiability.
if strcmp(obj,'x') || strcmp(obj,'BC')
    % For the state components, we order based on order of
    % differentiability as well.
    comp_tab = PDE.([obj,'_tab']);
    dep_tab = comp_tab(:,3:2+nvars);
    diff_tab = comp_tab(:,3+nvars:2*nvars+2);
    comp_tab_alt = [dep_tab(:,end:-1:1), diff_tab(:,end:-1:1), comp_tab(:,1), comp_tab(:,2)];
    comp_tab_alt = sortrows_integerTable(comp_tab_alt);
    obj_tab_new = [comp_tab_alt(:,end-1:end), comp_tab_alt(:,nvars:-1:1), comp_tab_alt(:,2*nvars:-1:nvars+1)];
else
    % For inputs, we have no order of differentiability, but we do need to
    % adjust the terms involving inputs to refer to the appropriate input.
    comp_tab = PDE.([obj,'_tab']);
    dep_tab = comp_tab(:,3:2+nvars);
    comp_tab_alt = [dep_tab(:,end:-1:1), comp_tab(:,1), comp_tab(:,2)];
    comp_tab_alt = sortrows_integerTable(comp_tab_alt);
    obj_tab_new = [comp_tab_alt(:,end-1:end), comp_tab_alt(:,nvars:-1:1)];
end

if strcmp(obj,'x') || strcmp(obj,'u') || strcmp(obj,'w')
    % Adjust the terms in each equation to match the new order of the
    % components.
    [~,new_order] = sort(obj_tab_new(:,1));    % comp_tab_new(new_order,:) = comp_tab_old
    for ii=1:numel(PDE.x)
        for jj=1:numel(PDE.x{ii}.term)
            if isfield(PDE.x{ii}.term{jj},obj)
                PDE.x{ii}.term{jj}.(obj) = new_order(PDE.x{ii}.term{jj}.(obj));
            end
        end
    end
    for ii=1:numel(PDE.y)
        for jj=1:numel(PDE.y{ii}.term)
            if isfield(PDE.y{ii}.term{jj},obj)
                PDE.y{ii}.term{jj}.(obj) = new_order(PDE.y{ii}.term{jj}.(obj));
            end
        end
    end
    for ii=1:numel(PDE.z)
        for jj=1:numel(PDE.z{ii}.term)
            if isfield(PDE.z{ii}.term{jj},obj)
                PDE.z{ii}.term{jj}.(obj) = new_order(PDE.z{ii}.term{jj}.(obj));
            end
        end
    end
    for ii=1:numel(PDE.BC)
        for jj=1:numel(PDE.BC{ii}.term)
            if isfield(PDE.BC{ii}.term{jj},obj)
                PDE.BC{ii}.term{jj}.(obj) = new_order(PDE.BC{ii}.term{jj}.(obj));
            end
        end
    end
end
% Establish the new order of the components.
comp_order = obj_tab_new(:,1);
% Re-arrange the equations to match the new order.
PDE.(obj) = PDE.(obj)(comp_order);
% Assign new indices 1:ncomps to the components.
obj_tab_new(:,1) = 1:size(obj_tab_new,1);
PDE.([obj,'_tab']) = obj_tab_new;
% All fields of the returned PDE should still be appropriately specified.
PDE.is_initialized = true;

% Print a summary, if desired.
if ~suppress_summary
    print_reorder_summary(PDE,obj,comp_order)
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_reorder_summary(PDE,obj,comp_order)
% print_reorder_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning the new order of
% the components "obj" in the PDE structure.
%
% INPUTS:
% - PDE:    A "pde_struct" class object defining a PDE.
% - obj:    Char 'x', 'u', 'w', 'y', 'z', or 'BC', indicating for which
%           object to display the new order of the variables.
% - comp_order: comp_order(j) provides the index of the component in the
%               original PDE associated to component j in the new PDE.
%
% OUTPUTS:
% Displays information in the command window concerning the new order of
% the components in PDE.obj, compared to the original PDE.
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
ncomps = numel(PDE.(obj));
if all(comp_order == (1:ncomps)')
    % The order of the components has not changed.
    fprintf(['\n','The order of the ',object_name,'s ',obj,' has not changed.\n']);
    return
else
    % Otherwise, we list the new order of the components.
    fprintf(['\n','The ',object_name,'s have been reordered as:\n']);
end

% Use UNICODE to add subscript indices to different components.
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
% Determine how many subscripts will be needed.
n_digits = 1+floor(log(ncomps)/log(10));


% % For the purpose of clearer display, we make an estimate of how long the
% % display will be. 
% First, determine which variables appear in the different state variables.
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
% Then, estimate the size of the string of characters denoting the components
% "PDE.x{ii}".
nvars_max = max(sum(PDE.([obj,'_tab'])(:,3:2+nvars),2));
lngth_varnames_mean = ceil(length(global_varnames)*nvars_max/nvars);
LHS_length_max = 1+1 + n_digits + lngth_varnames_mean+3; % e.g. ' x13(t,s1,s2,s3)', 


% % For each of the components, display its size, and which variables it
% % depends on.
for ii=1:ncomps
    old_idx = comp_order(ii);
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
    if numel(PDE.x)==1
        % There is only one state component --> no need to give index
        Lcomp_idx = '';
    elseif numel(PDE.x)<=9
        % The state number consists of a single decimal.
        Lcomp_idx = sub_num{old_idx+1};
        LHS_length = LHS_length + 1;
    else
        % The state number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(old_idx)')+1)');
        LHS_length = LHS_length + length(num2str(old_idx));
    end
    % Set the name of the component, including its depdence on spatial
    % variables.
    LHS_name = [' x',Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = 1 + LHS_length + 3;
        
    % For the added state components, indicate of which state component
    % they are the temporal derivative.
    new_idx = ii;
    Rcomp_idx = cell2mat(sub_num(str2num(num2str(new_idx)')+1)');
    RHS = [' -->   x',Rcomp_idx,'(',varnames_ii_t,')'];
%    RHS = '';
    
    % % % Finally, display:
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
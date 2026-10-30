function [PIE,tdiff_tab] = expand_tderivatives(PIE,tdiff_list,suppress_summary)
% PIE_OUT = EXPAND_TDERIVATIVES(PIE_IN,TDIFF_LIST,SUPPRESS_SUMMARY)
% takes a PIE involving higher-order temporal derivatives, and introduces
% auxiliary fundamental state components to model these higher-order
% temporal derivatives, representing the PIE as a first-order in time
% system.
%
% INPUTS
% - PIE:        'pie_struct' object representing a PIE with nx state
%               evolution equations, some of which involve higher-order
%               temporal derivatives.
% - tdiff_list: nx x 1 array of type 'double', specifying for each of the
%               nx fundamental state components the order of the temporal
%               derivative in the evolution equation of this state
%               component.
% - suppress_summary: (optional) boolean variable indicating whether to
%                   suppress the display of a summary of how the state
%                   components in the output PIE relate to those in the
%                   input PIE.
%
% OUTPUTS
% - PIE:    'pie_struct' object representing the same PIE but now involving
%           only first-order in time derivatives, introducing auxiliary
%           variables to model the higher-order temporal derivatives. See
%           also the notes.
% - tdiff_tab:  nx x 2 array specifying for each of the fundamental state
%               components in the new PIE how they relate to the state
%               components in the original PIE. In particular, for each
%               state component i in the output PIE, we have
%                  x_{i} = (d/dt)^{tdiff_tab(i,2)} x_{old,tdiff_tab(i,1)}.
%
% NOTES
% For each state component in the PIE, the function checks whether any of
% the associated evolution equations involves a higher-order temporal
% derivative of this component. If so, auxiliary variables are added to
% represent the higher-order temporal derivatives of this state component,
% expressing the dynamics in terms of only first-order temporal
% derivatives. For example, if equation ii involves the pth-order
% derivative of two state components,
%   (d/dt)^p Top(i,:)*v = Top(i,i)*(d/dt)^p v_{i} + Top(i,j)*(d/dt)^p v_{j}
%            Aop(i,:)*v = Aop(i,i)*v_{i} + Aop(i,k)*v_{k}
% then we introduce v_{i,l} = d/dt^(l-1) v_{i} and 
% v_{j,l} = (d/dt)^(l-1) v_{j} for l ranging from 1 to p 
% (so that v_{i,1}=v_{i}), and expand the dynamics as
%   (d/dt)Top_new(i,:) v = d/dt v_{i,l} = v_{i,l+1} = Aop_new(i,:)*v
% for l from 1 to p-1, and
%   (d/dt)Top_new(i+p,:) v = Top_new(i+k,i+p)*(d/dt)v_{i,k} + Top(i+p,j+q)*(d/dt)v_{j,p}
%                          = Aop_new(i+p,i)*v_{i} + Aop_new(i+p,k)*v_{k}
% NOTE that this means that Top*v may no longer return the PDE state (if it
% did so previously), as now the dynamics for state variable v_{i}=v_{i,1}
% are given by    d/dt v_{i,1} = v_{i,2},   so that Top_new(i,i) = I.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2025 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, 06/18/2025;


% % % Process the inputs
% Check that sufficient arguments are passed.
if nargin<2
    error("Insufficient input arguments.")
elseif nargin==2
    % Print summary by default.
    suppress_summary = false;
elseif nargin>3
    error("Too many input arguments.")
end
% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error("First argument must be of type 'pie_struct.")
else
    PIE = initialize(PIE);
end
% Check if the orders of the temporal derivatives are properly specified.
if ~isnumeric(tdiff_list) || any(tdiff_list<=0) || any(tdiff_list~=round(tdiff_list))
    error("Orders of temporal derivatives should be specified as nx1 array of positive integers.")
end
% Make sure we have a temporal derivative for each state variable.
tdiff_list = tdiff_list(:);
nx = size(PIE.T,1);
if isempty(PIE.x_tab)
    % Assume each state component is scalar-valued
    ncomps = nx;
    comp_IDs = (1:nx)';
    comp_IDs_full = comp_IDs;
    n_vars = ones(nx,1);
    nn_vars = (0:nx)';
    var_idcs = mat2cell([(1:nx)',(1:nx)'],ones(nx,1),2);
else
    % Get size of each state component from second column of x_tab
    ncomps = size(PIE.x_tab,1);
    comp_IDs = (1:ncomps)';
    n_vars = PIE.x_tab(:,2);
    comp_IDs_full = repelem(comp_IDs,n_vars);
    nn_vars = [0;cumsum(PIE.x_tab(:,2))];
    % Establish what variable indices in the state correspond to each
    % component
    var_idcs = cellfun(@(a)(a(1):a(2))',mat2cell([nn_vars(1:end-1)+1,nn_vars(2:end)],ones(ncomps,1),2),'UniformOutput',false);
end
if numel(tdiff_list)==ncomps
    tdiff_list = tdiff_list(comp_IDs_full);
elseif numel(tdiff_list)~=nx
    error("Number of temporal derivatives should match number of evolution equations in the PIE.")
end



% % % Determine what temporal derivatives must be taken, and of which state
% % % variables.

% Check which evolution equations involve higher-order temporal derivatives
eq_expnd = tdiff_list>1;
tdiff_list = tdiff_list(eq_expnd)-1;
if ~any(eq_expnd)
    % All evolution equations are first-order in time
    % --> nothing to do here.
    tdiff_tab = [(1:ncomps)',zeros(nx,1)];
    if ~suppress_summary
        fprintf(['\n','No higher-order temporal derivatives were encountered; returning original PIE'])
    end
    return
elseif any(~(PIE.Tw(eq_expnd,:)==0)) || any(~(PIE.Tu(eq_expnd,:)==0))
    error("Expansion of higher-order temporal derivatives is not supported for PDE states with boundary inputs.")
end

% Check which state variables require higher order-temporal derivatives.
Top_tdiff = PIE.T(eq_expnd,:);
comp_tdiffs = zeros(ncomps,1);
for tt = unique(tdiff_list)'
    r_idcs = tdiff_list==tt;
    Top_tt = Top_tdiff(r_idcs,:);
    for jj=1:ncomps
        c_idcs = var_idcs{jj};
        if any(~(Top_tt(:,c_idcs)==0))
            comp_tdiffs(jj) = tt;
        end
    end
end
comp_expnd = comp_tdiffs>0;
comp_tdiffs = comp_tdiffs(comp_expnd);



% % % Augment the fundamental, introducing auxiliary variables representing
% % % temporal derivatives of state variables where necessary.

% Add higher-order temporal derivatives.
n_aux = sum(comp_tdiffs);                          % number of auxiliary vars
aux_IDs = repelem(comp_IDs(comp_expnd),comp_tdiffs);
n_vars_aux = n_vars(aux_IDs);
nn_vars_aux = nn_vars(end) + cumsum(n_vars_aux);
% Assign to each auxiliary component a set of new state variables.
var_idcs_aux = cellfun(@(a)(a(1):a(2))',mat2cell([[nn_vars(end)+1;nn_vars_aux(1:end-1)+1],nn_vars_aux(1:end)],ones(n_aux,1),2),'UniformOutput',false);
% Also keep track of which of the original state variables are associated
% with each auxiliary component
var_idcs_aux_old = repelem(var_idcs(comp_expnd),comp_tdiffs);

% Update the state information table to incorporate the new variables
if ~isempty(PIE.x_tab)
    aux_tab = PIE.x_tab(aux_IDs,:);
    %aux_tab(:,2+PIE.dim+(1:PIE.dim)) = 0;      % should we set the order of regularity to zero?
    %aux_tab(:,1) = ncomps+(1:n_aux)';
    PIE.x_tab = [PIE.x_tab; aux_tab];
end

% Augmenting the state, for each state component, keep track of what order
% temporal derivative of which original state component it represents
tot_IDs = [comp_IDs; aux_IDs];                  % associated original state component
aux_tdiffs = cell2mat(cellfun(@(a)(1:a)',num2cell(comp_tdiffs),'UniformOutput',false));
comp_tdiffs = [zeros(ncomps,1); aux_tdiffs];    % order of temporal derivative
n_vars_tot = [n_vars; n_vars_aux];
var_idcs_tot = [var_idcs; var_idcs_aux];
var_idcs_old = [var_idcs; var_idcs_aux_old];
ncomps_tot = ncomps + n_aux;

% We reorder the state components so that the auxiliary variables
% immediately follow their corresponding original state
[tot_IDs_new,old2new_idcs] = sortrows_integerTable(tot_IDs);    % tot_IDs_new = tot_IDs(old2new_idcs);
n_vars_new = n_vars_tot(old2new_idcs);
nn_vars_new = [0;cumsum(n_vars_new)];
var_idcs_new = cellfun(@(a)(a(1):a(2))',mat2cell([nn_vars_new(1:end-1)+1,nn_vars_new(2:end)],ones(ncomps_tot,1),2),'UniformOutput',false);
var_idcs_old = var_idcs_old(old2new_idcs);
% We update the order of the state variables (comprising the state
% components) accordingly
old2new_idcs_full = cell2mat(var_idcs_tot(old2new_idcs));   % x_new = x_old(old2new_idcs);
[~,new2old_idcs_full] = sort(old2new_idcs_full);            % x_old = x_new(new2old_idcs);  

% Compute updated dimensions of the opvar object
nx_op = PIE.T.dim(:,1);
nnx_op = [0;cumsum(nx_op)];
nx_op_new = zeros(size(nx_op));
for jj=1:numel(nx_op)
    is_jj_comps = nn_vars(2:end)>nnx_op(jj) & nn_vars(2:end)<=nnx_op(jj+1);
    comp_IDs_jj = comp_IDs(is_jj_comps);
    nx_op_new(jj) = sum(n_vars_new(ismember(tot_IDs_new,comp_IDs_jj)));
end
comp_tdiffs_new = comp_tdiffs(old2new_idcs);
if ~isempty(PIE.x_tab)
    PIE.x_tab = PIE.x_tab(old2new_idcs,:);
end
    


% % % Update the PIE structure to express dynamics in terms of augmented
% % % state

% % Update the columns of the various PI operators to account for the new
% % order of the state components
fnames = {'A','C1','C2'};
for kk=1:numel(fnames)
    op_name = fnames{kk};
    % Declare zero operator with column dimension matching new state size.
    nr_op = PIE.(op_name).dim(:,1);
    if PIE.dim==1
        Pop_new = opvar();
    else
        Pop_new = opvar2d();
    end
    Pop_new.dim = [nr_op,nx_op_new];
    Pop_new.I = PIE.dom;
    Pop_new.var1 = PIE.vars(:,1);       Pop_new.var2 = PIE.vars(:,2);
    % Set terms of original operator in appropriate columns of Pop_new.
    if any(nr_op)
        Pop_new(:,new2old_idcs_full(1:nx)) = PIE.(op_name);
    end
    PIE.(op_name) = Pop_new;
end


% % Express higher-order temporal derivative in old evolution equations
% % in terms of auxiliary variables,
% %   Top(i,i)*(d/dt)^{k}v_{i} + Top(i,j)*(d/dt)^{k}v_{j} = Aop(i,:)*v
% %   => Top(i,i)*(d/dt)v_{i,k-1} + Top(i,j)*(d/dt)v_{j,k-1} = Aop(i,:)*v
% First, declare a zero operator with appropriate column dimension.
if PIE.dim==1
    Top_new = opvar();
else
    Top_new = opvar2d();
end
Top_new.dim = [nx_op,nx_op_new];
Top_new.I = PIE.dom;
Top_new.var1 = PIE.vars(:,1);       Top_new.var2 = PIE.vars(:,2);
% Set terms associated to first-order in time evolution equations.
Top_new(~eq_expnd,new2old_idcs_full(1:nx)) = PIE.T(~eq_expnd,:);
% Set terms associated to higher-order in time evolution equations, using
% the auxiliary variables to represent higher-order derivatives.
Top_new_tdiff = Top_new(eq_expnd,:);
for tt = unique(tdiff_list')
    r_idcs = tdiff_list==tt;
    is_tt_ID = comp_tdiffs_new==tt;
    c_idcs_new = cell2mat(var_idcs_new(is_tt_ID));
    c_idcs_old = cell2mat(var_idcs_old(is_tt_ID));
    Top_new_tdiff(r_idcs,c_idcs_new) = Top_tdiff(r_idcs,c_idcs_old);
end
Top_new(eq_expnd,:) = Top_new_tdiff;
PIE.T = Top_new;


% % Declare evolution equations for the auxiliary state variables,
% %   (d/dt) v_{i,k} = v_{i,k+1}
fnames = {'T','Tw','Tu','A','B1','B2'};
[~,strt_idcs] = unique(tot_IDs_new);            % idcs at which new component starts
end_idcs = [strt_idcs(2:end)-1; ncomps_tot];
end_idcs_full = cell2mat(var_idcs_new(end_idcs));
for kk=1:numel(fnames)
    op_name = fnames{kk};
    % Declare zero operators of row dimensions matching new state size.
    nc_op = PIE.(op_name).dim(:,2);
    if PIE.dim==1
        Pop_new = opvar();
    else
        Pop_new = opvar2d();
    end
    Pop_new.dim = [nx_op_new,nc_op];
    Pop_new.I = PIE.dom;
    Pop_new.var1 = PIE.vars(:,1);       Pop_new.var2 = PIE.vars(:,2);
    % Set original PIE equations in appropriate rows of Pop_new
    if any(nc_op)
        Pop_new(end_idcs_full,:) = PIE.(op_name);
    end
    PIE.(op_name) = Pop_new;
end
% Declare the equations of the form (d/dt) v_{i,k} = v_{i,k+1}
for ii=1:length(strt_idcs)
    r_idcs = cell2mat(var_idcs_new(strt_idcs(ii):end_idcs(ii)-1));
    c_idcs = cell2mat(var_idcs_new(strt_idcs(ii)+1:end_idcs(ii)));
    PIE.T(r_idcs,r_idcs) = eye(numel(r_idcs));
    PIE.A(r_idcs,c_idcs) = eye(numel(r_idcs));
end



% % % Process the outputs, and print a summary of the work, if desired.
tdiff_tab = [tot_IDs_new, comp_tdiffs_new];
if ~suppress_summary
    print_tdiff_expansion_summary(PIE,tdiff_tab);
end

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_tdiff_expansion_summary(PIE,tdiff_tab)
% print_tdifff_expanion_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning how the state
% variables in the PIE structure 'PIE' relate to the state variabels in the
% original PIE specified by the user in the call to "expand_tderivatives".
%
% INPUTS:
% - PIE:    A "pie_struct" class object defining a PIE.
% - tdiff_tab:  A nx_new x 3 array of integers indicating for each
%               of the j=1:nx state components what component i in 
%               the original PIE they correspond to, and what order
%               temporal derivative of that variable:
%                   x_{tdiff_tab(j,1)} = d/dt^{tdiff_tab(j,3)} x_{tdiff_tab(j.2)}.
%
% OUTPUTS:
% Displays information in the command window concerning how the state
% variables in the "PIE" structure relate to the fundamental state
% variables in the original PIE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncomps = size(tdiff_tab,1);
fprintf(['\n','Introduced the following (vector-valued) fundamental state components:\n'])

% Use UNICODE to add subscript indices to different components.
thin_space = char(8201);
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
partial = '\x2202';
sub_t = '\x209C';
%int = '\x222B';
%xdot = '\x1E8B';
%xdot = [partial,sub_t,' x'];
% Also use UNICODE for superscript orders of the temporal derivative
sup_num = cell(10,1);    % Superscripts
sup_num{1} = '\x2070';  % Superscript 0
sup_num{2} = '\xB9';
sup_num{3} = '\xB2';
sup_num{4} = '\xB3';
sup_num{5} = '\x2074';
sup_num{6} = '\x2075';
sup_num{7} = '\x2076';
sup_num{8} = '\x2077';
sup_num{9} = '\x2078';
sup_num{10} = '\x2079';

% Determine how many subscripts will be needed.
n_digits = 1+floor(log(ncomps)/log(10));


% % For the purpose of clearer display, we make an estimate of how long the
% % display will be. 
% First, determine which variables appear in the different state variables.
nvars = size(PIE.vars,1);
global_vars_obj = PIE.vars;
lngth_varnames = 3;      % just three characters: (t)
if ~isempty(global_vars_obj)
    varname_cell = cell(nvars,1);
    if nvars==1
        varname_cell{1} = 's';
        lngth_varnames = lngth_varnames+2;    % add two characters: ,s
    else
        varname_cell{1} = ['s',sub_num{2}];
        lngth_varnames = lngth_varnames+2;    % add three characters: ,s1
        for kk=2:nvars
            varname_cell{kk} = ['s',cell2mat(sub_num(str2num(num2str(kk)')+1)')];
            lngth_varnames = lngth_varnames+2;    % add three characters: ,sj
        end
    end
end
% Then, estimate the size of the string of characters denoting the components
% "PDE.x{ii}".
LHS_length_max = 1+1 + n_digits + lngth_varnames+2; % e.g. ' x13(t,s1,s2,s3)', 


% % For each of the state variables in the PIE structure, display
% % which of the state variables in the original PIE it corresponds to,
% % and what order temporal derivative.
for ii=1:ncomps
    Lstate_num = ii;   % State index associated to the new state.
    
    % Establish the names of the variables on which the component depends.
    varnames_ii_t = 't';
    var_idcs = find(PIE.x_tab(ii,3:2+PIE.dim));
    for kk=var_idcs
        varnames_ii_t = [varnames_ii_t,',',varname_cell{kk}];
    end
    
    % Establish the (subscript) index for the component.
    LHS_length = length(varnames_ii_t);
    if ncomps<=9
        % The state number consists of a single decimal.
        Lcomp_idx = sub_num{Lstate_num+1};
        LHS_length = LHS_length + 1;
    else
        % The state number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(Lstate_num)')+1)');
        LHS_length = LHS_length + length(num2str(Lstate_num));
    end
    % Set the name of the component, including its depdence on spatial
    % variables.
    LHS_name = [' x',Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = 1 + LHS_length + 3;
        
    % For the added state components, indicate of which state component
    % they are the temporal derivative.
    Rstate_num = tdiff_tab(ii,1);
    Rcomp_idx = cell2mat(sub_num(str2num(num2str(Rstate_num)')+1)');

    % Also indicate the order of the derivative
    if tdiff_tab(ii,2)==0
        tdiff_str = ['   ',thin_space];
    elseif tdiff_tab(ii,2)==1
        tdiff_str = [' ',partial,sub_t,' '];
    else
        t_sup = cell2mat(sup_num(str2num(num2str(tdiff_tab(ii,2))')+1)');
        tdiff_str = [partial,sub_t,t_sup,' '];
    end
    RHS = [' <-- ',tdiff_str,'x',Rcomp_idx,'(',varnames_ii_t,')'];
    
    % % % Finally, display
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
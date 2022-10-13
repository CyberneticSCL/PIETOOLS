function [PDE] = expand_delays(PDE)


PDE = initialize(PDE,true);


tau_vars = [PDE.tau(:,1),zeros(size(PDE.tau,1),1)];
for vv=1:size(tau_vars,1)
    tau_varname_vv = tau_vars(vv,1).varname{1};
    tau_varname_vv2 = [tau_varname_vv,'_dum'];
    tau_vars(vv,2) = polynomial({tau_varname_vv2});
end
tau_dom = [-abs(double(PDE.tau(:,2))),zeros(size(PDE.tau,1),1)];



nx = numel(PDE.x);
nw = numel(PDE.w);
nu = numel(PDE.u);
nz = numel(PDE.z);
ny = numel(PDE.y);
nBC = numel(PDE.BC);

% Expand the tables to include the delay variables
nvars = size(PDE.vars,1);
ndelays = size(tau_vars,1);

PDE.vars = [PDE.vars; tau_vars];
PDE.dom = [PDE.dom; tau_dom];

PDE.x_tab = [PDE.x_tab(:,1:2+nvars),zeros(nx,ndelays),PDE.x_tab(:,3+nvars:2+2*nvars),zeros(nx,ndelays)];
PDE.w_tab = [PDE.w_tab(:,1:2+nvars),zeros(nw,ndelays)];
PDE.u_tab = [PDE.u_tab(:,1:2+nvars),zeros(nu,ndelays)];
PDE.z_tab = [PDE.z_tab(:,1:2+nvars),zeros(nz,ndelays)];
PDE.y_tab = [PDE.y_tab(:,1:2+nvars),zeros(ny,ndelays)];
PDE.BC_tab = [PDE.BC_tab(:,1:2+nvars),zeros(nBC,ndelays),PDE.BC_tab(:,3+nvars:2+2*nvars),zeros(nBC,ndelays)];


% Loop over the terms, 
state_dep_tab = zeros(0,2);
[PDE,state_dep_tab] = expand_terms(PDE,state_dep_tab,nvars,'x',nx);
[PDE,state_dep_tab] = expand_terms(PDE,state_dep_tab,nvars,'z',nz);
[PDE,state_dep_tab] = expand_terms(PDE,state_dep_tab,nvars,'y',ny);
[PDE,state_dep_tab] = expand_terms(PDE,state_dep_tab,nvars,'BC',nBC);

% All delays have been replaced by spatial variables 
%   --> remove the delays.
PDE.tau = zeros(0,2);

% Display summary of the results.
ndelay_states = size(state_dep_tab,1);
ndelay_vars = size(PDE.vars,1)-nvars;
if ndelay_states==0
    fprintf('\n   Encountered no delays.\n')
else
    %print_summary_expand_delay(PDE,nx)
    if ndelay_vars==1
        fprintf(['\nAdded ',num2str(ndelay_vars),' spatial variable to represent the delays.\n'])
    else
        fprintf(['\nAdded ',num2str(ndelay_vars),' spatial variables to represent the delays.\n'])
    end
    if ndelay_states==1
        fprintf(['\nAdded ',num2str(ndelay_states),' state component.\n'])
    else
        fprintf(['\nAdded ',num2str(ndelay_states),' state components.\n'])
    end
end


end



function [PDE,state_dep_tab] = expand_terms(PDE,state_dep_tab,nvars_old,Lobj,ncomps)


ndelays = size(PDE.vars,1) - nvars_old;
nvars = size(PDE.vars,1);
nx = numel(PDE.x);

% Extract variable names
tau_name = cell(ndelays,1);
for vv=1:ndelays
    % varnames are extracted 1-by-1, as "polynomial" stores them
    % alphabetically.
    tau_name{vv} = PDE.tau(vv,1).varname{1};
end

for ii=1:ncomps
    nterms = numel(PDE.(Lobj){ii}.term);
    % % Loop over all terms and replace delays with new state.
    for jj=1:nterms
        term_jj = PDE.(Lobj){ii}.term{jj};
        
        % % Extract the delay.
        delay = term_jj.delay;
        if (isa(delay,'double') || isdouble(delay)) && double(delay)==0
            % If there is no delay, there is nothing to do.
            continue
        end
        
        % % Establish whether the term involves a state or input.
        if isfield(term_jj,'x')
            Robj = 'x';
            is_xcomp = true;
            Rindx = term_jj.x;
            state_idx = Rindx;
        elseif isfield(term_jj,'w')
            Robj = 'w';
            is_xcomp = false;
            Rindx = term_jj.w;
            state_idx = Rindx + nx;
        else
            Robj = 'u';
            is_xcomp = false;
            Rindx = term_jj.u;
            state_idx = Rindx + nx + numel(PDE.w);
        end
        
        % Distinguish case where delay is a fixed real value, and where
        % delay is distributed.
        if isa(delay,'double') || isdouble(delay)
            % Fixed delay.
            delay = abs(double(delay));
            % Check if this delay is already included in our table of
            % delays.
            tau_idx = find(delay==abs(double(PDE.tau(:,2))),1,'first');
            if ~isempty(tau_idx)
                % The considered delay already appears in our table.
                % We still have to check if this delay already appears in
                % combination with our considered RHS object Robj.
                tau_idx = tau_idx + nvars_old;
                state_tau_idx = find(all(state_dep_tab==[state_idx,tau_idx],2),1,'first');
                if ~isempty(state_tau_idx)
                    % The combination of Robj and delay is already present
                    % in state_dep_tab
                    % --> Adjust the term to replace delay by new state:
                    %   Robj(t-delay) = xnew(t,-delay)
                    nvars = size(PDE.vars,1);
                    new_state_idx = state_tau_idx + nx;
                    has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                    has_vars_Robj = has_vars_delay_state;
                    has_vars_Robj(tau_idx) = false;
                    new_D = zeros(1,nvars);
                    new_loc = PDE.vars(:,1)';
                    if is_xcomp
                        new_D(has_vars_Robj) = term_jj.D;
                        new_loc(has_vars_Robj) = term_jj.loc;                        
                    end
                    new_loc(tau_idx) = -delay;

                    term_jj = rmfield(term_jj,Robj);
                    term_jj.x = new_state_idx;
                    term_jj.D = new_D(has_vars_delay_state);
                    term_jj.loc = new_loc(has_vars_delay_state);
                    term_jj.delay = 0;
                else
                    % The combination of Robj and delay does not occur yet,
                    % but the delay is already included as one of the
                    % variables.
                    % Add a new state xnew(t,s) = Robj(t-s);
                    PDE = add_delay_state(PDE,Robj,Rindx,tau_idx);
                    nvars = size(PDE.vars,1);
                    new_state_idx = size(PDE.x_tab,1);

                    % Adjust the term to replace delay by new state var:
                    %   Robj(t-delay) = xnew(t,-delay)
                    has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                    has_vars_Robj = has_vars_delay_state;
                    has_vars_Robj(tau_idx) = false;
                    new_D = zeros(1,nvars);
                    new_loc = PDE.vars(:,1)';
                    if is_xcomp
                        new_D(has_vars_Robj) = term_jj.D;
                        new_loc(has_vars_Robj) = term_jj.loc;                        
                    end
                    new_loc(tau_idx) = -delay;

                    term_jj = rmfield(term_jj,Robj);
                    term_jj.x = new_state_idx;
                    term_jj.D = new_D(has_vars_delay_state);
                    term_jj.loc = new_loc(has_vars_delay_state);
                    term_jj.delay = 0;
                end
            else
                % The delay has not been included in our delay table yet
                % --> we have to add it
                % Add the combination of delay and Robj to the table.
                tau_idx = size(PDE.vars,1) + 1;
                state_dep_tab = [state_dep_tab; [state_idx,tau_idx]];

                % % % Add the delay to the PDE as a new spatial variable.
                % Define new spatial variables (s,th) to represent the delay.
                var1 = polynomial({['ntau_',num2str(tau_idx)]});
                var2 = polynomial({['tau_dum',num2str(tau_idx)]});
                PDE.vars = [PDE.vars; [var1,var2]];
                % Add the corresponding dom
                PDE.dom = [PDE.dom; [-delay,0]];
                
                % Add the new spatial variable to all the tables
                PDE.x_tab = [PDE.x_tab(:,1:2+nvars),zeros(size(PDE.x_tab,1),1),PDE.x_tab(:,3+nvars:2+2*nvars),zeros(size(PDE.x_tab,1),1)];
                PDE.w_tab = [PDE.w_tab,zeros(size(PDE.w_tab,1),1)];
                PDE.u_tab = [PDE.u_tab,zeros(size(PDE.u_tab,1),1)];
                PDE.z_tab = [PDE.z_tab,zeros(size(PDE.z_tab,1),1)];
                PDE.y_tab = [PDE.y_tab,zeros(size(PDE.y_tab,1),1)];
                PDE.BC_tab = [PDE.BC_tab(:,1:2+nvars),zeros(size(PDE.BC_tab,1),1),PDE.BC_tab(:,3+nvars:2+2*nvars),zeros(size(PDE.BC_tab,1),1)];
                
                % Add a new state xnew(t,s) = Robj(t-s);
                PDE = add_delay_state(PDE,Robj,Rindx,tau_idx);
                nvars = size(PDE.vars,1);
                new_state_idx = size(PDE.x_tab,1);

                % Adjust the term to replace delay by new state var:
                %   Robj(t-delay) = xnew(t,-delay)
                has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                has_vars_Robj = has_vars_delay_state;
                has_vars_Robj(tau_idx) = false;
                new_D = zeros(1,nvars);
                new_loc = PDE.vars(:,1)';
                if is_xcomp
                    new_D(has_vars_Robj) = term_jj.D;
                    new_loc(has_vars_Robj) = term_jj.loc;
                end
                new_loc(tau_idx) = -delay;

                term_jj = rmfield(term_jj,Robj);
                term_jj.x = new_state_idx;
                term_jj.D = new_D(has_vars_delay_state);
                term_jj.loc = new_loc(has_vars_delay_state);
                term_jj.delay = 0;
            end
        else
            % The delay is distributed: int_{-delay}^{0} Robj(t-tau) dtau.
            %tau_var1 = delay;
            tau_idx = find(ismember(tau_name,delay.varname{1}),1,'first');
            if isempty(tau_idx)
                error('The delay variable does not appear in the list PDE.tau...')
            end
            delay = double(PDE.tau(tau_idx,2));
            tau_idx = tau_idx + nvars_old;
            state_tau_idx = find(all(state_dep_tab==[state_idx,tau_idx],2),1,'first');
            if ~isempty(state_tau_idx)
                % The combination of Robj and delay is already present in
                % state_dep_tab
                % --> Adjust the term to replace delay by new state:
                %   int_{-delay}^{0} Robj(t-tau) dtau = int_{-delay}^{0} xnew(t,tau) dtau
                new_state_idx = nx + state_tau_idx;
                
                % Determine which spatial variables Robj and the delayed
                % state depend on.
                has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                has_vars_Robj = has_vars_delay_state;
                has_vars_Robj(tau_idx) = false;
                
                % Define the new derivative, position, and integral.
                new_D = zeros(1,nvars);
                new_loc = PDE.vars(:,1)';
                new_int = cell(nvars,1);
                if is_xcomp
                    new_D(has_vars_Robj) = term_jj.D;
                    new_loc(has_vars_Robj) = term_jj.loc;
                end
                new_int{tau_idx} = [-delay,0];
                new_int(has_vars_Robj) = term_jj.I;
                
                % Set the new term: 
                % int_{-delay}^{0} Robj(t-tau) dtau = int_{-delay}^{0} xnew(t,tau) dtau
                term_jj = rmfield(term_jj,Robj);
                term_jj.x = new_state_idx;
                term_jj.D = new_D(has_vars_delay_state);
                term_jj.loc = new_loc(has_vars_delay_state);
                term_jj.I = new_int(has_vars_delay_state);
                term_jj.delay = 0;
            else
                % The combination of Robj and delay does not occur yet
                % --> Add the combination to the table and the PDE
                state_dep_tab = [state_dep_tab; [state_idx,tau_idx]];

                % Add a new state xnew(t,s) = Robj(t-s);
                PDE = add_delay_state(PDE,Robj,Rindx,tau_idx);
                nvars = size(PDE.vars,1);
                new_state_idx = size(PDE.x_tab,1);
                
                % Determine which spatial variables Robj and the delayed
                % state depend on.
                has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                has_vars_Robj = has_vars_delay_state;
                has_vars_Robj(tau_idx) = false;
                
                % Define the new derivative, position, and integral.
                new_D = zeros(1,nvars);
                new_loc = PDE.vars(:,1)';
                new_int = cell(nvars,1);
                if is_xcomp
                    new_D(has_vars_Robj) = term_jj.D;
                    new_loc(has_vars_Robj) = term_jj.loc;
                end
                new_int{tau_idx} = [-delay,0];
                new_int(has_vars_Robj) = term_jj.I;
                
                % Set the new term: 
                % int_{-delay}^{0} Robj(t-tau) dtau = int_{-delay}^{0} xnew(t,tau) dtau
                term_jj = rmfield(term_jj,Robj);
                term_jj.x = new_state_idx;
                term_jj.D = new_D(has_vars_delay_state);
                term_jj.loc = new_loc(has_vars_delay_state);
                term_jj.I = new_int(has_vars_delay_state);
                term_jj.delay = 0;
            end
        end
        PDE.(Lobj){ii}.term{jj} = term_jj;
    end
end


end



function PDE = add_delay_state(PDE,Robj,Rindx,tau_var_idx)

full_vars = PDE.vars;
full_dom = PDE.dom;
nvars = size(full_vars,1);
Robj_tab = [PDE.([Robj,'_tab'])(Rindx,1:2+nvars), zeros(1,nvars)];
Robj_size = Robj_tab(1,2);

% % % Add a new state component to x_tab.
% The new state should depend on the same variables as Robj.
PDE.x_tab = [PDE.x_tab; Robj_tab];
has_vars_Robj = logical(PDE.x_tab(end,3:2+nvars));
% The new state should also depend on the new delay variable.
PDE.x_tab(end,[2+tau_var_idx,2+nvars+tau_var_idx]) = 1;
has_vars_delay_state = logical(PDE.x_tab(end,3:2+nvars));



% % % Add the new state xnew(t,s) = Robj(t-s).
% Initialize the new state.
PDE.x = [PDE.x; struct()];
PDE.x{end}.size = Robj_size;
PDE.x{end}.vars = full_vars(has_vars_delay_state,:);
PDE.x{end}.dom = full_dom(has_vars_delay_state,:);
PDE.x{end}.diff = [zeros(1,tau_var_idx-1),1,zeros(1,nvars-tau_var_idx)]; % Only explicitly differentiable wrt the delay variable
PDE.x{end}.diff = PDE.x{end}.diff(has_vars_delay_state);     
PDE.x{end}.tdiff = 1;

% Set up transport equation 
%   d.dt xnew(t,s) = d/ds xnew(t,s)
PDE.x{end}.term{1}.x = numel(PDE.x);
PDE.x{end}.term{1}.D = PDE.x{end}.diff;
PDE.x{end}.term{1}.loc = full_vars(has_vars_delay_state,1)';
PDE.x{end}.term{1}.C = eye(Robj_size);
PDE.x{end}.term{1}.I = cell(sum(has_vars_delay_state),1);
PDE.x{end}.term{1}.delay = 0;


% % % Add the boundary condition xnew(t,s=0) = Robj(t)
% Initialize the BC
PDE.BC = [PDE.BC; struct()];
PDE.BC{end}.size = Robj_size;
PDE.BC{end}.vars = full_vars(has_vars_Robj,:);
PDE.BC{end}.dom = full_dom(has_vars_Robj,:);

% Set up the BC:   Robj(t) = xnew(t,s=0)
% First term: Robj(t)
if strcmp(Robj,'x')
    PDE.BC{end}.term{1}.x = Rindx;
    PDE.BC{end}.term{1}.D = zeros(1,sum(has_vars_Robj));
    PDE.BC{end}.term{1}.loc = full_vars(has_vars_Robj,1)';
elseif strcmp(Robj,'w')
    PDE.BC{end}.term{1}.w = Rindx;
else
    PDE.BC{end}.term{1}.u = Rindx;
end
PDE.BC{end}.term{1}.C = eye(Robj_size);
PDE.BC{end}.term{1}.I = cell(sum(has_vars_Robj),1);
PDE.BC{end}.term{1}.delay = 0;

% Second term: -xnew(t,s=0) 
PDE.BC{end}.term{2}.x = numel(PDE.x);
PDE.BC{end}.term{2}.D = zeros(1,sum(has_vars_delay_state));
PDE.BC{end}.term{2}.loc = full_vars(:,1)';
PDE.BC{end}.term{2}.loc(tau_var_idx) = 0;
PDE.BC{end}.term{2}.loc = PDE.BC{end}.term{2}.loc(has_vars_delay_state);
PDE.BC{end}.term{2}.C = -eye(Robj_size);
PDE.BC{end}.term{2}.I = cell(sum(has_vars_delay_state),1);
PDE.BC{end}.term{2}.delay = 0;

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_summary_expand_delay(PDE,nx_old)
% print_initialization_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning how many
% state components have been added to the system.
%
% INPUTS:
% - PDE:        A "struct" or "pde_struct" class object defining a PDE.
% - ncomps_x:   Integer scalar specifying how many state components were
%               present in the PDE before expanding delays. 
%
% OUTPUTS:
% Displays information in the command window concerning the number of
% added components of type obj, along with the size of each of these
% components, what variables they depend on, and (if obj=='x') to what
% order they are differentiable in each of these variables.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


object_name = 'state component';
% If there are multiple components of the object, we add an 's' at the end.
num_comps = numel(PDE.x);
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
global_vars_obj = PDE.vars(any(PDE.x_tab(:,3:2+nvars),1),:);
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
LHS_length_max = 1+length('x')+n_digits+length(global_varnames) + 3; % e.g. ' x13(t,s1,s2,s3)', 
if num_comps==1
    LHS_length_max = LHS_length_max - 1;    % No subscript index.
end

% For each of the components, display its size, and which variables it
% depends on.
for ii=nx_old+1:num_comps
    comp_ii = PDE.x{ii};
    
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
        Lcomp_idx = sub_num{ii+1};
        LHS_length = LHS_length + 1;
    else
        % The state number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(ii)')+1)');
        LHS_length = LHS_length + length(num2str(ii));
    end
    % Set the name of the component, including its depdence on spatial
    % variables.
    LHS_name = [' x',Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = length('x') + LHS_length + 3;
    
    % If the state component is one that has been added to impose the
    % higher order temporal derivative, indicate to which component the
    % added state corresponds.
    if (numel(PDE.x)-nx_old)~=1
        add_s2 = 's';
    else
        add_s2 = '';
    end
    % Indicate that the remaining state components were all added
    % during initialization.
    if ii==nx_old+1
        fprintf(['\n','Added ',num2str(numel(PDE.x)-nx_old),' ',object_name,add_s2,': \n']);
    end
    % For the added state components, indicate of which state component
    % they are the temporal derivative.
    Rcomp_idx = cell2mat(sub_num(str2num(num2str(PDE.x{ii}.term{1}.x)')+1)');
    xdot = '\x1E8B';    % UNICODE \dot{x}
    %RHS = [' := ',xdot,Rcomp_idx,'(',varnames_ii_t,')'];
    RHS = '';
    
    % Determine the order of differentiability in each spatial variable for
    % the state component
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
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

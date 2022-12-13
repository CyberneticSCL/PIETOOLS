function [PDE,del_state_tab] = expand_delays(PDE,suppress_summary)
% PDE = expand_delays(PDE) takes a pde_struct and retruns a PDE with no
% time delays, using new state variables to represent the delayed terms.
%
% INPUTS:
% PDE:  A pde_struct class object, defining a PDE in a manner as
%       outlined in "@pde_struct/initialize".
% suppress_summary: logical value indicating whether to suppress the 
%                   output summary. Defaults to false.
% 
% OUTPUTS:
% PDE:  A pde_struct object describing the same system as the input
%       PDE, but with no temporal delay in any of the terms. For a delayed
%       object (d/dr)^{D} v(t-s,r1=a,r2) where v=x_k, w_k or u_k, a new 
%       spatial variable s is added to the PDE, as well as a new state x_n
%       such that x_n(t,s,r2) = (d/dr)^{D} v(t-s,r1=a,r2). This is done by
%       adding a transport equation \dot{x}_n(t,s,r2) = (d/ds) x(t,s,r_2),
%       and a BC x_n(t,0,r2) = (d/dr)^{D} v(t,r1=a,r2).
% del_state_tab:    A q x (2+2*nvars_old) table of integers, matching
%                   each of the q newly added state components with a
%                   unique combination of state/input, delay variable,
%                   spatial position, and derivative. In particular:
%                   The first are integer values indicating for each of the
%                       1:q new states which component obj_k it corresponds
%                       to, where obj_k = x_k if k<nx_old, 
%                       obj_k = w_{k-nx_old} if k<nx_old + nw,
%                       and obj_k = u_{k-nx_old-nw} else,
%                       if del_state_tab(j,1) = k.
%                   The second column are integer values indicating which
%                       of the spatial variables is used to represent the
%                       delay, so that x_j(t,s_l) = obj_k(t-s_l) if
%                       del_state_tab(j,2) = l.
%                   Columns l=3:2+nvars_old indicate for each of the
%                       spatial variables s_l in the original PDE whether
%                       the state obj_k is evaluated at s_l=s_l (0),
%                       s_l=a_l (1), or s_l=b_l (0).
%                   Columns l=3+nvars_old:end indicate for each of the
%                       spatial variables s_l in the original PDE to what
%                       order the state obj_k is differentiated with
%                       respect to this spatial variable.
%
% NOTES:
% No boundary conditions are enforced upon the newly added state
% components. For the ODE state components, on which no BCs are imposed
% anyway, this is no problem, but for PDE state components with delay, this 
% may introduce conservatism in e.g. stability tests. An updated version
% adding BCs may be incorporated in a later version of PIETOOLS.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 10/13/2022

if nargin==1
    suppress_summary = false;
end

% % % Initialize the PDE, and extract some necessary parameters
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
if ~PDE.has_delay
    fprintf(['\n','No delayed states or inputs were encountered.'])
    return
end

% Extract number of components of each type.
nx = numel(PDE.x);
nw = numel(PDE.w);
nu = numel(PDE.u);
nz = numel(PDE.z);
ny = numel(PDE.y);
nBC = numel(PDE.BC);

% Determine number of spatial variables and delay variables.
nvars = size(PDE.vars,1);
ndelays = size(PDE.tau,1);


% % % Next: include the delay variables as spatial variables in the PDE
% Assign a dummy variable to each delay variable.
tau_vars = [PDE.tau(:,1),zeros(size(PDE.tau,1),1)];
for vv=1:size(tau_vars,1)
    tau_varname_vv = tau_vars(vv,1).varname{1};
    tau_varname_vv2 = [tau_varname_vv,'_dum'];
    tau_vars(vv,2) = polynomial({tau_varname_vv2});
end
% Delay variables will exist on domain [-max_delay,0].
tau_dom = [-abs(double(PDE.tau(:,2))),zeros(size(PDE.tau,1),1)];

% Add the delay variables to the list of spatial variables.
PDE.vars = [PDE.vars; tau_vars];
PDE.dom = [PDE.dom; tau_dom];

% Expand the tables to include the delay variables
PDE.x_tab = [PDE.x_tab(:,1:2+nvars),zeros(nx,ndelays),PDE.x_tab(:,3+nvars:2+2*nvars),zeros(nx,ndelays)];
PDE.w_tab = [PDE.w_tab(:,1:2+nvars),zeros(nw,ndelays)];
PDE.u_tab = [PDE.u_tab(:,1:2+nvars),zeros(nu,ndelays)];
PDE.z_tab = [PDE.z_tab(:,1:2+nvars),zeros(nz,ndelays)];
PDE.y_tab = [PDE.y_tab(:,1:2+nvars),zeros(ny,ndelays)];
PDE.BC_tab = [PDE.BC_tab(:,1:2+nvars),zeros(nBC,ndelays),PDE.BC_tab(:,3+nvars:2+2*nvars),zeros(nBC,ndelays)];


% % % Loop over the terms in the PDE, output equations, and BC, replacing
% % % any instance of a delayed state (d/dr)^{D} x(t-s,r), input w(t-s,r),
% % % and input u(t-s,r), with a new unique state variable x_n(t,s,r).

% % We keep track of which new state variable x_n(t,s,r) is linked to which
% % derivative of which object x, w or r, and wrt which delay s, using the
% % array "del_state_tab"
% % First column provides an index j, uniquely identifying
% %     a state x_i, input w_{i-nx}, or input u_{i-nx-nw}
% % Second column provides an index j, uniquely identifying a delay var s_j
% % Next nvars columns indicate the position at which to evaluate state.
% % Remaining columns indicate the desired order of the
% %     derivative of x_i wrt each of the original spatial vars r.
del_state_tab = zeros(0,2+2*nvars);
[PDE,del_state_tab,is_PDE_delay_x] = expand_delays_terms(PDE,del_state_tab,'x',nx);
[PDE,del_state_tab,is_PDE_delay_z] = expand_delays_terms(PDE,del_state_tab,'z',nz);
[PDE,del_state_tab,is_PDE_delay_y] = expand_delays_terms(PDE,del_state_tab,'y',ny);
[PDE,del_state_tab,is_PDE_delay_b] = expand_delays_terms(PDE,del_state_tab,'BC',nBC);


% % % Finally, clean up the structure and return.
% Remove the delay variables (they have all been replaced by spatial vars).
PDE.tau = zeros(0,2);
PDE.has_delay = false;
PDE = initialize(PDE,true);

% Display summary of the results.
if ~suppress_summary
    ndelay_states = size(del_state_tab,1);
    if ndelay_states==0
        fprintf(['\n','No delays were encountered.\n'])
    else
        print_expand_delay_summary(PDE,del_state_tab)
    end
end
if any([is_PDE_delay_x, is_PDE_delay_z, is_PDE_delay_y, is_PDE_delay_b])
    fprintf(2,['\n Warning: No BCs have been imposed on the newly added state components representing the delayed PDE states. \n',...
                 '          Results of e.g. stability tests may be very conservative.\n'])
end

% Normalize domain to reduce the number of spatial variables.
PDE = combine_vars(PDE);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [PDE,del_state_tab,is_PDE_delay] = expand_delays_terms(PDE,del_state_tab,Lobj,ncomps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over all the terms in the equation for "Lobj", replacing
% any instance of a delayed state (d/dr)^{D} x(t-s,r), input w(t-s,r),
% and input u(t-s,r), with a new unique state variable x_j(t,s,r).
%
% INPUTS:
% - PDE:    A "pde_struct" class object defining a PDE.
% - del_state_tab:  A n x (2+nvars) array matching each new state j=1,...n             
%                   with a component x_{k}, w_{k-nx}, u_{k-nx-nw},
%                   a new delay spatial variable s_l,
%                   indices indicating position a to evaluate x_{k}(r=a)
%                   and a desired order of differentiation (d/dr)^D x_{k},
%                   as del_state_tab(j,:) = [k, l, a, D].
% - Lobj:   A char 'x', 'y', 'z' or 'BC', indicating which equations to
%           loop over.
% - ncomps: Original number of components of type Lobj (must be specified,
%           as the number of BCs and state components will increase).
% 
% OUTPUTS:
% - PDE:    A "pde_struct" representing the same PDE as input, but with all
%           delays in the equations of "Lobj" replaced by new states.
% - del_state_tab:  Updated array of size (n+m) x (2+nvars), with 
%                   information on the m newly added state components.
% - is_PDE_delay:   Logical value indicating whether delay occurs in any
%                   infinite-dimensional state component.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvars_old = (size(del_state_tab,2) - 2)/2;  % Original number of spatial vars.
ndelays = size(PDE.vars,1) - nvars_old;     % Number of spatial vars representing delays
nvars = size(PDE.vars,1);                   % Total number of spatial vars
nx = numel(PDE.x) - size(del_state_tab,1);  % Original number of state components

is_PDE_delay = false;   % Keep track of whether any infinite-dimensional state is delayed  

% Extract variable names
tau_name = cell(ndelays,1);
for vv=1:ndelays
    % varnames are extracted 1-by-1, as "polynomial" stores them
    % alphabetically.
    tau_name{vv} = PDE.vars(nvars_old+vv,1).varname{1};
end

% % Loop over all equations for the considered obj, replacing delays with 
% % new states.
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
        % We assign an index to each RHS component, depending on whether it
        % is a state, exogenous input, or actuator input, and which index
        % it already has.
        loc_Robj = zeros(1,nvars);
        diff_Robj = zeros(1,nvars);
        if isfield(term_jj,'x')
            Robj = 'x';
            Rindx = term_jj.x;
            state_idx = Rindx;

            % Establish what derivative is taken of the state component.
            has_vars_Robj = logical(PDE.x_tab(Rindx,3:2+nvars));
            is_PDE_delay = is_PDE_delay || any(has_vars_Robj);
            diff_Robj(has_vars_Robj) = term_jj.D;   
            
            % Use index 0, 1, 2 to indicate if x is evaluated at interior,
            % lower boundary, or upper boundary for each spatial variable.
            dom_Robj = PDE.dom(has_vars_Robj,:);
            loc_Robj_partial = zeros(1,sum(has_vars_Robj));
            for kk=1:numel(loc_Robj_partial)
                try loc_kk = double(term_jj.loc(kk));
                    loc_Robj_partial(kk) = find(dom_Robj(kk,:)==loc_kk,1,'first');
                catch
                    continue
                end
            end
            loc_Robj(has_vars_Robj) = loc_Robj_partial;
        elseif isfield(term_jj,'w')
            Robj = 'w';
            Rindx = term_jj.w;
            state_idx = Rindx + nx; % w inputs are assigned indices nx+1 through nx+nw
        else
            Robj = 'u';
            Rindx = term_jj.u;
            state_idx = Rindx + nx + numel(PDE.w);  % u inputs are assigned indices nx+nw+1 through nx+nw+nu
        end
       
        % Distinguish case where delay is a fixed real value, and where
        % delay is distributed.
        if isa(delay,'double') || isdouble(delay)
            % % % The delay is fixed.
            delay = abs(double(delay));
            % Check if the delay is already included in our list.
            tau_idx = find(delay==abs(double(PDE.tau(:,2))),1,'first');
            if ~isempty(tau_idx)
                % The considered delay already appears in our table.
                % We still have to check if this delay already appears in
                % combination with our particular derivative of the RHS 
                % object Robj, d/dr Robj.
                tau_idx = tau_idx + nvars_old;
                state_tau_idx = find(all(del_state_tab==[state_idx,tau_idx,loc_Robj(1:nvars_old),diff_Robj(1:nvars_old)],2));
                if ~isempty(state_tau_idx)
                    % The combination of Robj and delay is already present
                    % in state_dep_tab
                    % --> Adjust the term to replace delay by new state:
                    %   d/dr Robj(t-s,r1=a,r2) = xnew(t,s,r2)
                    new_state_idx = state_tau_idx + nx;
                    has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                    term_jj = rmfield(term_jj,Robj);
                    term_jj.x = new_state_idx;
                    term_jj.D = zeros(1,sum(has_vars_delay_state));
                    new_loc = PDE.vars(:,1)';
                    new_loc(tau_idx) = -delay;  % Evaluate at s=-delay
                    term_jj.loc = new_loc(has_vars_delay_state);
                    term_jj.delay = 0;
                else
                    % The combination of Robj and delay does not occur yet,
                    % but the delay is already included as one of the
                    % variables.
                    % Add a new state xnew(t,s) = d/dr Robj(t-s);
                    PDE = add_delay_state(PDE,Robj,Rindx,tau_idx,loc_Robj,diff_Robj);
                    nvars = size(PDE.vars,1);
                    new_state_idx = size(PDE.x_tab,1);
                    del_state_tab = [del_state_tab; [state_idx,tau_idx,loc_Robj(1:nvars_old),diff_Robj(1:nvars_old)]];

                    % Adjust the term to replace delay by new state var:
                    %   d/dr Robj(t-s,r1=a,r2) = xnew(t,s,r2)
                    has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                    term_jj = rmfield(term_jj,Robj);
                    term_jj.x = new_state_idx;
                    term_jj.D = zeros(1,sum(has_vars_delay_state));
                    new_loc = PDE.vars(:,1)';
                    new_loc(tau_idx) = -delay;  % Evaluate at s=-delay
                    term_jj.loc = new_loc(has_vars_delay_state);
                    term_jj.delay = 0;
                end
            else
                % The delay has not been included in our delay table yet
                % --> we have to add it.
                % % % Add the delay to the PDE as a new spatial variable.
                % Define new spatial variables (s,th) to represent the delay.
                tau_idx = size(PDE.vars,1) + 1;
                var1 = polynomial({['ntau_',num2str(tau_idx)]});
                var2 = polynomial({['tau_dum',num2str(tau_idx)]});
                PDE.vars = [PDE.vars; [var1,var2]];
                PDE.dom = [PDE.dom; [-delay,0]];
                
                % Add the new spatial variable to all the tables.
                PDE.x_tab = [PDE.x_tab(:,1:2+nvars),zeros(size(PDE.x_tab,1),1),PDE.x_tab(:,3+nvars:2+2*nvars),zeros(size(PDE.x_tab,1),1)];
                PDE.w_tab = [PDE.w_tab,zeros(size(PDE.w_tab,1),1)];
                PDE.u_tab = [PDE.u_tab,zeros(size(PDE.u_tab,1),1)];
                PDE.z_tab = [PDE.z_tab,zeros(size(PDE.z_tab,1),1)];
                PDE.y_tab = [PDE.y_tab,zeros(size(PDE.y_tab,1),1)];
                PDE.BC_tab = [PDE.BC_tab(:,1:2+nvars),zeros(size(PDE.BC_tab,1),1),PDE.BC_tab(:,3+nvars:2+2*nvars),zeros(size(PDE.BC_tab,1),1)];
                
                % Add the combination of delay, Robj, and diff to table.
                del_state_tab = [del_state_tab; [state_idx,tau_idx,loc_Robj(1:nvars_old),diff_Robj(1:nvars_old)]];
                loc_Robj = [loc_Robj,0];
                diff_Robj = [diff_Robj,0];


                % % % Add a new state xnew(t,s) = Robj(t-s);
                PDE = add_delay_state(PDE,Robj,Rindx,tau_idx,loc_Robj,diff_Robj);
                nvars = size(PDE.vars,1);
                new_state_idx = size(PDE.x_tab,1);

                % % % Adjust the term to replace delay by new state var:
                %   d/dr Robj(t-s,r1=a,r2) = xnew(t,s,r2)
                has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                term_jj = rmfield(term_jj,Robj);
                term_jj.x = new_state_idx;
                term_jj.D = zeros(1,sum(has_vars_delay_state));
                new_loc = PDE.vars(:,1)';
                new_loc(tau_idx) = -delay;  % Evaluate at s=-delay
                term_jj.loc = new_loc(has_vars_delay_state);
                term_jj.delay = 0;
            end
        else
            % % % The delay is distributed: int_{-delay}^{0} d/dr Robj(t-s,r) ds.
            tau_idx = find(ismember(tau_name,delay.varname{1}),1,'first');
            if isempty(tau_idx)
                error('The delay variable does not appear in the list PDE.tau...')
            end
            delay = double(PDE.tau(tau_idx,2));
            tau_idx = tau_idx + nvars_old;
            state_tau_idx = find(all(del_state_tab==[state_idx,tau_idx,loc_Robj(1:nvars_old),diff_Robj(1:nvars_old)],2),1,'first');
            if ~isempty(state_tau_idx)
                % The combination of Robj and delay is already present in
                % state_dep_tab
                % --> Adjust the term to replace delay by new state:
                %   int_{-delay}^{0} d/dr Robj(t-s,r) ds = int_{-delay}^{0} xnew(t,s,r) ds
                new_state_idx = nx + state_tau_idx;
                
                % Determine which spatial variables Robj and the delayed
                % state depend on.
                has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                has_vars_Robj = has_vars_delay_state;
                has_vars_Robj(tau_idx) = false;
                
                % Define the new integral.
                new_int = cell(nvars,1);
                new_int{tau_idx} = [-delay,0];
                new_int(has_vars_Robj) = term_jj.I;
                
                % Set the new term: 
                % int_{-delay}^{0} d/dr Robj(t-s) ds = int_{-delay}^{0} xnew(t,s) ds
                term_jj = rmfield(term_jj,Robj);
                term_jj.x = new_state_idx;
                term_jj.D = zeros(1,sum(has_vars_delay_state));
                term_jj.loc = PDE.vars(has_vars_delay_state,1)';
                term_jj.I = new_int(has_vars_delay_state);
                term_jj.delay = 0;
            else
                % The combination of Robj and delay does not occur yet
                % --> Add the combination to the table and the PDE
                del_state_tab = [del_state_tab; [state_idx,tau_idx,loc_Robj(1:nvars_old),diff_Robj(1:nvars_old)]];

                % Add a new state xnew(t,s) = Robj(t-s);
                PDE = add_delay_state(PDE,Robj,Rindx,tau_idx,loc_Robj,diff_Robj);
                nvars = size(PDE.vars,1);
                new_state_idx = size(PDE.x_tab,1);
                
                % Determine which spatial variables Robj and the delayed
                % state depend on.
                has_vars_delay_state = logical(PDE.x_tab(new_state_idx,3:2+nvars));
                has_vars_Robj = has_vars_delay_state;
                has_vars_Robj(tau_idx) = false;
                
                % Define the new integral.
                new_int = cell(nvars,1);
                new_int{tau_idx} = [-delay,0];
                new_int(has_vars_Robj) = term_jj.I;
                
                % Set the new term: 
                % int_{-delay}^{0} d/dr Robj(t-s,r) ds = int_{-delay}^{0} xnew(t,s,r) ds
                term_jj = rmfield(term_jj,Robj);
                term_jj.x = new_state_idx;
                term_jj.D = zeros(1,sum(has_vars_delay_state));
                term_jj.loc = PDE.vars(has_vars_delay_state,1)';
                term_jj.I = new_int(has_vars_delay_state);
                term_jj.delay = 0;
            end
        end
        PDE.(Lobj){ii}.term{jj} = term_jj;
    end
end


end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function PDE = add_delay_state(PDE,Robj,Rindx,tau_var_idx,loc_Robj,diff_Robj)
% % % This subroutine adds a state x(t,s,r2) to the PDE to represent a
% % % delayed object, x(t,s,r2) = d/dr^{diff} Robj(t-s,r1=a,r2)

% % % Extract the variables, domain, and info on the object Robj.
full_vars = PDE.vars;
full_dom = PDE.dom;
nvars = size(full_vars,1);
Robj_tab = [PDE.([Robj,'_tab'])(Rindx,1:2+nvars), zeros(1,nvars)];
Robj_size = Robj_tab(1,2);


% % % Add a new state component to x_tab.
% The new state should depend on the same variables as Robj.
PDE.x_tab = [PDE.x_tab; Robj_tab];
PDE.x_tab(end,1) = size(PDE.x_tab,1);
has_vars_Robj = logical(PDE.x_tab(end,3:2+nvars));
% The new state should also depend on the new delay variable.
PDE.x_tab(end,[2+tau_var_idx,2+nvars+tau_var_idx]) = 1;
% But the new state should not depend on variables wrt which Robj is
% evalauted at a boundary
PDE.x_tab(end,3:2+nvars) = double(PDE.x_tab(end,3:2+nvars) & ~loc_Robj);
has_vars_delay_state = logical(PDE.x_tab(end,3:2+nvars));


% % % Add the new state xnew(t,s) = d/dr Robj(t-s).
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


% % % Add the boundary condition xnew(t,s=0) = d/dr Robj(t)
% Initialize the BC
PDE.BC = [PDE.BC; struct()];
PDE.BC{end}.size = Robj_size;
PDE.BC{end}.vars = full_vars(has_vars_Robj,:);
PDE.BC{end}.dom = full_dom(has_vars_Robj,:);
PDE.BC_tab = [PDE.BC_tab;
              [size(PDE.BC_tab,1)+1,Robj_size,double(has_vars_Robj),zeros(1,nvars)]];
PDE.BC_tab(2+nvars+tau_var_idx) = 1;    % BC involves first order diff wrt delay variable.

% Set up the BC:   (d/dr)^{diff} Robj(t,r1=a,r2) = xnew(t,s=0,r2)
% First term: (d/dr)^{diff} Robj(t,r1=a,r2)
if strcmp(Robj,'x')
    PDE.BC{end}.term{1}.x = Rindx;
    PDE.BC{end}.term{1}.D = diff_Robj(has_vars_Robj);
    new_loc = [full_vars(:,1),full_dom]';
    new_loc = new_loc((loc_Robj+1) + 3*(0:nvars-1));
    PDE.BC{end}.term{1}.loc = new_loc(has_vars_Robj);
elseif strcmp(Robj,'w')
    PDE.BC{end}.term{1}.w = Rindx;
else
    PDE.BC{end}.term{1}.u = Rindx;
end
PDE.BC{end}.term{1}.C = eye(Robj_size);
PDE.BC{end}.term{1}.I = cell(sum(has_vars_Robj),1);
PDE.BC{end}.term{1}.delay = 0;

% Second term: -xnew(t,s=0,r2) 
PDE.BC{end}.term{2}.x = numel(PDE.x);
PDE.BC{end}.term{2}.D = zeros(1,sum(has_vars_delay_state));
PDE.BC{end}.term{2}.loc = full_vars(:,1)';
PDE.BC{end}.term{2}.loc(tau_var_idx) = 0;
PDE.BC{end}.term{2}.loc = PDE.BC{end}.term{2}.loc(has_vars_delay_state);
PDE.BC{end}.term{2}.C = -eye(Robj_size);
PDE.BC{end}.term{2}.I = cell(sum(has_vars_delay_state),1);
PDE.BC{end}.term{2}.delay = 0;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_expand_delay_summary(PDE,del_state_tab)
% print_expand_delay_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning how many
% state components have been added to the system, and which delayed
% states/inputs these represent.
%
% INPUTS:
% - PDE:        A "struct" or "pde_struct" class object defining a PDE.
% - del_state_tab:  A q x (2+2*nvars_old) array linking each of the q new
%                   state components to a particular combination of
%                   state or input, indicated by first column;
%                   delay variable, indicated by second column;
%                   boundary at which to evaluate state, indicated by
%                   columns 3:2+nvars_old;
%                   derivative wrt each of the original spatial vars,
%                   indicated by columns 3+nvars_old:2+2*nvars_old.
%
% OUTPUTS:
% Displays information in the command window concerning how many
% state components have been added to the system, and which delayed
% states/inputs these represent.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx_old = numel(PDE.x) - size(del_state_tab,1);  % Number of state components in original PDE.
nx_new = size(del_state_tab,1);                 % number of added state components
nvars_old = (size(del_state_tab,2)-2)/2;        % Number of spatial variables in original PDE.
nvars = size(PDE.vars,1);                       % total number of spatial variables in new PDE.

% Check if any state components have actually been added.
if nx_new==0
    fprintf(['\n','No delays were encountered.\n']);
    return
elseif nx_new==1
    fprintf(['\n','Added ',num2str(numel(PDE.x)-nx_old),' state components: \n']);
elseif nx_new>1
    fprintf(['\n','Added ',num2str(numel(PDE.x)-nx_old),' state components: \n']);
end


% % % Set up a library of char objects to use in the display
% UNICODE superscript digits.
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
% UNICODE subscript indices.
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
% Symbol for partial integration.
partial = '\x2202';
% Determine how many subscripts will be needed.
n_digits = 1+floor(log(numel(PDE.x))/log(10));


% % For the purpose of clearer display, we make an estimate of how long the
% % display will be. 
% First, determine which variables appear in the different state variables.
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
% Then, estimate the size of the string of characters denoting the components
% "PDE.x{ii}".
nvars_max = max(sum(PDE.x_tab(:,3:2+nvars),2));
lngth_varnames_mean = ceil(length(global_varnames)*nvars_max/nvars);
LHS_length_max = 1+1 + n_digits + lngth_varnames_mean+2; % e.g. ' x13(t,s1,s2,s3)', 


% % Loop over the new components, displaying the variables it depends on,
% % and which component in the original PDE it represents.
for ii=1:nx_new

    % % % Build a str to represent the new delay variable x_n(t,s,r2)
    Lidx = ii+nx_old;
    comp_ii = PDE.x{Lidx};
    
    % Establish the names of the variables on which the component depends.
    % Set a comma separated list.
    varnames_ii = comp_ii.vars(1,1).varname{1};
    for kk=2:size(comp_ii.vars,1)
        varnames_ii = [varnames_ii,',',comp_ii.vars(kk,1).varname{1}];
    end
    varnames_ii_t = ['t,',varnames_ii]; % All components may vary in time
    
    % Establish the (subscript) index for the component.
    LHS_length = length(varnames_ii_t);
    if numel(PDE.x)==1
        % There is only one state component --> no need to give index
        Lcomp_idx = '';
    elseif numel(PDE.x)<=9
        % The state number consists of a single decimal.
        Lcomp_idx = sub_num{Lidx+1};
        LHS_length = LHS_length + 1;
    else
        % The state number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(Lidx)')+1)');
        LHS_length = LHS_length + length(num2str(Lidx));
    end
    % Set the name of the component, including its depdence on spatial
    % variables.
    LHS_name = [' x',Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = 1 + LHS_length + 3;
    
    
    % % % Build a str to represent the component which the new state
    % % % represents: x_n(t,s,r2) = (d/dr)^D x_o(t-s,r1=a,r2)
    % The component can be a state x, or input w or u.
    Robj_full_idx = del_state_tab(ii,1);
    if Robj_full_idx <= nx_old
        Robj = 'x';
        Ridx = Robj_full_idx;
        has_vars_Rcomp = logical(PDE.x_tab(Ridx,3:2+nvars));
    elseif Robj_full_idx <= nx_old+numel(PDE.w)
        Robj = 'w';
        Ridx = Robj_full_idx - nx_old;
        has_vars_Rcomp = logical(PDE.w_tab(Ridx,3:2+nvars));
    else
        Robj = 'u';
        Ridx = Robj_full_idx - nx_old - numel(PDE.w);
        has_vars_Rcomp = logical(PDE.u_tab(Ridx,3:2+nvars));
    end
    Ridx_str = cell2mat(sub_num(str2num(num2str(Ridx)')+1)');
    Robj = [Robj,Ridx_str];

    % Establish on which variables the RHS object depends
    tau_idx = del_state_tab(ii,2);
    Rvar1_str = ['(t-',PDE.vars(tau_idx,1).varname{1}]; % Display delay
    Rvar1 = PDE.vars(has_vars_Rcomp,1);
    Rvar1_list = cell(size(Rvar1));
    for kk=1:numel(Rvar1_list)
        Rvar1_list{kk} = Rvar1(kk).varname{1};
    end

    % The new state L may correspond to the state R evaluated at a boundary
    % --> display this evaluation.
    Robj_loc = del_state_tab(ii,3:2+nvars_old);
    for kk=1:nvars_old
        if Robj_loc(kk)==0
            continue
        else
            Rvar1_list{kk} = num2str(PDE.dom(kk,Robj_loc(kk)));
        end
    end
    for idx=1:numel(Rvar1_list)
        Rvar1_str = [Rvar1_str,',',Rvar1_list{idx}];
    end
    Rvar1_str = [Rvar1_str,')'];

    % The new state may represent a derivative of the original one
    % --> include this derivative in the expression
    Robj_diff = del_state_tab(ii,3+nvars_old:2+2*nvars_old);
    Robj_diff = Robj_diff(has_vars_Rcomp(1:nvars_old));
    D_trm = '';
    for kk=1:numel(Robj_diff)
        D_kk = Robj_diff(kk);
        if D_kk==0
            % Don't display partial symbol to 0th power
            continue
        elseif D_kk==1
            % Don't display superscript 1 in derivative
            D_trm = [D_trm,'(',partial,'/',partial,Rvar1(kk).varname{1},') '];
        elseif D_kk<=9
            D_trm = [D_trm,'(',partial,sup_num{D_kk+1},'/',partial,Rvar1(kk).varname{1},sup_num{D_kk+1},') '];
        else
            D_sup = cell2mat(sup_num(str2num(num2str(D_kk)')+1)');
            D_trm = [D_trm,'(',partial,D_sup,'/',partial,Rvar1(kk).varname{1},D_sup,') '];
        end
    end

    % Define the RHS object
    RHS = [' := ',D_trm,Robj,Rvar1_str];

    % % % Finally, display
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,';\n'])
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
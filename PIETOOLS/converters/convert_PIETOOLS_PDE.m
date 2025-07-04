function PIE = convert_PIETOOLS_PDE(PDE,comp_order,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert a coupled ODE-1D_PDE-2D system to an equivalent PIE.
% 
% INPUTS:
% - PDE:    A struct or pde_struct object defining a PDE in the terms
%           format (see also "@pde_struct/initialize").
% - comp_order:     Optional argument specifying the order of the state
%                   variables, inputs, and outputs from the PDE as they 
%                   appear in the PIE, specified as a struct. That is, e.g.
%                   PIE.x{i} = PDE.x{comp_order.x(1)};
%                   If not specified, the converter will reorder the
%                   variables itself, and determine the new order.
% - flag:           (optional) 'char' object which can be set to 'silent'
%                   to suppress display window messages informing the user
%                   on the progress, or 'Top' to stop conversion after the
%                   T operators (from PIE to PDE state) have been computed,
%                   returning only those.
%
% OUTPUTS:
% - PIE:    A struct with fields:
%
%           Tx, Tw, Tu: opvar or opvar2d objects describing the map from
%           fundamental state x_f, exogenous input w, and actuator input u 
%           to the PDE state x.
%
%           A, Bw, Bu: opvar or opvar objects describing the PIE dynamics
%           for fundamental state x_f.
%
%           Cz, Dzw, Dzu: opvar or opvar2d objects describing the PIE for
%           the regulated output z.
%
%           Cy, Dyw, Dyu: opvar or opvar2d objects describing the PIE for
%           the observed output y.
%
%           dim, dom, vars,: The dimension of the spatial domain for the
%           PIE, the domain of the spatial variables, and the spatial
%           variables that appear in the PIE. dom should be a dimx2 array,
%           with each row defining the interval associated to a spatial
%           variable. vars should be a dimx2 pvar object, with the first
%           column defining the primary variables, and the second the
%           associated dummy variables as used in the opvar objects.
%
%           x_tab, y_tab, z_tab, u_tab, w_tab: Arrays describing for each
%           fundamental state, input, and output component, which PDE
%           state, input, or output component it corresponds to, what size
%           this component is, and which of the dim variables it depends
%           on.
%
% NOTES:
% - For a state variables x(s1,s2), differentiable up to order i wrt s1 and
%   order j wrt s2, the associated fundamental state variable xf(s1,s2) is
%   defined as
%       xf(s1,s2) = \partial_{s1}^{i} \partial_{s2}^{j} x(s1,s2);
%   For any well-posed PDE, with well-posed BCs of the form
%       0 = F(x,u,w);
%   there exist unique PI operator T, Tu and Tw such that
%       x(s1,s2) = (T * x)(s1,s2) + (Tu * u)(s1,s2) + (Tw * w)(s1,s2);
%   This function constructs these operators {T, Tx Tw}. Using these 
%   operators, it then constructs the PI operators A through Dzw such that 
%   the PDE may be equivalently represented by the PIE
%       Tw*w(t) + Tu*u + T*xf(t) = A  * xf(t) + Bu  * u(t) + Bw  * w(t);
%                           y(t) = Cy * xf(t) + Dyu * u(t) + Dyw * w(t);
%                           z(t) = Cz * xf(t) + Dzu * u(t) + Dzw * w(t);
%
% - The order of the state, input, and output components in the PIE
%   structure may not be the same as that in the PDE structure. Check e.g.
%   PIE.x_tab(:,1) to see how the state components of the PDE have been
%   re-arranged.
%
% - If the PDE is 2D and involves boundary conditions with inputs u or
% disturbances w that vary in space, the inputs and disturbance may be
% differentiated with respect to the spatial variables as well in the PIE
% representation. The function will warn the user of this fact, and the
% order of the derivative taken with respect to each spatial is stored in
% the outputs PIE.u_tab and PIE.w_tab, with e.g. the element
% PIE.u_tab(i,5) indicating the order of derivative taken of the ith input
% with respect to the first spatial variable, and
% PIE.u_tab(i,6) indicating the order of derivative taken of the ith input
% with respect to the second spatial variable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - convert_PIETOOLS_PDE
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/09/2022
% Add 1D conversion, DJ - 10/16/2024
% DJ, 12/16/2024: Bugfix for batch and terms_legacy conversion;
% DJ, 01/12/2025: Add check if sufficient boundary conditions are specified
%                   in 1D case;
% DJ, 03/24/2025: Add check if too many boundary conditions are specified
%                   in 1D case;
% DJ, 06/01/2025: Pre-process PDE to get rid of higher-order temporal 
%                   derivatives and delays;
% DJ, 06/18/2025: Expand higher-order temporal derivatives only after
%                   conversion. Also, display summary of how PIE variables
%                   relate to PDE variables, and add an optional argument
%                   to suppress this summary.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Initialization                                          % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pass PDEs in legacy format to legacy converters.
if isfield(PDE,'n0') || isfield(PDE,'n1') || isfield(PDE,'n2')
    PIE = convert_PIETOOLS_PDE_batch(PDE);                                  % DJ, 12/16/2024
    return
elseif isfield(PDE,'n') || isfield(PDE,'ODE')  || isfield(PDE,'PDE')
    PIE = convert_PIETOOLS_PDE_terms_legacy(PDE);                           % DJ, 12/16/2024
    return
elseif (isa(PDE,'struct') && isfield(PDE,'x'))
    % Convert 'struct' to 'pde_struct'.
    PDE = pde_struct(PDE);
    PDE = initialize(PDE,true);
elseif isa(PDE,'sys')
    % Extract PDE terms from 'sys' structure.
    PDE = PDE.params;
elseif ~isa(PDE,'pde_struct')
    error('The input PDE is not appropriately specified. Please define your PDE as a "pde_struct" class object, and consult the manual and examples for illustration of the structure.')
end

% % Check the optional input arguments                                      % DJ, 06/18/2025
suppress_summary = false;       % set to true to suppress display of conversion summary
return_Top = false;             % set to true to ONLY return the fundamental to PDE state map
if nargin==1
    comp_order = [];
    flag = {};
elseif nargin==2
    % Allow flag to be specified as second argument as well.
    if ischar(comp_order)
        flag = {comp_order};
        comp_order = [];
    elseif iscellstr(comp_order)
        flag = comp_order;
        comp_order = [];
    else
        flag = {};
    end
elseif nargin>3
    error("Too many input arguments.")
end 
if ismember('silent',flag)
    suppress_summary = true;
end
if ismember('Top',flag) || ismember('return_Top',flag)
    return_Top = true;
end

% Coefficients smaller than tol in PI operators will be discarded.
op_clean_tol = 1e-8;   

% Initialize PDE in case this has not been done.
PDE = pde_struct(PDE);
if ~PDE.is_initialized 
    PDE = initialize(PDE,true);
end

% % % Get rid of any higher-order temporal derivatives and delays in the
% % % system, and try to reduce the number of spatial variables.
% if PDE.has_hotd                                                             % DJ, 06/18/2025
%     fprintf(['\n','Higher-order temporal derivatives were encounterd:\n']);
%     fprintf([' --- Running "expand_tderivatives" to reduce to first-order temporal derivatives ---\n']);
%     PDE = expand_tderivatives(PDE);
% end
if PDE.has_delay
    fprintf(['\n','Delayed states on inputs were encounterd:\n']);
    fprintf([' --- Running "expand_delays" to remove the delays --- \n']);
    PDE = expand_delays(PDE);
end
if PDE.dim>2
    fprintf(['\n','The PDE involves more than 2 spatial variables:\n']);
    fprintf([' --- Attempting to reduce the dimensionality using "combine_vars" ---\n']);
    PDE = combine_vars(PDE);
    if PDE.dim>2
        error('The number of spatial variables could not be reduced to 2 or fewer; conversion to PIE is not currently supported.')
    end
end

% Extract spatial variables, and their domain (these will be used in
% multiple subroutines).
global dom;         dom = PDE.dom;
global vars;        vars = PDE.vars;
global nvars;       nvars = size(vars,1);

% Since opvar2d objects are constructed to strictly map states x that can
% be decomposed as
%       [x0]    [ R^n0          ]
%   x = [x1] in [ L2^n1 [s1]    ]
%       [x2]    [ L2^n2 [s2]    ]
%       [x3]    [ L2^n3 [s1,s2] ]
% we have to re-arrange the state components, inputs and outputs to also
% occur in this order. That is, if e.g.
%   x = [x_i; x_j; x_k] in [ L2^n_i[s2]; R^n_j; L2^n_k[s2] ], 
% we re-order the state components x_i, x_j and x_k such that 
%   x_new = [x_j; x_i; x_k] in [ R^n_j; L2^n_i[s2]; L2^n_k[s2] ];
% We do this using the "reorder_comps" subroutine.
if nargin==1 || ~isa(comp_order,'struct')
    if ~suppress_summary
        fprintf(['\n',' --- Reordering the state components to allow for representation as PIE ---\n']);
    end
    [PDE,xcomp_order] = reorder_comps(PDE,'x',true);
    [PDE,ycomp_order] = reorder_comps(PDE,'y',true);
    [PDE,zcomp_order] = reorder_comps(PDE,'z',true);
    [PDE,wcomp_order] = reorder_comps(PDE,'w',true);
    [PDE,ucomp_order] = reorder_comps(PDE,'u',true);
    
    comp_order = struct();
    comp_order.x = xcomp_order;
    comp_order.y = ycomp_order;
    comp_order.z = zcomp_order;
    comp_order.w = wcomp_order;
    comp_order.u = ucomp_order;
    comp_order.BC = (1:numel(PDE.BC))';
else
    xcomp_order = comp_order.x;
    ycomp_order = comp_order.y;
    zcomp_order = comp_order.z;
    wcomp_order = comp_order.w;
    ucomp_order = comp_order.u;
    if ~isfield(comp_order,'BC')
        comp_order.BC = (1:numel(PDE.BC))';
    end
end
% If the PDE is not 2D, we artificially augment.
is2D = nvars==2;
if nvars==0
    % Define a 1D domain with variables.
    pvar ss1 tt1
    vars = [ss1,tt1];
    dom = [0,1;0,1];
    nvars = 1;
    
    % Add a temporary state variable that exists on the 1D domain.
    ncomps = numel(PDE.x);
    PDE.x{ncomps+1}.dom = dom;
    PDE.x{ncomps+1}.vars = vars;
    PDE.x{ncomps+1}.term{1}.x = ncomps+1;
    % Initialize the augmented system, and get rid of the temporary state.
    PDE = initialize_PIETOOLS_PDE(PDE,true);
    PDE.x = PDE.x((1:ncomps)');
    PDE.x_tab = PDE.x_tab((1:ncomps)',:);

% elseif nvars==1
%     is2D = false;
%      % Define a 2D domain with variables.
%     xtra_var = polynomial({[vars(1).varname{1},'_2']; [vars(2).varname{1},'_2']});
%     vars = [vars; xtra_var'];
%     dom = [dom; [0,1]];
%     nvars = 2;
% 
%     % Add a temporary state variable that exists on the 2D domain.
%     ncomps = numel(PDE.x);
%     PDE.x{ncomps+1}.dom = dom;
%     PDE.x{ncomps+1}.vars = vars;
%     PDE.x{ncomps+1}.term{1}.x = ncomps+1;
%     % Initialize the augmented system, and get rid of the temporary state.
%     PDE = initialize_PIETOOLS_PDE(PDE,true);
%     PDE.x = PDE.x((1:ncomps)');
%     PDE.x_tab = PDE.x_tab((1:ncomps)',:);
%     PDE.is_initialized = true;  % No need to re-initialize
end

% Extract the tables providing information on the spatial dependence (and
% order of differentiability) of the state, input, and output components.
x_tab = PDE.x_tab;
u_tab = PDE.u_tab;
w_tab = PDE.w_tab;
y_tab = PDE.y_tab;
z_tab = PDE.z_tab;

% Extract dimensions of state components, inputs and outputs.
global np_op;               np_op = struct();
np_op.x = get_opdim(x_tab);
np_op.u = get_opdim(u_tab);
np_op.w = get_opdim(w_tab);
np_op.y = get_opdim(y_tab);
np_op.z = get_opdim(z_tab);


if ~suppress_summary
    fprintf('\n --- Converting the PDE to an equivalent PIE --- \n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Deriving a map from fundamental to PDE state            % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We first derive a map from the fundamental state xf to the PDE state x,
% such that
%    x = Tx*xf + Tu*u + Tw*w;
% where Tx, Tu and Tw are PI operators. This map is uniquely defined by the
% BCs (of the form 0 = F(x,u,w)), and we have two options for computing it:
%
% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% % % OPTION 1:
% % % We expand the state separately along each dimension as
% % %   x = Tx_1*D_1*x +Tw_1*w +Tu_1*u;
% % %   x = Tx_2*D_2*x +Tw_2*w +Tu_2*u;
% % % where D_1 and D_2 are diagonal differential operators taking the
% % % highest-order admissible derivative of each state component with
% % % respect to the first and second spatial variable, respectively.
% % % Then, we combine these expansions as either
% % %   D_2*x = Tx_1*D_2*D_1*x +Tw_1*D_2*w +Tu_1*D_2*u
% % %   --> x = Tx_2*Tx_1*D_12*x +Tx_2*Tw_1*D_2*w +Tw_2*w 
% % %               +Tx_2*Tu_1*D_2*u +Tu_2*u;
% % %   --> x = Tx*xf + Tw*D_w*w + Tu*D_u*u;
% % % or as
% % %   D_1*x = Tx_2*D_1*D_2*x +Tw_2*D_1*w +Tu_2*D_1*u
% % %   --> x = Tx_1*Tx_2*D_12*x +Tx_1*Tw_2*D_1*w +Tw_1*w 
% % %               +Tx_1*Tu_2*D_1*u +Tu_1*u;
% % %   --> x = Tx*xf + Tw*D_w*w + Tu*D_u*u;
% % % where xf = D_2*D_1*x.
% % % This approach only works for "separable" boundary conditions, so that
% % % we can separately enforce the boundary conditions along each spatial
% % % direction. Note that this approach may also require derivatives of
% % % the inputs to be taken.
if (is2D && (any(any(w_tab(:,3:2+nvars))) || any(any(u_tab(:,3:2+nvars)))))
    old_dom = dom;      old_vars = vars;        np_op_old = np_op;
    try [Top_x,Top_w,Top_u,Dvals_w,Dvals_u] = Compute_Tmap_PDE_2D_Separable(PDE,comp_order);
        use_Tmap_opt2 = false;
        Top = struct();
        Top.x = Top_x;      Top.u = Top_u;      Top.w = Top_w;
    catch
        use_Tmap_opt2 = true;
    end
    % Correct values of global variables.
    dom = old_dom;      vars = old_vars;        np_op = np_op_old;
    nvars = size(vars,1);
else
    use_Tmap_opt2 = true;
end


% % % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
if use_Tmap_opt2
% % % OPTION 2:
% % % We use a more "traditional" approach, which can handle more general
% % % boundary conditions, but is not supported for function-valued forcing
% % % at the boundaries. This approach works in three steps:
% % % 
% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % STEP 1: 
% % % Express the PDE state x in terms of the fundamental state xf and a
% % % set of "core boundary" states xb as
% % %   x = P_b2x * xb + P_f2x * xf,
% % % The values of the PI operators P_b2x and Pf2x are obtained using the 
% % % fundamental theorem of calculus, in the following function:
[Pop_f2x, Pop_b2x, bc_state_tab] = PIETOOLS_FTC_Expansion_2D(PDE);
np_op.b = get_opdim(bc_state_tab);
% Here, the bc_state_tab is a table similar to x_tab, providing for each
% component of the core boundary state what spatial variables it depends
% on, and its order of differentiability.


% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % STEP 2: 
% % % Express the BCs 0 = F(x,u,w) in terms of the fundamental state xf and
% % % core boundary state xb as
% % %   0 = P_b2BC * xb + P_f2BC * xf + P_u2BC * u + P_w2BC * w.
% % % To this end, we first expand the boundary conditions to be expressed
% % % in terms of xb and xf, without explicitly constructing PI operators
% % % P_b2BC through P_w2BC just yet.
PDE = expand_PIETOOLS_PDE_BCs(PDE);
np_op.BC = get_opdim(PDE.BC_tab);

% Then, we construct opvar2d objects to represent the operators P_b2BC 
% through P_w2BC such that 
%   0 = P_b2BC * xb + P_f2BC * xf + P_u2BC * u + P_w2BC * w,
% using a separate subroutine.
[Pop_b2BC, Pop_f2BC, Pop_u2BC, Pop_w2BC] = construct_PI_ops_BCs(PDE,bc_state_tab);


% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % STEP 3: 
% % % Invert the operator P_b2BC, to express xb in terms of xf, u and w as
% % %   xb = -(P_b2BC\P_f2BC)* xf - (P_b2BC\P_u2BC)* u - (P_b2BC\P_w2BC)* w
% % %      = P_f2b * xf + P_u2b * u + P_w2b * w
% % % It follows then that
% % %   x = P_b2x * xb + P_f2x * xf
% % %     = P_b2x * (P_f2b * xf + P_u2b * u + P_w2b * w) + P_f2x * xf
% % %     = (P_b2x*Pf2b + P_f2x) * xf + (P_b2x*P_u2b) * u + (P_b2x*P_w2b) * w
% % %     = Tx * xf + Tu * u + Tw * w;
% Start with Tx.
if Pop_f2BC==0
    % If boundary conditions are imposed only on lower boundaries, the core
    % boundary states must all be zero.
    % --> x = P_b2x * 0 + P_f2x * xf = P_f2x * xf
    Top_x = Pop_f2x;
else
    % Otherwise, compute Pop_f2b = - P_b2BC\P_f2BC;
    if nvars==1
        % In 1D case, boundary values are all finite-dimensional, so
        % Pop_b2BC should be just a matrix
        if ~all(Pop_b2BC.dim(2,:)==0)
            error('The boundary conditions appear to be infinite-dimensional, something is going wrong...')
        else
            % Test for well-posedness of boundary conditions
            if size(Pop_b2BC.P,1)<size(Pop_b2BC.P,2)                        % DJ, 01/12/2025
                error('It appears insufficient boundary conditions have been specified; conversion to PIE is not supported.')
            elseif size(Pop_b2BC.P,1)>size(Pop_b2BC.P,2)                    % DJ, 03/24/2025
                error('It appears too many boundary conditions have been specified; either declare a higher order of differentiability of the PDE state, or reduce the number of boundary conditions.')
            end    
            if abs(det(double(Pop_b2BC.P)))>eps
                Pop_b2BC_inv = pinv(double(Pop_b2BC.P));
                Pop_f2b = -Pop_b2BC_inv*Pop_f2BC;
            else
                error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                    ' Please make sure that your system is well-posed.'])
            end
        end
    else
        % Boundary conditions may involve infinite-dimensional parameters
        % --> inversion becomes more tedious...
        try Pop_f2b = -mldivide(Pop_b2BC,Pop_f2BC,1,1e-10,5);
        catch
            error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                    ' Please make sure that your system is well-posed.'])
        end
    end
    % Check that the inverse is sufficiently accurate.
    if eq(Pop_b2BC * Pop_f2b + Pop_f2BC,0,op_clean_tol)
        Top_x = Pop_b2x*Pop_f2b + Pop_f2x;
        Top_x = clean_opvar(Top_x,op_clean_tol);
    else
        error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                ' Please make sure that your system is well-posed.'])
    end
end
% Then Tu.
if ~any(np_op.u) || Pop_u2BC==0
    % Avoid computations if inputs are not present/do not contribute to BCs.
    if nvars<=1
        opvar Top_u;
        Top_u.dim = [np_op.x,np_op.u];  Top_u.I = dom;
        Top_u.var1 = vars(1,1);         Top_u.var2 = vars(1,2);
    else
        Top_u = opvar2d([],[np_op.x,np_op.u],dom,vars);
    end
else
    % Pop_f2b = - P_b2BC\P_u2BC;
    if nvars<=1
        Pop_u2b = -Pop_b2BC_inv*Pop_u2BC;
    else
        try Pop_u2b = -mldivide(Pop_b2BC,Pop_u2BC,1,1e-10,5);
        catch
            error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                    ' Please make sure that your system is well-posed.'])
        end
    end
    if eq(Pop_b2BC * Pop_u2b + Pop_u2BC,0,op_clean_tol)
        Top_u = Pop_b2x*Pop_u2b;
        Top_u = clean_opvar(Top_u,op_clean_tol);
    else
        error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                ' Please make sure that your system is well-posed.'])
    end
end
% Then Tw.
if ~any(np_op.w) || Pop_w2BC==0
    % Avoid computations if inputs are not present/do not contribute to BCs.
    if nvars<=1
        opvar Top_w;
        Top_w.dim = [np_op.x,np_op.w];  Top_w.I = dom;
        Top_w.var1 = vars(1,1);         Top_w.var2 = vars(1,2);
    else
        Top_w = opvar2d([],[np_op.x,np_op.w],dom,vars);
    end
else
    % Pop_w2b = - P_b2BC\P_w2BC;
    if nvars<=1
        Pop_w2b = -Pop_b2BC_inv*Pop_w2BC;
    else
        try Pop_w2b = -mldivide(Pop_b2BC,Pop_w2BC,1,1e-10,5);
        catch
            error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                    ' Please make sure that your system is well-posed.'])
        end
    end
    if eq(Pop_b2BC * Pop_w2b + Pop_w2BC,0,op_clean_tol)
        Top_w = Pop_b2x*Pop_w2b;
        Top_w = clean_opvar(Top_w,op_clean_tol);
    else
        error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                ' Please make sure that your system is well-posed.'])
    end
end

% % % With that, we have a map from the fundamental state space to the PDE
% % % state space as
% % %   x = Tx * xf + Tu * u + Tw * w;
Top = struct();
Top.x = Top_x;      Top.u = Top_u;      Top.w = Top_w;

% The obtained map does not require any derivative to be taken of the
% inputs.
Dvals_w = zeros(numel(PDE.w),2);
Dvals_u = zeros(numel(PDE.u),2);

end

if return_Top
    % If only the T operators are requested, we stop here.
    PIE = pie_struct();
    PIE.dom = dom;  PIE.vars = vars;
    PIE.T = Top_x;  PIE.Tw = Top_w;     PIE.Tu = Top_u;
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Deriving an equivalent PIE representation               % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now derive an equivalent PIE representation of the PDE, describing the
% dynamics of our fundamental state xf, which is free of the BCs and
% continuity constraints imposed upon the PDE state x. We do this in two
% steps:

% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % STEP 1: 
% % % Write the PDE using PI operators, so that
% % % (d/dt)^k x = sum_j=1^n1 (P_x2x{j} * Delta_j * D_j)* x + P_u2x* u + P_w2x* w; 
% % %          y = sum_j=1^n2 (P_x2y{j} * Delta_j * D_j)* x + P_u2y* u + P_w2y* w; 
% % %          z = sum_j=1^n3 (P_x2z{j} * Delta_j * D_j)* x + P_u2z* u + P_w2z* w; 
% % % where the P_ are all PI operators, D_j are differential operators,
% % % Delta_j are delta operators, evaluating the state x(s) at s=s, or
% % % s=a or s=b for s\in[a,b].

% We construct the different operators using a subroutine:
[Pop_x2x_cell, Pop_u2x, Pop_w2x, loc_diff_tab_x, retain_xvars_x] = construct_PI_ops_PDE(PDE,'x');
[Pop_x2y_cell, Pop_u2y, Pop_w2y, loc_diff_tab_y, retain_xvars_y] = construct_PI_ops_PDE(PDE,'y');
[Pop_x2z_cell, Pop_u2z, Pop_w2z, loc_diff_tab_z, retain_xvars_z] = construct_PI_ops_PDE(PDE,'z');
% Here, loc_diff_tab is a table in which each row has 2*nvars columns, the
% first nvars of which indicate whether the state is evaluated at the lower
% boundary (-1), interior (0), or upper boundary (1) for each spatial
% variable, and the remaining nvars columns indicate the order up to which
% the state component is differentiated wrt each variable.
% For each row j of the table, Pop_._cell{j} provides the associated
% opvar2d object.
% Each element of the cell retain_xvars provides indices of the state
% variables that can be differentiated up to the desired order.


% % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % STEP 2:
% % % Impose the relation
% % %   x = Tx * xf + Tu * u + Tw * w;
% % % so that e.g.
% % %  (d/dt)^k x = sum_j=1^n1 (P_x2x{j} * Delta_j * D_j)* x + P_u2x*u + P_w2x*w; 
% % %             = sum_j=1^n1 (P_x2x{j} * Delta_j * D_j * Tx)* xf
% % %                 + (P_u2x + sum_j=1^n1 (P_x2x{j} * Delta_j * D_j * Tu)* u
% % %                   + (P_w2x + sum_j=1^n1 (P_x2x{j} * Delta_j * D_j * Tw)* w
% % %             = A * xf + Bu * u + Bw * w;
% % % where the composition of the differential operators and the Delta
% % % operators with the PI operators T is a PI operator.
Pop = struct();
Pop.x = Pop_x2x_cell;   Pop.u = Pop_u2x;    Pop.w = Pop_w2x;
[Aop_x2x, Bop_u2x, Bop_w2x] = impose_fundamental_map(Top,Pop,loc_diff_tab_x,retain_xvars_x);
Pop.x = Pop_x2y_cell;   Pop.u = Pop_u2y;    Pop.w = Pop_w2y;
[Cop_x2y, Dop_u2y, Dop_w2y] = impose_fundamental_map(Top,Pop,loc_diff_tab_y,retain_xvars_y);
Pop.x = Pop_x2z_cell;   Pop.u = Pop_u2z;    Pop.w = Pop_w2z;
[Cop_x2z, Dop_u2z, Dop_w2z] = impose_fundamental_map(Top,Pop,loc_diff_tab_z,retain_xvars_z);


% % % With that, we can equivalently represent the PDE as a PIE:
% % %   D_{t} * (Tw * w + Tu * u + Tx * xf) = A  * xf + Bu  * u + Bw  * w ;
% % %                                    y  = Cy * xf + Dyu * u + Dyw * w ;
% % %                                    z  = Cz * xf + Dzu * u + Dzw * w ;
% % % where D_{t} is a diagonal temporal differential operator,
% % %   D_{t} = diag([(d/dt)^tdiff(1) , ... , (d/dt)^tdiff(nx)])
% % % We determine the orders of the temporal derivatives here
tdiff_list = ones(numel(PDE.x),1);
for ii=1:numel(PDE.x)
    if isfield(PDE.x{ii},'tdiff')
        tdiff_list(ii) = PDE.x{ii}.tdiff;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Defining the PIE structure                              % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, we define the PIE structure, collecting the PI operators that
% define our PIE.

% Keep track of the variables and domain of the system.
PIE = pie_struct();
PIE.dom = dom;
PIE.vars = vars;

% Store the PI operators defining the fundamental to PDE state map.
PIE.T = Top_x;
PIE.Tw = Top_w;
PIE.Tu = Top_u;

% Store the PI operators defining the PIE.
PIE.A = clean_opvar(Aop_x2x,op_clean_tol);
PIE.B1 = clean_opvar(Bop_w2x,op_clean_tol);
PIE.B2 = clean_opvar(Bop_u2x,op_clean_tol);
PIE.C1 = clean_opvar(Cop_x2z,op_clean_tol);
PIE.C2 = clean_opvar(Cop_x2y,op_clean_tol);
PIE.D11 = clean_opvar(Dop_w2z,op_clean_tol);
PIE.D12 = clean_opvar(Dop_u2z,op_clean_tol);
PIE.D21 = clean_opvar(Dop_w2y,op_clean_tol);
PIE.D22 = clean_opvar(Dop_u2y,op_clean_tol);

% Keep track of how the state components in the original PDE are
% now ordered in the PIE.
x_tab(:,1) = xcomp_order;       PIE.x_tab = x_tab;
y_tab(:,1) = ycomp_order;       PIE.y_tab = y_tab;
z_tab(:,1) = zcomp_order;       PIE.z_tab = z_tab;
u_tab(:,1) = ucomp_order;       PIE.u_tab = u_tab;
w_tab(:,1) = wcomp_order;       PIE.w_tab = w_tab;

% Also keep track of which derivatives of the inputs must be taken in the
% PIE representation.
if size(w_tab,2)==2+2*nvars
    PIE.w_tab(:,end-nvars+1:end) = Dvals_w;
else
    PIE.w_tab = [PIE.w_tab,Dvals_w];
end
if size(u_tab,2)==2+2*nvars
    PIE.u_tab(:,end-nvars+1:end) = Dvals_u;
else
    PIE.u_tab = [PIE.u_tab,Dvals_u];
end

% Get rid of any higher-order temporal derivatives in the PIE.              % DJ, 06/18/2025
if any(tdiff_list>1)
    [PIE,tdiff_tab] = expand_tderivatives(PIE,tdiff_list,true);
else
    tdiff_tab = [(1:numel(PDE.x))',zeros(numel(PDE.x),1)];
end

% Finally, provide an overview of how the new state variables relate to the
% PDE variables.                                                            % DJ, 06/18/2025
if ~suppress_summary
    objs = {'x';'u';'w';'y';'z'};
    for idx = 1:5
        obj = objs{idx};
        if size(PIE.([obj,'_tab']),1)>=1
            print_convert_summary(PIE,obj,tdiff_tab(:,2));
        end
    end
end

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function np_op = get_opdim(obj_tab)
% Establish dimensions of the space R^n0 x L2^n1[s1] x L2^n2[s2] x
% L2^n3[s1,s2] in which the full state, input, or output associated to
% table "obj_tab" exists. 
%
% INPUTS:
% - obj_tab:    A Nx(2+nvars) array. For j=1,...N, element
%               (j,2) should provide the size of this component, and
%               (j,3:2+nvars) should be binary indices, indicating for each
%               of the nvars variables whether the component j depends on
%               this variable.
%               NOTE: nvars should be 1 or 2
%
% OUTPUTS:
% - np_op:      A 4x1 array, indicating the total number of state, input,
%               or output variables presented in obj_tab map to
%               [ R^n0; L2^n1[s1]; L2^n2[s2]; L2^n2[s1,s2] ];
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/09/2022
%

global nvars

if nvars==0
    np_op = [sum(obj_tab(:,2)); 0; 0; 0];
elseif nvars==1
    dep_tab = obj_tab(:,3:2+nvars);
    rindcs_0 = ~any(dep_tab,2);                 % Row indices in comp_tab associated to finite dimensional states
    rindcs_1 = logical(dep_tab(:,1));           % Row indices in comp_tab associated to states that vary just in s1
    
    np_op = [sum(obj_tab(rindcs_0,2));
            sum(obj_tab(rindcs_1,2))];
elseif nvars==2
    dep_tab = obj_tab(:,3:2+nvars);
    rindcs_00 = ~any(dep_tab,2);                % Row indices in comp_tab associated to finite dimensional states
    rindcs_10 = dep_tab(:,1) & ~dep_tab(:,2);   % Row indices in comp_tab associated to states that vary just in s1
    rindcs_01 = ~dep_tab(:,1) & dep_tab(:,2);   % Row indices in comp_tab associated to states that vary just in s2
    rindcs_11 = dep_tab(:,1) & dep_tab(:,2);    % Row indices in comp_tab associated to states that vary in s1 and s2
    
    np_op = [sum(obj_tab(rindcs_00,2));
            sum(obj_tab(rindcs_10,2));
            sum(obj_tab(rindcs_01,2));
            sum(obj_tab(rindcs_11,2))];
else
    error('At most 2 spatial variables are currently supported.')
end
    
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function [Pop_b2BC, Pop_f2BC, Pop_u2BC, Pop_w2BC] = construct_PI_ops_BCs(PDE,bc_state_tab)
% Define opvar2d objects to represent the BCs of the PDE in terms of the
% fundamental state xf and a set of core boundary states xb (as defined by
% bc_state_tab) as
%   0 = Pop_b2BC * xb + Pop_f2BC * xf + Pop_u2BC * u + Pop_w2BC * w;
%
% INPUTS:
% - PDE:            A struct or pde_struct type object, defining a PDE in 
%                   the terms format (see also the "@pde_struct/initialize"
%                   function).
% - bc_state_tab:   A Nx(2+2*nvars) array, defining the core boundary state. 
%               For j=1,...N, element
%               (j,2) should provide the size of the state component;
%               (j,3:2+nvars) should be binary indices, indicating for
%               each of the nvars variables whether the component j depends
%               on this variable. If the associated fundamental state
%               component depends on a variable s, but the boundary state
%               component does not, the boundary state component is
%               evaluated at the lower boundary s=a (for s\in[a,b]).
%               (j,3+nvars:end) should be integer values, indicating for
%               each of the nvars variables to what degree the state is
%               differentiated wrt this variable.
%
% OUTPUTS:
% - Pop_b2BC, Pop_f2BC, Pop_u2BC, Pop_w2BC:
%               opvar2d objects such that, for the specified PDE structure,
%               the BCs may be represented as
%   0 = Pop_b2BC * xb + Pop_f2BC * xf + Pop_u2BC * u + Pop_w2BC * w;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/09/2022
%


% Extract necessary PDE information.
global dom vars nvars np_op
x_tab = PDE.x_tab;

% Establish indices associated to the variables in each state component,
% input component, and BC.    
nnx_arr = cumsum([0;x_tab(:,2)]);
nnu_arr = cumsum([0;PDE.u_tab(:,2)]);
nnw_arr = cumsum([0;PDE.w_tab(:,2)]);
nnBC_arr = cumsum([0;PDE.BC_tab(:,2)]);

% Establish dimensions of the core boundary state.
nb_arr = bc_state_tab(:,2);     nnb_arr = cumsum([0;nb_arr]);
% Establish input and output dimensions of the operators.
nb_op = np_op.b;                nnb_op = cumsum([0;nb_op]);
nx_op = np_op.x;                nnx_op = cumsum([0;nx_op]);
nu_op = np_op.u;                nnu_op = cumsum([0;nu_op]);
nw_op = np_op.w;                nnw_op = cumsum([0;nw_op]);
nBC_op = np_op.BC;              nnBC_op = cumsum([0;nBC_op]);

% Initialize empty operators.
if nvars<=1
    opvar Pop_f2BC Pop_b2BC Pop_w2BC Pop_u2BC;
    Pop_f2BC.dim = [nBC_op,nx_op];
    Pop_f2BC.I = dom;
    Pop_f2BC.var1 = vars(1,1);  Pop_f2BC.var2 = vars(1,2);

    Pop_b2BC.dim = [nBC_op,nb_op];
    Pop_b2BC.I = dom;
    Pop_b2BC.var1 = vars(1,1);  Pop_b2BC.var2 = vars(1,2);

    Pop_w2BC.dim = [nBC_op,nw_op];
    Pop_w2BC.I = dom;
    Pop_w2BC.var1 = vars(1,1);  Pop_w2BC.var2 = vars(1,2);

    Pop_u2BC.dim = [nBC_op,nu_op];
    Pop_u2BC.I = dom;
    Pop_u2BC.var1 = vars(1,1);  Pop_u2BC.var2 = vars(1,2);

    % For each of the operators, we first set the parameters in a cell.
    Rparams = {'P', 'Q1';
               'Q2', 'R'};
    % Certain parameters pertain multiple elements, involving multipliers or
    % partial integrators along different spatial directions.
    iscell_Rparam = [false, false;
                     false, true];

    % We set the parameters first in a struct object.
    f2BC_params = struct();    
    b2BC_params = struct();
    w2BC_params = struct();
    u2BC_params = struct();

    for kk=1:numel(Rparams)-1
        f2BC_params.(Rparams{kk}) = polynomial(Pop_f2BC.(Rparams{kk}));
        b2BC_params.(Rparams{kk}) = polynomial(Pop_b2BC.(Rparams{kk}));
        w2BC_params.(Rparams{kk}) = polynomial(Pop_w2BC.(Rparams{kk}));
        u2BC_params.(Rparams{kk}) = polynomial(Pop_u2BC.(Rparams{kk}));
    end
    f2BC_params.R = {polynomial(Pop_f2BC.R.R0); polynomial(Pop_f2BC.R.R1); polynomial(Pop_f2BC.R.R2)};
    b2BC_params.R = {polynomial(Pop_b2BC.R.R0); polynomial(Pop_b2BC.R.R1); polynomial(Pop_b2BC.R.R2)};
    w2BC_params.R = {polynomial(Pop_w2BC.R.R0); polynomial(Pop_w2BC.R.R1); polynomial(Pop_w2BC.R.R2)};
    u2BC_params.R = {polynomial(Pop_u2BC.R.R0); polynomial(Pop_u2BC.R.R1); polynomial(Pop_u2BC.R.R2)};
else
    Pop_f2BC = opvar2d([],[nBC_op,nx_op],dom,vars);
    Pop_b2BC = opvar2d([],[nBC_op,nb_op],dom,vars);
    Pop_w2BC = opvar2d([],[nBC_op,nw_op],dom,vars);
    Pop_u2BC = opvar2d([],[nBC_op,nu_op],dom,vars);

    % For each of the operators, we first set the parameters in a cell.
    Rparams = {'R00', 'R0x', 'R0y', 'R02';
               'Rx0', 'Rxx', 'Rxy', 'Rx2';
               'Ry0', 'Ryx', 'Ryy', 'Ry2';
               'R20', 'R2x', 'R2y', 'R22'};
    % Certain parameters pertain multiple elements, involving multipliers or
    % partial integrators along different spatial directions.
    iscell_Rparam = [false, false, false, false;
                     false, true,  false, true;
                     false, false, true,  true;
                     false, true,  true,  true];
           
    f2BC_params = struct(Pop_f2BC);
    b2BC_params = struct(Pop_b2BC);
    w2BC_params = struct(Pop_w2BC);
    u2BC_params = struct(Pop_u2BC);
end




% % % Loop over all BCs.
% % % For each BC, loop over all the terms, and add the coefficients that
% % % appear in each term to the appropriate parameter in the appropriate
% % % operator to represent the BCs as
% % %   0 = b2BC_op * x_bc + f2BC_op * x_f + u2BC_op * u + w2BC_op * w;
for BCnum = 1:numel(PDE.BC)
    
    % Establish which row indices in the operators are linked to this BC.
    rindcs = nnBC_arr(BCnum)+1:nnBC_arr(BCnum+1);
    Pop_rnum = (rindcs(1)>nnBC_op(1:end-1) & rindcs(1)<=nnBC_op(2:end));
    rindcs = rindcs - nnBC_op(Pop_rnum);     % Row indices in the actual parameters, associated to the current component
    
    % % Loop over all terms in the current BC.
    for jj=1:numel(PDE.BC{BCnum}.term)
        term_jj = PDE.BC{BCnum}.term{jj};
        
        if isfield(term_jj,'w')
            % % The term involves an exogenous input.
            % Determine which input component is involved in the term.
            comp_indx = term_jj.w;
            % Establish the column numbers in the operator associated to
            % this input component.
            cindcs = nnw_arr(comp_indx)+1:nnw_arr(comp_indx+1);
            Pop_cnum = (cindcs(1)>nnw_op(1:end-1) & cindcs(1)<=nnw_op(2:end));
            cindcs = cindcs - nnw_op(Pop_cnum);
            % Set the appropriate rows and columns of the parameter.
            Rparam = Rparams{Pop_rnum,Pop_cnum};
            if iscell_Rparam(Pop_rnum,Pop_cnum)
            %w2BC_params{Pop_rnum,Pop_cnum} = polynomial(w2BC_params{Pop_rnum,Pop_cnum}{1});
                w2BC_params.(Rparam){1}(rindcs,cindcs) = w2BC_params.(Rparam){1}(rindcs,cindcs) + term_jj.C;
            else
                w2BC_params.(Rparam)(rindcs,cindcs) = w2BC_params.(Rparam)(rindcs,cindcs) + term_jj.C;
            end
        elseif isfield(term_jj,'u')
            % % The term involves an actuator input.
            % Determine which input component is involved in the term.
            comp_indx = term_jj.u;
            % Establish the column numbers in the operator associated to
            % this input component.
            cindcs = nnu_arr(comp_indx)+1:nnu_arr(comp_indx+1);
            Pop_cnum = (cindcs(1)>nnu_op(1:end-1) & cindcs(1)<=nnu_op(2:end));
            cindcs = cindcs - nnu_op(Pop_cnum);
            % Set the appropriate rows and columns of the parameter.
            Rparam = Rparams{Pop_rnum,Pop_cnum};
            if iscell_Rparam(Pop_rnum,Pop_cnum)
                u2BC_params.(Rparam){1}(rindcs,cindcs) = u2BC_params.(Rparam){1}(rindcs,cindcs) + term_jj.C;
            else
                u2BC_params.(Rparam)(rindcs,cindcs) = u2BC_params.(Rparam)(rindcs,cindcs) + term_jj.C;
            end
        else
            % % The term involves an state component.
            % % --> We have to check if a fundamental state component x_f,
            % %     or core boundary component x_bc is involved.
            % Determine which state component appears.
            comp_indx = term_jj.x;
            has_vars_xcomp = logical(x_tab(comp_indx,3:2+nvars));
            nvars_xcomp = sum(has_vars_xcomp);
            % Establish the maximal order of differentiability of the
            % considered state component.
            Dmax = x_tab(comp_indx,3+nvars:2+2*nvars);
            Dmax = Dmax(has_vars_xcomp);
            if all(term_jj.D==Dmax)
                % % A fundamental state component is considered:
                % Establish the column numbers in the operator associated to
                % this fundamental state component.
                cindcs = nnx_arr(comp_indx)+1:nnx_arr(comp_indx+1);
                Pop_cnum = (cindcs(1)>nnx_op(1:end-1) & cindcs(1)<=nnx_op(2:end));
                cindcs = cindcs - nnx_op(Pop_cnum);
                Rparam = Rparams{Pop_rnum,Pop_cnum};
                
                % % Establish which elements of the parameter must be set;
                % % are we performing integration?
                Cval = polynomial(term_jj.C);
                if ~iscell_Rparam(Pop_rnum,Pop_cnum)
                    % If the involved parameter does not involve partial
                    % integration, just set the term.
                    Cval = subs(Cval,vars(:,2),vars(:,1));
                    f2BC_params.(Rparam)(rindcs,cindcs) = f2BC_params.(Rparam)(rindcs,cindcs) + Cval;
                else
                    % Otherwise, we'll have to check whether multiplier or
                    % integration is desired.
                    Idoms = term_jj.I;

                    % First, estbalish along which directions integration is
                    % allowed, by checking the dimensions of the parameter.
                    param_sz_full = cell(1,nvars);
                    [param_sz_full{:}] = size(f2BC_params.(Rparam));
                    param_sz_full = cell2mat(param_sz_full);
                    param_sz = param_sz_full(has_vars_xcomp);
                    % For each direction along which integration is possible,
                    % establish whether an integral on _a^s, _s^b or _a^b is
                    % desired
                    Pop_int_indx = ones(1,nvars_xcomp);
                    for ll=1:nvars_xcomp
                        if param_sz(ll)==1
                            % An integral is taken over the full domain
                            % anyway --> replace dummy vars in the
                            % kernel with primary vars.
                            vars_ll = vars(has_vars_xcomp,:);
                            vars_ll = vars_ll(ll,:);
                            Cval = subs(Cval,vars_ll(2),vars_ll(1));
                        elseif ~isempty(term_jj.I{ll})
                            if isa(Idoms{ll},'double') || isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll})
                                % Integrating over full spatial domain.
                                % A full integral has to be constructed
                                % using two partial integrals.
                                Pop_int_indx1 = Pop_int_indx;
                                Pop_int_indx2 = Pop_int_indx;
                                Pop_int_indx1(:,ll) = 2;
                                Pop_int_indx2(:,ll) = 3;
                                Pop_int_indx = [Pop_int_indx1; Pop_int_indx2];
                            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(1))
                                % Integrating over lower "half" of domain.
                                Pop_int_indx(:,ll) = 2;
                            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(2))
                                % Integrating over upper "half" of domain.
                                Pop_int_indx(:,ll) = 3;
                            end
                        end
                    end
                    % Each colum in "Pop_int_indx_full" corresponds to one of
                    % the global variables in the PDE. Each row corresponds to
                    % an integral that the user has specified for this term. 
                    Pop_int_indx_full = ones(size(Pop_int_indx,1),nvars);
                    Pop_int_indx_full(:,has_vars_xcomp) = Pop_int_indx;

                    % Set the parameters for each of the desired partial
                    % integrals/multipliers.
                    param_linsz = cumprod([1,param_sz_full]);
                    param_linsz = param_linsz(1:end-1);
                    for ll=1:size(Pop_int_indx_full,1)
                        %int_lindx = (Pop_int_indx(ll,:)-1)*param_linsz' + 1;
                        int_lindx = (Pop_int_indx_full(ll,:)-1)*param_linsz' + 1;
                        f2BC_params.(Rparam){int_lindx}(rindcs,cindcs) = f2BC_params.(Rparam){int_lindx}(rindcs,cindcs) + Cval;
                    end
                end
            else
                % % A core boundary state component is considered:
                % To determine the column indices in the operator
                % associated to this state component, we have to check at
                % what position it is evaluated, and what derivative is
                % taken.
                % First check on which of the global variables the state
                % still depends.
                if isa(term_jj.loc,'double') || isdouble(term_jj.loc)
                    is_intrr_full = zeros(1,nvars);
                else
                    is_intrr_full = zeros(1,nvars);
                    is_intrr = zeros(1,nvars_xcomp);
                    for ll=1:nvars_xcomp
                        is_intrr(ll) = ~isdouble(term_jj.loc(ll));
                    end
                    is_intrr_full(has_vars_xcomp) = is_intrr;
                end
                % Further check to what order the state is differentiated
                % wrt each of the global variables.
                Dval_full = zeros(1,nvars);
                Dval_full(has_vars_xcomp) = term_jj.D;
                % Finally, also check the size of the component.
                sz_xcomp = x_tab(comp_indx,2);
                % Collect all the info on the considered state component,
                % and check which of the core boundary components this
                % state component must be.
                xcomp_info = [comp_indx,sz_xcomp,is_intrr_full,Dval_full];
                cindx1 = all(bc_state_tab==xcomp_info,2);
                cindcs = nnb_arr(cindx1)+1 : nnb_arr(cindx1)+sz_xcomp;
                Pop_cnum = (cindcs(1)>nnb_op(1:end-1) & cindcs(1)<=nnb_op(2:end));
                cindcs = cindcs - nnb_op(Pop_cnum);
                Rparam = Rparams{Pop_rnum,Pop_cnum};
                
                % % Establish which elements of the parameter must be set;
                % % are we performing integration?
                % % Establish which elements of the parameter must be set;
                % % are we performing integration?
                Cval = polynomial(term_jj.C);
                if ~iscell_Rparam(Pop_rnum,Pop_cnum)
                    % If the considered parameter does not involve partial
                    % integration, just set the term.
                    Cval = subs(Cval,vars(:,2),vars(:,1));
                    b2BC_params.(Rparam)(rindcs,cindcs) = b2BC_params.(Rparam)(rindcs,cindcs) + Cval;
                else
                    % Otherwise, we'll have to check whether multiplier or
                    % integration is desired.
                    Idoms = term_jj.I;
                    
                    % First, estbalish along which directions integration is
                    % allowed, by checking the dimensions of the parameter.
                    param_sz_full = cell(1,nvars);
                    [param_sz_full{:}] = size(b2BC_params.(Rparam));
                    param_sz_full = cell2mat(param_sz_full);
                    param_sz = param_sz_full(has_vars_xcomp);
                    % For each direction along which integration is possible,
                    % establish whether an integral on _a^s, _s^b or _a^b is
                    % desired
                    Pop_int_indx = ones(1,nvars_xcomp);
                    for ll=1:nvars_xcomp
                        if param_sz(ll)==1
                            % An integral is taken over the full domain
                            % anyway --> replace dummy vars in the
                            % kernel with primary vars.
                            vars_ll = vars(has_vars_xcomp,:);
                            vars_ll = vars_ll(ll,:);
                            Cval = subs(Cval,vars_ll(2),vars_ll(1));
                        elseif ~isempty(term_jj.I{ll})
                            if isa(Idoms{ll},'double') || isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll})
                                % Integrating over full spatial domain.
                                Pop_int_indx1 = Pop_int_indx;
                                Pop_int_indx2 = Pop_int_indx;
                                Pop_int_indx1(:,ll) = 2;
                                Pop_int_indx2(:,ll) = 3;
                                Pop_int_indx = [Pop_int_indx1; Pop_int_indx2];
                            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(1))
                                % Integrating over lower "half" of domain.
                                Pop_int_indx(:,ll) = 2;
                            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(2))
                                % Integrating over upper "half" of domain.
                                Pop_int_indx(:,ll) = 3;
                            end
                        end
                    end
                    % Each colum in "Pop_int_indx_full" corresponds to one of
                    % the global variables in the PDE. Each row corresponds to
                    % an integral that the user has specified for this term. 
                    Pop_int_indx_full = ones(size(Pop_int_indx,1),nvars);
                    Pop_int_indx_full(:,has_vars_xcomp) = Pop_int_indx;

                    % Set the parameters for each of the desired partial
                    % integrals/multipliers.
                    param_linsz = cumprod([1,param_sz_full]);
                    param_linsz = param_linsz(1:end-1);
                    for ll=1:size(Pop_int_indx_full,1)
                        int_lindx = (Pop_int_indx_full(ll,:)-1)*param_linsz' + 1;
                        %b2BC_params{Pop_rnum,Pop_cnum}{int_lindx} = polynomial(b2BC_params{Pop_rnum,Pop_cnum}{int_lindx});
                        b2BC_params.(Rparam){int_lindx}(rindcs,cindcs) = b2BC_params.(Rparam){int_lindx}(rindcs,cindcs) + Cval;
                    end
                end
            end
        end
    end
end

% Having determined values for all the parameters, now set the parameters 
% of the operators.
if nvars<=1
    for kk=1:numel(Rparams)-1
        Pop_b2BC.(Rparams{kk}) = b2BC_params.(Rparams{kk});
        Pop_f2BC.(Rparams{kk}) = f2BC_params.(Rparams{kk});
        Pop_w2BC.(Rparams{kk}) = w2BC_params.(Rparams{kk});
        Pop_u2BC.(Rparams{kk}) = u2BC_params.(Rparams{kk});
    end
    RRparams = {'R0','R1','R2'};
    for kk=1:numel(RRparams)
        Pop_b2BC.R.(RRparams{kk}) = b2BC_params.R{kk};
        Pop_f2BC.R.(RRparams{kk}) = f2BC_params.R{kk};
        Pop_w2BC.R.(RRparams{kk}) = w2BC_params.R{kk};
        Pop_u2BC.R.(RRparams{kk}) = u2BC_params.R{kk};
    end
else
    for kk=1:numel(Rparams)
        Pop_b2BC.(Rparams{kk}) = b2BC_params.(Rparams{kk});
        Pop_f2BC.(Rparams{kk}) = f2BC_params.(Rparams{kk});
        Pop_w2BC.(Rparams{kk}) = w2BC_params.(Rparams{kk});
        Pop_u2BC.(Rparams{kk}) = u2BC_params.(Rparams{kk});
    end
end

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function [Pop_x_cell, Pop_u, Pop_w, loc_diff_tab, retain_xvars_cell] = construct_PI_ops_PDE(PDE,obj)
% Construct PI operators such that the presented PDE for obj can be wirtten
% as
%   obj = sum_j=1^n1 (P_x2obj{j} * Delta_j * D_j)* x + P_u2obj* u + P_w2obj* w;
% where Delta_j is a Delta operator, evaluating the PDE state at a
% particular spatial position, and D_j is a differential operator.
%
% INPUTS:
% - PDE:    A struct or pde_struct type object, defining a PDE in the terms
%           format (see also the "@pde_struct/initialize" function).
% - obj:    'x', 'y', or 'z', indicating which equations to be represented
%           using PI operators.
%
% OUTPUTS:
% - Pop_x_cell: Nx1 cell of opvar2d objects. Each object corresponds to a
%               PI operator describing the contribution of a particular
%               spatial derivative of the PDE state, evaluated at a
%               particular spatial position, to the PDE of "obj".
% - Pop_u, Pop_w:   opvar2d objects describing the contributions of the
%                   inputs u and w to the PDE of "obj".
% - loc_diff_tab:   Nx2*nvars array, of which each row is linked to an 
%                   element of Pop_x_cell, where
%                   (j,1:nvars) are indices in {-1,0,1}, indicating whether
%                   the state is evaluated at the lower boundary s=a (-1),
%                   interior s=s (0) or upper boundary s=b (+1) of the
%                   domain, for each of the nvars variables s\in[a,b];
%                   (j,nvars+1:end) are integer values indicating the order
%                   of the derivative of the PDE state wrt each of the
%                   nvars variables.
% - retain_xvars_cell:  Nx1 cell of integer arrays, indicating for each of
%                       the operators Pop_x_cell{j} which state components
%                       they map i.e. which state components are
%                       differentiable up to the required degree as
%                       specified in loc_diff_tab(j,:).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/09/2022
%

% Extract some necessary information from the PDE.
global dom vars nvars np_op
x_tab = PDE.x_tab;

% Establish indices associated to the variables in each state component,
% input component, and output component.    
nnx_arr = cumsum([0;PDE.x_tab(:,2)]);
nnu_arr = cumsum([0;PDE.u_tab(:,2)]);
nnw_arr = cumsum([0;PDE.w_tab(:,2)]);
nnr_arr = cumsum([0;PDE.([obj,'_tab'])(:,2)]);

% Establish input and output dimensions of the operators.
%nx_op = np_op.x;        %nnx_op = cumsum([0;nx_op]);
nu_op = np_op.u;        nnu_op = cumsum([0;nu_op]);
nw_op = np_op.w;        nnw_op = cumsum([0;nw_op]);
nr_op = np_op.(obj);    nnr_op = cumsum([0;nr_op]);

% Initialize empty operators.
%Pop_x = opvar2d([],[nr_op,nx_op],dom,vars);
if nvars<=1
    opvar Pop_u Pop_w;
    Pop_u.dim = [nr_op,nu_op];
    Pop_u.I = dom;
    Pop_u.var1 = vars(1,1);     Pop_u.var2 = vars(1,2);
    Pop_w.dim = [nr_op,nw_op];
    Pop_w.I = dom;
    Pop_w.var1 = vars(1,1);     Pop_w.var2 = vars(1,2);
else
    Pop_u = opvar2d([],[nr_op,nu_op],dom,vars);
    Pop_w = opvar2d([],[nr_op,nw_op],dom,vars);
end
% For the state components, we consider a PI operator for each combination
% of position and derivative of the state.
loc_diff_tab = zeros(0,2*nvars);
Pop_x_cell = cell(0,1);
retain_xvars_cell = cell(0,1);

% If the "obj" has no components, we are done.
if nnr_op(end)==0    
    return
end

% Otherwise for each of the operators, we set the parameters in a cell.
if nvars<=1
    Rparams = {'P', 'Q1';
               'Q2', 'R'};
    % Certain parameters pertain multiple elements, involving multipliers or
    % partial integrators along different spatial directions.
    iscell_Rparam = [false, false;
                     false, true];

    params_u = struct();
    params_w = struct();
    for kk=1:numel(Rparams)-1
        params_u.(Rparams{kk}) = polynomial(Pop_u.(Rparams{kk}));
        params_w.(Rparams{kk}) = polynomial(Pop_w.(Rparams{kk}));
    end
    params_u.R = {polynomial(Pop_u.R.R0); polynomial(Pop_u.R.R1); polynomial(Pop_u.R.R2)};
    params_w.R = {polynomial(Pop_w.R.R0); polynomial(Pop_w.R.R1); polynomial(Pop_w.R.R2)};
else
    Rparams = {'R00', 'R0x', 'R0y', 'R02';
               'Rx0', 'Rxx', 'Rxy', 'Rx2';
               'Ry0', 'Ryx', 'Ryy', 'Ry2';
               'R20', 'R2x', 'R2y', 'R22'};
           
    % Certain parameters pertain multiple elements, involving multipliers or
    % partial integrators along different spatial directions.
    iscell_Rparam = [false, false, false, false;
                     false, true,  false, true;
                     false, false, true,  true;
                     false, true,  true,  true];
    
    params_u = struct(Pop_u);
    params_w = struct(Pop_w);
end

% For each PI operator for the state components, we keep track of which
% state components this operator maps, as certain state components may not
% be differentiable up to the desired degree.
params_x_cell = cell(0,1);
nnx_op_cell = cell(0,1);
nnx_arr_cell = cell(0,1);


% % % Loop over all equations of the considered "obj" (PDEs or outputs)
% % % For each equation, loop over all the terms, and add the coefficients
% % % that appear in each term to the appropriate parameter in the
% % % appropriate operator to represent the equation as
% % %   obj = sum_j=1^n1 (P_x2obj{j} * Delta_j * D_j)* x
% % %               + P_u2obj* u + P_w2obj* w;
% % % where Delta_j and D_j are differential and Delta operators,
% % % corresponding to the jth row of loc_diff_tab.
for eqnum=1:numel(PDE.(obj))
        
    % Establish which row indices in the operators are linked to this 
    % equation.
    rindcs = nnr_arr(eqnum)+1:nnr_arr(eqnum+1);
    Pop_rnum = (rindcs(1)>nnr_op(1:end-1) & rindcs(1)<=nnr_op(2:end));
    rindcs = rindcs - nnr_op(Pop_rnum);     % Row indices in the actual parameters, associated to the current component
    
    for jj=1:numel(PDE.(obj){eqnum}.term)
        term_jj = PDE.(obj){eqnum}.term{jj};
        
        if isfield(term_jj,'w')
            % Determine which input component is involved in the term.
            comp_indx = term_jj.w;
            % Establish the column numbers in the operator associated to
            % this input component.
            cindcs = nnw_arr(comp_indx)+1:nnw_arr(comp_indx+1);
            Pop_cnum = (cindcs(1)>nnw_op(1:end-1) & cindcs(1)<=nnw_op(2:end));
            cindcs = cindcs - nnw_op(Pop_cnum);
            
            % Add the coefficients to the appropriate parameter.
            Rparam = Rparams{Pop_rnum,Pop_cnum};
            if iscell_Rparam(Pop_rnum,Pop_cnum)
                params_w.(Rparam){1}(rindcs,cindcs) = params_w.(Rparam){1}(rindcs,cindcs) + term_jj.C;
            else
                params_w.(Rparam)(rindcs,cindcs) = params_w.(Rparam)(rindcs,cindcs) + term_jj.C;
            end
            
        elseif isfield(term_jj,'u')
            % Determine which input component is involved in the term.
            comp_indx = term_jj.u;
            % Establish the column numbers in the operator associated to
            % this input component.
            cindcs = nnu_arr(comp_indx)+1:nnu_arr(comp_indx+1);
            Pop_cnum = (cindcs(1)>nnu_op(1:end-1) & cindcs(1)<=nnu_op(2:end));
            cindcs = cindcs - nnu_op(Pop_cnum);
            
            % Add the coefficients to the appropriate parameter.
            Rparam = Rparams{Pop_rnum,Pop_cnum};
            if iscell_Rparam(Pop_rnum,Pop_cnum)
                params_u.(Rparam){1}(rindcs,cindcs) = params_u.(Rparam){1}(rindcs,cindcs) + term_jj.C;
            else
                params_u.(Rparam)(rindcs,cindcs) = params_u.(Rparam)(rindcs,cindcs) + term_jj.C;
            end
            
        else
            % Determine which PDE state component is involved in this term.
            comp_indx = term_jj.x;
            
            % Determine which variables the state component depends on.
            has_vars_xcomp = logical(PDE.x_tab(comp_indx,3:2+nvars));
            nvars_xcomp = sum(has_vars_xcomp);
            
            % Establish what derivative of the state is considered.
            Dval = term_jj.D;
            Dval_full = zeros(1,nvars);
            Dval_full(has_vars_xcomp) = Dval;
            
            % Establish at what position the state is evaluated.
            loc_val_full = zeros(1,nvars);
            loc_val = zeros(1,nvars_xcomp); % -1 for lower boundary, 0 for interior, 1 for upper boundary
            Rloc = term_jj.loc;
            Rdom = PDE.dom(has_vars_xcomp,:);
            if ~isa(Rloc,'double') && ~isdouble(Rloc)
                for kk=1:nvars_xcomp
                    if isdouble(Rloc(kk)) && double(Rloc(kk))==Rdom(kk,1)
                        loc_val(kk) = -1;
                    elseif isdouble(Rloc(kk)) && double(Rloc(kk))==Rdom(kk,2)
                        loc_val(kk) = 1;
                    end
                end
                loc_val_full(has_vars_xcomp) = loc_val;
            elseif ~isempty(Rloc)
                % If the state is evaluated at a corner, we can avoid a for
                % loop.
                is_lower_bndry = double(Rloc)==Rdom(:,1)';
                loc_val(is_lower_bndry) = -1;
                loc_val(~is_lower_bndry) = 1;
                loc_val_full(has_vars_xcomp) = loc_val;
            end
            
            % Check whether the desired combination of location and
            % derivative is already accounted for in our set of parameters.
            Pop_op_indx = ismember(loc_diff_tab,[loc_val_full,Dval_full],'rows');
            if ~any(Pop_op_indx)
                % Establish which state components are differentiable up to
                % the necessary degree.
                use_loc = abs(loc_val_full);
                retain_rows = all(x_tab(:,3+nvars:2+2*nvars)>=Dval_full+use_loc,2);
                nx_arr_new = x_tab(:,2);
                nx_arr_new(~retain_rows) = 0;
                nnx_arr_new = cumsum([0; nx_arr_new]);
                nc_op = get_opdim([x_tab(:,1),nx_arr_new,x_tab(:,3:end)]);
                
                % Determine for each combination of variables whether it
                % allows for the desired substitution to be performed. That
                % is, if substituting variable s, then variable s must be
                % present in the combination of variables.
                issub_op = false([2*ones(1,nvars),1]);
                for ll = 2:numel(issub_op)
                    use_var_list = cell(1,nvars);
                    [use_var_list{:}] = ind2sub(size(issub_op),ll);
                    use_var_list = logical(cell2mat(use_var_list)-1);
                    issub_op(ll) = any(use_loc(use_var_list));                    
                end
                % After substitution, states that depend on the substituted
                % variable will no longer depend on this variable. The size
                % of these states will therefore have to replace that
                % of the state components that do not depend on this
                % variable.
                nc_op = reshape(nc_op,[2*ones(1,nvars),1]);
                sub_indcs = find(issub_op)';
                for ll = sub_indcs(end:-1:1)
                    % Which combination of variables are we considering?
                    use_var_list = cell(1,nvars);
                    [use_var_list{:}] = ind2sub(size(issub_op),ll);
                    use_var_list = cell2mat(use_var_list)-1;
                    % Remove the last variable from the combination.
                    use_var_sub_list = use_var_list & use_loc;
                    last_var = find(use_var_sub_list,1,'last');
                    use_var_list(last_var) = 0;
                    % Add the size of state components that depend on the
                    % considered combination to that of components that
                    % depend on the combination excluding the last var.
                    new_idx = use_var_list*(2.^(0:nvars-1))' + 1;
                    nc_op(new_idx) = nc_op(ll); % Note that the old nc_op(ll) will be overwritten, as it does not depend on all required variables
                    nc_op(ll) = 0;
                end
                nc_op = nc_op(:);
                
                % Initalize a new operator that maps only the allowed state
                % components, evaluated at the desired boundaries.
                if nvars<=1
                    opvar Pop_new;
                    Pop_new.dim = [nr_op,nc_op];
                    Pop_new.I = dom;
                    Pop_new.var1 = vars(1,1);   Pop_new.var2 = vars(1,2);
                    
                    params_new = struct();
                    for kk=1:numel(Rparams)-1
                        params_new.(Rparams{kk}) = polynomial(Pop_new.(Rparams{kk}));
                    end
                    params_new.R = {polynomial(Pop_new.R.R0), polynomial(Pop_new.R.R1), polynomial(Pop_new.R.R2)};

                else
                    Pop_new = opvar2d([],[nr_op,nc_op],dom,vars);
                    params_new = struct(Pop_new);
                end
                
                % Keep track of which state variables the operator maps.
                indx_l = nnx_arr(1:end-1);      indx_u = nnx_arr(2:end);
                retain_xvars = mat2cell([indx_l(retain_rows),indx_u(retain_rows)],ones(sum(retain_rows),1));
                retain_xvars = cellfun(@(x) (x(1)+1:x(2))',retain_xvars,'UniformOutput',false);
                retain_xvars = cell2mat(retain_xvars);
                retain_xvars_cell = [retain_xvars_cell; retain_xvars];
                
                % Also keep track of the dimensions of the new operator.
                nnx_op_cell = [nnx_op_cell; cumsum([0;nc_op])];
                nnx_arr_cell = [nnx_arr_cell; nnx_arr_new];
                
                % Add the new combination to the list of combinations.
                loc_diff_tab = [loc_diff_tab;
                                [loc_val_full, Dval_full]];
                % Add a new operator to the cell.
                params_x_cell = [params_x_cell;
                                 {params_new}];
                Pop_x_cell = [Pop_x_cell;
                              {Pop_new}];
                
                % We will modify the parameters of the new element.
                Pop_op_indx = length(params_x_cell);
            end
            
            % Establish which parameter in the operator correspond to the
            % desired state component, and which column indices in this
            % parameter.
            cindcs = nnx_arr_cell{Pop_op_indx}(comp_indx)+1:nnx_arr_cell{Pop_op_indx}(comp_indx+1);
            Pop_cnum = (cindcs(1)>nnx_op_cell{Pop_op_indx}(1:end-1) & cindcs(1)<=nnx_op_cell{Pop_op_indx}(2:end));
            cindcs = cindcs - nnx_op_cell{Pop_op_indx}(Pop_cnum);
            Rparam = Rparams{Pop_rnum, Pop_cnum};

            % % Establish which elements of the parameter must be set;
            % % are we performing integration?
            Cval = polynomial(term_jj.C);
            if ~iscell_Rparam(Pop_rnum,Pop_cnum)
                % If the involved parameter does not involve partial
                % integration, just set the term.
                Cval = subs(Cval,vars(:,2),vars(:,1));
                params_x_cell{Pop_op_indx}.(Rparam)(rindcs,cindcs) = params_x_cell{Pop_op_indx}.(Rparam)(rindcs,cindcs) + Cval;
            else
                Idoms = term_jj.I;
                % First, estbalish along which directions integration is
                % allowed, by checking the dimensions of the parameter.
                param_sz_full = cell(1,nvars);
                [param_sz_full{:}] = size(params_x_cell{Pop_op_indx}.(Rparam));
                param_sz_full = cell2mat(param_sz_full);
                param_sz = param_sz_full(has_vars_xcomp);
                % For each direction along which integration is possible,
                % establish whether an integral on _a^s, _s^b or _a^b is
                % desired
                Pop_int_indx = ones(1,nvars_xcomp);
                for ll=1:nvars_xcomp
                    if param_sz(ll)==1
                        % An integral is taken over the full domain
                        % anyway --> replace dummy vars in the
                        % kernel with primary vars.
                        vars_ll = PDE.vars(has_vars_xcomp,:);
                        vars_ll = vars_ll(ll,:);
                        Cval = subs(Cval,vars_ll(2),vars_ll(1));
                    elseif ~isempty(term_jj.I{ll})
                        if isa(Idoms{ll},'double') || isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll})
                            % Integrating over full spatial domain.
                            % A full integral has to be constructed
                            % using two partial integrals.
                                Pop_int_indx1 = Pop_int_indx;
                                Pop_int_indx2 = Pop_int_indx;
                                Pop_int_indx1(:,ll) = 2;
                                Pop_int_indx2(:,ll) = 3;
                                Pop_int_indx = [Pop_int_indx1; Pop_int_indx2];
                        elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(1))
                            % Integrating over lower "half" of domain.
                            Pop_int_indx(:,ll) = 2;
                        elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(2))
                            % Integrating over upper "half" of domain.
                            Pop_int_indx(:,ll) = 3;
                        end
                    end
                end
                % Each colum in "Pop_int_indx_full" corresponds to one of
                % the global variables in the PDE. Each row corresponds to
                % an integral that the user has specified for this term. 
                Pop_int_indx_full = ones(size(Pop_int_indx,1),nvars);
                Pop_int_indx_full(:,has_vars_xcomp) = Pop_int_indx;
            
                % Set the parameters for each of the desired partial
                % integrals/multipliers.
                if nvars<=1
                    param_linsz = 1;
                else
                    param_linsz = cumprod([1,param_sz_full]);
                    param_linsz = param_linsz(1:end-1);
                end
                for ll=1:size(Pop_int_indx_full,1)
                    int_lindx = (Pop_int_indx_full(ll,:)-1)*param_linsz' + 1;
                    params_x_cell{Pop_op_indx}.(Rparam){int_lindx}(rindcs,cindcs) = params_x_cell{Pop_op_indx}.(Rparam){int_lindx}(rindcs,cindcs) + Cval;
                end
            end               
        end
    end
end


% Having determined values for all the parameters, now set the parameters 
% of the operators.
if nvars<=1
    RRparams = {'R0','R1','R2'};
    for ll=1:numel(params_x_cell)
        for kk=1:numel(Rparams)-1
            Pop_x_cell{ll}.(Rparams{kk}) = params_x_cell{ll}.(Rparams{kk});
        end
        for kk=1:numel(RRparams)
            Pop_x_cell{ll}.R.(RRparams{kk}) = params_x_cell{ll}.R{kk};
        end
    end
    for kk=1:numel(Rparams)-1
        Pop_w.(Rparams{kk}) = params_w.(Rparams{kk});
        Pop_u.(Rparams{kk}) = params_u.(Rparams{kk});
    end
    for kk=1:numel(RRparams)
        Pop_w.R.(RRparams{kk}) = params_w.R{kk};
        Pop_u.R.(RRparams{kk}) = params_u.R{kk};
    end
else
    for kk=1:numel(Rparams)
        for ll=1:numel(params_x_cell)
            Pop_x_cell{ll}.(Rparams{kk}) = params_x_cell{ll}.(Rparams{kk});
        end
        Pop_w.(Rparams{kk}) = params_w.(Rparams{kk});
        Pop_u.(Rparams{kk}) = params_u.(Rparams{kk});
    end
end

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function [Pop_x, Pop_u, Pop_w] = impose_fundamental_map(Top,Pop,loc_diff_tab,retain_xvars_cell)
% Construct PI operators such that the PDE corresponding to operators Pop,
%   obj = sum_j=1^n1 (P_x2obj{j} * Delta_j * D_j)* x + P_u2obj* u + P_w2obj* w;
% can be represented in terms of the fundamental state as
%   obj = Pop_x * xf + Pop_u * u + Pop_w;
%
% INPUTS:
% - Top:    Struct with field 'x', 'u' and 'w', containing PI operators Tx,
%           Tu and Tw such that the PDE state x can be expressed in terms
%           of the fundamental state xf as
%               x = Tx * xf + Tu * u + Tw * w;
% - Pop:    Struct with field 'x', 'u' and 'w'
%           Pop.x should be an Nx1 cell of opvar2d objects, and Pop_u and
%           Pop_w should be opvar2d objects, defining a PDE as in
%           "construct_PI_ops_PDE".
% - loc_diff_tab:   Nx2*nvars array that comes with Pop.x, see the outputs
%                   of "construct_PI_ops_PDE".
% - retain_xvars_cell:  Nx1 cell providing indices of the state components
%                       that each element Pop.x{j} maps, see the outputs
%                       of "construct_PI_ops_PDE".
%
% OUTPUTS:
% - Pop_x, Pop_u, Pop_w:    opvar2d objects describing the PIE associated
%                           to the input operators:
%       Pop_x = sum_j=1^n1 (Pop.x{j} * Delta_j * D_j * Top.x);
%       Pop_u = sum_j=1^n1 (Pop.x{j} * Delta_j * D_j * Top.u);
%       Pop_w = sum_j=1^n1 (Pop.x{j} * Delta_j * D_j * Top.w);
% where the composition of the differential and Delta operators with the
% operators Top are PI operators, guaranteeing that Pop_x, Pop_u and Pop_w
% are also PI operators.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/09/2022
%

% Construct an array indicating for each spatial direction the lower
% boundary, spatial variable, and upper boundary.
global dom vars nvars
all_locs = [dom(:,1)'; vars(:,1)'; dom(:,2)'];

% Extract the pre-defined PI operators.
Top_x = Top.x;      Pop_x_cell = Pop.x;
Top_u = Top.u;      Pop_u = Pop.u;
Top_w = Top.w;      Pop_w = Pop.w;


% Initialize an operator to map the fundamental state.
if nvars<=1
    opvar Pop_x;
    Pop_x.dim = [Pop_u.dim(:,1),Top_x.dim(:,2)];
    Pop_x.I = dom;
    Pop_x.var1 = vars(1,1);     Pop_x.var2 = vars(1,2);
else
    Pop_x = opvar2d([],[Pop_u.dim(:,1),Top_x.dim(:,2)],dom,vars);
end

for jj = 1:numel(Pop_x_cell)
    % Establish what derivative of the state component is considered.
    Dval = loc_diff_tab(jj,nvars+1:2*nvars);
    
    % Establish at what position the state is evaluated.
    loc_val = loc_diff_tab(jj,1:nvars);
    loc_val_lin = loc_val+2+(0:nvars-1)*3;
    Rloc = all_locs(loc_val_lin);

    % If the state is differentiated wrt some variable s, the state must be
    % differentiable wrt this variable up to the desired order. We extract
    % only the PDE state variables (rows of Tx) satisfying this condition.
    retain_xvars = retain_xvars_cell{jj};
                 
    % Take the composition of the desired differential and Delta operator
    % with the PI operators T.
    Top_x_jj = Top_x(retain_xvars,:);
    Top_x_jj = diff(Top_x_jj,vars(:,1),Dval','pure');
    Top_x_jj = subs(Top_x_jj,vars(:,1),Rloc','pure');
    
    % Take the composition of the resulting operator with the coefficient
    % operator, and add to the PIE.
    Pop_jj = Pop_x_cell{jj};
    PTop = Pop_jj * Top_x_jj;
    Pop_x = Pop_x + PTop;
    
    % Repeat for the input signals.
    if any(Top_u.dim(:,2))
        Top_u_jj = Top_u(retain_xvars,:);
        Top_u_jj = diff(Top_u_jj,vars(:,1),Dval','pure');
        Top_u_jj = subs(Top_u_jj,vars(:,1),Rloc','pure');
        PTop = Pop_jj * Top_u_jj;
        Pop_u = Pop_u + PTop;
    end
    if any(Top_w.dim(:,2))
        Top_w_jj = Top_w(retain_xvars,:);
        Top_w_jj = diff(Top_w_jj,vars(:,1),Dval','pure');
        Top_w_jj = subs(Top_w_jj,vars(:,1),Rloc','pure');
        PTop = Pop_jj * Top_w_jj;
        Pop_w = Pop_w + PTop;
    end    
end

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_convert_summary(PIE,obj,tdiff_list)
% print_reorder_summary(PIE,obj,ncomps_x)
% prints information in the command window on how the state variables,
% inputs, and outputs in the computed PIE structure relate to those in the
% original PDE.
%
% INPUTS:
% - PIE:    A "pie_struct" class object defining a PIE.
% - obj:    Char 'x', 'u', 'w', 'y', or 'z', indicating for which
%           object to display how they relate to the associated object in 
%           the PDE.
% - tdiff_list: If 'obj' = 'x', this should be an nx x 1 array indicating
%               for each of the nx (vector-valued_ fundamental state 
%               components what temporal derivative of the corresponding
%               PDE state component they correspond to.
%
% OUTPUTS:
% Displays information in the command window on how the state components,
% inputs, and outputs in the PIE relate to those in the PDE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set a name associated to each object.
if strcmp(obj,'x')
    object_name = 'fundamental state component';
elseif strcmp(obj,'y')
    object_name = 'observed output';
elseif strcmp(obj,'z')
    object_name = 'regulated output';
elseif strcmp(obj,'u')
    object_name = 'actuator input';
elseif strcmp(obj,'w')
    object_name = 'exogenous input';
end
obj_tab = PIE.([obj,'_tab']);
ncomps = size(obj_tab,1);
nvars = size(PIE.vars,1);
has_diffs = false;
if all(obj_tab(:,1)==(1:ncomps)') && (size(obj_tab,2)<=2+nvars || ~any(any(obj_tab(:,3+nvars:2+2*nvars))))
    % The components in the PIE correspond exactly to the components in the
    % PDE.
    return
elseif size(obj_tab,2)<=2+nvars || ~any(any(obj_tab(:,3+nvars:2+2*nvars)))
    % Otherwise, we list the new order of the components.
    fprintf(['\n','The ',object_name,'s have been reindexed as:\n']);
else
    % Indicate how the new components relate to those in the PDE
    has_diffs = true;
    fprintf(['\n','The following ',object_name,'s have been introduced:\n']);  
end

% Use UNICODE to add subscript indices to different components.
thin_space = char(8201);
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
partial = '\x2202';
sub_t = '\x209C';
%sub_s = '\x209B';
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
global_vars_obj = PIE.vars(any(PIE.([obj,'_tab'])(:,3:2+nvars),1),:);
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
nvars_max = max(sum(PIE.([obj,'_tab'])(:,3:2+nvars),2));
lngth_varnames_mean = ceil(lngth_varnames*nvars_max/nvars);
LHS_length_max = 1+1 + n_digits + lngth_varnames_mean+3; % e.g. ' x13(t,s1,s2,s3)', 


% % For each of the components, display its size, and which variables it
% % depends on.
for ii=1:ncomps
    old_ID = obj_tab(ii,1);

    % Establish the names of the variables on which the component depends.
    varnames_ii_t = 't';
    LHS_length = 1;
    var_idcs = find(obj_tab(ii,3:2+PIE.dim));
    for kk=var_idcs
        varnames_ii_t = [varnames_ii_t,',',varname_cell{kk}];
        LHS_length = LHS_length+3;  % add 3 characters: ,s1
    end
    
    % Establish the (subscript) index for the new component.
    new_idx = ii;
    if ncomps==1
        Lcomp_idx = '';
    elseif ii<=9
        % The component number consists of a single decimal.
        Lcomp_idx = sub_num{new_idx+1};
        LHS_length = LHS_length + 1;
    else
        % The component number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(new_idx)')+1)');
        LHS_length = LHS_length + length(num2str(new_idx));
    end
    % Set the name of the component, including its dependence on spatial
    % variables.
    LHS_name = [' ',obj,Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = 1 + LHS_length + 3;
        
    % Establish the index for the old component
    if ncomps==1
        Rcomp_idx = '';
    elseif old_ID<=9
        % The component number consists of a single decimal.
        Rcomp_idx = sub_num{old_ID+1};
    else
        % The component number consists of multiple decimals.
        Rcomp_idx = cell2mat(sub_num(str2num(num2str(old_ID)')+1)');
    end
    % Set the name of the component, including its dependence on spatial
    % variables.
    RHS_name = [thin_space,obj,Rcomp_idx,'(',varnames_ii_t,')'];
    % For state components, also indicate what spatial/temporal derivative
    % of the PDE state they correspond to.
    if strcmp(obj,'x') && any(tdiff_list)
        tdiff = tdiff_list(ii);
        if tdiff==0
            diff_str = ['       ',thin_space];
        elseif tdiff==1
            diff_str = [' (',partial,'/',partial,'t)',thin_space];
        else
            t_sup = cell2mat(sup_num(str2num(num2str(tdiff)')+1)');
            diff_str = ['(',partial,'/',partial,'t)',t_sup,thin_space];
        end
    else
        diff_str = '';
    end
    if has_diffs
        Dvals = obj_tab(ii,2+nvars+1:2+2*nvars);
        for jj=1:nvars
            sdiff = Dvals(jj);
            if sdiff==0
                diff_str = ['        ',diff_str,thin_space];
            elseif sdiff==1
                diff_str = [diff_str,'(',partial,'/',partial,varname_cell{jj},')',thin_space];
            else
                s_sup = cell2mat(sup_num(str2num(num2str(sdiff)')+1)');
                diff_str = [diff_str,'(',partial,'/',partial,varname_cell{jj},')',s_sup,thin_space];
            end
        end
        RHS_name = [diff_str,RHS_name];
    end
    
    % % % Finally, display:
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),' <--  ',RHS_name,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
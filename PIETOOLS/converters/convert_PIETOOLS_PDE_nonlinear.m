function PIE = convert_PIETOOLS_PDE_nonlinear(PDE,comp_order,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_2D.m     PIETOOLS 2026
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
%           f, g, h: 'polyopvar' objects representing distributed
%           polynomials defining the right-hand side of the PIE for the
%           state, observed output, and regulated output, respectively.
%
%           Bu, Dyu, Dzu: opvar or opvar2d objects representing the
%           contribution of actuator inputs to the PIE and outputs.
%
%           Bw, Dyw, Dzw: opvar or opvar2d objects representing the
%           contribution of exogenous inputs to the PIE and outputs.
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
%       Tw*w(t) + Tu*u + T*xf(t) = f(xf(t)) + Bu  * u(t) + Bw  * w(t);
%                           y(t) = g(xf(t)) + Dyu * u(t) + Dyw * w(t);
%                           z(t) = h(xf(t)) + Dzu * u(t) + Dzw * w(t);
%
% - The order of the state, input, and output components in the PIE
%   structure may not be the same as that in the PDE structure. Check e.g.
%   PIE.x_tab(:,1) to see how the state components of the PDE have been
%   re-arranged.
%
% - If the PDE is 2D and involves boundary conditions with inputs u or
%   disturbances w that vary in space, the inputs and disturbance may be
%   differentiated with respect to the spatial variables as well in the PIE
%   representation. The function will warn the user of this fact, and the
%   order of the derivative taken with respect to each spatial is stored in
%   the outputs PIE.u_tab and PIE.w_tab, with e.g. the element
%   PIE.u_tab(i,5) indicating the order of derivative taken of the ith
%   input with respect to the first spatial variable, and
%   PIE.u_tab(i,6) indicating the order of derivative taken of the ith
%   input with respect to the second spatial variable.
%   It is also possible that boundary values of the input are used as
%   inputs in the PIE. In this case, multiple elements of PIE.w_tab(:,1)
%   (or PIE.u_tab(:,1)) will have the same value, indicating that they
%   correspond to the same input variable. To check if the ith input in
%   the PIE corresponds to a boundary value, compute
%       idx = find(PIE.w_tab(:,1)==PIE.w_tab(i,1),1,'last');
%   If idx~=i, then row i in PIE.w_tab corresponds to a boundary value of
%   the same input as specified by PDE.w_tab(idx,:). To check at what
%   boundary this input is evaluated, compute
%       is_bval = PIE.w_tab(idx,3:2+nvars) - PIE.w_tab(i,3:2+nvars);
%   where nvars = PIE.dim. For each j in 1 to nvars, if is_bval(j)=1, then
%   the input is evaluated at the lower-boundary PIE.dom(j,1) in the
%   spatial variable PIE.vars(j,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - convert_PIETOOLS_PDE_nonlinear
%
% Copyright (C)2026 PIETOOLS Team
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
% DJ, 02/24/2026: Initial coding;
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Initialization                                          % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pass PDEs in legacy format to legacy converters.
if ~isa(PDE,'pde_struct')
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
% If the PDE is just an ODE, we artificially augment.
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
% % % direction. Note that this approach may also require derivatives and
% % % boundary values of the inputs to be taken.
if (is2D && (any(any(w_tab(:,3:2+nvars))) || any(any(u_tab(:,3:2+nvars)))))
    old_dom = dom;      old_vars = vars;        np_op_old = np_op;
    [Top,Dvals] = Compute_Tmap_PDE_2D_Separable(PDE,comp_order,true);
    use_Tmap_opt2 = false;
    Dvals_w = Dvals.w;      Dvals_u = Dvals.u;
    Wop = Top.wmap;         Uop = Top.umap;                                 % DJ, 02/15/2026
    w_idcs = Dvals.w_idcs;  u_idcs = Dvals.u_idcs;

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
if Pop_f2BC==0 && (~any(np_op.u) || Pop_u2BC==0) && (~any(np_op.w) || Pop_w2BC==0)      % DJ, 01/25/2026
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
% inputs, nor does it take any boundary values of the inputs
Dvals_w = zeros(numel(PDE.w),2);
Dvals_u = zeros(numel(PDE.u),2);
Wop = 1;         
Uop = 1;
w_idcs = (1:numel(PDE.w));  
u_idcs = (1:numel(PDE.u));

end

if return_Top
    % If only the T operators are requested, we stop here.
    PIE = pie_struct();
    PIE.dom = dom;  PIE.vars = vars;
    PIE.T = Top.x;  PIE.Tw = Top.w;     PIE.Tu = Top.u;
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Building the PIE                                        % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Having derived the map from fundamental to PDE state,
% % %   x = Top.x*xf + Top.w*wf + Top.u*uf
% % % we now substitute this relation into the PDE,
% % % (d/dt)^k x = sum_j=1^n1 (P_x2x{j} * Delta_j * D_j)* x + P_u2x* u + P_w2x* w; 
% % %          y = sum_j=1^n2 (P_x2y{j} * Delta_j * D_j)* x + P_u2y* u + P_w2y* w; 
% % %          z = sum_j=1^n3 (P_x2z{j} * Delta_j * D_j)* x + P_u2z* u + P_w2z* w; 
% % % to obtain an equivalent PIE representation,
% % % (d/dt)^k x = fx(xf) + Bwop*w + Buop*u;
% % %          y = fy(xf) + Dywop*w + Dyuop*u;
% % %          z = fz(xf) + Dzwop*w + Dzuop*u;

% First, check which input signals appear as spatial derivatives in the PIE
% --> those cannot also appear without derivative
diff_idcs_w = find(any(Dvals_w,2));
is_diff_w = false(size(w_tab,1),1);
for j=1:size(w_tab,1)
    % Check if the input appears ONLY as derivative
    idx = find(w_idcs==j);
    if isscalar(idx) && any(Dvals_w(idx,:))
        is_diff_w(j) = true;
    end
end
diff_idcs_u = find(any(Dvals_u,2));
is_diff_u = false(size(u_tab,1),1);
for j=1:size(u_tab,1)
    % Check if the input appears ONLY as derivative
    idx = find(u_idcs==j);
    if isscalar(idx) && any(Dvals_u(idx,:))
        is_diff_u(j) = true;
    end
end

% Substitute the map from the fundamental to PDE state into the PDE and
% output equations
[fx, Bop_u2x, Bop_w2x] = compute_RHS_PIE(PDE,'x',Top,Wop,Uop,is_diff_w,is_diff_u,op_clean_tol);
[fy, Dop_u2y, Dop_w2y] = compute_RHS_PIE(PDE,'y',Top,Wop,Uop,is_diff_w,is_diff_u,op_clean_tol);
[fz, Dop_u2z, Dop_w2z] = compute_RHS_PIE(PDE,'z',Top,Wop,Uop,is_diff_w,is_diff_u,op_clean_tol);

% % % With that, we can equivalently represent the PDE as a PIE:
% % %   D_{t} * (Tw * w + Tu * u + Tx * xf) = fx(xf) + Bu  * u + Bw  * w ;
% % %                                    y  = fy(xf) + Dyu * u + Dyw * w ;
% % %                                    z  = fx(xf) + Dzu * u + Dzw * w ;
% % % where D_{t} is a diagonal temporal differential operator,
% % %   D_{t} = diag([(d/dt)^tdiff(1) , ... , (d/dt)^tdiff(nx)])
% % % We determine the orders of the temporal derivatives here
tdiff_list = ones(numel(PDE.x),1);
for ii=1:numel(PDE.x)
    if isfield(PDE.x{ii},'tdiff')
        tdiff_list(ii) = PDE.x{ii}.tdiff;
    end
end

% Finally, we define the PIE structure, collecting the PI operators that
% define our PIE.

% Keep track of the variables and domain of the system.
PIE = struct();
PIE.dom = dom;
PIE.vars = vars;
PIE.dim = size(PIE.vars,1);

% Store the PI operators defining the fundamental to PDE state map.
PIE.T = clean_opvar(Top.x,op_clean_tol);
PIE.Tw = clean_opvar(Top.w,op_clean_tol);
PIE.Tu = clean_opvar(Top.u,op_clean_tol);

% Store the PI operators defining the PIE.
PIE.f = fx;     PIE.B1 = Bop_w2x;       PIE.B2 = Bop_u2x;
PIE.g = fy;     PIE.D11 = Dop_w2z;      PIE.D12 = Dop_u2z;
PIE.h = fz;     PIE.D21 = Dop_w2y;      PIE.D22 = Dop_u2y;

% Keep track of how the state components in the original PDE are
% now ordered in the PIE.
x_tab(:,1) = xcomp_order;       PIE.x_tab = x_tab;
y_tab(:,1) = ycomp_order;       PIE.y_tab = y_tab;
z_tab(:,1) = zcomp_order;       PIE.z_tab = z_tab;
u_tab(:,1) = ucomp_order;       PIE.u_tab = u_tab;
w_tab(:,1) = wcomp_order;       PIE.w_tab = w_tab;

% Account for the introduction of boundary values of the inputs
PIE.w_tab = PIE.w_tab(w_idcs,:);
PIE.u_tab = PIE.u_tab(u_idcs,:);
% Indicate that the introduced boundary values do not vary in space
PIE.w_tab(1:numel(w_idcs)-numel(wcomp_order),3:2+nvars) = 0;
PIE.u_tab(1:numel(u_idcs)-numel(ucomp_order),3:2+nvars) = 0;

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

% Get rid of any higher-order temporal derivatives in the PIE.
if any(tdiff_list>1)
    error("Higher-order termporal derivatives in nonlinear PDEs are currently not supported.")
    %[PIE,tdiff_tab] = expand_tderivatives(PIE,tdiff_list,true);
else
    tdiff_tab = [(1:numel(PDE.x))',zeros(numel(PDE.x),1)];
end

% Finally, provide an overview of how the new state variables relate to the
% PDE variables.
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [fx, Pop_u, Pop_w] = compute_RHS_PIE(PDE,obj,Top,Wop,Uop,is_diff_w,is_diff_u,ztol)
% Construct distributed polynomial FX and PI operators POP_U, POP_W such 
% that the PDE corresponding to operators Pop,
%   obj = sum_j=1^n1 (Cop_x{j} * Delta_j * D_j)* x + Cop_u{j}* u + Cop_w{j}* w;
% can be represented in terms of the fundamental state as
%   obj = FX(xf) + POP_U * u + POP_W * w;
%
% INPUTS:
% - PDE:    A struct or pde_struct type object, defining a PDE in the terms
%           format (see also the "@pde_struct/initialize" function).
% - obj:    'x', 'y', or 'z', indicating which equations to be represented
%           using PI operators.
% - Top:    Struct with field 'x', 'u' and 'w', containing PI operators Tx,
%           Tu and Tw such that the PDE state x can be expressed in terms
%           of the fundamental state xf and inputs u and w as
%               x = Tx * xf + Tu * uf + Tw * w;
% - Wop:    'opvar' or 'opvar2d' object representing the map from the
%           fundamental input to the fundamental input to the PDE input;
% - Uop:    'opvar' or 'opvar2d' object representing the map from the
%           fundamental input to the fundamental input to the PDE input;
% - is_diff_w:  nw x 1 boolean array specifying for each of the exogenous
%               inputs in the PDE whether the PIE is expressed in terms of
%               a spatial derivative of this input. Currently, we do not
%               support inputs that appear both as the PDE input and as a
%               derivative in the PIE;
% - is_diff_u:  nu x 1 boolean array, specifying the same information as
%               'is_diff_w' but for actuator inputs;
%
% OUTPUTS:
% - fx: 'polyopvar' object representing the equations for "obj" in the PDE
%       in terms of the fundamental state variables (in the absence of
%       input signals);
% - Pop_u, Pop_w:   opvar2d objects describing the contributions of the
%                   inputs u and w to the PDE of "obj".
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 02/24/2026
%

% Extract some necessary information from the PDE.
global dom vars nvars np_op
x_tab = PDE.x_tab;
obj_tab = PDE.([obj,'_tab']);

% Extract the pre-defined PI operators.
Top_x = Top.x;
Top_u = Top.u;
Top_w = Top.w;

% Establish indices associated to the variables in each state component,
% input component, and output component.    
nx_arr = PDE.x_tab(:,2);    nnx_arr = cumsum([0;nx_arr]);
nnu_arr = cumsum([0;PDE.u_tab(:,2)]);
nnw_arr = cumsum([0;PDE.w_tab(:,2)]);
nnr_arr = cumsum([0;obj_tab(:,2)]);

% Establish input and output dimensions of the operators.
nx_op = np_op.x;        %nnx_op = cumsum([0;nx_op]);
nu_op = np_op.u;        nnu_op = cumsum([0;nu_op]);
nw_op = np_op.w;        nnw_op = cumsum([0;nw_op]);
nr_op = np_op.(obj);    nnr_op = cumsum([0;nr_op]);

% Initialize an empty distributed polynomial
fx = polyopvar();
varname = cell(size(x_tab,1),1);
pvarname = cell(nvars,1);
for j=1:size(x_tab,1)
    varname{j} = ['x',num2str(j)];
end
for j=1:nvars
    pvarname(j) = vars(j,1).varname;
end
fx.varname = varname;
fx.pvarname = pvarname;
fx.varmat = logical(x_tab(:,3:2+nvars));
fx.dom = dom;

% Initialize empty operators.
if nvars<=1
    opvar Pop_tmp Pop_u Pop_w;
    Pop_tmp.dim = [nr_op,nx_op];
    Pop_tmp.I = dom;
    Pop_tmp.var1 = vars(1,1);     Pop_tmp.var2 = vars(1,2);
    Pop_u.dim = [nr_op,nu_op];
    Pop_u.I = dom;
    Pop_u.var1 = vars(1,1);     Pop_u.var2 = vars(1,2);
    Pop_w.dim = [nr_op,nw_op];
    Pop_w.I = dom;
    Pop_w.var1 = vars(1,1);     Pop_w.var2 = vars(1,2);
else
    Pop_tmp = opvar2d([],[nr_op,nx_op],dom,vars);
    Pop_u = opvar2d([],[nr_op,nu_op],dom,vars);
    Pop_w = opvar2d([],[nr_op,nw_op],dom,vars);
end
if isa(Uop,'double')
    Pop_u2x = Pop_u;
else
    Pop_u2x = Pop_u;    Pop_u2x.dim(:,2) = Uop.dim(:,2);
end
if isa(Wop,'double')
    Pop_w2x = Pop_w;
else
    Pop_w2x = Pop_w;    Pop_w2x.dim(:,2) = Wop.dim(:,2);
end

% If the "obj" has no components, there is nothing to do
if nnr_op(end)==0    
    return
end

% Declare a matrix of degrees of monomials that appear in the PIE, and
% associated cell of operators
degmat = zeros(0,size(x_tab,1));
op_cell_full = cell(1,0);

% % % Loop over all equations of the considered "obj" (PDEs or outputs)
% % % For each equation, loop over all the terms, and add the coefficients
% % % that appear in each term to the appropriate parameter in the
% % % appropriate operator
for eqnum=1:numel(PDE.(obj))
        
    % Establish which row indices in the operators are linked to this 
    % equation.
    rindcs = nnr_arr(eqnum)+1:nnr_arr(eqnum+1);
    Cop_rnum = find((rindcs(1)>nnr_op(1:end-1) & rindcs(1)<=nnr_op(2:end)),1);

    % For each monomial, keep track of which term involving this monomial
    % we are considering for this equation.
    trm_cntr = ones(size(degmat,1),1);
    
    for jj=1:numel(PDE.(obj){eqnum}.term)
        % Extract the jth term
        term_jj = PDE.(obj){eqnum}.term{jj};
        
        if isfield(term_jj,'w')
            % Determine which input component is involved in the term.
            comp_indx = term_jj.w;
            if is_diff_w(comp_indx)
                error("An exogenous input will appear as both itself and as its own derivative in the PIE; this is currently not supported.")
            end
            % Establish the column numbers in the operator associated to
            % this input component.
            cindcs = nnw_arr(comp_indx)+1:nnw_arr(comp_indx+1);
            Cop_cnum = find((cindcs(1)>nnw_op(1:end-1) & cindcs(1)<=nnw_op(2:end)));
            
            % Build the coefficient operator acting on this term
            Idoms = term_jj.I;
            Cval = polynomial(term_jj.C);
            parnum = [Cop_rnum,Cop_cnum];
            Cop_jj = coeff2op(Cval,Idoms,parnum,vars,dom);
            Cop_jj = clean_opvar(Cop_jj,ztol);
            
            % Add the coefficients to the appropriate parameter.
            Pop_w(rindcs,cindcs) = Pop_w(rindcs,cindcs) + Cop_jj;
            
        elseif isfield(term_jj,'u')
            % Determine which input component is involved in the term.
            comp_indx = term_jj.u;
            if is_diff_u(comp_indx)
                error("An actuator input will appear as both itself and as its own derivative in the PIE; this is currently not supported.")
            end
            % Establish the column numbers in the operator associated to
            % this input component.
            cindcs = nnu_arr(comp_indx)+1:nnu_arr(comp_indx+1);
            Cop_cnum = find((cindcs(1)>nnu_op(1:end-1) & cindcs(1)<=nnu_op(2:end)));

            % Build the coefficient operator acting on this term
            Idoms = term_jj.I;
            Cval = polynomial(term_jj.C);
            parnum = [Cop_rnum,Cop_cnum];
            Cop_jj = coeff2op(Cval,Idoms,parnum,vars,dom);
            Cop_jj = clean_opvar(Cop_jj,ztol);
            
            % Add the coefficients to the appropriate parameter.
            Pop_u(rindcs,cindcs) = Pop_u(rindcs,cindcs) + Cop_jj;
            
        else
            % % The component involves a state variable, which may appear
            % % in a nonlinear manner

            % Determine the number of factors in the term
            nfctrs = numel(term_jj);
            state_nums_jj = zeros(1,0);
            op_cell_jj = cell(1,0);

            for fctr_num=1:nfctrs
                % Extract the factor from the term
                fctr_jj = term_jj(fctr_num);

                % Determine which PDE state component is involved in this 
                % factor, and associated rows of the T operator
                comp_indx = fctr_jj.x;
                r_idcs_tmp = nnx_arr(comp_indx)+1:nnx_arr(comp_indx+1);

                % Determine which variables the state component depends on.
                has_vars_xcomp = logical(PDE.x_tab(comp_indx,3:2+nvars));
                nvars_xcomp = sum(has_vars_xcomp);
                xvars = vars(has_vars_xcomp,1);
               
                % Establish what derivative of the state is considered.
                Dval = fctr_jj.D;
                
                % Establish at what position the state is evaluated.
                is_sub = false(1,nvars_xcomp); % 1 for substitution at boundary, 0 else
                Rloc = fctr_jj.loc;
                if ~isa(Rloc,'double') && ~isdouble(Rloc)
                    for kk=1:nvars_xcomp
                        if isdouble(Rloc(kk))
                            is_sub(kk) = true;
                        end
                    end
                elseif ~isempty(Rloc)
                    % If the state is evaluated at a corner, we can avoid a for
                    % loop.
                    is_sub = true(1,nvars_xcomp);
                end
                
                % Perform the desired differentiation/substitution of the
                % T operators mapping the fundamental state to PDE state
                Top_x_jj = Top_x(r_idcs_tmp,:);
                Top_w_jj = Top_w(r_idcs_tmp,:);
                Top_u_jj = Top_u(r_idcs_tmp,:);

                Top_x_jj = diff(Top_x_jj,xvars,Dval','pure');
                Top_w_jj = diff(Top_w_jj,xvars,Dval','pure');
                Top_u_jj = diff(Top_u_jj,xvars,Dval','pure');

                if any(is_sub)
                    Top_x_jj = subs(Top_x_jj,xvars(is_sub),Rloc(is_sub)','pure');
                    Top_w_jj = subs(Top_w_jj,xvars(is_sub),Rloc(is_sub)','pure');                
                    Top_u_jj = subs(Top_u_jj,xvars(is_sub),Rloc(is_sub)','pure');
                end

                % Declare the PI operator that is acting on this PDE state
                % variable
                Cop_cnum = find(Top_x_jj.dim(:,1)>0,1);
                Idoms = fctr_jj.I(~is_sub);
                Cval = polynomial(fctr_jj.C);
                parnum = [Cop_rnum,Cop_cnum];
                Cop_jj = coeff2op(Cval,Idoms,parnum,vars,dom);

                % Multiply coefficients with map from fundamental to PDE
                % state
                Pop_x2x_ll = clean_opvar(Cop_jj*Top_x_jj,ztol);
                Pop_u2x_ll = clean_opvar(Cop_jj*Top_u_jj,ztol);
                Pop_w2x_ll = clean_opvar(Cop_jj*Top_w_jj,ztol);

                % Add the input operators to the full operator
                if nfctrs>1 && (~all(all(Pop_u2x_ll==0)) || ~all(all(Pop_w2x_ll==0)))
                    error("An input signal appears nonlinearly in the PIE; this is not supported.")
                end
                Pop_u2x(rindcs,:) = Pop_u2x(rindcs,:) + Pop_u2x_ll;
                Pop_w2x(rindcs,:) = Pop_w2x(rindcs,:) + Pop_w2x_ll;

                % Augment x operator to map to all equations
                Pop_x2x_ll = [Pop_tmp(1:rindcs(1)-1,:); Pop_x2x_ll; Pop_tmp(rindcs(end)+1:end,:)];

                % Split the operator into columns associated with each
                % fundamental state variable
                state_nums_tmp = state_nums_jj;
                op_cell_tmp = op_cell_jj;
                nreps = size(state_nums_tmp,1);
                state_nums_jj = zeros(0,size(state_nums_tmp,2)+1);
                op_cell_jj = cell(0,size(op_cell_tmp,2)+1);
                for state_num=1:numel(nx_arr)
                    % Extract columns of the operator associated with
                    % fundamental state variable state_num
                    c_idcs = nnx_arr(state_num)+1:nnx_arr(state_num+1);
                    Pop_x_tmp = Pop_x2x_ll(:,c_idcs);
                    %Pop_x_tmp = dopvar2ndopvar(Pop_x2x_tmp(:,c_idcs));
                    if all(all(Pop_x_tmp==0))
                        % The current fundamental state variable does not
                        % contribute to this equation
                        continue
                    end
                    % Add the state number to the list of variables that
                    % contribute to this particular PDE term
                    state_nums_jj = [state_nums_jj; [state_nums_tmp,state_num*ones(nreps,1)]];
                    op_cell_jj = [op_cell_jj; [op_cell_tmp,repmat({Pop_x_tmp},nreps,1)]];
                end
            end
            % Add the terms to the PIE
            for ll=1:size(state_nums_jj,1)
                % Reorder factors based on order of states
                [state_nums_ll,new_order] = sort(state_nums_jj(ll,:));
                op_cell_ll = op_cell_jj(ll,new_order);
                % Determine the degree of the monomial in each state
                [~,strt_idx] = unique(state_nums_ll);
                degs_part = [strt_idx(2:end); numel(state_nums_ll)+1] - strt_idx;
                degs_ll = zeros(1,size(x_tab,1));
                degs_ll(state_nums_ll) = degs_part;
                % Determine whether the monomial already appears
                deg_idx = find(ismember(degmat,degs_ll,'rows'));
                if isempty(deg_idx)
                    % The monomial does not appear yet
                    % --> add it to the list
                    degmat = [degmat; degs_ll];
                    if sum(degs_ll)==1
                        op_cell_full = [op_cell_full,op_cell_ll];
                    else
                        op_cell_full = [op_cell_full,{op_cell_ll}];
                    end
                    trm_cntr = [trm_cntr; 2];
                else
                    % The monomial already appears
                    % --> add the new operator to the old one
                    if sum(degs_ll)==1
                        op_cell_full{deg_idx} = op_cell_full{deg_idx}+op_cell_ll{1};
                    else
                        trm_cntr_ll = trm_cntr(deg_idx);
                        nterms_ll = size(op_cell_full{deg_idx},1);
                        if trm_cntr_ll<=nterms_ll
                            % For terms belonging to a new equation, we can
                            % just add to terms already declared for
                            % previous equations
                            for fctr_num=1:size(op_cell_ll,2)
                                op_cell_full{deg_idx}{trm_cntr_ll,fctr_num} = op_cell_full{deg_idx}{trm_cntr_ll,fctr_num}+op_cell_ll{fctr_num};
                            end
                        else
                            % Otherwise, we concatenate to represent the
                            % newly added term
                            op_cell_full{deg_idx} = [op_cell_full{deg_idx};op_cell_ll];
                        end
                        trm_cntr(deg_idx) = trm_cntr_ll + 1;
                    end
                end
            end
        end
    end
end

% % Finally, declare the distributed polynomial based on the encountered
% % monomials and associated operators
% Declare the degrees of the monomials appearing in f
[degmat_f,deg_order] = sortrows_integerTable([sum(degmat,2),degmat(:,end:-1:1)]);
fx.degmat = degmat_f(:,end:-1:2);
% Declare the coefficients acting on these monomials
Cop = fx.C;
Cop.ops = op_cell_full(1,deg_order);
fx.C = Cop;

% Add the contribution of the inputs
Pop_u = Pop_u2x + Pop_u*Uop;
Pop_w = Pop_w2x + Pop_w*Uop;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function Cop = coeff2op(Cval,Idoms,parnum,vars,dom)
% COP = COEFF2OP(CVAL,IDOMS,PARNUM,VARS,DOM) constructs an operator 
% COP representing the integral over the domains IDOMS with kernel CVAL
%
% INPUTS
% - Cval:   m x n 'double' or polynomial object representing a multiplier
%           or a kernel in an integral operator acting on x;
% - Idoms:  q x 1 cell specifying the domains over which to integrate x
%           with respect to each of the q variables on which x depends. If
%           Idoms{i} = [], then no integration is performed along the ith
%           variable;
% - parnum: 1 x 2 array of integers specifying which parameter in the
%           operator Cop will need to be set. This depends on what function
%           space Cop maps between: 
%               parnum(1) = 1 implies a map to R
%               parnum(1) = 2 implies a map to L2[s1]
%               parnum(1) = 3 implies a map to L2[s2]
%               parnum(1) = 4 implies a map to L2[s1,s2]
%           and similarly, parnum(2) specifies what space Cop maps from;
% - vars:   N x 2 'pvar' array specifying the primary and dummy variables
%           in the operator;
% - dom:    N x 2 array specifying for each of the N primary variables on
%           what interval it is defined;
%
% OUTPUTS
% - Cop:    'opvar' or 'opvar2d' object specifying the map defined by the
%           coefficients 'Cval' and domains of integration 'Idoms';

% Initialize an operator of desired dimensions
nvars = size(vars,1);
if nvars==1
    Cop = opvar();
    Rparams = {'P', 'Q1';
               'Q2', 'R'};
else
    Cop = opvar2d();
    Rparams = {'R00', 'R0x', 'R0y', 'R02';
               'Rx0', 'Rxx', 'Rxy', 'Rx2';
               'Ry0', 'Ryx', 'Ryy', 'Ry2';
               'R20', 'R2x', 'R2y', 'R22'};
end
Cop.dim(parnum(1),1) = size(Cval,1);
Cop.dim(parnum(2),2) = size(Cval,2);
Cop.var1 = vars(:,1);
Cop.var2 = vars(:,2);
Cop.I = dom;

% Establish between what function spaces the operator maps
Rname = Rparams{parnum(1),parnum(2)};
[has_vars_Lcomp{1:nvars+1}] = ind2sub([2*ones(1,nvars),1],parnum(1));
has_vars_Lcomp = logical(cell2mat(has_vars_Lcomp(1:nvars))-1);
[has_vars_Rcomp{1:nvars+1}] = ind2sub([2*ones(1,nvars),1],parnum(2));
has_vars_Rcomp = logical(cell2mat(has_vars_Rcomp(1:nvars))-1);
nvars_Rcomp = sum(has_vars_Rcomp);

% Check whether the coefficients act as multiplier or as
% kernel of some integral
Cop_int_indx = ones(1,nvars);
if isempty(Idoms)
    % The coefficients act as just a multiplier operator
    Cval = subs(Cval,vars(:,2),vars(:,1));
else
    % The coefficients may represent kernels
    % Check first along which directions we can have
    % partial integrals
    param_sz_full = 3.^(has_vars_Lcomp & has_vars_Rcomp);
    param_sz = param_sz_full(has_vars_Rcomp);
    % For each direction along which integration is possible,
    % establish whether an integral over _a^s, _s^b or _a^b is taken
    Cop_int_indx_tmp = ones(1,nvars_Rcomp);
    for ll=1:numel(Idoms)
        if param_sz(ll)==1
            % An integral is taken over the full domain
            % anyway --> replace dummy vars in the
            % kernel with primary vars.
            vars_ll = vars(has_vars_Rcomp,:);
            vars_ll = vars_ll(ll,:);
            Cval = subs(Cval,vars_ll(2),vars_ll(1));
        elseif ~isempty(Idoms{ll})
            if isa(Idoms{ll},'double') || isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll})
                % Integrating over full spatial domain.
                % A full integral has to be constructed
                % using two partial integrals.
                Cop_int_indx1 = Cop_int_indx_tmp;
                Cop_int_indx2 = Cop_int_indx_tmp;
                Cop_int_indx1(:,ll) = 2;
                Cop_int_indx2(:,ll) = 3;
                Cop_int_indx_tmp = [Cop_int_indx1; Cop_int_indx2];
            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(1))
                % Integrating over lower "half" of domain.
                Cop_int_indx_tmp(:,ll) = 2;
            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(2))
                % Integrating over upper "half" of domain.
                Cop_int_indx_tmp(:,ll) = 3;
            end
        end
    end
    Cop_int_indx(:,has_vars_Rcomp) = Cop_int_indx_tmp;
end

% Declare coefficients in appropriate element of Cop
if nvars==1
    % Set elements of 4-PI operator.
    if parnum(1)==1 || parnum(2)==1
        % Map to or from R, no need to account for partial integrals
        Cop.(Rname) = Cval;
    else
        % Map L2-->L2, set appropriate element depending on integral
        for ll=1:size(Cop_int_indx,1)
            switch Cop_int_indx(ll)
                case 1
                    Cop.R.R0 = Cval;
                case 2
                    Cop.R.R1 = Cval;
                case 3
                    Cop.R.R2 = Cval;
            end
        end
    end
else
    % Set elements of 2D PI operator
    Cparam = Cop.(Rname);
    if ~isa(Cparam,'cell')
        Cparam = polynomial(Cval);
    else
        sz_param = size(Cparam);
        for ll=1:size(Cop_int_indx,1)
            lindx = cumprod([1,sz_param(1:end-1)])*(Cop_int_indx(ll,:)-1)' + 1;
            Cparam{lindx} = polynomial(Cval);
        end
    end
    Cop.(Rname) = Cparam;
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



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
    % Check if this component is derived from some higher-dimensional
    % component in the original PDE (e.g. a boundary value of some input)
    comp_idx = find(obj_tab(:,1)==old_ID,1,'last');

    % Establish the names of the variables on which the component depends.
    varnames_LHS_t = 't';
    LHS_length = 1;
    var_idcs = find(obj_tab(ii,3:2+PIE.dim));
    for kk=var_idcs
        varnames_LHS_t = [varnames_LHS_t,',',varname_cell{kk}];
        LHS_length = LHS_length+3;  % add 3 characters: ,s1
    end
    if ii~=comp_idx
        % In case the current component is derived from some other
        % component, it may not depend on the same variables
        var_idcs_RHS = find(obj_tab(comp_idx,3:2+PIE.dim));
        varnames_RHS_t = 't';
        for kk=var_idcs_RHS
            is_var = ismember(kk,var_idcs);
            if is_var
                % The new component still depends on this variable
                varnames_RHS_t = [varnames_RHS_t,',',varname_cell{kk}];
            else
                % The component is evaluated at the lower boundary w/r to
                % this spatial variable
                varnames_RHS_t = [varnames_RHS_t,',',num2str(PIE.dom(kk,1))];
            end
        end
    else
        varnames_RHS_t = varnames_LHS_t;
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
    LHS_name = [' ',obj,Lcomp_idx,'(',varnames_LHS_t,')'];
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
    RHS_name = [thin_space,obj,Rcomp_idx,'(',varnames_RHS_t,')'];
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
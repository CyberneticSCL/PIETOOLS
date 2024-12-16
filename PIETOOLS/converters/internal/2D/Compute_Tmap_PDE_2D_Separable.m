function [Top,Twop,Tuop,Dvals_w,Dvals_u] = Compute_Tmap_PDE_2D_Separable(PDE,comp_order)
% Compute opvar2d objects Top, Twop, Tuop defining the map from the
% fundamental state to the PDE state
%   v = Top*Dop*v + Twop*Dwop*w + Twop*Duop*u
% where v is the PDE state, Dop*v the fundamental state, w the exogenous
% input, and u the actuator input. Operators Dwop and Duop are diagonal
% differential operators, with the order of the derivative specified by
% Dvals_w and Dvals_u.
% 
% INPUTS:
% - PDE:        A struct or pde_struct object defining a PDE in the terms
%               format (see also "@pde_struct/initialize").
% - comp_order: Optional input specifying the adjusted order of components
%               after running "reorder_comps" for the PDE. Since opvar2d
%               objects map RxL2[x]xL2[y]xL2[x,y], and cannot map e.g.
%               RxL2[y]xL2[x,y]xL2[x], the components of the PDE must be
%               reordered to ensure the opvar2d object can act on them. If
%               no argument "comp_order" is passed, the function assumes no
%               re-ordering has been perform, and will perform re-ordering
%               itself.
%
% OUTPUTS:
% - Top:        opvar2d object mapping fundamental state to PDE state.
% - Twop:       opvar2d object mapping exogenous inputs to PDE state.
% - Tuop:       opvar2d object mapping actuator inputs to PDE state.
% - Dvals_w:    nw x 2 array of type 'double' specifying for each of the nw
%               exogenous inputs that appear in the (reordered) PDE to what
%               order they are differentiated with respect to each spatial
%               variable in the map back to the PDE state.
% - Dvals_u:    nu x 2 array of type 'double' specifying for each of the nu
%               actuator inputs that appear in the (reordered) PDE to what
%               order they are differentiated with respect to each spatial
%               variable in the map back to the PDE state.
%
% NOTES:
% - This routine only works for PDEs with "separable" boundary conditions.
%   That is, there must be no ambiguity whether a boundary condiiton
%   corresponds to a condition along the x-direction or the y-direction,
%   and we should be able to (easily) separate these conditions. 
%   For example, standard Robin-type boundary conditions are accepted:
%       v(a,y)=0, v_{x}(b,y)=0, v(x,c)=w1(x), v_{y}(x,d)+v(x,d)=w2(x)
%
% - The routine DOES NOT warn for potential conflicts in the  boundary
%   conditions the produced relation 
%       v = Top*D_xy*v + Twop*D_w*w + Tuop*D_u*u
%   holds only if the boundary conditions do not conflict. For example,
%   if we enforce
%       v(a,y)=0,   v(x,c)=w1(x),
%   then we must have w1(a)=0 to avoid a conflict v(a,c)=0\neq w1(a).
%
% - For a state variable v differentiable up to order i in x and order j in 
%   y, we must have i boundary conditions along the x-direction, and j
%   boundary conditions along the y-direction.
%
% - The 2D map is constructed by first expanding only along each spatial
%   dimension separately, as
%       v = Top_x*D_x*v +Twop_x*w +Tuop_x*u;
%       v = Top_y*D_y*v +Twop_y*w +Tuop_y*u;
%   Then, we combine these expressions to expand v in terms of its
%   highest-order derivative with respect to both spatial variables as
%       D_y*v = Top_x*D_xy*v +Twop_x*D_y*w +Tuop_x*D_y*u
%   --> v = Top_y*Top_x*D_xy*v +Top_y*Twop_x*D_y*w +Twop_y*w 
%               +Top_y*Tuop_x*D_y*u +Tuop_y*u;
%   --> v = Top*D_xy*v + Twop*D_w*w + Tuop*D_u*u;
%
% - Certain inputs w and u may appear as derivatives in the map from the
%   fundamental state to the PDE state. The routine should provide a
%   warning of which inputs have been differentiated and reordered, and the
%   outputs "Dvals_w" and "Dvals_u" specify the orders of the derivatives
%   taken of each input IN THE NEW ORDER with respect to each spatial
%   variable.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  PIETOOLS Team
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
% Initial coding DJ - 04/15/2024
% DJ, 12/16/2024: Fix call to 2D converter;
%

% % % Extract the inputs
% Make sure the specified PDE is initialized.
PDE = pde_struct(PDE);
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end

if nargin==1
    % If no re-ordering of the components in the PDE has been performed
    % yet, perform this re-ordering and keep track of the new order.
    [PDE,comp_order] = reorder_comps(PDE,'all',false);
end

% % % Check that the BCs of the PDE can be separated
% % % Limited tests so far, further tests may be necessary

% % Make sure that no 2D state variable appears in a finite-dimensional BC
% Determine indices of 2D state components
is_2D_x = PDE.x_tab(:,3) & PDE.x_tab(:,4);
indcs_2D_x = find(is_2D_x);

% Determine indices of finite-dimensional BCs
is_0D_BC = ~PDE.BC_tab(:,3) & ~PDE.BC_tab(:,4);

% Check for appearance of 2D state in finite-dimensional BC
exit_function = false;
if any(is_0D_BC)
for ii=find(is_0D_BC)'
    for jj=1:numel(PDE.BC{ii}.term)
        if isfield(PDE.BC{ii}.term{jj},'x') && ismember(PDE.BC{ii}.term{jj}.x,indcs_2D_x)
            exit_function = true;
            break
        end
    end
    if exit_function
        break
    end
end
end
if exit_function
    error('The input PDE has mixed boundary conditions, which is not supported.')
end


% % % Build the operators

% % First, we expand the state v only along 1 spatial direction:

% Compute Top along x-direction such that
%   v = Top_x*D_x*v +Twop_x*w +Tuop_x*u;
% for a suitable differential operator D_x along x-direction
[Top_x,Twop_x,Tuop_x] = construct_1D_Tops(PDE,1);

% Compute Top along y-direction such that
%   v = Top_y*D_y*v +Twop_y*w +Tuop_y*u;
% for a suitable differential operator D_y along y-direction
[Top_y,Twop_y,Tuop_y] = construct_1D_Tops(PDE,2);

% % Then, we combine our expansions, using one of two options
% %
% % Dy option:
% %     Dop_y*v = Dop_y*Top_x*Dop_x*v +Twop_x*Dop_y*w   +...
% % --> v = Top_y*Tnop_x*Dop_xy*v +Top_y*Twop_x*Dop_y*w +Twop_y*w    +...;
% % --> v = Top*Dop_xy*v + Twop*Dop_w*w    +...;
% %
% % Dx option:
% %     Dop_x*v = Dop_x*Top_y*Dop_y*v +Twop_y*Dop_x*w   +....
% % --> v = Top_x*Tnop_y*Dop_xy*v +Top_x*Twop_y*Dop_x*w +Twop_x*w    +...;
% % --> v = Top*Dop_xy*v + Twop*Dop_w*w    +...;
% %
% % for suitable integral operators Tnop_x and Tnop_y and
% % for suitable differential operators D_xy=D_x*D_y and D_w, D_u;
% % We check if either of these options is viable, and determine which
% % derivatives must be taken.

% Determine which derivatives are taken of the state components
Dvals_xy = PDE.x_tab(:,end-1:end);

% Determine whether we can express
%   Dop_y*(Top_x*Dop_x*v) = Tnop_x*Dop_xy*v;
% or
%   Dop_x*(Top_y*Dop_y*v) = Tnop_y*Dop_xy*v;
% for some PI operator Tnop_x or Tnop_y. In particular, note that e.g.
%   Dop_y*(Top_x*v) = (Dop_y*Top_x)*v + Top_x*(Dop_y*v);
% We need the operator Dop_y*Top_x to be zero, so that the expression is
% solely in terms of Dop_y*v;
% Since e.g. Top_x acts only along the x-direction, we may assume that
% Dop_y*Top_x = 0. However, we can check if we want to be sure:
check_Dv = false;
if check_Dv
    [Tnop_x,Tnop_y,use_Dx_v,use_Dy_v] = check_state_diffs(PDE,Top_x,Top_y,Dvals_xy)
else
    Tnop_x = Top_x;     Tnop_y = Top_y;
    use_Dx_v = false;   use_Dy_v = false;
end

% Determine which derivatives would need to be taken of the inputs w
if numel(PDE.w)>=1
    % Compute contribtuion of exogenous inputs.
    [Dvals_w,use_Dx_w,use_Dy_w] = compute_input_diffs(Twop_x,Twop_y,Dvals_xy,PDE.x_tab(:,2),PDE.w_tab(:,2));
else
    Dvals_w = zeros(0,2);
    use_Dx_w = false;       use_Dy_w = false;
end

% Determine which derivatives would need to be taken of the inputs u
if numel(PDE.u)>=1
    % Compute contribtuion of actuator inputs.
    [Dvals_u,use_Dx_u,use_Dy_u] = compute_input_diffs(Tuop_x,Tuop_y,Dvals_xy,PDE.x_tab(:,2),PDE.u_tab(:,2));
else
    Dvals_u = zeros(0,2);
    use_Dx_u = false;       use_Dy_u = false;
end

% Check if we are required to use the Dx option or Dy option
use_Dx = use_Dx_v || use_Dx_w || use_Dx_u;
use_Dy = use_Dy_v || use_Dy_w || use_Dy_u;

% We can only choose 1!
if use_Dx && use_Dy
    error('The provided boundary conditions require multiple different derivative of certain inputs to be taken; this is not currently supported.')
end

% If both Dx and Dy option are viable, choose option that requires most
% derivatives to be taken
if ~use_Dx && ~use_Dy
    if nnz(Dvals_w(:,1))+nnz(Dvals_u(:,1))<nnz(Dvals_w(:,2))+nnz(Dvals_u(:,2))
        use_Dx = false;
        Dvals_w = [0,1].*Dvals_w;
        Dvals_u = [0,1].*Dvals_u;
    else
        use_Dx = true;
        Dvals_w = [1,0].*Dvals_w;
        Dvals_u = [1,0].*Dvals_u;
    end
end

% % Having determined which option to proceed with, we finally compute
% % the actual operators Top, Twop, and Tuop
if use_Dx
    % % Proceed with the expansion
    % % v = Top_x*Tnop_y*Dop_xy*v +Top_x*Twop_y*Dop_x*w +Twop_x*w
    Top = Top_x*Tnop_y;
    if numel(PDE.w)>=1
        Twop = Top_x*Twop_y +Twop_x;
        print_reorder_diff_summary(PDE,'w',comp_order.w,Dvals_w);
    else
        Twop = Twop_x;
    end
    if numel(PDE.u)>=1
        Tuop = Top_x*Tuop_y +Tuop_x;
        print_reorder_diff_summary(PDE,'u',comp_order.u,Dvals_u);
    else
        Tuop = Tuop_x;
    end
else
    % % Proceed with the expansion
    % % v = Top_y*Tnop_x*Dop_xy*v +Top_y*Twop_x*Dop_x*w +Twop_y*w
    Top = Top_y*Tnop_x;
    if numel(PDE.w)>=1
        Twop = Top_y*Twop_x +Twop_y;
        print_reorder_diff_summary(PDE,'w',comp_order.w,Dvals_w);
    else
        Twop = Twop_y;
    end
    if numel(PDE.u)>=1
        Tuop = Top_y*Tuop_x +Tuop_y;
        print_reorder_diff_summary(PDE,'u',comp_order.u,Dvals_u);
    else
        Tuop = Tuop_y;
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [Top,Twop,Tuop] = construct_1D_Tops(PDE,var_idx)
% % Build a 1D PDE along just the spatial direction specified by var_idx
% % satisfying the same BCs as specified in the PDE
alt_idx = setdiff([1,2],var_idx);


% % Check which states and inputs depend on the alternative variable
isvar_state = PDE.x_tab(:,2+var_idx);
is2D_state = PDE.x_tab(:,3) & PDE.x_tab(:,4);
isalt_state = PDE.x_tab(:,2+alt_idx) & ~PDE.x_tab(:,2+var_idx);
is2D_w = PDE.w_tab(:,3) & PDE.w_tab(:,4);
isalt_w = PDE.w_tab(:,2+alt_idx) & ~PDE.w_tab(:,2+var_idx);
is2D_u = PDE.u_tab(:,3) & PDE.u_tab(:,4);
isalt_u = PDE.u_tab(:,2+alt_idx) & ~PDE.u_tab(:,2+var_idx);

% % Initialize the 1D PDE with the same state components and inputs as
% % the 2D PDE, but get rid of dependence on alt_var
pde_struct PDE_1D;
PDE_1D.x = PDE.x;
for ii=1:length(PDE.x)
    if is2D_state(ii)
        % Treat 2D state components as 1D components along just var_idx
        PDE_1D.x{ii}.dom = PDE.x{ii}.dom(var_idx,:);
        PDE_1D.x{ii}.vars = PDE.x{ii}.vars(var_idx,:);
        PDE_1D.x{ii}.diff = PDE.x{ii}.diff(var_idx);
    elseif isalt_state(ii)
        % Treat state components along alt_idx as ODE state components
        PDE_1D.x{ii}.dom = [];
        PDE_1D.x{ii}.vars = [];
        PDE_1D.x{ii}.diff = [];
    end
end
PDE_1D.x_tab = PDE.x_tab(:,[1,2,2+var_idx,4+var_idx]);

% Disturbances
PDE_1D.w = PDE.w;
for ii=1:length(PDE.w)
    if is2D_w(ii)
        % Treat 2D disturbances as 1D disturbances along just var_idx
        PDE_1D.w{ii}.dom = PDE_1D.w{ii}.dom(var_idx,:);
        PDE_1D.w{ii}.vars = PDE_1D.w{ii}.vars(var_idx,:);
    elseif isalt_w(ii)
        % Treat 1D disturbances as finite-dimensional disturbances
        PDE_1D.w{ii}.dom = [];
        PDE_1D.w{ii}.vars = [];
    end
end
PDE_1D.w_tab = PDE.w_tab(:,[1,2,2+var_idx]);

% Controlled inputs
PDE_1D.u = PDE.u;
for ii=1:length(PDE.u)
    if is2D_u(ii)
        % Treat 2D inputs as 1D inputs along just var_idx
        PDE_1D.u{ii}.dom = PDE_1D.u{ii}.dom(var_idx,:);
        PDE_1D.u{ii}.vars = PDE_1D.u{ii}.vars(var_idx,:);
    elseif isalt_u(ii)
        % Treat 1D inputs as finite-dimensional disturbances
        PDE_1D.u{ii}.dom = [];
        PDE_1D.u{ii}.vars = [];
    end
end
PDE_1D.u_tab = PDE.u_tab(:,[1,2,2+var_idx]);

% % Declare the 1D PDE along the direction of var_idx
% We don't care about the actual PDE dynamics, just make sure the state is
% differentiated up to the desired value.
for ii=1:length(PDE_1D.x)
    PDE_1D.x{ii} = rmfield(PDE_1D.x{ii},'term');    % Get rid of old terms
    PDE_1D.x{ii}.term{1}.x = ii;                    % Introduce new term
    PDE_1D.x{ii}.term{1}.D = PDE_1D.x{ii}.diff;     % Introduce new term
end

% Reduce the BCs to just BCs along var_idx.
BC_idcs = ~PDE.BC_tab(:,2+var_idx);
PDE_1D.BC = PDE.BC(BC_idcs);
% Keep track of which BCs are only on states that depend on alt_var;
is_altvar_BC = true(length(PDE_1D.BC),1);
for ii=1:length(PDE_1D.BC)
    % Boundary conditions in 1D cannot vary in space.
    PDE_1D.BC{ii}.vars = [];     PDE_1D.BC{ii}.dom = [];
    for jj=1:length(PDE_1D.BC{ii}.term)
        if isfield(PDE_1D.BC{ii}.term{jj},'x')
            if is2D_state(PDE_1D.BC{ii}.term{jj}.x)
                if PDE_1D.BC{ii}.term{jj}.D(alt_idx)
                    error('The boundary conditions of the specified PDE are not supported.')
                end
                % Treat 2D state as 1D state along just var_idx
                PDE_1D.BC{ii}.term{jj}.D = PDE_1D.BC{ii}.term{jj}.D(var_idx);
                PDE_1D.BC{ii}.term{jj}.loc = PDE_1D.BC{ii}.term{jj}.loc(var_idx);
                PDE_1D.BC{ii}.term{jj}.I = PDE_1D.BC{ii}.term{jj}.I(var_idx);
            elseif isalt_state(PDE_1D.BC{ii}.term{jj}.x)
                % Treat 1D state along alt_idx as ODE state
                PDE_1D.BC{ii}.term{jj}.D = [];
                PDE_1D.BC{ii}.term{jj}.loc = [];
                PDE_1D.BC{ii}.term{jj}.I = {};
            end
            % Check if BC is indeed a BC on some state that depends on the
            % considered variable.
            is_altvar_BC(ii) = is_altvar_BC(ii) && ~isvar_state(PDE_1D.BC{ii}.term{jj}.x);
        elseif isfield(PDE_1D.BC{ii}.term{jj},'w')
            if is2D_w(PDE_1D.BC{ii}.term{jj}.w)
                % Treat 2D disturbance as 1D disturbance along just var_idx
                PDE_1D.BC{ii}.term{jj}.loc = PDE_1D.BC{ii}.term{jj}.loc(var_idx);
                PDE_1D.BC{ii}.term{jj}.I = PDE_1D.BC{ii}.term{jj}.I(var_idx);
            elseif isalt_w(PDE_1D.BC{ii}.term{jj}.w)
                % Treat 1D diturbance along alt_idx as finite-dimensional
                PDE_1D.BC{ii}.term{jj}.loc = [];
                PDE_1D.BC{ii}.term{jj}.I = {};
            end
        else
            if is2D_u(PDE_1D.BC{ii}.term{jj}.u)
                % Treat 2D disturbance as 1D disturbance along just var_idx
                PDE_1D.BC{ii}.term{jj}.loc = PDE_1D.BC{ii}.term{jj}.loc(var_idx);
                PDE_1D.BC{ii}.term{jj}.I = PDE_1D.BC{ii}.term{jj}.I(var_idx);
            elseif isalt_u(PDE_1D.BC{ii}.term{jj}.u)
                % Treat 1D diturbance along alt_idx as finite-dimensional
                PDE_1D.BC{ii}.term{jj}.loc = [];
                PDE_1D.BC{ii}.term{jj}.I = {};
            end
        end
    end
end
% Get rid of BCs that are only on states that depend on alt_var, since we
% cannot enforce them.
PDE_1D.BC = PDE_1D.BC(~is_altvar_BC);

% % Compute the 1D PIE representation associated to the 1D PDE
PDE_1D = initialize(PDE_1D,true);
[PDE_1D,comp_order] = reorder_comps(PDE_1D,'all',true);
try PIE_1D = convert_PIETOOLS_PDE(PDE_1D);                                  % DJ, 12/16/2024
catch
    error('The expansion of the PDE state in terms of the PIE state fails; perhaps the PDE is ill-posed, or not currently supported.')
end

% Extract the operators Top such that
% v(x,t) = Top*d_{var_idx}^{diff} v +Twop*w +Tuop*u;
Top_1D = PIE_1D.T;        % Note that these operators are only 1D opvars!
Twop_1D = PIE_1D.Tw;
Tuop_1D = PIE_1D.Tu;



% % % Convert the opvar object back to an opvar2d object
% % Initialize the operators
opvar2d Top Twop Tuop;
Top.I = PDE.dom;    Twop.I = PDE.dom;       Tuop.I = PDE.dom;
Top.var1 = PDE.vars(:,1);   Top.var2 = PDE.vars(:,2);
Twop.var1 = PDE.vars(:,1);  Twop.var2 = PDE.vars(:,2);
Tuop.var1 = PDE.vars(:,1);  Tuop.var2 = PDE.vars(:,2);

vars_new = PDE.vars(var_idx,:)';        vars_old = PDE_1D.vars(1,:)';

% Set the dimensions of the operators
nv_op_2D = get_opdim(PDE.x_tab);        Top.dim = [nv_op_2D,nv_op_2D];
nw_op_2D = get_opdim(PDE.w_tab);        Twop.dim = [nv_op_2D,nw_op_2D];
nu_op_2D = get_opdim(PDE.u_tab);        Tuop.dim = [nv_op_2D,nu_op_2D];

% Determine which state and input variables appear in which parameter
nnv_op_2D = [0;cumsum(nv_op_2D)];
nnw_op_2D = [0;cumsum(nw_op_2D)];
nnu_op_2D = [0;cumsum(nu_op_2D)];


% Determine which component in the 1D PDE corresponds to each component in
% the 2D PIE, e.g.
%   v_2D(k) = v_1D(idcs_v(k))   and  v_1D(k) = v_2D(comp_order.x(k))
[~,comp_idcs_v] = sort(comp_order.x);
[~,comp_idcs_w] = sort(comp_order.w);
[~,comp_idcs_u] = sort(comp_order.u);

% Since each state and input component may consist of multiple variables,
% we have to account for the sizes of these components in relating 
% different rows and columns of the 2D operators to those of the 1D 
% operators.
% We start with v,
%   column k in Top_2D should correspond to column idcs_v(k) in Top_1D
csz_v = cumsum([0;PDE_1D.x_tab(:,2)]);
strt_idcs_v = csz_v(1:end-1)+1;     end_idcs_v = csz_v(2:end);
idcs_v = zeros(csz_v(end),1);       strt_k = 0;
for k=1:length(comp_idcs_v)
    v_idx = comp_idcs_v(k);
    idcs_v(strt_k+(1:PDE_1D.x_tab(v_idx,2))') = (strt_idcs_v(v_idx):end_idcs_v(v_idx))';
    strt_k = strt_k +PDE_1D.x_tab(v_idx,2);
end
% Next for w,
%   column k in Twop_2D should correspond to column idcs_v(k) in Twop_1D
csz_w = cumsum([0;PDE_1D.w_tab(:,2)]);
strt_idcs_w = csz_w(1:end-1)+1;     end_idcs_w = csz_w(2:end);
idcs_w = zeros(csz_w(end),1);       strt_k = 0;
for k=1:length(comp_idcs_w)
    w_idx = comp_idcs_w(k);
    idcs_w(strt_k+(1:PDE_1D.w_tab(w_idx,2))') = (strt_idcs_w(w_idx):end_idcs_w(w_idx))';
    strt_k = strt_k +PDE_1D.w_tab(w_idx,2);
end
% Finally for u,
%   column k in Tuop_2D should correspond to column idcs_v(k) in Tuop_1D
csz_u = cumsum([0;PDE_1D.u_tab(:,2)]);
strt_idcs_u = csz_u(1:end-1)+1;     end_idcs_u = csz_u(2:end);
idcs_u = zeros(csz_u(end),1);       strt_k = 0;
for k=1:length(comp_idcs_u)
    u_idx = comp_idcs_u(k);
    idcs_u(strt_k+(1:PDE_1D.u_tab(u_idx,2))') = (strt_idcs_u(u_idx):end_idcs_u(u_idx))';
    strt_k = strt_k +PDE_1D.u_tab(u_idx,2);
end

% [~,idcs_v] = sort(comp_order.x);        nv_op_1D = Top_1D.dim(:,2);
% [~,idcs_w] = sort(comp_order.w);        nw_op_1D = Twop_1D.dim(:,2);
% [~,idcs_u] = sort(comp_order.u);        nu_op_1D = Tuop_1D.dim(:,2);

% Since rows/columns of the operators are divided over different
% parameters in the opvar structure, we adjust the indices to account for 
% which parameter they appear in
nv_op_1D = Top_1D.dim(:,2);     
nw_op_1D = Twop_1D.dim(:,2);
nu_op_1D = Tuop_1D.dim(:,2);

idcs_v(idcs_v>nv_op_1D(1)) = idcs_v(idcs_v>nv_op_1D(1)) - nv_op_1D(1);
idcs_w(idcs_w>nw_op_1D(1)) = idcs_w(idcs_w>nw_op_1D(1)) - nw_op_1D(1);
idcs_u(idcs_u>nu_op_1D(1)) = idcs_u(idcs_u>nu_op_1D(1)) - nu_op_1D(1);

% Determine which parameter in the 1D opvar corresponds to each parameter
% corresponds to each parameter in the opvar2d
fnames_2D = {'R00','R0x','R0y','R02';
             'Rx0','Rxx','Rxy','Rx2';
             'Ry0','Ryx','Ryy','Ry2';
             'R20','R2x','R2y','R22'};
if var_idx==1
    fnames_1D = {'P','Q1','P','Q1';
                'Q2','R','Q2','R';
                'P','Q1','P','Q1';
                'Q2','R','Q2','R'};
else
    fnames_1D = {'P','P','Q1','Q1';
                'P','P','Q1','Q1';
                'Q2','Q2','R','R';
                'Q2','Q2','R','R'};
end

% Set the value of each of the parameters in Top, Twop, and Tuop
for k=1:numel(fnames_2D)
    % Determine which rows and columns of Top the parameter corresponds to
    [rnum,cnum] = ind2sub(size(fnames_2D),k);

    % Extract parameter in 1D operators corresponding to parameter in 2D
    v_param = Top_1D.(fnames_1D{k});
    w_param = Twop_1D.(fnames_1D{k});
    u_param = Tuop_1D.(fnames_1D{k});

    % Determine row and column indices of parameters in the 1D objects that
    % correspond to the current parameter in the 2D object
    rv_idcs = idcs_v(nnv_op_2D(rnum)+1:nnv_op_2D(rnum+1));     
    cv_idcs = idcs_v(nnv_op_2D(cnum)+1:nnv_op_2D(cnum+1));
    cw_idcs = idcs_w(nnw_op_2D(cnum)+1:nnw_op_2D(cnum+1));
    cu_idcs = idcs_u(nnu_op_2D(cnum)+1:nnu_op_2D(cnum+1));

    if ~any(rv_idcs)
        % If the current parameter is empty, move on to the next one
        continue
    end
    
    % Set the value of the current parameter in Top based on value of
    % associated elements in Top_1D
    if any(cv_idcs)
    if ~strcmp(fnames_1D{k},'R') && ~isa(Top.(fnames_2D{k}),'cell')
        Top.(fnames_2D{k}) = subs(polynomial(v_param(rv_idcs,cv_idcs)),vars_old,vars_new);
    elseif ~strcmp(fnames_1D{k},'R')
        Top.(fnames_2D{k}){1} = subs(polynomial(v_param(rv_idcs,cv_idcs)),vars_old,vars_new);
    else
        idx = [1,1];
        for ii=1:3
            idx(var_idx) = ii;
            Top.(fnames_2D{k}){idx(1),idx(2)} = subs(v_param.(['R',num2str(ii-1)])(rv_idcs,cv_idcs),vars_old,vars_new);
        end
    end
    end

    % Set the value of the current parameter in Twop based on value of
    % associated elements in Twop_1D
    if any(cw_idcs)
    if ~strcmp(fnames_1D{k},'R') && ~isa(Top.(fnames_2D{k}),'cell')
        Twop.(fnames_2D{k}) = subs(polynomial(w_param(rv_idcs,cw_idcs)),vars_old,vars_new);
    elseif ~strcmp(fnames_1D{k},'R')
        Twop.(fnames_2D{k}){1} = subs(polynomial(w_param(rv_idcs,cw_idcs)),vars_old,vars_new);
    else
        idx = [1,1];
        for ii=1:3
            idx(var_idx) = ii;
            Twop.(fnames_2D{k}){idx(1),idx(2)} = subs(w_param.(['R',num2str(ii-1)])(rv_idcs,cw_idcs),vars_old,vars_new);
        end
    end
    end

    % Set the value of the current parameter in Tuop based on value of
    % associated elements in Tuop_1D
    if any(cu_idcs)
    if ~strcmp(fnames_1D{k},'R') && ~isa(Top.(fnames_2D{k}),'cell')
        Tuop.(fnames_2D{k}) = subs(polynomial(u_param(rv_idcs,cu_idcs)),vars_old,vars_new);
    elseif ~strcmp(fnames_1D{k},'R')
        Tuop.(fnames_2D{k}){1} = subs(polynomial(u_param(rv_idcs,cu_idcs)),vars_old,vars_new);
    else
        idx = [1,1];
        for ii=1:3
            idx(var_idx) = ii;
            Tuop.(fnames_2D{k}){idx(1),idx(2)} = subs(u_param.(['R',num2str(ii-1)])(rv_idcs,cu_idcs),vars_old,vars_new);
        end
    end
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [Tnop_x,Tnop_y,use_Dx,use_Dy] = check_state_diffs(PDE,Top_x,Top_y,Dvals_xy)
% Determine which option to use in the expansion of the PDE state, and
% determine associated operators. In particular, note that we can expand
% using one of two options:
% 
% Dy option:
%     Dop_y*v = Dop_y*Top_x*Dop_x*v +Twop_x*Dop_y*w   +...
% --> v = Top_y*Tnop_x*Dop_xy*v +Top_y*Twop_x*Dop_y*w +Twop_y*w    +...;
% --> v = Top*Dop_xy*v + Twop*Dop_w*w    +...;
%
% Dx option:
%     Dop_x*v = Dop_x*Top_y*Dop_y*v +Twop_y*Dop_x*w   +....
% --> v = Top_x*Tnop_y*Dop_xy*v +Top_x*Twop_y*Dop_x*w +Twop_x*w    +...;
% --> v = Top*Dop_xy*v + Twop*Dop_w*w    +...;
%
% However, if we use e.g. the Dy option, we need to be able to express
%   Dop_y*(Top_x*Dop_x*v) = Tnop_x*(Dop_y*Dop_x*v) = Tnop_x*Dop_xy*v;
% In particular, note that e.g.
%   Dop_y*(Top_x*v) = (Dop_y*Top_x)*v + Top_x*(Dop_y*v);
% We need the operator Dop_y*Top_x to be zero, so that the expression is
% solely in terms of Dop_y*v;
% Since e.g. Top_x acts only along the x-direction, we may assume that
% Dop_y*Top_x = 0. However, to be sure, this function checks if this is 
% indeed the case.

% We assume for now neither option is strictly necessary
use_Dx = false;   use_Dy = false;

% We will build the new operators one row at a time.
Tnop_x = 0*Top_x;               Tnop_y = 0*Top_y;
Tnop_x.dim(:,1) = zeros(4,1);   Tnop_y.dim(:,1) = zeros(4,1);
xvar = PDE.vars(1,1);   yvar = PDE.vars(2,1);
nv_op = Top_x.dim(:,2);
ridcs = 0;
for ii=1:numel(PDE.x)
    % Extract iith row of each operator
    ridcs = ridcs(end)+(1:PDE.x{ii}.size);
    Tnop_x_ii = Top_x(ii,:);    Tnop_y_ii = Top_y(ii,:);

    % Determine what order derivative must be taken of this row
    dval_x = Dvals_xy(ii,1);      dval_y = Dvals_xy(ii,2);
    if dval_y
        % Take derivative of the operator Top_x
        Tnop_x_ii = clean_opvar(diff(Tnop_x_ii,yvar,dval_y),1e-12);
        
        % Determine which columns in Tnop_x correspond to the highest-order
        % derivative of the state
        fctr_x = 2^dval_y;
        diff_idcs_x = 1:nv_op(1);
        diff_idcs_x = [diff_idcs_x, diff_idcs_x(end)+(1:nv_op(2))];
        diff_idcs_x = [diff_idcs_x, diff_idcs_x(end)+fctr_x*nv_op(3)+(-nv_op(3)+1:0)];
        diff_idcs_x = [diff_idcs_x, diff_idcs_x(end)+fctr_x*nv_op(4)+(-nv_op(4)+1:0)];

        nodiff_idcs_x = setdiff((1:sum(Tnop_x_ii.dim(:,2))),diff_idcs_x);
        if ~(Tnop_x_ii(1,nodiff_idcs_x)==0)
            use_Dx = true;
        end
        Tnop_x_ii = Tnop_x_ii(1,diff_idcs_x);
    end
    Tnop_x = [Tnop_x;Tnop_x_ii];

    if dval_x
        % Take derivative of the operator Top_y
        Tnop_y_ii = clean_opvar(diff(Tnop_y_ii,xvar,dval_x),1e-12);
        
        % Determine which columns in Tnop_y correspond to the highest-order
        % derivative of the state
        fctr_y = 2^dval_x;
        diff_idcs_y = 1:nv_op(1);
        diff_idcs_y = [diff_idcs_y, diff_idcs_y(end)+fctr_y*nv_op(2)+(-nv_op(2)+1:0)];
        diff_idcs_y = [diff_idcs_y, diff_idcs_y(end)+(1:nv_op(3))];        
        diff_idcs_y = [diff_idcs_y, diff_idcs_y(end)+fctr_y*nv_op(4)+(-nv_op(4)+1:0)];

        nodiff_idcs_y = setdiff((1:sum(Tnop_y_ii.dim(:,2))),diff_idcs_y);
        if ~(Tnop_y_ii(1,nodiff_idcs_y)==0)
            if use_Dx
                error('The input PDE does not allow boundary conditions to be enforced separately along each dimension.')
            else
                use_Dy = true;
            end
        end
        Tnop_y_ii = Tnop_y_ii(:,diff_idcs_y);
    end
    Tnop_y = [Tnop_y;Tnop_y_ii];
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [Dvals_w,use_Dx,use_Dy] = compute_input_diffs(Twop_x,Twop_y,Dvals_v,size_x,size_w)
% Determine what derivatives of the input w must be taken in the expansion
% of the state v in terms of its derivatives as specified by Dvals_v
% Dy option:
%     Dop_y*v = Top_x*Dop_xy*v +Twop_x*Dop_y*w   +...
% --> v = Top_y*Top_x*Dop_xy*v +Top_y*Twop_x*Dop_y*w +Twop_y*w    +...;
% --> v = Top*Dop_xy*v + Twop*Dop_w*w    +...;
% 
% Dx option:
%     Dop_x*v = Top_y*Dop_xy*v +Twop_y*Dop_x*w   +....
% --> v = Top_x*Top_y*Dop_xy*v +Top_x*Twop_y*Dop_x*w +Twop_x*w    +...;
% --> v = Top*Dop_xy*v + Twop*Dop_w*w    +...;
%
% INPUTS
% - Twop_x, Twop_y: opvar2d objects such that
%                       v = Top_x*Dop_x*v +Twop_x*w +gx
%                       v = Top_y*Dop_y*v +Twop_y*w +gy
%                   where Dop_x and Dop_y are diagonal differential
%                   operators such that Dop_x(k,k)=(d/dx)^Dvals_vy(k,1) and
%                   Dop_y(k,k) =(d/dy)^Dvals_v(k,2), and where gx and gy
%                   are potential additional forcings;
% - Dvals_xy:       nv x 2 array of type 'double' providing the order of
%                   the derivatives of the PDE state as used in the PIE
%                   representation, so that v_PIE = Dop_xy*v_PDE, where
%                   Dop_xy(k,k)=(d/dx)^Dvals_v(k,1)*(d/dy)^Dvals_v(k,2);
% - size_w:         nw x 1 array specifying for each (vector-valued) input 
%                   component of how many elements it consists; 
% 
% OUTPUTS
% - Dvals_w:        nw x 2 array of type'double' providing the order of the
%                   derivatives of the input w taken in the proposed
%                   expansion, so that
%                   Dop_w(k,k)=(d/dx)^Dvals_w(k,1)*(d/dy)^Dvals_w(k,2);
% - use_Dx:         Logical object set to true if it strictly necessary to
%                   use the Dx option in the expansion.
% - use_Dy:         Logical object set to true if it strictly necessary to
%                   use the Dy option in the expansion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                   
% We have
% v(x,y) = Top_x*Dop_x*v +Qop_x*w   +gx
% v(x,y) = Top_y*Dop_y*v +Qop_y*w   +gy
% It follows that we can express
%
% Dop_y*v = Top_x*Dop_xy*v +Qop_x*Dop_y*w   +...
% --> v = Top_y*Top_x*Dop_xy*v +Top_y*Qop_x*Dop_y*w +Qop_y*w    +...;
% --> v = Top*Dop_xy*v + Qop*Dop_w*w    +...;
%
% or we can express
%
% Dop_x*v = Top_y*Dop_xy*v +Qop_y*Dop_x*w   +....
% --> v = Top_x*Top_y*Dop_xy*v +Top_x*Qop_y*Dop_x*w +Qop_x*w    +...;
% --> v = Top*Dop_xy*v + Qop*Dop_w*w    +...;
%
% In this code, we check which derivatives of w would need to be taken
% using either of these options.
%

% % Extract dimensions of the different operators
nw_op_x = Twop_x.dim(:,2);       nw_op_y = Twop_y.dim(:,2);
ncomps = length(size_w);  % should be equal to sum(nw_op_y);

% Check which derivatives we take of each state component
Dop_x = Dvals_v(:,1);        Dop_x_w = zeros(ncomps,1);
Dop_y = Dvals_v(:,2);        Dop_y_w = zeros(ncomps,1);

% Check which derivatives we take of each state variable
Dop_x = repelem(Dop_x,size_x);
Dop_y = repelem(Dop_y,size_x);

% % We have two options to expand v:
% % Dx option: v = Top_x*Top_y*Dop_xy*v +Top_x*Qop_y*Dop_x*w +Qop_x*w
% % Dy option: v = Top_y*Top_x*Dop_xy*v +Top_y*Qop_x*Dop_y*w +Qop_y*w;
% % We have to choose which option to go for, and what the operator Dop_x
% % or Dop_y would need to look like
use_Dx = false;
use_Dy = false;

% Loop over all inputs w(k)
c_idcs = 0;
for k=1:ncomps
    % % Check what derivative of w(k) would need to be taken if we use Dx
    % % option
    % Check to which states the input w(k) maps through Qop_y
    c_idcs = c_idcs(end) +(1:size_w(k));
    if c_idcs(end)<=nw_op_y(1)
        idcs = c_idcs;
        tmp_log = [~isequal(polynomial(Twop_y.R00(:,idcs)),0); 
                   ~isequal(polynomial(Twop_y.Rx0(:,idcs)),0);
                   ~isequal(polynomial(Twop_y.Ry0(:,idcs)),0);
                   ~isequal(polynomial(Twop_y.R20(:,idcs)),0)];
    elseif c_idcs(end)<=sum(nw_op_y(1:2))
        idcs = c_idcs - nw_op_y(1);
        tmp_log1 = ~isequal(polynomial(Twop_y.R0x(:,idcs)),0);
        tmp_log2 = ~isequal(polynomial(Twop_y.Rxx{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_y.Rxx{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_y.Rxx{3}(:,idcs)),0);
        tmp_log3 = ~isequal(polynomial(Twop_y.Ryx(:,idcs)),0);
        tmp_log4 = ~isequal(polynomial(Twop_y.R2x{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_y.R2x{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_y.R2x{3}(:,idcs)),0);
        tmp_log = [tmp_log1; tmp_log2; tmp_log3; tmp_log4];
    elseif c_idcs(end)<=sum(nw_op_y(1:3))
        idcs = c_idcs - sum(nw_op_y(1:2));
        tmp_log1 = ~isequal(polynomial(Twop_y.R0y(:,idcs)),0);
        tmp_log2 = ~isequal(polynomial(Twop_y.Rxy(:,idcs)),0);
        tmp_log3 = ~isequal(polynomial(Twop_y.Ryy{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_y.Ryy{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_y.Ryy{3}(:,idcs)),0);
        tmp_log4 = ~isequal(polynomial(Twop_y.R2y{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_y.R2y{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_y.R2y{3}(:,idcs)),0);
        tmp_log = [tmp_log1; tmp_log2; tmp_log3; tmp_log4];
    else
        idcs = c_idcs - sum(nw_op_y(1:3));
        tmp_log1 = ~isequal(polynomial(Twop_y.R02(:,idcs)),0);
        tmp_log2 = ~isequal(polynomial(Twop_y.Rx2{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_y.Rx2{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_y.Rx2{3}(:,idcs)),0);
        tmp_log3 = ~isequal(polynomial(Twop_y.Ry2{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_y.Ry2{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_y.Ry2{3}(:,idcs)),0);
        tmp_log4 = ~isequal(polynomial(Twop_y.R22{1,1}(:,idcs)),0);
        for ll=2:numel(Twop_y.R22)
            tmp_log4 = tmp_log4 | ~isequal(polynomial(Twop_y.R22{ll}(:,idcs)),0);
        end
        tmp_log = [tmp_log1; tmp_log2; tmp_log3; tmp_log4];
    end
    % Each state is differentiated up to certain order Dop_x
    % --> w(k) would need to be differentiated up to same order
    tmp_log = any(tmp_log,2);
    xdiffs_k = unique(Dop_x(tmp_log));

    % % Check what derivative of w(k) would need to be taken if we use Dy
    % % option
    % Check to which states the input w(k) maps through Qop_x
    if c_idcs(end)<=nw_op_x(1)
        idcs = c_idcs;
        tmp_log = [~isequal(polynomial(Twop_x.R00(:,idcs)),0); 
                   ~isequal(polynomial(Twop_x.Rx0(:,idcs)),0);
                   ~isequal(polynomial(Twop_x.Ry0(:,idcs)),0);
                   ~isequal(polynomial(Twop_x.R20(:,idcs)),0)];
    elseif c_idcs(end)<=sum(nw_op_x(1:2))
        idcs = c_idcs - nw_op_x(1);
        tmp_log1 = ~isequal(polynomial(Twop_x.R0x(:,idcs)),0);
        tmp_log2 = ~isequal(polynomial(Twop_x.Rxx{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_x.Rxx{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_x.Rxx{3}(:,idcs)),0);
        tmp_log3 = ~isequal(polynomial(Twop_x.Ryx(:,idcs)),0);
        tmp_log4 = ~isequal(polynomial(Twop_x.R2x{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_x.R2x{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_x.R2x{3}(:,idcs)),0);
        tmp_log = [tmp_log1; tmp_log2; tmp_log3; tmp_log4];
    elseif c_idcs(end)<=sum(nw_op_x(1:3))
        idcs = c_idcs - sum(nw_op_x(1:2));
        tmp_log1 = ~isequal(polynomial(Twop_x.R0y(:,idcs)),0);
        tmp_log2 = ~isequal(polynomial(Twop_x.Rxy(:,idcs)),0);
        tmp_log3 = ~isequal(polynomial(Twop_x.Ryy{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_x.Ryy{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_x.Ryy{3}(:,idcs)),0);
        tmp_log4 = ~isequal(polynomial(Twop_x.R2y{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_x.R2y{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_x.R2y{3}(:,idcs)),0);
        tmp_log = [tmp_log1; tmp_log2; tmp_log3; tmp_log4];
    else
        idcs = c_idcs - sum(nw_op_x(1:3));
        tmp_log1 = ~isequal(polynomial(Twop_x.R02(:,idcs)),0);
        tmp_log2 = ~isequal(polynomial(Twop_x.Rx2{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_x.Rx2{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_x.Rx2{3}(:,idcs)),0);
        tmp_log3 = ~isequal(polynomial(Twop_x.Ry2{1}(:,idcs)),0) | ...
                    ~isequal(polynomial(Twop_x.Ry2{2}(:,idcs)),0) | ...
                     ~isequal(polynomial(Twop_x.Ry2{3}(:,idcs)),0);
        tmp_log4 = ~isequal(polynomial(Twop_x.R22{1,1}(:,idcs)),0);
        for ll=2:numel(Twop_x.R22)
            tmp_log4 = tmp_log4 | ~isequal(polynomial(Twop_x.R22{ll}(:,idcs)),0);
        end
        tmp_log = [tmp_log1; tmp_log2; tmp_log3; tmp_log4];        
    end
    % Each state is differentiated up to certain order Dop_y
    % --> w(k) would need to be differentiated up to same order
    tmp_log = any(tmp_log,2);
    ydiffs_k = unique(Dop_y(tmp_log));

    % Check what derivative we would need to take if we use Dx or Dy option
    if ~isempty(xdiffs_k) && ~isempty(ydiffs_k)
        % Both Dx and Dy option require a derivative to be taken
        if any(xdiffs_k) && any(ydiffs_k)
            % Whatever option we choose, w(k) would have to appear as both
            % a derivative wrt x and a derivative wrt y
            % --> not supported
            error(['It seems some input appears as both w(k) (or u(k)) and a derivative of w(k) (or u(k)) in the PIE; this case is not yet supported.'])
        else
            % Both options would require derivative of order 0 to be taken
            % --> acceptable, take no derivative of w(k)
            Dop_x_w(k) = 0;     Dop_y_w(k) = 0;
        end
    elseif ~isempty(xdiffs_k) && isempty(ydiffs_k)
        % The input does not contribute to Qop_x
        % --> if we use Dy option, the input need not be differentiated
        Dop_y_w(k) = 0;
        if numel(xdiffs_k)>=2
            % If we use Dx option, then w(k) must appear as two different
            % derivatives, which is not supported
            if use_Dx
                error(['It seems some input appears as both w(k) (or u(k)) and a derivative of w(k) (or u(k)) in the PIE; this case is not yet supported.'])
            else
                use_Dy = true;
            end
        else
            % If we use Dx option, then w(k) must appear as derivative of
            % order xdiffs_k
            Dop_x_w(k) = xdiffs_k;
        end
    elseif isempty(xdiffs_k) && ~isempty(ydiffs_k)
        % The input does not contribute to Qop_y
        % --> if we use Dx option, the input need not be differentiated
        Dop_x_w(k) = 0;
        if numel(ydiffs_k)>=2
            % If we use Dy option, then w(k) must appear as two different
            % derivatives, which is not supported
            if use_Dy
                error(['It seems some input appears as both w(k) (or u(k)) and a derivative of w(k) (or u(k)) in the PIE; this case is not yet supported.'])
            else
                use_Dx = true;
            end
        else
            % If we use Dy option, then w(k) must appear as derivative of
            % order ydiffs_k
            Dop_y_w(k) = ydiffs_k;
        end
    else
        % The input does not contribute to either Qop_x or Qop_y
        % --> whatever option we choose, we don't need to take a derivative
        Dop_x_w(k) = 0;         Dop_y_w(k) = 0;
    end
end
Dvals_w = [Dop_x_w,Dop_y_w];

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

% Assume PDE is 2D
nvars = 2;

if nvars==0
    np_op = [sum(obj_tab(:,2)); 0; 0; 0];
elseif nvars==1
    dep_tab = obj_tab(:,3:2+nvars);
    rindcs_0 = ~any(dep_tab,2);                 % Row indices in comp_tab associated to finite dimensional states
    rindcs_1 = logical(dep_tab(:,1));           % Row indices in comp_tab associated to states that vary just in s1
    
    np_op = [sum(obj_tab(rindcs_0,2));
            sum(obj_tab(rindcs_1,2));
            0;
            0];
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
function print_reorder_diff_summary(PDE,obj,comp_order,Dvals)
% print_reorder_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning the new order of
% the components "obj" from the PDE structure as they will appear in the
% PIE structure.
%
% INPUTS:
% - PDE:    A "pde_struct" class object defining a PDE.
% - obj:    Char 'x', 'u', 'w', 'y', 'z', or 'BC', indicating for which
%           object to display the new order of the variables.
% - comp_order: comp_order(j) provides the index of the component in the
%               original PDE associated to component j in the new PDE.
% - Dvals:  A nobj x nvars object of 'type' double specifying for each
%           element PDE.obj(ii) the order of the derivative taken with
%           respect to each of the nvars global variables in PDE.vars, as
%           the object appears in the PIE representation
%
% OUTPUTS:
% Displays information in the command window concerning the new order of
% the components in PDE.obj, compared to the original PDE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set a name associated to each object.
if strcmp(obj,'x')
    obj_name = 'state components';
elseif strcmp(obj,'y')
    obj_name = 'observed outputs';
elseif strcmp(obj,'z')
    obj_name = 'regulated outputs';
elseif strcmp(obj,'u')
    obj_name = 'actuator inputs';
elseif strcmp(obj,'w')
    obj_name = 'exogenous inputs';
elseif strcmp(obj,'BC')
    obj_name = 'boundary conditions';
end
ncomps = numel(PDE.(obj));
if ~any(any(Dvals)) %all(comp_order == (1:ncomps)')
    % No derivative of the components is takenThe order of the components has not changed.
    %fprintf(['\n','The order of the ',object_name,'s ',obj,' has not changed.\n']);
    return
else
    % Otherwise, we list the new order of the components.
    %fprintf(['\n','The ',object_name,'s have been reordered as:\n']);
    fprintf(['\n','WARNING: Some ',obj_name,' from the PDE will appear as derivatives in the PIE representation:\n']);
end

% Use UNICODE to add subscript indices to different components.
sub_s = '\x209B';
partial = '\x2202';
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
sup_num = cell(10,1);       
sup_num{1} = '\x2070';      % Superscript 0
sup_num{2} = '\xB9';        % Superscript 1
sup_num{3} = '\xB2';        % etc.
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
% Estimate maximal size of derivatives applied to the component
diff_length = floor(log(max(Dvals,1))./log(10)).*logical(Dvals);    % determine length of power of diff
diff_length = diff_length + 7*logical(Dvals);     % add length of "(\partial/\partial s)"
diff_length = max(sum(diff_length,2));          % determine maximal total length
LHS_length_max = diff_length + 1+1 + n_digits + lngth_varnames_mean+3; % e.g. 'd_s1^2 x13(t,s1,s2,s3)', 


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
    if numel(PDE.(obj))==1
        % There is only one component --> no need to give index
        Lcomp_idx = '';
    elseif numel(PDE.(obj))<=9
        % The component number consists of a single decimal.
        Lcomp_idx = sub_num{old_idx+1};
        LHS_length = LHS_length + 1;
    else
        % The component number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(old_idx)')+1)');
        LHS_length = LHS_length + length(num2str(old_idx));
    end

    % Establish which derivaitve of the component is taken.
    dval = Dvals(ii,:);
    %dtot = sum(dval);
    % Write d^{k}obj/
    % if dtot==0
    %     Ldiff = [obj,Lcomp_idx];
    % elseif dtot==1
    %     Ldiff = [partial,obj,Lcomp_idx,'/'];
    %     LHS_length = LHS_length+1;
    % elseif dtot<=9
    %     sup_diff = sup_num{dtot+1};
    %     Ldiff = [partial,sup_diff,obj,Lcomp_idx,'/'];
    %     LHS_length = LHS_length+2;
    % else
    %     sup_diff = cell2mat(sup_num(str2num(num2str(dtot)')+1)');
    %     Ldiff = [partial,sup_diff,obj,Lcomp_idx,'/'];
    %     LHS_length = LHS_length+1+length(sup_diff);
    % end
    Ldiff = [];
    for jj=1:nvars
        if dval(jj)==0
            % No derivative is taken wrt this variable --> move on to next one
            continue
        end
        Ldiff = ['(',partial,'/',partial,PDE.vars(jj,1).varname{1},')'];
        LHS_length = LHS_length +6+length(PDE.vars(jj,1).varname{1});
        if dval(jj)>=2
            if dval(jj)<=9
                sup_diff = sup_num{dval(jj)+1};
            else
                sup_diff = cell2mat(sup_num(str2num(num2str(dval(jj))')+1)');
            end
            Ldiff = [Ldiff,sup_diff];
            LHS_length = LHS_length + ceil(log(dval(jj)/log(10)));
        end
        % if dval(jj)==1
        %     % The order of the derivative is 1
        %     Ldiff = ['(',partial,'/',partial,PDE.vars(jj,1).varname{1},')'];
        %     %Ldiff = [Ldiff,diff_str];
        %     LHS_length = LHS_length +5+length(PDE.vars(jj,1).varname{1});
        % elseif dval(jj)<=9
        %     % The order of the derivative consists of a single decimal
        %     sup_diff = sup_num{dval(jj)+1};
        %     Ldiff = ['(',partial,'/',partial,PDE.vars(jj,1).varname{1},')',sup_diff];
        %     %Ldiff = [Ldiff,diff_str];
        %     LHS_length = LHS_length +6+length(PDE.vars(jj,1).varname{1});
        % else
        %     % The order of the derivative consists of muliple decimals
        %     sup_diff = cell2mat(sup_num(str2num(num2str(dval(jj))')+1)');
        %     diff_str = ['(',partial,'/',partial,PDE.vars(jj,1).varname{1},')',sup_diff];
        %     Ldiff = [Ldiff,diff_str];
        %     LHS_length = LHS_length +1+length(PDE.vars(jj,1).varname{1})+length(sup_diff);
        % end
    end

    % Set the name of the component, including its depdence on spatial
    % variables.
    LHS_name = [' ',Ldiff,obj,Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = 1 + LHS_length + 3;
        
    % For the added state components, indicate of which state component
    % they are the temporal derivative.
    new_idx = ii;
    Rcomp_idx = cell2mat(sub_num(str2num(num2str(new_idx)')+1)');
    RHS = [' -->   ',obj,Rcomp_idx,'(',varnames_ii_t,')'];
%    RHS = '';
    
    % % % Finally, display:
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_input_diff_summary(obj,comp_order,var_num,Dvals)
% print_reorder_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning the new order of
% the components "obj" in the PDE structure.
%
% INPUTS:
% - obj:    'u' or 'w', indicating for which object to display how it
%           appears in the PIE.
% - comp_order: comp_order(j) provides the index of the component in the
%               original PDE associated to component j in the new PDE.
% - Dvals:  ncompsx1 'double' array of integer values, specifying what
%           order derivative is taken of component j w/r tp variable varnum
%
% OUTPUTS:
% Displays information in the command window concerning the new order of
% the components in PDE.obj, compared to the original PDE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set a name associated to each object.
if strcmp(obj,'x')
    obj_name = 'state components';
elseif strcmp(obj,'y')
    obj_name = 'observed outputs';
elseif strcmp(obj,'z')
    obj_name = 'regulated outputs';
elseif strcmp(obj,'u')
    obj_name = 'actuator inputs';
elseif strcmp(obj,'w')
    obj_name = 'exogenous inputs';
elseif strcmp(obj,'BC')
    obj_name = 'boundary conditions';
end
ncomps = length(comp_order);
if all(comp_order == (1:ncomps)')
    % The order of the components has not changed.
    fprintf(['\n','The order of the ',obj_name,'s ',obj,' has not changed.\n']);
    return
else
    % Otherwise, we list the new order of the components.
    fprintf(['\n','The ',obj_name,'s have been reordered as:\n']);
end

% Use UNICODE to add subscript indices to different components.
sub_s = '\x209B';
partial = '\x2202';
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
sup_num = cell(10,1);       
sup_num{1} = '\x2070';      % Superscript 0
sup_num{2} = '\xB9';        % Superscript 1
sup_num{3} = '\xB2';        % etc.
sup_num{4} = '\xB3';
sup_num{5} = '\x2074';
sup_num{6} = '\x2075';
sup_num{7} = '\x2076';
sup_num{8} = '\x2077';
sup_num{9} = '\x2078';
sup_num{10} = '\x2079';
sub_var_num = sub_num{var_num+1};
% Determine how many subscripts will be needed.
n_digits = 1+floor(log(ncomps)/log(10));


fprintf(['\n','WARNING: The following ',obj_name,' from the original PDE will appear as derivatives in the PIE representation:\n']);
for ii=1:ncomps
    if Dvals(ii)==0
        % No need to list components that are unchanged
        continue
    end
    old_idx = comp_order(ii);   new_idx = ii;

    % % % Display the component as it appears in the PDE
    %LHS_length = 0;
    if ncomps==1
        % There is only one component --> no need to give index
        Lcomp_idx = '';
    elseif ncomps<=9
        % The component number consists of a single decimal.
        Lcomp_idx = sub_num{old_idx+1};
        %LHS_length = LHS_length + 1;
    else
        % The component number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(old_idx)')+1)');
        %LHS_length = LHS_length + length(num2str(old_idx));
    end
    % Set the name of the component.
    LHS_name = [' ',obj,Lcomp_idx];
    %LHS_length = 2 + LHS_length;

    % % % Display the component as it appears in the PIE
    % Establish the (subscript) index for the component.
    %RHS_length = 0;
    if ncomps==1
        % There is only one component --> no need to give index
        Rcomp_idx = '';
    elseif ncomps<=9
        % The component number consists of a single decimal.
        Rcomp_idx = sub_num{new_idx+1};
        %RHS_length = RHS_length + 1;
    else
        % The component number consists of multiple decimals.
        Rcomp_idx = cell2mat(sub_num(str2num(num2str(new_idx)')+1)');
        %RHS_length = RHS_length + length(num2str(new_idx));
    end
    % Set the derivative
    if Dvals(ii)==1
        % The order of the derivative is 1
        RHS = [partial,sub_s,sub_var_num,' ',obj,Rcomp_idx];
        %RHS_length = RHS_length + 6;
    elseif Dvals(ii)<=9
        % The order of the derivative consists of a single decimal
        sup_diff = sup_num{Dvals(ii)+1};
        RHS = [partial,sub_s,sub_var_num,sup_diff,' ',obj,Rcomp_idx];
    else
        % The order of the derivative consists of muliple decimals
        sup_diff = cell2mat(sup_num(str2num(num2str(Dvals(ii))')+1)');
        RHS = [partial,sub_s,sub_var_num,sup_diff,' ',obj,Rcomp_idx];
        %RHS_length = RHS_length + 5 + length(num2str(tdiff));
    end
    RHS = [' -->   ',RHS];
    
    % % % Finally, display:
    %MT_space = max(LHS_length_max-LHS_length,1);
%    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
    fprintf(['  ',LHS_name,RHS,'\n'])
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%      OLD STUFF    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% function Top = build_2D_Top(Top,Top_x,Top_y,comp_order_x,comp_order_y)
% 
% % % Extract the appropriate spatial variables
% xx_n = [Top.var1(1); Top.var2(1)];      xx_o = [Top_x.var1; Top_x.var2];
% yy_n = [Top.var1(2); Top.var2(2)];      yy_o = [Top_y.var1; Top_y.var2];
% 
% 
% % % Extract dimensions of the different operators
% % (should have sum(n_op)=sum(n_op_x)=sum(n_op_y))
% n_op = Top.dim(:,1);
% nn_op = [0;cumsum(n_op)];
% n_op_x = Top_x.dim(:,1);    n_op_y = Top_y.dim(:,1);
% 
% % Determine which elements of Top correspond to which element of Top_x and
% % Top_y: Row k in Top corresponds to row idcs_x(k) in Top_x
% [~,idcs_x] = sort(comp_order_x.x);
% [~,idcs_y] = sort(comp_order_y.x);
% 
% % Adjust for splitting of parameters in the opvar structure
% idcs_x(idcs_x>n_op_x(1)) = idcs_x(idcs_x>n_op_x(1)) - n_op_x(1);
% idcs_y(idcs_y>n_op_y(1)) = idcs_y(idcs_y>n_op_y(1)) - n_op_y(1);
% 
% % % Determine which state elements are treated as finite-dimensional in the
% % % operators Top_x and Top_y;
% % is_0D_x = idcs_x <= n_op_x(1);      is_1D_x = idcs_x > n_op_x(1);
% % is_0D_y = idcs_y <= n_op_y(1);      is_1D_y = idcs_y > n_op_y(1);
% % 
% % % Extract indices of Top_x and Top_y associated to ODE and 1D variables
% % idcs_x0 = idcs_x(is_0D_x);          idcs_x1 = idcs_x(is_1D_x);
% % idcs_y0 = idcs_y(is_0D_y);          idcs_y1 = idcs_y(is_1D_y);
% % %NO, have to separate more carefully into 0D, 1D, etc.
% 
% % Declare which parameters in Top correspond to which parameters in Top_x
% % and Top_y
% fnames_2D = {'R00','R0x','R0y','R02';
%              'Rx0','Rxx','Rxy','Rx2';
%              'Ry0','Ryx','Ryy','Ry2';
%              'R20','R2x','R2y','R22'};
% 
% fnames_x = {'P','Q1','P','Q1';
%             'Q2','R','Q2','R';
%             'P','Q1','P','Q1';
%             'Q2','R','Q2','R'};
% 
% fnames_y = {'P','P','Q1','Q1';
%             'P','P','Q1','Q1';
%             'Q2','Q2','R','R';
%             'Q2','Q2','R','R'};
% 
% % Set the value of each of the parameters in Top
% for k=1:numel(fnames_2D)
%     % Determine which rows and columns of Top the parameter corresponds to
%     [rnum,cnum] = ind2sub(size(fnames_2D),k);
%     r_idcs = nn_op(rnum)+1:nn_op(rnum+1);     
%     c_idcs = nn_op(cnum)+1:nn_op(cnum+1);
% 
%     % If the parameter is empty, move on to the next one
%     if ~any(r_idcs) || ~any(c_idcs)
%         continue
%     end
% 
%     % Extract parameter in Top_x associated to current parameter in Top
%     x_param = Top_x.(fnames_x{k});
%     % Determine row and column indices in parameter in Top_x corresponding 
%     % to current parameter in Top
%     r_idcs_x = idcs_x(r_idcs);      c_idcs_x = idcs_x(c_idcs);
% 
%     % Extract parameter in Top_y associated to current parameter in Top
%     y_param = Top_y.(fnames_y{k});
%     % Determine row and column indices in parameter in Top_x corresponding 
%     % to current parameter in Top
%     r_idcs_y = idcs_y(r_idcs);      c_idcs_y = idcs_y(c_idcs);
% 
%     % Compute the value of the current parameter in Top as product of
%     % parameters in Top_x and Top_y
%     if ~strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         T_param = subs(polynomial(x_param(r_idcs_x,c_idcs_x)),xx_o,xx_n).*...
%                     subs(polynomial(y_param(r_idcs_y,c_idcs_y)),yy_o,yy_n);
%     elseif strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         T_param = cell(3,1);
%         for ii=1:3
%             T_param{ii} = subs(x_param.(['R',num2str(ii-1)])(r_idcs_x,c_idcs_x),xx_o,xx_n).*...
%                         subs(polynomial(y_param(r_idcs_y,c_idcs_y)),yy_o,yy_n);
%         end
%     elseif ~strcmp(fnames_x{k},'R') && strcmp(fnames_y{k},'R')
%         T_param = cell(1,3);
%         for ii=1:3
%             T_param{ii} = subs(polynomial(x_param(r_idcs_x,c_idcs_x)),xx_o,xx_n).*...
%                         subs(y_param.(['R',num2str(ii-1)])(r_idcs_y,c_idcs_y),yy_o,yy_n);
%         end
%     else
%         T_param = cell(3,3);
%         for ii=1:3
%         for jj=1:3
%             T_param{ii,jj} = subs(x_param.(['R',num2str(ii-1)])(r_idcs_x,c_idcs_x),xx_o,xx_n).*...
%                                 subs(y_param.(['R',num2str(jj-1)])(r_idcs_y,c_idcs_y),yy_o,yy_n);
%         end
%         end
%     end
% 
%     % Set the parameter value
%     Top.(fnames_2D{k}) = T_param;
% end
% 
% 
% 
% % Top.R22{1,1} = subs(Top_x.R.R0,xx_o,xx_n).*subs(Top_y.R.R0,yy_o,yy_n);
% % Top.R22{2,1} = subs(Top_x.R.R1,xx_o,xx_n).*subs(Top_y.R.R0,yy_o,yy_n);
% % Top.R22{3,1} = subs(Top_x.R.R2,xx_o,xx_n).*subs(Top_y.R.R0,yy_o,yy_n);
% % Top.R22{1,2} = subs(Top_x.R.R0,xx_o,xx_n).*subs(Top_y.R.R1,yy_o,yy_n);
% % Top.R22{2,2} = subs(Top_x.R.R1,xx_o,xx_n).*subs(Top_y.R.R1,yy_o,yy_n);
% % Top.R22{3,2} = subs(Top_x.R.R2,xx_o,xx_n).*subs(Top_y.R.R1,yy_o,yy_n);
% % Top.R22{1,3} = subs(Top_x.R.R0,xx_o,xx_n).*subs(Top_y.R.R2,yy_o,yy_n);
% % Top.R22{2,3} = subs(Top_x.R.R1,xx_o,xx_n).*subs(Top_y.R.R2,yy_o,yy_n);
% % Top.R22{3,3} = subs(Top_x.R.R2,xx_o,xx_n).*subs(Top_y.R.R2,yy_o,yy_n);
% 
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Finally, set the values of the parameters defining Qop:
% fnames_2D = {'R00','R0x','R0y','R02';
%              'Rx0','Rxx','Rxy','Rx2';
%              'Ry0','Ryx','Ryy','Ry2';
%              'R20','R2x','R2y','R22'};
% 
% fnames_x = {'P','Q1','P','Q1';
%             'Q2','R','Q2','R';
%             'P','Q1','P','Q1';
%             'Q2','R','Q2','R'};
% 
% fnames_y = {'P','P','Q1','Q1';
%             'P','P','Q1','Q1';
%             'Q2','Q2','R','R';
%             'Q2','Q2','R','R'};
% 
% % Set the value of each of the parameters in Qop
% if use_Dx
% % % % Proceed with the expansion
% % % % v = Top_x*Top_y*Dop_xy*v +Top_x*Qop_y*Dop_x*w +Qop_x*w
% for k=1:numel(fnames_2D)
%     % Determine which rows and columns of Qop the parameter corresponds to
%     [rnum,cnum] = ind2sub(size(fnames_2D),k);
%     r_idcs = nnr_op(rnum)+1:nnr_op(rnum+1);     
%     c_idcs = nnc_op(cnum)+1:nnc_op(cnum+1);
% 
%     % If the parameter is empty, move on to the next one
%     if ~any(r_idcs) || ~any(c_idcs)
%         continue
%     end
% 
%     % Extract parameter in Qop_x associated to current parameter in Qop
%     x_param = Qop_x.(fnames_x{k});
%     % Determine row and column indices in parameter in Qop_x corresponding 
%     % to current parameter in Qop
%     r_idcs_x = idcs_vx(r_idcs);      c_idcs_x = idcs_wx_op(c_idcs);
% 
%     % Extract parameter in Top_x and Qop_y associated to parameter in Qop
%     vx_param = Top_x.(fnames_x{k});      wy_param = Qop_y.(fnames_y{k});
%     % Determine column indices in parameter in Qop_y associated to current
%     % parameter in Qop  
%     c_idcs_y = idcs_wy_op(c_idcs);
% 
%     % Compute the value of the current parameter in Qop as 
%     %   Top_x*Qop_y +Qop_x
%     if ~strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         Q_param = subs(polynomial(x_param(r_idcs_x,c_idcs_x)),xx_o,xx_n) +...
%                    subs(polynomial(vx_param(r_idcs_x,:)),xx_o,xx_n)*...
%                     subs(polynomial(wy_param(:,c_idcs_y)),yy_o,yy_n);
%     elseif strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         Q_param = cell(3,1);
%         for ii=1:3
%             Q_param{ii} = subs(x_param.(['R',num2str(ii-1)])(r_idcs_x,c_idcs_x),xx_o,xx_n) +...
%                            subs(vx_param.(['R',num2str(ii-1)])(r_idcs_x,:),xx_o,xx_n)*...
%                             subs(polynomial(wy_param(:,c_idcs_y)),yy_o,yy_n);
%         end
%     elseif ~strcmp(fnames_x{k},'R') && strcmp(fnames_y{k},'R')
%         Q_param = cell(1,3);
%         for ii=1:3
%             Q_param{ii} = subs(polynomial(x_param(r_idcs_x,c_idcs_x)),xx_o,xx_n) +...
%                            subs(polynomial(vx_param(r_idcs_x,:)),xx_o,xx_n)*...
%                             subs(wy_param.(['R',num2str(ii-1)])(:,c_idcs_y),yy_o,yy_n);
%         end
%     else
%         Q_param = cell(3,3);
%         for ii=1:3
%         for jj=1:3
%             Q_param{ii,jj} = subs(x_param.(['R',num2str(ii-1)])(r_idcs_x,c_idcs_x),xx_o,xx_n) +...
%                               subs(vx_param.(['R',num2str(ii-1)])(r_idcs_x,:),xx_o,xx_n)*...
%                                subs(wy_param.(['R',num2str(jj-1)])(:,c_idcs_y),yy_o,yy_n);
%         end
%         end
%     end
%     % Set the parameter value
%     Qop.(fnames_2D{k}) = Q_param;
% end
% else
% % % % Proceed with the expansion
% % % % v = Top_x*Top_y*Dop_xy*v +Top_y*Qop_x*Dop_y*w +Qop_y*w
% for k=1:numel(fnames_2D)
%     % Determine which rows and columns of Qop the parameter corresponds to
%     [rnum,cnum] = ind2sub(size(fnames_2D),k);
%     r_idcs = nnr_op(rnum)+1:nnr_op(rnum+1);     
%     c_idcs = nnc_op(cnum)+1:nnc_op(cnum+1);
% 
%     % If the parameter is empty, move on to the next one
%     if ~any(r_idcs) || ~any(c_idcs)
%         continue
%     end
% 
%     % Extract parameter in Qop_y associated to current parameter in Qop
%     y_param = Qop_y.(fnames_y{k});
%     % Determine row and column indices in parameter in Qop_y corresponding 
%     % to current parameter in Qop
%     r_idcs_y = idcs_vy(r_idcs);      c_idcs_y = idcs_wy_op(c_idcs);
%     % NO! still need to adjust r_idcs_y for the parameter!!!
% 
%     % Compute the contribution of Qop_y to the current parameter in Qop 
%     % Note that     Qop = Top_y*Qop_x +Qop_y
%     if ~strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         Q_param = subs(polynomial(y_param(r_idcs_y,c_idcs_y)),yy_o,yy_n);
%     elseif strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         Q_param = cell(3,1);
%         Q_param{1} = subs(polynomial(y_param(r_idcs_y,c_idcs_y)),yy_o,yy_n);
%         Q_param{2} = polynomial(zeros(size(Q_param{1})));
%         Q_param{3} = polynomial(zeros(size(Q_param{1})));
%     elseif ~strcmp(fnames_x{k},'R') && strcmp(fnames_y{k},'R')
%         Q_param = cell(1,3);
%         for ii=1:3
%             Q_param{ii} = subs(y_param.(['R',num2str(ii-1)])(r_idcs_y,c_idcs_y),yy_o,yy_n);
%         end
%     else
%         Q_param = cell(3,3);
%         for jj=1:3
%             Q_param{1,jj} = subs(y_param.(['R',num2str(jj-1)])(r_idcs_y,c_idcs_y),yy_o,yy_n);
%             Q_param{2,jj} = polynomial(zeros(size(Q_param{1,jj})));
%             Q_param{3,jj} = polynomial(zeros(size(Q_param{1,jj})));
%         end
%     end
% 
%     % % This leaves the contribution of Top_x*Qop_y, which is more
%     % complicated...
%     c_idcs_x = idcs_wx_op(c_idcs);
%     if rnum==1 || rnum==2
%         % State components do not vary in y
%         vy_param1 = Top_y.P(r_idcs_y,:);        % NO, must extract certain columns!!!    
%         vy_param2 = subs(Top_y.Q1(r_idcs_y,:),yy_o,yy_n);
%     else
%         vy_param1 = subs(Top_y.Q2(r_idcs_y,:),yy_o,yy_n);    
%         vy_param2 = cell(1,3);
%         vy_param2{1} = subs(Top_y.R.R0(r_idcs_y,:),yy_o,yy_n);
%         vy_param2{2} = subs(Top_y.R.R1(r_idcs_y,:),yy_o,yy_n);
%         vy_param2{3} = subs(Top_y.R.R2(r_idcs_y,:),yy_o,yy_n);
%     end
%     if cnum==1 || cnum==3
%         % Inputs do not vary in x
%         wx_param1 = Qop_x.P(:,c_idcs_x);    
%         wx_param2 = subs(Qop_x.Q2(:,c_idcs_x,:),xx_o,xx_n);
%     else
%         wx_param1 = subs(Qop_x.Q1(:,c_idcs_x),xx_o,xx_n);    
%         wx_param2 = cell(1,3);
%         wx_param2{1} = subs(Qop_x.R.R0(:,c_idcs_x),xx_o,xx_n);
%         wx_param2{2} = subs(Qop_x.R.R1(:,c_idcs_x),xx_o,xx_n);
%         wx_param2{3} = subs(Qop_x.R.R2(:,c_idcs_x),xx_o,xx_n);
%     end
%     if ~isa(vy_param2,'cell') && ~isa(wx_param2,'cell')
% 
%     end
% 
%     % Extract parameter in Top_y and Qop_x associated to parameter in Qop
%     vy_param = Top_y.(fnames_y{k});      wx_param = Qop_x.(fnames_x{k});
%     % Determine column indices in parameter in Qop_y associated to current
%     % parameter in Qop  
%     c_idcs_x = idcs_wx_op(c_idcs);
% 
% 
%     % Consider element R2x
%     % Suppose that Qop_x has columns Q1_w_x(:,c) mapping w(x) to [v0;v2(y)] 
%     % and R_w_x(:,c) mapping w(x) to [v1(x);v3(x,y)]
%     % Suppose that Top_y has rows Q2_v_y(r,:) mapping [v0;v1(x)] to v3(x,y)
%     % and rows R_v_y(r,:) mapping [v2(y);v3(x,y)] to v3(x,y);
%     % Then these combine into a contribution
%     % R0 = Q2_v_y(r,0)*Q1_w_x(0,c) +R0_v_y(r,2)*Q1_w_x(2,c)
%     %       +Q2_v_y(r,1)*R0_w_x(1,c) +R0_v_y(r,3)*R0_w_x(3,c)
%     %      +int_{c}^{y} R1_v_y(r,2)(y,nu)*Q1_w_x(2,c)(nu) dnu
%     %       +int_{c}^{y} R1_v_y(r,3)(y,nu)*R0_w_x(3,c)(x,nu) dnu
%     %      +int_{y}^{d} R2_v_y(r,2)(y,nu)*Q1_w_x(2,c)(nu) dnu
%     %       +int_{y}^{d} R2_v_y(r,3)(y,nu)*R0_w_x(3,c)(x,nu) dnu
%     % R1 = 
% 
%     % Compute the value of the current parameter in Qop as 
%     %   Top_x*Qop_y +Qop_x
%     if ~strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         Q_param = subs(polynomial(y_param(r_idcs_y,c_idcs_y)),yy_o,yy_n) +...
%                    subs(polynomial(vy_param(r_idcs_y,:)),yy_o,yy_n)*...
%                     subs(polynomial(wx_param(:,c_idcs_x)),xx_o,xx_n);
%     elseif strcmp(fnames_x{k},'R') && ~strcmp(fnames_y{k},'R')
%         Q_param = cell(3,1);
%         for ii=1:3
%             Q_param{ii} = subs(polynomial(y_param(r_idcs_y,c_idcs_y)),yy_o,yy_n) +...
%                            subs(polynomial(vy_param(r_idcs_y,:)),yy_o,yy_n)*...
%                             subs(wx_param.(['R',num2str(ii-1)])(:,c_idcs_x),xx_o,xx_n);
%         end
%     elseif ~strcmp(fnames_x{k},'R') && strcmp(fnames_y{k},'R')
%         Q_param = cell(1,3);
%         for ii=1:3
%             Q_param{ii} = subs(y_param.(['R',num2str(ii-1)])(r_idcs_y,c_idcs_y),yy_o,yy_n) +...
%                            subs(vy_param.(['R',num2str(ii-1)])(r_idcs_y,:),yy_o,yy_n)*...
%                             subs(polynomial(wx_param(:,c_idcs_x)),xx_o,xx_n);
%         end
%     else
%         Q_param = cell(3,3);
%         for ii=1:3
%         for jj=1:3
%             Q_param{ii,jj} = subs(y_param.(['R',num2str(ii-1)])(r_idcs_y,c_idcs_y),yy_o,yy_n) +...
%                               subs(vy_param.(['R',num2str(ii-1)])(r_idcs_y,:),yy_o,yy_n)*...
%                                subs(wx_param.(['R',num2str(jj-1)])(:,c_idcs_x),xx_o,xx_n);
%         end
%         end
%     end
%     % Set the parameter value
%     Qop.(fnames_2D{k}) = Q_param;
%end
%
%end



% % % % OLD Stuff
% 
% % % Check which of the inputs must be differentiated in the map
% xx_idcs = any(~isequal(Qop_y.Q2,0),1);      % Elements of sgnl on which Qop_y acts
% yy_idcs = any(~isequal(Qop_x.Q2,0),1);      % Elements of sgnl on which Qop_x acts
% if ~any(xx_idcs) && ~any(yy_idcs)
%     % Neither operator acts on sgnl --> contribution of sgnl is just 0
%     return
% elseif ~any(xx_idcs)
%     % No nonzero columns in Qop_y means we can expand
%     %   v(x,y) = Top_x*Top_y*v_{xxyy} +Qop_x*sgnl;
%     % not requiring any derivatives of sgnl
%     % --> only need to morph opvar Qop_x into opvar2d Qop
%     Qop.R2y{1} = subs(Qop_x.Q2(:,find(yy_idcs)),xx_o,xx_n);
%     return
% elseif ~any(yy_idcs)
%     % No nonzero columns in Qop_x means we can expand
%     %   v(x,y) = Top_y*Top_x*v_{xxyy} +Qop_y*sgnl;
%     % not requiring any derivatives of sgnl
%     % --> only need to morph opvar Qop_y into opvar2d Qop
%     Qop.R2x{1} = subs(Qop_y.Q2(:,find(xx_idcs)),yy_o,yy_n);
%     return
% end
% 
% % % If the input sgnl contributes to both the x and y BCs, default to case
% % %     v(x,y) = Top_y*Top_x*v_{xxyy} +Top_y*Qop_x*sgnl_{yy} +Qop_y*sgnl;
% 
% % If a particular element of "sgnl" is acted on by both Qop_x and Qop_y,
% % then both sgnl(k) and sgnl_{yy}(k) would need to appear in the expression
% % for v. This is not supported.
% if any(xx_idcs & yy_idcs)
%     error(['It seems some input appears as both ',obj,'(k) and ',obj,'_{yy}(k) in the PIE; this case is not yet supported.'])
% end
% xx_idcs = find(xx_idcs);    yy_idcs = find(yy_idcs);
% 
% fprintf(['\n','WARNING: The following inputs ',obj,' from the original PDE will appear as derivatives in the PIE representation:\n']);
% for j=1:length(yy_idcs)
%     idx = yy_idcs(j);
%     sgnl_idx_o = num2str(PDE.([obj,'_tab'])(idx,1));
%     sgnl_idx_n = num2str(idx);
% 
%     fprintf(['  ',obj,'(',sgnl_idx_o,') --> ',obj,'_{yy}(',sgnl_idx_n,')','\n']);
% end
% 
% % fprintf(['\n WARNING: Disturbances ',sgnl,'(k) with k in {']);
% % fprintf('%g, ', yy_idcs(1:end-1));
% % fprintf('%g}', yy_idcs(end));
% % fprintf([' will appear as ',sgnl,'(k) --> ',sgnl,'_{yy}(k) in the PIE representation.\n'])
% 
% % Set the parameters of the operator
% Qop.R2x{1} = subs(Qop_y.Q2(:,xx_idcs),yy_o,yy_n);
% Qop.R2y{1} = subs(Qop_x.Q2(:,yy_idcs),xx_o,xx_n).*subs(Top_y.R.R0,yy_o,yy_n);
% Qop.R2y{2} = subs(Qop_x.Q2(:,yy_idcs),xx_o,xx_n).*subs(Top_y.R.R1,yy_o,yy_n);
% Qop.R2y{3} = subs(Qop_x.Q2(:,yy_idcs),xx_o,xx_n).*subs(Top_y.R.R2,yy_o,yy_n);
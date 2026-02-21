function [Top_f,Dvals] = Compute_Tmap_PDE_2D_Separable(PDE,comp_order,suppress_summary)
% Compute opvar2d objects Top, Twop, Tuop defining the map from the
% fundamental state to the PDE state
%   v = Top*vf + Twop*wf + Twop*uf
% where v is the PDE state, vf the fundamental state, wf the fundamental
% exogenous input, and uf the fundamental actuator input. The fundamental
% state will be given by the highest-order mixed spatial derivative of the
% PDE state, and the fundamental input signals may be given by derivatives
% or boundary values of the PDE input signals.
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
% - suppress_summary:   (optional) set to true to suppress the display of
%               an overview of the introduced input signals as they will
%               appear in the PIE.
%
% OUTPUTS:
% - Top_f:      struct with fields:
%               - x:  opvar2d object representing the operator Top mapping 
%                       the fundamental state to PDE state;
%               - w:  opvar2d object representing the operator Twop mapping 
%                       the fundamental exogenous inputs to PDE state;
%               - u:  opvar2d object representing the operator Tuop mapping 
%                       the fundamental actuator inputs to PDE state;
%               - wmap: opvar2d object representing the operator Wop
%                       mapping the set of fundamental exogenous inputs to
%                       a set of derivatives of the PDE exogenous inputs,
%                       Dw*w = Wop*wf, where the order of the derivatives
%                       is specified by Dvals.w; 
%               - umap: opvar2d object representing the operator Uop
%                       mapping the set of fundamental actuator inputs to a
%                       set of derivatives of the PDE actuator inputs,
%                       Du*u = Uop*uf, where the order of the derivatives
%                       is specified by Dvals.u;
% - Dvals:      struct with fields:
%               - w_idcs: nwf x 1 arrray of integers, where nwf is the
%                   number of fundamental exogenous inputs, specifying for
%                   each of these inputs which input in the (reordered) PDE
%                   it corresponds to. The number of fundamental inputs,
%                   nwf, may exceed the number of PDE inputs, nw, in which
%                   case the first nwf-nw fundamental inputs correspond to
%                   boundary values of certain PDE inputs, with the first
%                   nwf-nw elements of w_idcs specifying which;
%               - u_idcs: nuf x 1 array of integers, specifying the same
%                   information as w_idcs, but for the actuator inputs;
%               - x:    nx 2 array of integers, specifying for each of
%                   the fundamental states what order derivative is taken
%                   with respect to x (column 1) and y (column 2) in
%                   computing this fundamental state from the PDE input;
%               - w:    nwf x 2 array of integers, specifying for each of
%                   the fundamental exogenous inputs what order derivative
%                   is taken with respect to x (column 1) and y (column 2)
%                   in computing this fundamental input from the PDE input;
%               - u:    nuf x 2 array of integers, specifying the same
%                   information as w but for the actuator inputs;
%
% NOTES:
% - This routine only works for PDEs with "separable" boundary conditions.
%   That is, there must be no ambiguity whether a boundary condiiton
%   corresponds to a condition along the x-direction or the y-direction,
%   and we should be able to (easily) separate these conditions. 
%   For example, standard Robin-type boundary conditions are accepted,
%       v(a,y)=0, v_{x}(b,y)=0, v(x,c)=w1(x), v_{y}(x,d)+v(x,d)=w2(x),
%   but something with mixed boundary values is not, e.g.
%       v(a,y)+v(x,c)=0
%
% - The routine DOES NOT warn for potential conflicts in the  boundary
%   conditions. The produced relation 
%       v = Top*vf + Twop*wf + Twop*uf
%   holds only if the boundary conditions do not conflict. For example,
%   if we enforce
%       v(a,y)=0,   v(x,c)=w1(x),
%   then we must have w1(a)=0 to avoid a conflict v(a,c) = 0 ~= w1(a).
%
% - For a state variable v differentiable up to order i in x and order j in 
%   y, we must have i boundary conditions along the x-direction, and j
%   boundary conditions along the y-direction.
%
% - The 2D map is constructed by first expanding only along each spatial
%   dimension separately, as
%       v = Top_x*D_x*v +Twop_x*w +Tuop_x*u;
%       v = Top_y*D_y*v +Twop_y*w +Tuop_y*u;
%   If Twop_x=0 and Tuop_x=0, we can then expand v in terms of its
%   highest-order derivative with respect to both spatial variableas as
%       v = Top_x*D_x*v
%         = Top_x*D_x*Top_y*D_y*v +Top_x*D_x*Twop_y*w +Top_x*D_x*Tuop_y*u
%         = Top_x*Top_y*D_xy*v + Top_x*Twop_y*D_x*w + Top_x*Tuop_y*D_x*u
%         = Top*vf + Twop*wf + Tuop*uf
%   where now vf=D_xy*v, wf=D_x*w, and uf=D_x*u. In this case,
%   Dvals.w_idcs=(1:nw), and Dvals.w will specify the order of the
%   derivative defining D_x. Similarly for the case that Twop_y=0 and
%   Tuop_y=0.
%   If Twop_x~=0 and Twop_y~=0, we also expand w in terms of boundary
%   values and the highest-order derivative using Taylor's theorem with
%   integral form of the remainder, so that if y in [c,d],
%       w(x,y) = w(x,c) + (y-c)*(d/dy)w(x,c) + ... 
%                   + (y-c)^n/n!*(d/dy)^n w(x,c) 
%                   + int_{c}^{y}(y-z)^n/n!*(d/dz)^{n+1} w(x,z) dz
%              = Wop*wb
%   where wb = [w(x,c); (d/dy)w(x,c); .; (d/dy)^n w(x,c); (d/dy)^n w(x,.)].
%   Then, (letting u=0 for simplicity) we expand
%       v = Top_x*D_x*v + Twop_x*w
%         = Top_x*D_x*Top_y*D_y*v +Top_x*D_x*Twop_y*w + Twop*Wop*wb
%         = Top*vf + Twop*wf,
%   where now wf is a re-ordered combination of [D_x*w; wb] that places the
%   boundary values from wb in the first nwf-nw positions.
%
% - The fundamental inputs wf and uf may corresponds to boundary values
%   or derivatives of the PDE input signals. Boundary values are always
%   evaluated at the lower boundary of the interval of the respective
%   spatial variable, and are always listed as the first nwf-nw elements of
%   w_idcs (similarly for u_idcs). The routine provides a warning of which
%   inputs have been differentiated and reordered, and what boundary values
%   have been taken.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Initial coding DJ - 04/15/2024
% DJ, 12/16/2024: Fix call to 2D converter;
% DJ, 02/15/2026: Take derivative of all boundary inputs, introducing
%                   boundary values of inputs if necessary;
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
if nargin<=2
    suppress_summary = false;
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
    else
        use_Dx = true;
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
        % Expand the inputs of which a y-derivative should be taken in      % DJ, 02/15/2026
        % terms of boundary values and higher-order y-derivative
        [Wop,Dvals_w,w_idcs] = expand_input(PDE,'w',Dvals_w,2);
        Twop = Twop*Wop;
        if ~suppress_summary
            print_reorder_diff_summary(PDE,'w',comp_order.w,Dvals_w,w_idcs);
        end
    else
        Twop = Twop_x;
        Wop = 1;
        w_idcs = 1:size(Dvals_w,1);
    end
    if numel(PDE.u)>=1
        Tuop = Top_x*Tuop_y +Tuop_x;
        % Expand the inputs of which a y-derivative should be taken in      % DJ, 02/15/2026
        % terms of boundary values and higher-order y-derivative
        [Uop,Dvals_u,u_idcs] = expand_input(PDE,'u',Dvals_u,2);
        Tuop = Tuop*Uop;
        if ~suppress_summary
            print_reorder_diff_summary(PDE,'u',comp_order.u,Dvals_u,u_idcs);
        end
    else
        Tuop = Tuop_x;
        Uop = 1;
        u_idcs = 1:size(Dvals_u,1);
    end
else
    % % Proceed with the expansion
    % % v = Top_y*Tnop_x*Dop_xy*v +Top_y*Twop_x*Dop_x*w +Twop_y*w
    Top = Top_y*Tnop_x;
    if numel(PDE.w)>=1
        Twop = Top_y*Twop_x +Twop_y;
        % Expand the inputs of which a x-derivative should be taken in      % DJ, 02/15/2026
        % terms of boundary values and higher-order x-derivative
        [Wop,Dvals_w,w_idcs] = expand_input(PDE,'w',Dvals_w,1);
        Twop = Twop*Wop;
        if ~suppress_summary
            print_reorder_diff_summary(PDE,'w',comp_order.w,Dvals_w,w_idcs);
        end
    else
        Twop = Twop_y;
        Wop = 1;
        w_idcs = 1:size(Dvals_w,1);
    end
    if numel(PDE.u)>=1
        Tuop = Top_y*Tuop_x +Tuop_y;
        % Expand the inputs of which a x-derivative should be taken in      % DJ, 02/01/2026
        % terms of boundary values and higher-order x-derivative
        [Uop,Dvals_u,u_idcs] = expand_input(PDE,'u',Dvals_u,1);
        Tuop = Tuop*Uop;
        if ~suppress_summary
            print_reorder_diff_summary(PDE,'u',comp_order.u,Dvals_u,u_idcs);
        end
    else
        Tuop = Tuop_y;
        Uop = 1;
        u_idcs = 1:size(Dvals_u,1);
    end
end

% Collect the operators defining the map from fundamental state to PDE
% state in a struct
Top_f = struct();
Top_f.x = Top;              Top_f.u = Tuop;      Top_f.w = Twop;
Top_f.umap = Uop;           Top_f.wmap = Wop;
% Collect the parameters defining the map from PDE state to fundamental
% state in a struct
Dvals = struct();
Dvals.x = Dvals_xy;         Dvals.u = Dvals_u;          Dvals.w = Dvals_w;
Dvals.u_idcs = u_idcs;      Dvals.w_idcs = w_idcs;

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
try PIE_1D = convert_PIETOOLS_PDE(PDE_1D,[],{'silent','Top'});                      % DJ, 12/16/2024
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
function [Pop,Dvals_full,obj_idcs] = expand_input(PDE,obj,Dvals,dir)
% Expand the input signal specified by 'obj' (exogenous or actuator input)
% in terms of its highest-order spatial derivative and lower-boundary
% values along direction "dir"
%
% INPUTS
% - PDE:    'pde_struct' object specifying the PDE for which to expand the
%           inputs;
% - obj:    one of 'w' or 'u', specifying whether to consider the exogenous
%           inputs or actuator inputs;
% - Dvals:  n_obj x 2 array specifying for each element of PDE.(obj) what
%           order derivative is taken of this element in computing the
%           associated input in the PIE, along the first (column 1) and
%           second (column 2) spatial variable;
% - dir:    one of {1,2}, specifying whether to expand the input along the
%           first spatial variable or the second;
%
% OUTPUTS
% - Pop:    n_obj x n_objf 'opvar2d' object, representing the 2D PI
%           operator mapping the set of fundamental input signals back to
%           the original input. If the input is not differentiated along
%           direction 'dir', a value of Pop = 1 is returned, representing 
%           the identity operator;
% - Dvals_full: n_objf x 2 array of integers, specifying for each of the
%           fundamental input variables what order derivative is taken
%           along each of the spatial directions in computing this
%           variable;
% - obj_idcs:   n_objf x 1 array of integers specifying for each of the
%           fundamental input variables which element of PDE.(obj) they
%           correspond to. If n_objf > n_obj, the first n_objf-n_obj
%           elements always correspond to lower-boundary values of the
%           input signals specified in PDE.(obj)(obj_idcs);
%
% EXAMPLE
% Suppose dir = 1, and we have PDE.vars(1,1)=x with PDE.dom(1,:)=[a,b].
% Suppose that Dvals(:,dir)=Dvals(:,1) = [0;0;2;0], and obj = 'w'. Then, we
% expand the third exogenous input in terms of its 2nd-order derivative as
%   w3(x) = w3(a) + (x-a)*w3_{x}(a) + int_{a}^{x} (x-z)*w3_{xx}(z) dz   (*)
% We introduce the fundamental state vector
%   wf = [w3(a); w3_{x}(a); w1; w2; w3; w4]
% and use relation (*) to define the operator Pop such that w = Pop*wf. The
% function returns
%   Dvals_full = [ [0;1;Dvals(:,1)], [0;0;Dvals(:,2)]];
%   obj_idcs = [3;3;1;2;3;4];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 02/01/2026
%

% Extract necessary information from the PDE
var1 = PDE.vars(dir,1);     var2 = PDE.vars(dir,2);
dom = PDE.dom(dir,:);       nvars = PDE.dim;
obj_tab = PDE.([obj,'_tab']);
nobj_op = get_opdim(obj_tab);

% Establish which inputs have to be expanded in terms of boundary values
% and highest-order derivative
expand_idcs = find(Dvals(:,dir));
if isempty(expand_idcs)
    Pop = 1;
    Dvals_full = Dvals;
    obj_idcs = 1:size(obj_tab,1);
    return
end

% Establish the number of boundary values that need to be defined
nobj_arr = obj_tab(:,2);
nobj_arr_cum = cumsum([0;nobj_arr]);
nbval_arr = Dvals(expand_idcs,dir).*nobj_arr(expand_idcs);
nbval_arr_cum = cumsum([0;nbval_arr]);
nbvals = nbval_arr_cum(end);

% Keep track of which boundary value corresponds to which input
obj_idcs = (1:size(Dvals,1))';
obj_idcs = [repelem(obj_idcs(expand_idcs),Dvals(expand_idcs,dir),1); obj_idcs];

% Keep track of what derivative is taken of the boundary value
Dvals_B = zeros(size(nbval_arr,1),nvars);

% Initialize operators Bop and Cop mapping the expanded input back to the 
% original:
%   w(x) = Bop*[w(a); w_{x}(a)] + Cop*w_{xx}(.)
%        = [1, (x-a)]*[w(a);w_{x}(a)] + int_{a}^{x}(x-z) w_{xx}(z) dz
tmp_opvar = opvar2d();
tmp_opvar.var1 = PDE.vars(:,1);      
tmp_opvar.var2 = PDE.vars(:,2);
tmp_opvar.I = PDE.dom;
Bop = tmp_opvar;
Bop.dim = [nobj_op, [nbvals;0;0;0]];
% Initialize Cop as identity operator
Cop = mat2opvar(eye(sum(nobj_arr)),[nobj_op,nobj_op],PDE.vars,PDE.dom);

% For each of the inputs to expand, set the appropriate elements of Bop and
% Cop to express this input in terms the full input vector
if dir==1
    param_name = 'Rxx';
else
    param_name = 'Ryy';
end
strt_idx = 0;
for i=1:numel(expand_idcs)
    % Set the factors multiplying each boundary value based on Taylor's Thm
    idx = expand_idcs(i);
    dval = Dvals(idx,dir);
    fctrs = (var1-dom(1)).^(0:dval-1)./factorial(0:dval-1);
    % Account for the size of the input vectors
    Ri = kron(fctrs,eye(nobj_arr(idx)));
    % Determine which element of the original vector belong to this input
    ridcs = nobj_arr_cum(idx)+1:nobj_arr_cum(idx+1);
    cidcs = nbval_arr_cum(i)+1:nbval_arr_cum(i+1);
    % Set the new parameter
    Bop(ridcs,cidcs) = Ri;
    Dvals_B(strt_idx+(1:dval),dir) = 0:dval-1;
    strt_idx = strt_idx + dval;

    % Also replace the identity operator in Pop by the appropriate integral
    Pop_i = tmp_opvar;
    Pop_i.(param_name){2} = kron((var1-var2)^(dval-1)./factorial(dval-1),eye(nobj_arr(idx)));
    Cop(ridcs,ridcs) = Pop_i;
end

% Return the full operator from fundamental input to original input
Pop = [Bop,Cop];
Dvals_full = [Dvals_B; Dvals];

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_reorder_diff_summary(PDE,obj,comp_order,Dvals,obj_idcs)
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
% - Dvals:  A m x nvars object of 'type' double specifying for each
%           element PDE.obj(i) the order of the derivative taken with
%           respect to each of the nvars global variables in PDE.vars, as
%           the object appears in the PIE representation. If the value of m
%           exceeds the number of elements of PDE.obj, the first elements
%           of Dvals are taken to correspond to boundary values of some
%           object;
% - obj_idcs:   A m x 1 array of integers specifying for each row of Dvals
%               to which element of PDE.obj it corresponds, which is
%               necessary in case boundary values of the object have been
%               introduced;
%
% OUTPUTS:
% Displays information in the command window concerning the new order of
% the components in PDE.obj, as well as the derivatives taken of those
% components.
%
% NOTES
% If any boundary values are taken, these are assumed to be specified by
% the first m-ncomps rows of Dvals and obj_idcs, where ncomps =
% numel(PDE.(obj)), and where m = size(Dvals,1);
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
obj_tab = PDE.([obj,'_tab']);
ncomps = size(Dvals,1);
if ~any(any(Dvals)) && ncomps==numel(PDE.(obj))
    % No derivative of the components is taken
    %fprintf(['\n','The order of the ',object_name,'s ',obj,' has not changed.\n']);
    return
else
    % Otherwise, we list the new order of the components.
    %fprintf(['\n','The ',object_name,'s have been reordered as:\n']);
    fprintf(['\n','WARNING: Some ',obj_name,' from the PDE will appear as derivatives (and boundary values) in the PIE representation:\n']);
end

% Use UNICODE to add subscript indices to different components.
%sub_s = '\x209B';
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
    new_idx = obj_idcs(ii);
    old_idx = comp_order(new_idx);
    
    % Check if the component corresponds to a boundary value of an input
    % signal
    is_bval = ii<=ncomps-numel(PDE.(obj));
    if is_bval
        % Check if the boundary value is along x- or y-direction
        dir = logical(obj_tab(new_idx,3:2+nvars));
        % Set the actual boundary value: always the lower boundary
        Lval = PDE.dom(dir,1);
    end
    comp_ii = PDE.(obj){new_idx};
    
    % Establish the names of the variables on which the component depends.
    if is_bval
        varnames_ii = num2str(Lval);
        varnames_LHS_t = ['t,',varnames_ii]; % the boundary values still vary in time
        varnames_RHS_t = 't';
    elseif isempty(comp_ii.vars)
        varnames_LHS_t = 't';    % All components may vary in time
        varnames_RHS_t = varname_LHS_t;
    else
        % Set a comma separated list.
        varnames_ii = comp_ii.vars(1,1).varname{1};
        for kk=2:size(comp_ii.vars,1)
            varnames_ii = [varnames_ii,',',comp_ii.vars(kk,1).varname{1}];
        end
        varnames_ii_t = ['t,',varnames_ii]; % All components may vary in time
        varnames_LHS_t = varnames_ii_t;
        varnames_RHS_t = varnames_ii_t;
    end
    
    % Establish the (subscript) index for the component.
    LHS_length = length(varnames_LHS_t);
    if isscalar(PDE.(obj))
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

    % Establish which derivative of the component is taken.
    dval = Dvals(ii,:);
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
    end

    % Set the name of the component, including its depdence on spatial
    % variables.
    LHS_name = [' ',Ldiff,obj,Lcomp_idx,'(',varnames_LHS_t,')'];
    LHS_length = 1 + LHS_length + 3;
        
    % For the added state components, indicate of which state component
    % they are the temporal derivative.
    Rcomp_idx = cell2mat(sub_num(str2num(num2str(ii)')+1)');
    RHS = [' -->   ',obj,Rcomp_idx,'(',varnames_RHS_t,')'];
%    RHS = '';
    
    % % % Finally, display:
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
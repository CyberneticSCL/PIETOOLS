classdef npie_struct
% npie_struct is a class of objects used to describe nonlinear PIE systems
% in PIETOOLS, taking the form
%
%   (d/dt)* (Tu*u + Tw*w + T*x)(t) = f(x(t))  + Bw*w(t)  + Bu*u(t);
%                             z(t) = gz(x(t)) + Dzw*w(t) + Dzu*u(t);
%                             y(t) = gy(x(t)) + Dyw*w(t) + Dyu*u(t);
%       Ku*u(t) + Kw*w(t) + K*x(t) = 0
%
% for distributed polynomials f, gz, and gy in the fundamental state x,
% where w is the exogenous input, u the actuator input, z the regulated
% output, and y the observed output. For linear PIEs, see `pie_struct`.
%
%
% === The Equations ===
%
% - PIE.T, PIE.Tw, PIE.Tu, PIE.Bw, PIE.Bu, PIE.Dzw, PIE.Dzu, PIE.Dyw,
%   PIE.Dyu:
%               tensopmat objects representing the linear operators in the
%               equations above. Each operator maps to the same function
%               space as the other operators in its row, e.g. T and Bw map
%               to the same space, and Dyu maps u to y. Operators that have
%               not been specified are initialized as zero operators;
% - PIE.f, PIE.gz, PIE.gy:
%               polyopvar objects (never tensopmat objects) representing
%               the distributed polynomials in the state x, replacing the
%               operators A, Cz, and Cy of a linear PIE, in terms of the
%               state variables named by PIE.space.x.varname. A polynomial
%               need not depend on all components of x. Only f is
%               mandatory, gz and gy are initialized as zero polynomials;
% - PIE.K, PIE.Kw, PIE.Ku:
%               tensopmat objects representing the algebraic constraint on
%               the fundamental state. After initialize, these have no rows
%               if the PIE has no constraint.
%
%
% === The Function Spaces ===
%
% The states, inputs, and outputs may consist of multiple, vector-valued
% components, each in a function space of the form
%       L2^{ni}[[a1,b1]^{di1} x ... x [aN,bN]^{diN}],
% with dij in {0,1} and L2^{n}[\emptyset] = R^n. These spaces are described
% by
%
% - PIE.space:  struct with fields
%        vars:  1 x N cellstr with the names of the N spatial variables
%               (dummy variables are not stored);
%         dom:  N x 2 array with the interval of each spatial variable;
%    x, y, z, w, u:
%               structs describing the state, observed output, regulated
%               output, exogenous input, and actuator input components,
%               each with fields
%         dims: ncomps x 1 array with the vector size ni of each component;
%       depmat: ncomps x N logical array with element (i,j) true if
%               component i depends on variable j. A row of all-false
%               indicates a finite-dimensional component;
%      varname: ncomps x 1 cellstr with the (unique) name of each
%               component, used as the name of the state variable in
%               polynomials, e.g. 'x1_pie'. If not specified, `initialize`
%               assigns the names [type,num,'_pie'], with zero-padded
%               numbers;
%
% The state variables, as a `polyopvar` object, are available as
%
% - PIE.x:      (dependent) 'polyopvar' with variables named
%               PIE.space.x.varname, of the sizes PIE.space.x.dims,
%               depending on the variables PIE.space.vars as indicated by
%               PIE.space.x.depmat;
%
%
% === The PIE to PDE Bijection ===
%
% A PIE constructed from a PDE (see `convert_PIETOOLS_PDE_nonlinear`) keeps
% track of the equivalence between the PDE state space and the PIE state
% space,
%
%       x(t) = Dx*x_pde(t) + Dxw*w_pde(t) + Dxu*u_pde(t);
%       w(t) = Dw*w_pde(t);         w_pde(t) = W_pde*w(t);
%       u(t) = Du*u_pde(t);         u_pde(t) = U_pde*u(t);
%       x_pde(t) = X_pde*(T*x(t) + Tw*w(t) + Tu*u(t));
%       y_pde(t) = Y_pde*(gy(x(t)) + Dyw*w(t) + Dyu*u(t));
%       z_pde(t) = Z_pde*(gz(x(t)) + Dzw*w(t) + Dzu*u(t))
%
% where Dx, Dw, and Du are differential/trace operators, possibly composed
% with integral operators and delay maps, and X_pde, Y_pde, Z_pde, W_pde,
% and U_pde are PI operators, using the fields
%
% - PIE.space_pde:
%               struct of the same form as PIE.space for the PDE, whose
%               variables may differ from those of the PIE (e.g. delays).
%               Each of the fields x, y, z, w, u has a field varname with
%               the name of each PDE component (default [type,num],
%               assigned by `initialize`), and a field ID, identifying the
%               PDE components as in PDE.x{i}.ID, so that row i of ID,
%               dims, and depmat describes the same component;
% - PIE.pie2pde:
%               struct with fields x, y, z, u, w, with tensopmat objects
%               representing X_pde, Y_pde, Z_pde, U_pde and W_pde, with a
%               column for each PIE component and a row for each PDE
%               component. For the state, the columns correspond to the
%               components of T*x+Tw*w+Tu*u. Y_pde and Z_pde are
%               permutation matrices, accounting for a reordering of the
%               outputs, and U_pde and W_pde are often identity or
%               permutation matrices;
% - PIE.pde2pie:
%               struct with fields x, u, w, y, z, describing for each PIE
%               component how it is obtained from a PDE component. Each
%               field is a struct with fields (y and z have only idx)
%        idx:   ncomps x 3 (for x) or ncomps x 1 array, with the row in
%               PIE.space_pde.(type).ID of the associated PDE component.
%               For x, the nonzero column indicates whether this is a PDE
%               state (1), disturbance (2), or actuator input (3).
%               Different PIE components may use the same PDE component;
%       diff:   ncomps x N array with the derivative taken with respect to
%               each spatial variable;
%      tdiff:   ncomps x 1 array with the temporal derivative taken;
%        loc:   ncomps x N array in {-1,0,1}, with the location at which
%               the PDE component is evaluated, wrt each spatial variable:
%               lower boundary (-1), interior (0), or upper boundary (1);
%       delay:  (x only) ncomps x 1 array in 0:N. If delay(i)=k>0,
%               component i is the history of a PDE component, depending on
%               the k-th spatial variable on domain [-tau,0], so that
%                   x_i(t,s_k) = Dop x_{pde,idx(i)}(t+s_k);
%          C:   tensopmat object with a row for each PIE component, applied
%               after the derivatives and evaluation, e.g. multiplication
%               by a coefficient, or an integral operator. Otherwise the
%               identity.
%
% So PIE component i is obtained from PDE component idx(i) by (1) taking
% the derivatives diff(i,:) and tdiff(i), (2) evaluating at loc(i,:), and
% at time t+s_k if delay(i)=k>0, and (3) applying C.
%
% EXAMPLE: For the PDE u_t = u_ss on [0,1], with BCs u(0)=0 and u_s(1)=0,
% the fundamental state is the single component v = u_ss, with
%       pde2pie.x:  idx=[1,0,0],  diff=2,  tdiff=0,  loc=0,  delay=0,  C=1,
% i.e. v is PDE state 1, differentiated twice in s, evaluated in the
% interior, and left unchanged by C. Conversely, T is the integral operator
% with kernel -min(s,theta), and pie2pde.x = X_pde = 1, so that
%    u(t,s) = X_pde*(T*v)(t,s) = -int_0^1 min(s,theta) v(t,theta) dtheta.
% A component d/dt u would instead have tdiff=1, and the boundary value
% u(t,1) would have diff=0 and loc=1.
%
% The relation is a bijection only between admissible trajectories: x_pde
% must satisfy the BCs, i.e. be of the form X_pde*(T*x+Tw*w+Tu*u), and, if
% K is nontrivial, Ku*u + Kw*w + K*x = 0. For the components with tdiff=0
% and delay=0 this implies Dx*X_pde*T = I. If the PIE is not constructed
% from a PDE, the fields space_pde, pie2pde and pde2pie have no components.
%
%
% === Other Fields ===
%
% - PIE.misc:   struct in which any other objects relevant to the PIE can
%               be stored;
% - PIE.is_initialized:
%               logical value indicating whether the PIE has been
%               initialized using `initialize`, so that it is well-defined
%               and consistent. Set to false when a field of the PIE is
%               adjusted;
%
% For backward compatibility, the hidden fields B1, B2, D11, ..., D22
% (aliases for Bw, ..., Dyu), g1, g2 (aliases for gz, gy), vars, dom, dim,
% and x_tab, w_tab, u_tab, y_tab, z_tab can still be used to extract (and
% set) the associated information in the format of earlier versions,
% computed from PIE.space and PIE.pde2pie.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - npie_struct
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
% DJ, 10/29/2026: Initial coding;

properties
    % The dynamics
    T;      Tu;     Tw;
    f;      Bu;     Bw;
    gy;     Dyu;    Dyw;
    gz;     Dzu;    Dzw;
    K;      Ku;     Kw;
    % Function space of the PIE and of the PDE
    space;
    space_pde;
    % Maps between the PIE and the PDE (declared in the constructor)
    pie2pde;
    pde2pie;
    % "misc" field to add any other information relevant to the PIE.
    misc;
    % Has the PIE already been initialized (set by initialize, reset when
    % the PIE is adjusted).
    is_initialized = false;
end

properties (Hidden,Dependent)
    % Allow e.g. Bw to be extracted as B1, for backward compatibility.
    B1;     B2;
    D11;    D12;
    D21;    D22;
    g1;     g2;
    % The PIE state x, constructed from the function space of the PIE.
    x;
    % The spatial variables and domain, and the tables with standardized
    % information on each of the components, are also included for
    % backward compatibility.
    vars;
    dom;
    x_tab;  % [state_num, state_size, isvar1, isvar2, vardiff1, vardiff2]
    u_tab;
    w_tab;
    y_tab;
    z_tab;
end

properties (Dependent)
    % Spatial dimension of the PIE.
    dim;
end

methods
function obj = npie_struct(varargin)
    %NPIE_STRUCT Construct an instance of this class
    %   Use e.g. npie_struct PIE1 PIE2 to declare multiple empty PIEs, or
    %   convert a struct with fields named as the properties (including the
    %   hidden fields used in earlier versions) to an npie_struct.
    % Declare the fields of the function spaces and the maps between the
    % PIE and PDE, without any components.
    obj.space = empty_space(false);
    obj.space_pde = empty_space(true);
    obj.pie2pde = empty_pie2pde();
    obj.pde2pie = empty_pde2pie();
    obj.misc = struct();

    if nargin==1
        if ischar(varargin{1})
            if nargout==0
                assignin('caller', varargin{1}, npie_struct());
            end
        elseif isa(varargin{1},'npie_struct')
            obj = varargin{1};
        elseif isa(varargin{1},'struct')
            % If a struct is provided, convert to an "npie_struct".
            % Copy information from any field in the struct that matches a
            % property of the npie_struct. The spaces, variables and domain
            % are set first, as the tables depend on them.
            fnames = fieldnames(varargin{1})';
            first = {'space','space_pde','vars','dom'};
            last = {'x_tab','u_tab','w_tab','y_tab','z_tab'};
            order = [intersect(first,fnames,'stable'),...
                        setdiff(fnames,[first,last,{'dim','is_initialized'}],'stable'),...
                        intersect(last,fnames,'stable')];
            for ii=1:numel(order)
                if isprop(obj,order{ii})
                    obj.(order{ii}) = varargin{1}.(order{ii});
                end
            end
        else
            error("Input must be strings");
        end
    else
        for i=1:nargin
            if ischar(varargin{i})
                % Allow new npie_structs to be defined using
                % > npie_struct PIE1 PIE2 PIE3 ...
                if nargout==0
                    assignin('caller', varargin{i}, npie_struct());
                end
            else
                error("Input must be strings");
            end
        end
    end
end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% Operators named as in earlier versions
function val = get.B1(obj),         val = obj.Bw;           end
function obj = set.B1(obj,val),     obj.Bw = val;           end
function val = get.B2(obj),         val = obj.Bu;           end
function obj = set.B2(obj,val),     obj.Bu = val;           end
function val = get.g1(obj),         val = obj.gz;           end
function obj = set.g1(obj,val),     obj.gz = val;           end
function val = get.g2(obj),         val = obj.gy;           end
function obj = set.g2(obj,val),     obj.gy = val;           end
% The PIE state (polyopvar)
function val = get.x(obj)
    % Construct the PIE state from the function space of the PIE. If the
    % components are not yet named, they are assigned default names.
    sp = obj.space.x;
    if isempty(sp.dims)
        val = polyopvar();
        return
    end
    if numel(sp.varname)~=numel(sp.dims)
        obj = initialize(obj);
        sp = obj.space.x;
    end
    val = polyopvar(sp.varname,sp.dims,obj.space.vars,obj.space.dom,sp.depmat);
end

function val = get.D11(obj),        val = obj.Dzw;          end
function obj = set.D11(obj,val),    obj.Dzw = val;          end
function val = get.D12(obj),        val = obj.Dzu;          end
function obj = set.D12(obj,val),    obj.Dzu = val;          end
function val = get.D21(obj),        val = obj.Dyw;          end
function obj = set.D21(obj,val),    obj.Dyw = val;          end
function val = get.D22(obj),        val = obj.Dyu;          end
function obj = set.D22(obj,val),    obj.Dyu = val;          end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% Spatial variables, domain, and dimension
function val = get.vars(obj),       val = pie_compat_vars(obj);         end
function obj = set.vars(obj,val),   obj = pie_compat_vars(obj,val);     end
function val = get.dom(obj),        val = obj.space.dom;                end
function obj = set.dom(obj,val)
    if isempty(val)
        val = zeros(0,2);
    end
    obj.space.dom = val;
end
function dim_val = get.dim(obj)
    % Get spatial dimension of PIE.
    if size(obj.space.dom,1)==numel(obj.space.vars)
        % Each variable should be assigned a domain.
        dim_val = numel(obj.space.vars);
    else
        dim_val = nan;
    end
end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% Tables with information on the components
function val = get.x_tab(obj),      val = pie_compat_tab(obj,'x');      end
function obj = set.x_tab(obj,val),  obj = pie_compat_tab(obj,'x',val);  end
function val = get.u_tab(obj),      val = pie_compat_tab(obj,'u');      end
function obj = set.u_tab(obj,val),  obj = pie_compat_tab(obj,'u',val);  end
function val = get.w_tab(obj),      val = pie_compat_tab(obj,'w');      end
function obj = set.w_tab(obj,val),  obj = pie_compat_tab(obj,'w',val);  end
function val = get.y_tab(obj),      val = pie_compat_tab(obj,'y');      end
function obj = set.y_tab(obj,val),  obj = pie_compat_tab(obj,'y',val);  end
function val = get.z_tab(obj),      val = pie_compat_tab(obj,'z');      end
function obj = set.z_tab(obj,val),  obj = pie_compat_tab(obj,'z',val);  end

end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function sp = empty_space(with_ID)
% SP = EMPTY_SPACE(WITH_ID) returns the struct describing the function
% spaces of a PIE or PDE without any variables or components.
%
% INPUT
% - with_ID:    logical, set to true to add a field ID to each of the
%               components, as for the PDE function spaces. Each of the
%               components has a field varname, for the names of the
%               components;
%
% OUTPUT
% - sp:         struct with fields vars, dom, x, y, z, w, u, as described
%               in the header of npie_struct;
sp = struct('vars',{cell(1,0)},'dom',zeros(0,2));
comps = {'x','y','z','w','u'};
for ii=1:numel(comps)
    sp.(comps{ii}) = struct('dims',zeros(0,1),'depmat',false(0,0),'varname',{cell(0,1)});
    if with_ID
        sp.(comps{ii}).ID = zeros(0,1);
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function m = empty_pie2pde()
% M = EMPTY_PIE2PDE() returns the struct describing the map from PIE to PDE
% states, inputs and outputs, without any components, with an empty
% tensopmat object for each of the fields x, y, z, u, w (see the header of
% npie_struct).                                                           
m = struct('x',tensopmat(),'y',tensopmat(),'z',tensopmat(),...
           'u',tensopmat(),'w',tensopmat());
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function m = empty_pde2pie()
% M = EMPTY_PDE2PIE() returns the struct describing the map from PDE to PIE
% states, inputs and outputs, without any components. Each of the fields x,
% u, w, y, z is a struct with the fields specified in the header of
% npie_struct: idx for all components, diff, tdiff, loc and C for x, u and
% w, and delay for x only.                                                 
full_map = struct('idx',zeros(0,1),'diff',zeros(0,0),'tdiff',zeros(0,1),...
                  'loc',zeros(0,0),'C',tensopmat());
m = struct();
m.x = full_map;
m.x.idx = zeros(0,3);                   % state, exogenous input, or actuator input
m.x.delay = zeros(0,1);
m.u = full_map;
m.w = full_map;
m.y = struct('idx',zeros(0,1));
m.z = struct('idx',zeros(0,1));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

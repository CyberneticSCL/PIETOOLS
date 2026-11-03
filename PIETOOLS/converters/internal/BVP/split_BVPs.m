function [PDEs,info] = split_BVPs(PDE)
% [PDEs,info] = SPLIT_BVPS(PDE) decomposes the boundary conditions (BCs) of
% an ND PDE along each of its N spatial directions, returning a 1D
% boundary value problem (BVP) for each direction s_i, so that they can be
% solved separately. For example, the conditions
%   u(0,y)=u(1,y)=0,        u(x,0)=u_{y}(x,1)=0
% on u(x,y) are split into the 1D problems
%   u(0)=u(1)=0     along x,        u(0)=u_{y}(1)=0     along y.
% Only the BCs are considered; the dynamics are ignored.
%
% INPUT
% - PDE:    'pde_struct' object in the terms format, with N spatial
%           variables. Each BC must be a condition along exactly one
%           direction, see NOTES;
%
% OUTPUT
% - PDEs:   1 x N cell, with PDEs{i} a 'pde_struct' defining the 1D BVP
%           along PDE.vars(i,1): its states are the components of PDE that
%           depend on s_i, as functions of s_i only, and its BCs those of
%           PDE along s_i, restricted to a slice (fixed s_{~=i}), with
%           only the terms of states that depend on s_i. The other terms
%           are moved to the forcing. PDEs{i} has no equations or inputs;
% - info:   struct with fields
%           - BC_dir:   nBC x 1 array, the direction of each BC of PDE;
%           - BC_idcs:  1 x N cell, with PDEs{i}.BC{k} the BC
%                       PDE.BC{BC_idcs{i}(k)} restricted to the slice;
%           - x_idcs:   1 x N cell, with PDEs{i}.x{k} the component
%                       PDE.x{x_idcs{i}(k)};
%           - forcing:  1 x N cell, with forcing{i}{k} a cell of terms (in
%                       the format of PDE.BC{*}.term) such that
%                       (state terms of BC k) = sum(forcing{i}{k}).
%                       These are the terms of the original BC that do not
%                       involve a state depending on s_i (inputs, and
%                       states of fewer variables, e.g. ODE states), with
%                       C negated. They give the nonhomogeneous term w of
%                       compute_Tmap_1D;
%
% NOTES
% - A BC is a condition along s_i if it fixes s_i, by evaluating a state at
%   a boundary or integrating over s_i, and varies in the other variables
%   of the states it involves, as in
%       x(a,y)=0,       x_{x}(b,y)+int_{c}^{d} x(b,z)dz = w(y).
%   Not supported (error): conditions that fix several variables (corners
%   or edges), that fix s_i but also vary in s_i, that involve derivatives
%   of a state along directions it is not fixed in, with coefficients
%   varying along the boundary for states depending on s_i (coefficients
%   of the forcing are unrestricted), that couple states depending on s_i
%   on different sets of variables, or that involve delays.
% - States of fewer variables, including ODE states, may enter the BCs of
%   other states through the forcing. The reverse coupling only enters the
%   dynamics, which are ignored here.
% - This routine does not check that the 1D problems, or their combination,
%   are admissible (e.g. the number of conditions, or well-posedness);
%   see solve_BVP_1D.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - split_BVPs
%
% Copyright (C) 2026 PIETOOLS Team
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
% DJ (CC), 11/02/2026: Initial coding;

% % % Make sure the PDE is initialized
PDE = pde_struct(PDE);
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
N = size(PDE.vars,1);
if N==0
    error("The PDE does not depend on any spatial variables.")
end
nBC = numel(PDE.BC);

% Which variables does each component depend on, and what are the orders
x_vars = logical(PDE.x_tab(:,3:2+N));
w_vars = logical(PDE.w_tab(:,3:2+N));
u_vars = logical(PDE.u_tab(:,3:2+N));
x_diff = PDE.x_tab(:,2+N+(1:N));
dep_vars = struct('x',x_vars,'w',w_vars,'u',u_vars);
var_names = cell(N,2);
for j=1:N
    var_names{j,1} = PDE.vars(j,1).varname{1};
    var_names{j,2} = PDE.vars(j,2).varname{1};
end


% % % Determine to which direction each condition belongs, and check that it
% % % is a condition that can be imposed along this direction (see NOTES)
types = {'x','w','u'};
BC_dir = zeros(nBC,1);
for ii=1:nBC
    BC = PDE.BC{ii};
    BC_name = ['BC{',num2str(ii),'}'];
    varies_in = logical(PDE.BC_tab(ii,3:2+N));      % Variables in which the BC function varies

    % % Determine which variables are fixed by the terms in the condition
    fixed = false(1,N);
    for jj=1:numel(BC.term)
        term = BC.term{jj};
        obj_type = types{find(isfield(term,types),1)};
        comp_vars = dep_vars.(obj_type)(term.(obj_type),:);
        pos = cumsum(comp_vars);
        for j=find(comp_vars)
            kk = pos(j);
            % A variable is fixed if the term is evaluated at a boundary,
            % or integrated over the full domain.
            if strcmp(obj_type,'x') && (isa(term.loc(kk),'double') || isdouble(polynomial(term.loc(kk))))
                fixed(j) = true;
            elseif numel(term.I)>=kk && ~isempty(term.I{kk})
                I_lim = polynomial(term.I{kk});
                if isdouble(I_lim)
                    I_lim = double(I_lim);
                    fixed(j) = fixed(j) || (numel(I_lim)==2 && isequal(sort(I_lim(:))',sort(PDE.dom(j,:))));
                end
            end
        end
    end
    if ~any(fixed)
        error('split_BVPs:noBoundary',['The condition ',BC_name,' does not evaluate any state at a boundary, and is not supported as a boundary condition along any spatial direction.'])
    elseif nnz(fixed)>1
        error('split_BVPs:multipleDirections',['The condition ',BC_name,' fixes the variables ',strjoin(var_names(fixed,1)',', '),' simultaneously (e.g. at a corner or edge); it cannot be separated into conditions along a single spatial direction.'])
    end
    i = find(fixed);
    if varies_in(i)
        error('split_BVPs:mixedDirection',['The condition ',BC_name,' involves a term that is evaluated at a boundary in ',var_names{i,1},', but also a term that varies in ',var_names{i,1},'; it cannot be separated into conditions along a single spatial direction.'])
    end

    % % Check the terms
    % Determine which variables the states that depend on s_i vary in
    vary_vars = [];
    for jj=1:numel(BC.term)
        term = BC.term{jj};
        obj_type = types{find(isfield(term,types),1)};
        comp_vars = dep_vars.(obj_type)(term.(obj_type),:);
        if strcmp(obj_type,'x') && comp_vars(i)
            if isfield(term,'delay') && any(term.delay(:)~=0)
                error('split_BVPs:delay',['The condition ',BC_name,' involves a delay of a state; delays are not supported.'])
            end
            if isempty(vary_vars)
                vary_vars = comp_vars;
            elseif ~isequal(vary_vars,comp_vars)
                error('split_BVPs:mixedDomains',['The condition ',BC_name,' couples states that depend on ',var_names{i,1},' but vary in different sets of variables; this would couple the problems along ',var_names{i,1},' for different values of the other variables.'])
            end
        end
    end
    if isempty(vary_vars)
        error('split_BVPs:noState',['The condition ',BC_name,' does not involve any state that depends on ',var_names{i,1},', and cannot be a condition along this direction.'])
    end
    if any(varies_in & ~vary_vars)
        error('split_BVPs:extraVars',['The condition ',BC_name,' varies in variables ',strjoin(var_names(varies_in & ~vary_vars,1)',', '),', on which the states it constrains do not depend.'])
    end

    % Check the derivatives and coefficients of the states depending on s_i
    % (all other terms are moved to the forcing, which is not restricted)
    for jj=1:numel(BC.term)
        term = BC.term{jj};
        obj_type = types{find(isfield(term,types),1)};
        comp_vars = dep_vars.(obj_type)(term.(obj_type),:);
        pos = cumsum(comp_vars);
        if strcmp(obj_type,'x') && comp_vars(i)
            % Derivatives of boundary values along other directions are not
            % separable
            D_tangent = term.D(pos(comp_vars & ~((1:N)==i)));
            if any(D_tangent)
                error('split_BVPs:tangentialDerivative',['The condition ',BC_name,' involves a derivative of a state along a direction other than ',var_names{i,1},'; this is not supported.'])
            end
            % The coefficients cannot vary along the boundary, as they
            % define the 1D problem, except on the variable of integration
            % of an integral over s_i (s_i or its dummy variable)
            allowed = cell(1,0);
            if numel(term.I)>=pos(i) && ~isempty(term.I{pos(i)})
                allowed = var_names(i,:);
            end
            if isa(term.C,'polynomial') && ~isdouble(term.C) && any(~ismember(term.C.varname,allowed))
                error('split_BVPs:varyingCoefficient',['The condition ',BC_name,' involves a state depending on ',var_names{i,1},' with coefficients that vary along the boundary; this is not supported.'])
            end
        end
    end
    BC_dir(ii) = i;
end


% % % Build the 1D problem along each direction, from the conditions along
% % % it, and collect the forcing of each condition (the terms moved to its
% % % right-hand side)
PDEs = cell(1,N);
info = struct('BC_dir',BC_dir,'BC_idcs',{cell(1,N)},'x_idcs',{cell(1,N)},'forcing',{cell(1,N)});
for i=1:N
    var_i = PDE.vars(i,:);          dom_i = PDE.dom(i,:);
    x_idcs = find(dep_vars.x(:,i))';
    BC_idcs = find(BC_dir==i)';
    nx_i = numel(x_idcs);
    x_new_idx = zeros(numel(PDE.x),1);   x_new_idx(x_idcs) = 1:nx_i;

    % Declare the states and conditions
    PDE_struct = struct();
    PDE_struct.vars = var_i;         PDE_struct.dom = dom_i;
    PDE_struct.x = cell(nx_i,1);
    for k=1:nx_i
        PDE_struct.x{k} = struct('size',PDE.x{x_idcs(k)}.size,'vars',var_i(1),'dom',dom_i,'diff',x_diff(x_idcs(k),i));
    end
    PDE_struct.BC = cell(numel(BC_idcs),1);
    forcing = cell(1,numel(BC_idcs));
    for kk=1:numel(BC_idcs)
        ii = BC_idcs(kk);
        terms = cell(1,0);
        forcing{kk} = cell(1,0);
        for jj=1:numel(PDE.BC{ii}.term)
            term = PDE.BC{ii}.term{jj};
            obj_type = types{find(isfield(term,types),1)};
            obj_idx = term.(obj_type);
            comp_vars = dep_vars.(obj_type)(obj_idx,:);
            pos = cumsum(comp_vars);
            if strcmp(obj_type,'x') && comp_vars(i)
                % State along s_i: keep only the derivative, position, and
                % integral along s_i.
                new_term = struct();
                new_term.x = x_new_idx(obj_idx);
                new_term.D = term.D(pos(i));
                new_term.loc = term.loc(pos(i));
                new_term.I = term.I(pos(i));
                new_term.C = term.C;
                new_term.delay = 0;
                terms{end+1} = new_term;                        %#ok<AGROW>
            else
                % Any other term goes to the right-hand side
                term.C = -term.C;
                forcing{kk}{end+1} = term;
            end
        end
        PDE_struct.BC{kk} = struct('size',PDE.BC{ii}.size,'term',{terms});
    end
    PDEs{i} = initialize(pde_struct(PDE_struct),true);
    info.BC_idcs{i} = BC_idcs;      
    info.x_idcs{i} = x_idcs;        
    info.forcing{i} = forcing;
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

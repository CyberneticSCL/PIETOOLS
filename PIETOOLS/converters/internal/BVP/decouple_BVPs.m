function [PDEs,info] = decouple_BVPs(PDE)
% [PDEs,info] = DECOUPLE_BVPS(PDE) splits the boundary conditions (BCs)
% of a PDE into groups that share no state or input, and returns a
% boundary value problem (BVP) for each group, so that they can be
% solved separately. Two BCs are in the same group if they involve the
% same state or input, or are linked by a chain of BCs that do. The
% dynamics are ignored.
%
% INPUT
% - PDE:    initialized 'pde_struct' object in the terms format;
%
% OUTPUT
% - PDEs:   1 x nG cell, with PDEs{g} the 'pde_struct' of the BVP of
%           group g: the states and inputs of PDE in its BCs (with the
%           same sizes, domains and orders, and no equations), and these
%           BCs (in their original order, with the indices of the states
%           and inputs of the group), on the variables that these
%           components depend on, with the delays of PDE. The groups with
%           BCs come first, ordered by their first BC, followed by the
%           states in no BC, each a BVP without BCs (unconstrained).
%           Inputs in no BC are in no group;
% - info:   struct with fields BC_idcs, x_idcs, w_idcs, u_idcs, 1 x nG
%           cells, with e.g. PDEs{g}.x{k} = PDE.x{x_idcs{g}(k)};
%
% NOTES
% - Components are grouped as a whole, so that e.g. a state of size 2 is
%   in a single group, even if one BC only involves its first element.
% - The BCs of a group need its states to depend on the same variables.
%   A group that would have states on different variables because their
%   BCs share an input (e.g. x1(s1)=w and x2(s2)=w) is split into the
%   groups that share a state, each including the input.
% - If no component of a group depends on a spatial variable (e.g. a BC
%   of an ODE state only), the BVP keeps all variables of PDE, as
%   initialize requires BCs to be constant in at least one variable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - decouple_BVPs
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
% DJ (CC), 11/13/2026: Require an initialized 'pde_struct' (no check);

N = size(PDE.vars,1);
nBC = numel(PDE.BC);
nx = numel(PDE.x);      nw = numel(PDE.w);      nu = numel(PDE.u);


% % % Link each BC to the components in its terms, as edges in a graph with
% % % a node for each BC, state, input w, and input u (in that order)
obj_off = struct('x',nBC,'w',nBC+nx,'u',nBC+nx+nw);
types = {'x','w','u'};
edge_BC = zeros(1,0);       edge_comp = zeros(1,0);
for ii=1:nBC
    for jj=1:numel(PDE.BC{ii}.term)
        term = PDE.BC{ii}.term{jj};
        obj_type = types{find(isfield(term,types),1)};        % x, w or u
        edge_BC(end+1) = ii;                            %#ok<AGROW>
        edge_comp(end+1) = obj_off.(obj_type) + term.(obj_type);          %#ok<AGROW>
    end
end
n_nodes = nBC+nx+nw+nu;
comp_labels = conncomp(graph(edge_BC,edge_comp,[],n_nodes));

% Inputs that are not in any BC are not part of any BVP
keep = true(1,n_nodes);
keep(nBC+nx+1:end) = ismember(nBC+nx+1:n_nodes,edge_comp);

% % % Order the groups by their first BC, followed by the states without
% % % BCs by index, and build the BVP of each group
labels = unique(comp_labels(keep));
nodes = arrayfun(@(b) find(comp_labels==b & keep),labels,'UniformOutput',false);
% A group with states on different variables, linked through inputs shared
% by BCs of these states, is split into the groups that share a state, as
% the BCs along the directions of one group need states on the same
% variables. The inputs are then part of each group whose BCs involve them
isvar_x = logical(PDE.x_tab(:,3:2+N));
is_state_edge = edge_comp<=nBC+nx;
comp_state = conncomp(graph(edge_BC(is_state_edge),edge_comp(is_state_edge),[],n_nodes));
nodes_split = cell(1,0);
for g=1:numel(nodes)
    BC_g = nodes{g}(nodes{g}<=nBC);
    x_g = nodes{g}(nodes{g}>nBC & nodes{g}<=nBC+nx) - nBC;
    if size(unique(isvar_x(x_g,:),'rows'),1)<=1
        nodes_split{end+1} = nodes{g};                  %#ok<AGROW>
        continue
    end
    for b=unique(comp_state(BC_g))
        BC_b = BC_g(comp_state(BC_g)==b);
        nodes_split{end+1} = unique([BC_b,edge_comp(ismember(edge_BC,BC_b))]);   %#ok<AGROW>
    end
end
nodes = nodes_split;
nG = numel(nodes);
% Sort key of each group: its first BC, then its first node
sort_key = cell2mat(cellfun(@(n) [min([n(n<=nBC),inf]),min(n)],nodes(:),'UniformOutput',false));
[~,order] = sortrows(sort_key);
nodes = nodes(order);
info = struct('BC_idcs',{cell(1,nG)},'x_idcs',{cell(1,nG)},'w_idcs',{cell(1,nG)},'u_idcs',{cell(1,nG)});
PDEs = cell(1,nG);
for g=1:nG
    info.BC_idcs{g} = nodes{g}(nodes{g}<=nBC);
    info.x_idcs{g} = nodes{g}(nodes{g}>nBC & nodes{g}<=nBC+nx) - nBC;
    info.w_idcs{g} = nodes{g}(nodes{g}>nBC+nx & nodes{g}<=nBC+nx+nw) - (nBC+nx);
    info.u_idcs{g} = nodes{g}(nodes{g}>nBC+nx+nw) - (nBC+nx+nw);

    % The BVP of the group: the BCs and components in the group,
    % renumbered, with no equations, and on the variables that these
    % components depend on (all variables if none does)
    idcs = struct('x',info.x_idcs{g},'w',info.w_idcs{g},'u',info.u_idcs{g});
    new_idx = struct('x',zeros(nx,1),'w',zeros(nw,1),'u',zeros(nu,1));
    used = false(1,N);
    for obj_type = 'xwu'
        new_idx.(obj_type)(idcs.(obj_type)) = 1:numel(idcs.(obj_type));
        used = used | any(PDE.([obj_type,'_tab'])(idcs.(obj_type),3:2+N),1);
    end
    if ~any(used)
        used = true(1,N);
    end
    PDE_struct = struct();
    PDE_struct.vars = PDE.vars(used,:);      PDE_struct.dom = PDE.dom(used,:);
    if PDE.has_delay
        PDE_struct.tau = PDE.tau;
    end
    PDE_struct.x = cell(numel(idcs.x),1);
    for k=1:numel(idcs.x)
        x_k = PDE.x{idcs.x(k)};
        PDE_struct.x{k} = struct('size',x_k.size,'vars',x_k.vars,'dom',x_k.dom,'diff',x_k.diff);
    end
    for obj_type = 'wu'
        PDE_struct.(obj_type) = cell(numel(idcs.(obj_type)),1);
        for k=1:numel(idcs.(obj_type))
            obj_k = PDE.(obj_type){idcs.(obj_type)(k)};
            PDE_struct.(obj_type){k} = struct('size',obj_k.size,'vars',obj_k.vars,'dom',obj_k.dom);
        end
    end
    PDE_struct.BC = cell(numel(info.BC_idcs{g}),1);
    for kk=1:numel(info.BC_idcs{g})
        terms = PDE.BC{info.BC_idcs{g}(kk)}.term;
        for jj=1:numel(terms)
            obj_type = types{find(isfield(terms{jj},types),1)};
            terms{jj}.(obj_type) = new_idx.(obj_type)(terms{jj}.(obj_type));
        end
        PDE_struct.BC{kk} = struct('size',PDE.BC{info.BC_idcs{g}(kk)}.size,'term',{terms});
    end
    PDEs{g} = initialize(pde_struct(PDE_struct),true);
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

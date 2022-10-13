%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function [PDE,comp_order] = reorder_comps(PDE,obj)
% This function reorders the rows of tables PDE.x_tab, PDE.w_tab,
% PDE.u_tab, PDE.y_tab and PDE.z_tab, in preparation of converting the PDE
% to a PIE.
%
% INPUTS:
% - PDE:    A struct or pde_struct type object, defining a PDE in the terms
%           format (see also the "@pde_struct/initialize" function).
% - obj:    'x', 'y', 'z', 'u', 'w', indicating for what object to reorder
%           the components.
%
% OUTPUTS:
% - PDE:    A structure of the same type as the input, with the table
%           PDE.([obj'_tab']) and the fields PDE.(obj) adjusted to match
%           the new order of the components of obj. If obj='x', obj='u' or
%           obj='w', the terms in each of the PDE equations are also
%           adjusted to make sure each term still refers to the appropriate
%           component of PDE.x, PDE.u, or PDE.w.
% - comp_order: An nx1 array indicating the new order of the components of
%               obj, so that PDE_new.(obj){j} = PDE_old.(obj){comp_order(j)}.
%
%
% NOTES:
% The components are seperated first based on which variables they depend
% on, combining e.g. the state components x_{j} into a full state
% x = [x0; x1(s1); x2(s2); x3(s1,s2)];
% For each function space, the components x_{j} within x1, x2 and x3 are
% then ordered based on their order of differentiability wrt each of the
% variables they depend on. For example, letting x_{i,j}(s1,s2) within x2
% be differentiable up to degree i wrt s1 and degree j wrt s2, they will be
% ordered as
% x3 = [x_{0,0}; x_{1,0}; ...; x_{N1,0}; x_{0,1}; x_{1,1}; ...; x_{N1,N2}];
% Finally, if different state components are differentiable up to the same
% degree, they will be ordered based on their original order in the PDE
% structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/09/2022
%

% Extract PDE information.
nvars = PDE.dim;

% Re-order the rows of x_tab such that:
% - the index of the state components increases fastest;
% - the order of differentiability wrt each of the nvars variables
%   increase next, with the order of differentiability wrt the first
%   variable increasing fastest and the order of differentiability wrt the
%   last variable increases slowest;
% - the dependence of the componetn on each of the nvars variables
%   increase last, with the dependence on the first variable increasing
%   fastest and the dependence on the last variable increasing slowest.
% For inputs, there is no order of differentiability.
if strcmp(obj,'x')
    % For the state components, we order based on order of
    % differentiability as well.
    comp_tab = PDE.x_tab;
    dep_tab = comp_tab(:,3:2+nvars);
    diff_tab = comp_tab(:,3+nvars:2*nvars+2);
    comp_tab_alt = [dep_tab(:,end:-1:1), diff_tab(:,end:-1:1), comp_tab(:,1), comp_tab(:,2)];
    comp_tab_alt = sortrows_integerTable(comp_tab_alt);
    obj_tab_new = [comp_tab_alt(:,end-1:end), comp_tab_alt(:,nvars:-1:1), comp_tab_alt(:,2*nvars:-1:nvars+1)];
else
    % For inputs, we have no order of differentiability, but we do need to
    % adjust the terms involving inputs to refer to the appropriate input.
    comp_tab = PDE.([obj,'_tab']);
    dep_tab = comp_tab(:,3:2+nvars);
    comp_tab_alt = [dep_tab(:,end:-1:1), comp_tab(:,1), comp_tab(:,2)];
    comp_tab_alt = sortrows_integerTable(comp_tab_alt);
    obj_tab_new = [comp_tab_alt(:,end-1:end), comp_tab_alt(:,nvars:-1:1)];
end

if strcmp(obj,'x') || strcmp(obj,'u') || strcmp(obj,'w')
    % Adjust the terms in each equation to match the new order of the
    % components.
    [~,new_order] = sort(obj_tab_new(:,1));    % comp_tab_new(new_order,:) = comp_tab_old
    for ii=1:numel(PDE.x)
        for jj=1:numel(PDE.x{ii}.term)
            if isfield(PDE.x{ii}.term{jj},obj)
                PDE.x{ii}.term{jj}.(obj) = new_order(PDE.x{ii}.term{jj}.(obj));
            end
        end
    end
    for ii=1:numel(PDE.y)
        for jj=1:numel(PDE.y{ii}.term)
            if isfield(PDE.y{ii}.term{jj},obj)
                PDE.y{ii}.term{jj}.(obj) = new_order(PDE.y{ii}.term{jj}.(obj));
            end
        end
    end
    for ii=1:numel(PDE.z)
        for jj=1:numel(PDE.z{ii}.term)
            if isfield(PDE.z{ii}.term{jj},obj)
                PDE.z{ii}.term{jj}.(obj) = new_order(PDE.z{ii}.term{jj}.(obj));
            end
        end
    end
    for ii=1:numel(PDE.BC)
        for jj=1:numel(PDE.BC{ii}.term)
            if isfield(PDE.BC{ii}.term{jj},obj)
                PDE.BC{ii}.term{jj}.(obj) = new_order(PDE.BC{ii}.term{jj}.(obj));
            end
        end
    end
end

comp_order = obj_tab_new(:,1);
%if strcmp(obj,'x') || strcmp(obj,'y') || strcmp(obj,'z')
% Re-arrange the equations to match the new order
PDE.(obj) = PDE.(obj)(comp_order);
%end

obj_tab_new(:,1) = 1:size(obj_tab_new,1);
PDE.([obj,'_tab']) = obj_tab_new;

end
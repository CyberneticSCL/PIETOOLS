function [PDE_1_out,PDE_2_out,new_obnums_1,new_obnums_2] = pde_common_basis(PDE_1,PDE_2)
% [PDE_1_OUT,PDE_2_OUT] = PDE_COMMON_BASIS(PDE_1,PDE_2) takes two PDE
% structures, and expresses them in terms of the same state components,
% inputs and outputs.
%
% INPUT
% - PDE_1, PDE_2:   pde_struct objects representing two PDEs or loose terms
%                   to build PDEs, expressed in terms of potentially
%                   different states, input, and output variables.
%
% OUTPUT
% - PDE_1_out, PDE_2_out:
%                   pde_struct objects representing the same PDEs or loose
%                   terms as the inputs, but now expressed in terms of a
%                   common set of state, input, and output variables. In
%                   particular, the fields PDE.x, PDE.y, PDE.z, PDE.u and
%                   PDE.w will have the same number of elements and define
%                   the same objects, though equations specified for these
%                   objects (the fields PDE.x{ii}.term) may be diffferent.
% - new_obnums_1, new_obnums_2:
%                   1x5 cell arrays specifying for each of the objects
%                   {'x','y','z','u','w'} how the new object indices
%                   correspond to the old ones. That is, e.g.
%                   PDE_1_out.x_tab(new_obnums_1{1},:) = PDE_1.x_tab;
%                   PDE_1_out.y_tab(new_obnums_1{2},:) = PDE_1.y_tab;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 06/23/2024

% % Since state variables and inputs may appear both in PDE_1 and PDE_2, we
% % will have to be careful to combine them in the structure.
PDE_1_out = PDE_1;      PDE_2_out = PDE_2;
objs = {'x';'y';'z';'u';'w'};
n_objs = length(objs);
new_obnums_1 = cell(1,n_objs);   new_obnums_2 = cell(1,n_objs);
for kk=1:n_objs
    % % Combine the information tables of the object from PDE_1 and PDE_2.
    obj = objs{kk};
    tab_1 = PDE_1.([obj,'_tab']);     tab_2 = PDE_2.([obj,'_tab']);
    tab_full = [tab_1;tab_2];
    [new_IDs,idcs1,idcs2] = unique(tab_full(:,1),'stable');   % new_IDs = tab_full(idcs1,1);      tab_full(:,1) = new_IDs(idcs2);
    % Output PDE structures both have the combined table
    PDE_1_out.([obj,'_tab']) = tab_full(idcs1,:);
    PDE_2_out.([obj,'_tab']) = tab_full(idcs1,:);
    % Keep track of combined number of components.
    n_obj_kk_1 = size(tab_1,1);      n_obj = length(new_IDs);
    
    % % Copy components from both PDEs into the output PDE_1_out.
    is_PDE1_comp = idcs1<=n_obj_kk_1;
    PDE_1_out.(obj) = cell(n_obj,1);
    PDE_1_out.(obj)(is_PDE1_comp) = PDE_1.(obj)(idcs1(is_PDE1_comp));
    PDE_1_out.(obj)(~is_PDE1_comp) = PDE_2.(obj)(idcs1(~is_PDE1_comp)-n_obj_kk_1);
    % Remove any equations for the state or output object that might have
    % been copied from PDE_2.
    if isfield(PDE_1_out.(obj)(~is_PDE1_comp),'term')
        PDE_1_out.(obj)(~is_PDE1_comp).term = {};
    end

    % % Copy components from both PDEs into the output PDE_2_out.
    [is_PDE2_comp,PDE2_comp_idcs] = ismember(new_IDs,tab_2(:,1));   % new_IDs(is_PDE2_comp) = tab_2(PDE2_comp_idcs(is_PDE2_comp),1)
    PDE_2_out.(obj) = cell(n_obj,1);
    PDE_2_out.(obj)(is_PDE2_comp) = PDE_2.(obj)(PDE2_comp_idcs(is_PDE2_comp));
    PDE_2_out.(obj)(~is_PDE2_comp) = PDE_1.(obj)(idcs1(~is_PDE2_comp));

    % % Check that shared components in the two PDEs are indeed the same,
    % % and remove equations from the output PDEs that don't belong to them.
    for ii=1:n_obj
        if ismember(new_IDs(ii),tab_1(:,1)) && ismember(new_IDs(ii),tab_2(:,1))
            % Make sure the objects in both equations indeed match
            if PDE_1_out.(obj){ii}.size~=PDE_2_out.(obj){ii}.size ||...
                 ~isequal(sort(PDE_1_out.(obj){ii}.vars.varname),sort(PDE_2_out.(obj){ii}.vars.varname))
                error("Components that appear in both PDEs should have the same size and spatial variables.")
            end
            continue
        end
        if ~ismember(new_IDs(ii),tab_1(:,1)) && isfield(PDE_1_out.(obj){ii},'term')
            % Remove equations from PDE_1 that belong to PDE_1.
            PDE_1_out.(obj){ii}.term = {};
        end
        if ~ismember(new_IDs(ii),tab_2(:,1)) && isfield(PDE_2_out.(obj){ii},'term')
            % Remove equations from PDE_2 that belong to PDE_1.
            PDE_2_out.(obj){ii}.term = {};
        end
    end

    % Keep track of how old component numbers relate to new ones
    n_obj_kk_1 = size(tab_1,1);
    new_obnums_1{kk} = idcs2(1:n_obj_kk_1);
    new_obnums_2{kk} = idcs2(n_obj_kk_1+1:end);   % new_idx = new_obnums_2(old_idx);
end

% % Update the indices of the object in each term, in each equation, in
% % each PDE, to match the combined list of indices.
PDE_1_out = update_obnums(PDE_1_out,new_obnums_1);
PDE_2_out = update_obnums(PDE_2_out,new_obnums_2);

end

function PDE = update_obnums(PDE,new_obnums)
% Update the indices associated to the objects 'x', 'y', 'z', 'u', and 'w'
% in the PDE structure, based on the new list of indices "new_obnums".
%
% INPUT
% - PDE:            pde_struct object
% - new_obnums:     1x5 cell providing updated indices for the objects
%                   {'x','y','z','u','w'}, in that order! Each element
%                   new_obnums{j} should be a nx1 array where n is the
%                   number of element of the field PDE.(obj), with obj
%                   being either 'x', 'y', etc.. Updated indices should be
%                   such that   new_idx = new_obnums{j}(old_idx);
% OUTPUT
% - PDE:            pde_struct object representing the same PDE, but now
%                   expressed in terms of the updated indices associated to
%                   each object.


% % Loop over all types of equations that may appear
eq_types = {'free';'x';'y';'z';'BC'};
for kk=1:numel(eq_types)
    eq_type_kk = eq_types{kk};
    % % Loop over all equations in each object
    for ii=1:numel(PDE.(eq_type_kk))
        if ~isfield(PDE.(eq_type_kk){ii},'term') || isempty(PDE.(eq_type_kk){ii}.term)
            % No adjusting of terms if there are no terms...
            continue
        end
        % Deal with the case that the first term corresponds to the temporal 
        % derivative of a state, or to an output (only for "free" type
        % equation)
        if strcmp(eq_type_kk,'free')
            [is_LHS_kk,obj_kk] = is_LHS_term(PDE.free{ii}.term{1});
            if is_LHS_kk
                % Update the index, depending on what type of object is involved.
                if strcmp(obj_kk,'x')
                    PDE.free{ii}.term{1}.x = new_obnums{1}(PDE.free{ii}.term{1}.x);
                elseif strcmp(obj_kk,'y')
                    PDE.free{ii}.term{1}.y = new_obnums{2}(PDE.free{ii}.term{1}.y);
                elseif strcmp(obj_kk,'z')
                    PDE.free{ii}.term{1}.z = new_obnums{3}(PDE.free{ii}.term{1}.z);
                % elseif strcmp(obj_kk,'u')
                %     PDE.free{ii}.term{1}.u = new_obnums{4}(PDE.free{ii}.term{1}.u);
                % elseif strcmp(obj_kk,'w')
                %     PDE.free{ii}.term{1}.w = new_obnums{5}(PDE.free{ii}.term{1}.w);
                else
                    error(["Field 'free{",num2str(ii),"}.term{1}' of one of the input PDEs refers to something other than 'x', 'y', or 'z'; how did this happen?"])
                end
                start_idx = 2;      % Start with second term
            else
                start_idx = 1;
            end
        else
            start_idx = 1;
        end
    
        % % For the remaining terms, update the index for either the state
        % % or input component.
        for jj=start_idx:numel(PDE.(eq_type_kk){ii}.term)
            if isfield(PDE.(eq_type_kk){ii}.term{jj},'x')
                new_obnums_jj = new_obnums{1};
                PDE.(eq_type_kk){ii}.term{jj}.x = new_obnums_jj(PDE.(eq_type_kk){ii}.term{jj}.x);
            elseif isfield(PDE.(eq_type_kk){ii}.term{jj},'u')
                new_obnums_jj = new_obnums{4};
                PDE.(eq_type_kk){ii}.term{jj}.u = new_obnums_jj(PDE.(eq_type_kk){ii}.term{jj}.u);
            else 
                new_obnums_jj = new_obnums{5};
                PDE.(eq_type_kk){ii}.term{jj}.w = new_obnums_jj(PDE.(eq_type_kk){ii}.term{jj}.w);
            end
        end
    end
end

end
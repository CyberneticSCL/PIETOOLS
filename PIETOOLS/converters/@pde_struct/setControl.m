function [PDE_out] = setControl(PDE_in,control_IDs)
% PDE_OUT = SETCONTROL(PDE_IN,CONTROL_IDS) changes the exogenous inputs
% with IDs "control_IDs" in the structure "PDE_in" to controlled inputs.
%
% INPUT
% - PDE_in:         'pde_struct' object specifying a PDE with possible
%                   exogenous inputs.
% - control_IDs:    nx1 integer array specifying which exogenous input
%                   signals should be converted to controlled inputs. Each
%                   element should correspond to "PDE_in.w{i}.ID" for some
%                   index i.
%
% OUTPUT
% - PDE_out:        'pde_struct' object specifying the same PDE as
%                   "PDE_in", but now with the desired exogenous inputs
%                   changed to controlled inputs.
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
% Initial coding DJ - 06/30/2024

% % Check that the input indices are properly specified.
if ~isnumeric(control_IDs) || any(control_IDs~=round(control_IDs),'all') || any(control_IDs<1,'all')
    error("The input names to convert to controlled inputs should be specified as nx1 array of positive integers.")
end
% Get rid of possible duplicates.
control_IDs = unique(control_IDs(:));

% % For each of the specified input IDs, replace the associated exogenous
% % input with a controlled input.
PDE_out = PDE_in;
for ll=1:length(control_IDs)
    ID_ll = control_IDs(ll);
    % Determine which row in PDE_out.w corresponds to the ID.
    is_wrow_ll = PDE_out.w_tab(:,1)==ID_ll;
    if ~any(is_wrow_ll)
        continue
    end
    % Move the object of interest to the controlled inputs.
    PDE_out.u = [PDE_out.u; PDE_out.w(is_wrow_ll)];
    PDE_out.w = PDE_out.w(~is_wrow_ll);
    % Move the associated table elements to the controlled inputs table.
    PDE_out.u_tab = [PDE_out.u_tab; PDE_out.w_tab(is_wrow_ll,:)];
    PDE_out.w_tab = PDE_out.w_tab(~is_wrow_ll,:);
    

    % % Loop over all equations, adjusting the terms to match the fact that
    % % the exogenous has been changed to a controlled input.
    % Establish the old and new index associated to the object.
    old_idx = find(PDE_in.w_tab(:,1)==ID_ll);
    new_idx = numel(PDE_out.u);
    % Loop over all types of equations.
    eq_types = {'x';'y';'z';'BC'};
    for kk=1:numel(eq_types)
        % Loop over all equations of the considered type.
        for ii=1:numel(PDE_out.(eq_types{kk}))
            if ~isfield(PDE_out.(eq_types{kk}){ii},'term') || isempty(PDE_out.(eq_types{kk}){ii}.term)
                continue
            end
            % Loop over all terms in the equation.
            for jj=1:numel(PDE_out.(eq_types{kk}){ii}.term)
                if isfield(PDE_out.(eq_types{kk}){ii}.term{jj},'w') &&...
                        PDE_out.(eq_types{kk}){ii}.term{jj}.w==old_idx
                    % Replace the exogenous input with a controlled one.
                    PDE_out.(eq_types{kk}){ii}.term{jj} = rmfield(PDE_out.(eq_types{kk}){ii}.term{jj},'w');
                    PDE_out.(eq_types{kk}){ii}.term{jj}.u = new_idx;
                end
            end
        end
    end
end
% This change should not require re-initialization.
if PDE_in.is_initialized
    PDE_out.is_initialized = true;
end

fprintf('%d inputs were designated as controlled inputs\n', numel(PDE_out.u)-numel(PDE_in.u));

end
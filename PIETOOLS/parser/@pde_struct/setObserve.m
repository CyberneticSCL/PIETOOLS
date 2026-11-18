function [PDE_out] = setObserve(PDE_in,observe_IDs)
% PDE_OUT = SETOBSERVE(PDE_IN,OBSERVE_IDS) changes the regulated outputs
% with IDs "observe_IDs" in the structure "PDE_in" to observed outputs.
%
% INPUT
% - PDE_in:         'pde_struct' object specifying a PDE with possible
%                   regulated outputs.
% - observe_IDs:    nx1 integer array specifying which regulated output
%                   signals should be converted to observed outputs. Each
%                   element should correspond to "PDE_in.z{i}.ID" for some
%                   index i.
%
% OUTPUT
% - PDE_out:        'pde_struct' object specifying the same PDE as
%                   "PDE_in", but now with the desired regulated outputs
%                   changed to observed outputs.
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
% DJ, 12/28/2024: Allow output to be declared as 'pde_var' object.

% % Check that the output indices are properly specified.
if isa(observe_IDs,'pde_struct')
    [isvar,obj] = is_pde_var(observe_IDs);
    if isvar
        if ~strcmp(obj,'z')
            error("Only conversion of regulated outputs to observed outputs is supported.")
        else
            observe_IDs = observe_IDs.z{1}.ID;
        end
    end
if ~isnumeric(observe_IDs) || any(observe_IDs~=round(observe_IDs),'all') || any(observe_IDs<1,'all')
    error("The output names to convert to observed outputs should be specified as nx1 array of positive integers.")
end
% Get rid of possible duplicates.
observe_IDs = unique(observe_IDs(:));

% % For each of the specified output IDs, replace the associated regulated
% % output with an observed output.
PDE_out = PDE_in;
for ii=1:length(observe_IDs)
    ID_ii = observe_IDs(ii);
    % Determine which row in PDE_out.z corresponds to the ID.
    is_zrow_ii = ismember(PDE_out.z_tab(:,1),ID_ii);
    if ~any(is_zrow_ii)
        continue
    end
    % Move the equation of interest to the observed outputs.
    PDE_out.y = [PDE_out.y; PDE_out.z(is_zrow_ii)];
    PDE_out.z = PDE_out.z(~is_zrow_ii);
    % Move the associated table elements to the observed outputs table.
    PDE_out.y_tab = [PDE_out.y_tab; PDE_out.z_tab(is_zrow_ii,:)];
    PDE_out.z_tab = PDE_out.z_tab(~is_zrow_ii,:);
end
% This change should not require re-initialization.
if PDE_in.is_initialized
    PDE_out.is_initialized = true;
end

fprintf('%d outputs were designated as observed outputs\n', numel(PDE_out.y)-numel(PDE_in.y));

end
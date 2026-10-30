function Pop = pie_rename_dummy(Pop)
% POP = PIE_RENAME_DUMMY(POP) renames the dummy variables of an 'opvar' or
% 'opvar2d' object to the convention of the PIE structure, adding '_dum' to
% the name of the associated spatial variable, e.g. theta --> s_dum if the
% spatial variable of the operator is s.
%
% INPUT
% - Pop:    'opvar' or 'opvar2d' object;
%
% OUTPUT
% - Pop:    'opvar' or 'opvar2d' object representing the same operator,
%           with the dummy variables Pop.var2 renamed, both in the property
%           var2 and in all the polynomial kernels of the operator.
%
% NOTES
% The conversion of operators to 'copvar' objects requires the dummy
% variables to follow this convention, as the dummy variables are not
% recorded in the 'sopvar' blocks.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie_rename_dummy
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

old_vars = Pop.var2;
new_vars = polynomial(strcat(reshape(Pop.var1.varname,[],1),'_dum'));
if all(all(isequal(old_vars,new_vars)))
    % The dummy variables already follow the convention.
    return
end

% Substitute the new dummy variables in each of the kernels.
if isa(Pop,'opvar')
    fnames = {'P','Q1','Q2'};
    for ii=1:numel(fnames)
        Pop.(fnames{ii}) = rename_in(Pop.(fnames{ii}),old_vars,new_vars);
    end
    Pop.R = rename_in(Pop.R,old_vars,new_vars);
else
    fnames = {'R00','R0x','R0y','R02';
              'Rx0','Rxx','Rxy','Rx2';
              'Ry0','Ryx','Ryy','Ry2';
              'R20','R2x','R2y','R22'};
    for ii=1:numel(fnames)
        Pop.(fnames{ii}) = rename_in(Pop.(fnames{ii}),old_vars,new_vars);
    end
end
Pop.var2 = new_vars;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function val = rename_in(val,old_vars,new_vars)
% VAL = RENAME_IN(VAL,OLD_VARS,NEW_VARS) substitutes the variables OLD_VARS
% by NEW_VARS in the polynomial VAL, or in each of the elements of a cell
% or struct VAL.
if iscell(val)
    for ii=1:numel(val)
        val{ii} = rename_in(val{ii},old_vars,new_vars);
    end
elseif isstruct(val)
    fnames = fieldnames(val);
    for ii=1:numel(fnames)
        val.(fnames{ii}) = rename_in(val.(fnames{ii}),old_vars,new_vars);
    end
elseif isa(val,'polynomial') || isa(val,'dpvar')
    val = subs(val,old_vars,new_vars);
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
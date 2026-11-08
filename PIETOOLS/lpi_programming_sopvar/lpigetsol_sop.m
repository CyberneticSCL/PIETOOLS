function dvar_sol = lpigetsol_sop(prog,dvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DVAR_SOL = LPIGETSOL_SOP(PROG,DVAR) returns DVAR with each of its
% decision variables replaced by its value in the solved LPI program PROG:
% an 'sdopvar' or 'cdopvar' DVAR by GETSOL_LPIVAR_SOP; a cell that holds,
% at any depth, an object of class 'sopvar', 'sdopvar', 'copvar' or
% 'cdopvar', element by element; and a 'dpvar', char or other cell by
% SOSGETSOL. A 'sopvar', 'copvar', 'double' or 'polynomial' DVAR, or [], is
% returned unchanged.
%
% INPUT
% - prog:     'struct', an LPI program solved by LPISOLVE, with the fields
%             solinfo.info, decvartable and solinfo.RRx;
% - dvar:     'sdopvar' or 'cdopvar' object, a decision operator; a
%             'sopvar' or 'copvar' object; a cell of objects of the classes
%             listed here; or [], a 'double', 'polynomial', 'dpvar', char
%             or 'cellstr' object.
%
% OUTPUT
% - dvar_sol: DVAR at the solution: for an 'sdopvar' ('cdopvar') DVAR, a
%             'sopvar' ('copvar') object, or an 'sdopvar' ('cdopvar')
%             object on the names that prog.decvartable does not hold; for
%             a cell, a cell of the same size; for [], 'double' and
%             'polynomial', DVAR itself; for 'dpvar', char and a cell
%             without a 'sopvar', 'sdopvar', 'copvar' or 'cdopvar' element,
%             the value SOSGETSOL returns.
%
% NOTES
% A DVAR of another class, a PROG that is not a 'struct' and a PROG whose
% solinfo.info is empty are errors.
%
% See also SOSGETSOL, GETSOL_LPIVAR_SOP, SUBS_DVAR_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpigetsol_sop
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
% Initial coding MMP, 09/29/2026. Dispatcher over both operator families,
%                in a folder parallel to lpi_programming and suffixed _sop
%                so that it shadows nothing (every folder is on the path).
%                Mirrors lpigetsol (DJ, 10/19/2024), which it calls for
%                every legacy class.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was the class table with its comparison to
%                lpigetsol and the reason a cell is taken apart here.
% MMP, 11/04/2026: [], 'double' and 'polynomial' are returned unchanged and
%                a 'dpvar', char or cell goes to SOSGETSOL here (what
%                LPIGETSOL did for these classes, in the same order), so the
%                container path no longer calls lpi_programming.
% MMP, 11/08/2026: The branches of the opvar classes are gone, an 'opvar'
%                or 'opvar2d' with 'dpvar' fields (local SUBS_FIELDS) and the
%                legacy classes (LPIGETSOL): both are retired (cleaning/
%                private, 11/08/2026), and such a DVAR now raises
%                lpigetsol_sop:badClass. The local functions HASDPVAR and
%                SUBS_FIELDS, used only there, are deleted. Help text to
%                match.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): the 11/04/2026
%                entry without its superseded sentence, two body notes that
%                named the retired lpigetsol reworded; code unchanged.

% The program must be a 'struct' with a solver report (the checks LPIGETSOL
% made)
if ~isa(prog,'struct')
    error('LPI program structure should be specified as object of type ''struct''.')
end
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'info') || isempty(prog.solinfo.info)
    error('The specified LPI optimization program has not been solved; either ''lpisolve'' has not been run, or no solution was produced.')
end

if isa(dvar,'sdopvar') || isa(dvar,'cdopvar')
    % Decision operator of the container family: getsol_lpivar_sop
    dvar_sol = getsol_lpivar_sop(prog,dvar);
elseif isa(dvar,'sopvar') || isa(dvar,'copvar')
    dvar_sol = dvar;                            % fixed: nothing to substitute
elseif iscell(dvar) && ~iscellstr(dvar) && any(cellfun(@hasnew,dvar(:)))
    % A cell holding a container-family object: one element at a time,
    % since sosgetsol cannot read them
    dvar_sol = cell(size(dvar));                % nested cells recurse here too
    for k = 1:numel(dvar)
        dvar_sol{k} = lpigetsol_sop(prog,dvar{k});
    end
% elseif (isa(dvar,'opvar') || isa(dvar,'opvar2d')) && hasdpvar(dvar)       % MMP, 11/08/2026 (was)
%     % Legacy [dpvar, opvar] concatenation leaves dpvar-valued fields in an  % MMP, 11/08/2026 (was)
%     % opvar, which lpigetsol returns UNSUBSTITUTED; substitute each field.  % MMP, 11/08/2026 (was)
%     dvar_sol = subs_fields(prog,dvar);                                    % MMP, 11/08/2026 (was)
elseif isempty(dvar) || isa(dvar,'double') || isa(dvar,'polynomial')        % MMP, 11/04/2026
    dvar_sol = dvar;    % fixed: nothing to substitute, as in LPIGETSOL     % MMP, 11/04/2026
elseif isa(dvar,'dpvar') || ischar(dvar) || iscell(dvar)                    % MMP, 11/04/2026
    dvar_sol = sosgetsol(prog,dvar);    % as LPIGETSOL does for these       % MMP, 11/04/2026
% elseif isempty(dvar) || isa(dvar,'double') || isa(dvar,'polynomial') ...
%         || isa(dvar,'opvar') || isa(dvar,'opvar2d') || isa(dvar,'dpvar') ...
%         || ischar(dvar) || iscell(dvar) || isa(dvar,'dopvar') || isa(dvar,'dopvar2d')
%     dvar_sol = lpigetsol(prog,dvar);            % legacy classes, unchanged % MMP, 11/08/2026 (was)
else
    % Any other class
%   error('lpigetsol_sop:badClass',...
%       ['Variable of which to extract the solution should be specified as object of type '...
%        '''dpvar'', ''dopvar'', ''dopvar2d'', ''sdopvar'' or ''cdopvar'' (or a fixed '...
%        'opvar, opvar2d, sopvar, copvar, polynomial or double); got ''%s''.'],class(dvar)) % MMP, 11/08/2026 (was)
    error('lpigetsol_sop:badClass',...
        ['Variable of which to extract the solution should be specified as object of type '...
         '''dpvar'', ''sdopvar'' or ''cdopvar'' (or a fixed sopvar, copvar, '...
         'polynomial or double); got ''%s''.'],class(dvar))                 % MMP, 11/08/2026
end

end


function tf = hasnew(x)
% TF = HASNEW(X) returns true when X is, or is a cell holding at any depth,
% a 'sopvar', 'sdopvar', 'copvar' or 'cdopvar' object.
% INPUT  x:  any value.
% OUTPUT tf: logical.
tf = isa(x,'sopvar') || isa(x,'sdopvar') || isa(x,'copvar') || isa(x,'cdopvar') ...
    || (iscell(x) && any(cellfun(@hasnew,x(:))));
end


% The local functions HASDPVAR and SUBS_FIELDS (MMP, 09/29/2026), which     % MMP, 11/08/2026
% served only the 'opvar' branch above, are deleted.                        % MMP, 11/08/2026

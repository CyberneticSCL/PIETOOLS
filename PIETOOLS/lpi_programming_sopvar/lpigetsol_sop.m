function dvar_sol = lpigetsol_sop(prog,dvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DVAR_SOL = LPIGETSOL_SOP(PROG,DVAR) returns DVAR with each of its
% decision variables replaced by its value in the solved LPI program PROG:
% an 'sdopvar' or 'cdopvar' DVAR with every parameter evaluated at the
% solution; a cell that holds, at any depth, an object of class 'sopvar',
% 'sdopvar', 'copvar' or 'cdopvar', element by element; and a 'dpvar',
% char or other cell as a variable of the SOS program. A 'sopvar',
% 'copvar', 'double' or 'polynomial' DVAR, or [], is returned unchanged.
%
% INPUT
% - prog:     'struct', a solved LPI program, with the fields
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
%             its value in the SOS solution.
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
% Initial coding MMP, 09/29/2026

% The program must be a 'struct' with a solver report
if ~isa(prog,'struct')
    error('LPI program structure should be specified as object of type ''struct''.')
end
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'info') || isempty(prog.solinfo.info)
    error('The specified LPI optimization program has not been solved; either ''lpisolve'' has not been run, or no solution was produced.')
end

if isa(dvar,'sdopvar') || isa(dvar,'cdopvar')
    % Decision operator: its parameters at the solution
    dvar_sol = getsol_lpivar_sop(prog,dvar);
elseif isa(dvar,'sopvar') || isa(dvar,'copvar')
    dvar_sol = dvar;                            % fixed: nothing to substitute
elseif iscell(dvar) && ~iscellstr(dvar) && any(cellfun(@hasnew,dvar(:)))
    % A cell holding an operator object: one element at a time, since the
    % SOS solution routine does not read these classes
    dvar_sol = cell(size(dvar));                % nested cells recurse here too
    for k = 1:numel(dvar)
        dvar_sol{k} = lpigetsol_sop(prog,dvar{k});
    end
elseif isempty(dvar) || isa(dvar,'double') || isa(dvar,'polynomial')
    dvar_sol = dvar;    % fixed: nothing to substitute
elseif isa(dvar,'dpvar') || ischar(dvar) || iscell(dvar)
    dvar_sol = sosgetsol(prog,dvar);    % the value in the SOS solution
else
    % Any other class
    error('lpigetsol_sop:badClass',...
        ['Variable of which to extract the solution should be specified as object of type '...
         '''dpvar'', ''sdopvar'' or ''cdopvar'' (or a fixed sopvar, copvar, '...
         'polynomial or double); got ''%s''.'],class(dvar))
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

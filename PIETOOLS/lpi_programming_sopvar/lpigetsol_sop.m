function dvar_sol = lpigetsol_sop(prog,dvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DVAR_SOL = LPIGETSOL_SOP(PROG,DVAR) takes a solved LPI program structure
% 'prog' and returns the solved value of the variable or function 'dvar',
% for the LPI objects of both operator families. It is 'lpigetsol' with the
% 'sopvar' family added, so that one executive body serves both:
%
%   class of DVAR                               handled by
%   'sdopvar', 'cdopvar'                        getsol_lpivar_sop (new)
%   'sopvar', 'copvar'                          returned as is (no variables)
%   cell holding any of the four classes above  this function, per element
%     (at any depth)
%   'opvar', 'opvar2d' with dpvar-valued fields this function, per field
%   [], 'double', 'polynomial', 'opvar',        lpigetsol, unchanged
%   'opvar2d', 'dpvar', char, other cell,
%   'dopvar', 'dopvar2d'
%
% Any other class is an error, as in 'lpigetsol'.
%
% INPUT
% - prog:   'struct' specifying a solved LPI optimization program structure
%           (see also 'lpiprogram'). The format should match that used by
%           SOSTOOLS;
% - dvar:   object of one of the classes above: an array of decision
%           variables, a function of decision variables, or a PI operator
%           decision variable of the LPI program.
%
% OUTPUT
% - dvar_sol:   solved value of 'dvar' as per the solved program. For
%               'sdopvar'/'cdopvar' a 'sopvar'/'copvar' (see
%               'getsol_lpivar_sop' for names absent from the program).
%
% NOTES
% A cell holding a new-class object is taken apart here, one element per
% call, because 'lpigetsol' hands cells to 'sosgetsol', which cannot read
% the new classes. A cell with no such element is handed to 'lpigetsol'
% whole, so its result is exactly the legacy one.
%
% See also LPIGETSOL, GETSOL_LPIVAR_SOP, SUBS_DVAR_SOP.
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

% Same program checks as lpigetsol, made here too so that the new-class
% branches, which do not pass through lpigetsol, apply them.
if ~isa(prog,'struct')
    error('LPI program structure should be specified as object of type ''struct''.')
end
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'info') || isempty(prog.solinfo.info)
    error('The specified LPI optimization program has not been solved; either ''lpisolve'' has not been run, or no solution was produced.')
end

if isa(dvar,'sdopvar') || isa(dvar,'cdopvar')
    dvar_sol = getsol_lpivar_sop(prog,dvar);
elseif isa(dvar,'sopvar') || isa(dvar,'copvar')
    dvar_sol = dvar;                            % fixed: nothing to substitute
elseif iscell(dvar) && ~iscellstr(dvar) && any(cellfun(@hasnew,dvar(:)))
    dvar_sol = cell(size(dvar));                % nested cells recurse here too
    for k = 1:numel(dvar)
        dvar_sol{k} = lpigetsol_sop(prog,dvar{k});
    end
elseif (isa(dvar,'opvar') || isa(dvar,'opvar2d')) && hasdpvar(dvar)
    % Legacy [dpvar, opvar] concatenation leaves dpvar-valued fields in an
    % opvar, which lpigetsol returns UNSUBSTITUTED; substitute each field.
    dvar_sol = subs_fields(prog,dvar);
elseif isempty(dvar) || isa(dvar,'double') || isa(dvar,'polynomial') ...
        || isa(dvar,'opvar') || isa(dvar,'opvar2d') || isa(dvar,'dpvar') ...
        || ischar(dvar) || iscell(dvar) || isa(dvar,'dopvar') || isa(dvar,'dopvar2d')
    dvar_sol = lpigetsol(prog,dvar);            % legacy classes, unchanged
else
    error('lpigetsol_sop:badClass',...
        ['Variable of which to extract the solution should be specified as object of type '...
         '''dpvar'', ''dopvar'', ''dopvar2d'', ''sdopvar'' or ''cdopvar'' (or a fixed '...
         'opvar, opvar2d, sopvar, copvar, polynomial or double); got ''%s''.'],class(dvar))
end

end


function tf = hasnew(x)
% x is, or is a cell holding at any depth, a sopvar-family object.
tf = isa(x,'sopvar') || isa(x,'sdopvar') || isa(x,'copvar') || isa(x,'cdopvar') ...
    || (iscell(x) && any(cellfun(@hasnew,x(:))));
end


function tf = hasdpvar(X)
% X (an opvar/opvar2d, one of its fields, a struct or a cell) holds a dpvar.
if isa(X,'dpvar')
    tf = true;
elseif isa(X,'opvar') || isa(X,'opvar2d') || isstruct(X)
    f = fieldnames(X);      tf = false;         % fieldnames lists an object's public properties
    for k = 1:numel(f)
        if hasdpvar(X.(f{k})),  tf = true;  return,  end
    end
elseif iscell(X)
    tf = any(cellfun(@hasdpvar,X(:)));
else
    tf = false;
end
end


function X = subs_fields(prog,X)
% Replace every dpvar held in X by its solved value (lpigetsol, i.e.
% sosgetsol). Only fields that hold a dpvar are assigned, so no other
% property (dim, I, var1, ...) passes through a set method.
if isa(X,'dpvar')
    X = lpigetsol(prog,X);
elseif isa(X,'opvar') || isa(X,'opvar2d') || isstruct(X)
    f = fieldnames(X);
    for k = 1:numel(f)
        if hasdpvar(X.(f{k})),  X.(f{k}) = subs_fields(prog,X.(f{k}));  end
    end
elseif iscell(X)
    for k = 1:numel(X)
        if hasdpvar(X{k}),  X{k} = subs_fields(prog,X{k});  end
    end
end
end

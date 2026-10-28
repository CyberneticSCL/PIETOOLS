function dvar_sol = lpigetsol_sop(prog,dvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DVAR_SOL = LPIGETSOL_SOP(PROG,DVAR) returns DVAR with each of its
% decision variables replaced by its value in the solved LPI program PROG:
% an 'sdopvar' or 'cdopvar' DVAR by GETSOL_LPIVAR_SOP; a cell that holds,
% at any depth, an object of class 'sopvar', 'sdopvar', 'copvar' or
% 'cdopvar', element by element; an 'opvar' or 'opvar2d' DVAR with
% 'dpvar'-valued fields, field by field, by LPIGETSOL; and any other DVAR
% by LPIGETSOL. A 'sopvar' or 'copvar' DVAR is returned unchanged.
%
% INPUT
% - prog:     'struct', an LPI program solved by LPISOLVE, with the fields
%             solinfo.info, decvartable and solinfo.RRx;
% - dvar:     'sdopvar' or 'cdopvar' object, a decision operator; a
%             'sopvar' or 'copvar' object; a cell of objects of the classes
%             listed here; an 'opvar' or 'opvar2d' object, possibly with
%             'dpvar'-valued fields; or [], a 'double', 'polynomial',
%             'dpvar', char, 'cellstr', 'dopvar' or 'dopvar2d' object.
%
% OUTPUT
% - dvar_sol: DVAR at the solution: for an 'sdopvar' ('cdopvar') DVAR, a
%             'sopvar' ('copvar') object, or an 'sdopvar' ('cdopvar')
%             object on the names that prog.decvartable does not hold; for
%             a cell, a cell of the same size; for an 'opvar' or 'opvar2d'
%             DVAR, the object with each 'dpvar' it holds replaced by its
%             value; for [], 'double' and 'polynomial', DVAR itself; for
%             'dpvar', char and a cell without a 'sopvar', 'sdopvar',
%             'copvar' or 'cdopvar' element, the value SOSGETSOL returns;
%             for 'dopvar' and 'dopvar2d', the 'opvar' or 'opvar2d' object
%             GETSOL_LPIVAR returns.
%
% NOTES
% A DVAR of another class, a PROG that is not a 'struct' and a PROG whose
% solinfo.info is empty are errors.
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was the class table with its comparison to
%                lpigetsol and the reason a cell is taken apart here.

% Same program checks as lpigetsol, made here too so that the new-class
% branches, which do not pass through lpigetsol, apply them.
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
    % since sosgetsol, which lpigetsol hands a cell to, cannot read them
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
    % Any other class
    error('lpigetsol_sop:badClass',...
        ['Variable of which to extract the solution should be specified as object of type '...
         '''dpvar'', ''dopvar'', ''dopvar2d'', ''sdopvar'' or ''cdopvar'' (or a fixed '...
         'opvar, opvar2d, sopvar, copvar, polynomial or double); got ''%s''.'],class(dvar))
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


function tf = hasdpvar(X)
% TF = HASDPVAR(X) returns true when X is a 'dpvar' or holds one at any
% depth, as a property of an 'opvar' or 'opvar2d', a field of a struct or
% an element of a cell.
% INPUT  X:  any value.
% OUTPUT tf: logical.
if isa(X,'dpvar')
    tf = true;
elseif isa(X,'opvar') || isa(X,'opvar2d') || isstruct(X)
    % An object or a struct: any of its fields
    f = fieldnames(X);      tf = false;         % fieldnames lists an object's public properties
    for k = 1:numel(f)
        if hasdpvar(X.(f{k})),  tf = true;  return,  end
    end
elseif iscell(X)
    % A cell: any of its elements
    tf = any(cellfun(@hasdpvar,X(:)));
else
    tf = false;
end
end


function X = subs_fields(prog,X)
% X = SUBS_FIELDS(PROG,X) returns X with every 'dpvar' it holds replaced by
% its solved value (lpigetsol, i.e. sosgetsol).
% INPUT  prog: 'struct', the solved LPI program;
%        X:    'dpvar', 'opvar' or 'opvar2d' object, struct or cell.
% OUTPUT X:    X with each 'dpvar' replaced by its value.
% Only fields that hold a dpvar are assigned, so no other property (dim,
% I, var1, ...) passes through a set method.
if isa(X,'dpvar')
    X = lpigetsol(prog,X);
elseif isa(X,'opvar') || isa(X,'opvar2d') || isstruct(X)
    % An object or a struct: each field that holds a dpvar
    f = fieldnames(X);
    for k = 1:numel(f)
        if hasdpvar(X.(f{k})),  X.(f{k}) = subs_fields(prog,X.(f{k}));  end
    end
elseif iscell(X)
    % A cell: each element that holds a dpvar
    for k = 1:numel(X)
        if hasdpvar(X{k}),  X{k} = subs_fields(prog,X{k});  end
    end
end
end

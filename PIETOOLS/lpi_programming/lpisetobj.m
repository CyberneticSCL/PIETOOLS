function prog = lpisetobj(prog,obj)
% LPISETOBJ takes an LPI program structure 'prog' and declares some
% function 'obj' to be the objective function to minimize in the
% optimization program. The function is mostly a shell for 'sossetobj'.
%
% INPUT
% - prog:       'struct' specifying an LPI optimization program structure
%               (see also 'lpiprogram'). The format should match that used
%               by SOSTOOLS.
% - obj:        1x1 object of type 'dpvar', specifying the objective
%               function to minimize in the optimization program. Must
%               correspond to an affine combination of decision variables,
%               cannot depend on any independent variables.
%
% OUTPUT
% - prog:       'struct' specifying the same optimization program as passed
%               by the user, but now with the specified function 'obj' set
%               as objective function to minimize.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpisetobj
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 10/19/2024

% Check that the inputs make sense.
if ~isa(prog,'struct')
    error('LPI program structure should be specified as object of type ''struct''.')
end
if ~isa(obj,'dpvar')
    error('Objective function in the LPI optimization program should be specified as object of type ''dpvar''.')
end

% Check that the objective function does not involve independent variables.
obj = combine(obj);
if ~isempty(obj.varname)
    error('Objective function in the LPI optimization program cannot depend on independent variables.')
end

% Add new decision variables to the optimization program.
dvarname = obj.dvarname;
isnew_dvar = ~ismember(dvarname,prog.decvartable);
if any(isnew_dvar)
    prog = sosdecvar(prog,dpvar(dvarname(isnew_dvar)));
end

% Finally, declare the objective function.
prog = sossetobj(prog,obj);

end
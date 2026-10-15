function prog = piesos_setobj(prog,obj)
% PROG = PIESOS_SETOBJ(PROG,OBJ) declares the function OBJ to be the
% objective function to minimize in the PIESOS program PROG. Decision
% variables in OBJ that are not yet in PROG are added to the program, and
% the objective is set using `sossetobj`.
%
% INPUTS
% - prog:       'struct' specifying a PIESOS optimization program
%               (see also `piesos_program`). The format largely matches
%               that used by SOSTOOLS;
% - obj:        1 x 1 'dpvar' object, specifying the objective function
%               to minimize in the optimization program. Must be an
%               affine combination of decision variables, and cannot
%               depend on any independent variables;
%
% OUTPUTS
% - prog:       'struct' specifying the same optimization program as the
%               input, but with OBJ set as the objective function to
%               minimize;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_setobj
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
% DJ, 10/15/2026: Initial coding;

% Check that the inputs make sense.
if ~isa(prog,'struct')
    error("PIESOS program structure should be specified as object of type 'struct'.")
end
if ~isa(obj,'dpvar')
    error("Objective function in the PIESOS optimization program should be specified as object of type 'dpvar'.")
end

% Check that the objective function does not involve independent variables.
obj = combine(obj);
if ~isempty(obj.varname)
    error("Objective function in the PIESOS optimization program cannot depend on independent variables.")
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
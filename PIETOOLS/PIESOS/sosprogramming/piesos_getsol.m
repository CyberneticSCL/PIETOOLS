function dvar_sol = piesos_getsol(prog,dvar)
% DVAR_SOL = PIESOS_GETSOL(PROG,DVAR) returns the value of the operator or
% polynomial decision variable DVAR in the solved PIESOS program PROG.
%
% INPUTS
% - prog:       'struct' specifying a solved PIESOS optimization program
%               (see also `piesos_program`, `piesos_solve`). The format
%               matches that used by SOSTOOLS;
% - dvar:       object of type 'polyopvar', 'tensopmat', 'tensopvar',
%               'intop', 'dopvar', 'dopvar2d', or 'dpvar', or a cell of
%               such objects, specifying a decision variable in the
%               optimization program, corresponding to a function,
%               operator, or just a matrix. A 'char' name of a decision
%               variable is also accepted. Objects of type 'double',
%               'polynomial', 'opvar', 'opvar2d', or 'nopvar', which do
%               not contain decision variables, are returned as is;
%
% OUTPUTS
% - dvar_sol:   object of the same format as DVAR, with the decision
%               variables replaced by their values in the solved PIESOS
%               program;
%
% NOTES
% - Decision variables of type 'ndopvar' are not supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_getsol
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

% Check that the program structure is properly specified.
if ~isa(prog,'struct')
    error('PIESOS program structure should be specified as object of type ''struct''.')
end
if isempty(prog.solinfo.info)
    error('The specified PIESOS optimization program has not been solved; either ''piesos_solve'' has not been run, or no solution was produced.')
end

% Extract the solution.
if isa(dvar,'polyopvar')
    % The variable is a distributed polynomial
    % --> only the coefficients will be variable
    dvar_sol = dvar;
    dvar_sol.C = piesos_getsol(prog,dvar.C);
elseif isa(dvar,'tensopmat')
    dvar_sol = dvar;
    for j=1:numel(dvar.ops)
        dvar_sol.ops{j} = piesos_getsol(prog,dvar.ops{j});
    end
elseif isa(dvar,'tensopvar')
    % The variable defines coefficients acting on distributed monomials
    dvar_sol = dvar;
    for j=1:numel(dvar.ops)
        dvar_sol.ops{j} = piesos_getsol(prog,dvar.ops{j});
    end
elseif isa(dvar,'cell')
    % A cell of variables is specified
    dvar_sol = cell(size(dvar));
    for k=1:numel(dvar)
        dvar_sol{k} = piesos_getsol(prog,dvar{k});
    end
elseif isempty(dvar) || isa(dvar,'double') || isa(dvar,'polynomial') ||...
        isa(dvar,'opvar') || isa(dvar,'opvar2d') || isa(dvar,'nopvar')
    % The variable is not really a variable at all;
    dvar_sol = dvar;
elseif isa(dvar,'dpvar') || ischar(dvar)
    % The variable is some decision variable, or linear combination of
    % decision variables.
    dvar_sol = sosgetsol(prog,dvar);
elseif isa(dvar,'dopvar') || isa(dvar,'dopvar2d')
    % The variable is a PI operator decision variable.
    dvar_sol = getsol_lpivar(prog,dvar);
elseif isa(dvar,'intop')
    % The variable is a functional operator
    % Solve the kernels, and declare the functional in terms of these
    dvar_sol = intop(piesos_getsol(prog,kernels(dvar)),dvar.intmat,dvar.vars,dvar.dom);
elseif isa(dvar,'ndopvar')
    error("Decision variables of type 'ndopvar' are currently not supported.")
else
    error("Variable of which to extract the solution should be specified as object of type 'polyopvar', 'dpvar' or 'dopvar'.")
end

end
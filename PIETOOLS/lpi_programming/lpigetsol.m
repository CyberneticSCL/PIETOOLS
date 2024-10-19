function dvar_sol = lpigetsol(prog,dvar)
% DVAR_SOL = LPIGETSOL(PROG,DVAR) takes a solved LPI program structure
% 'prog' and returns the corresponding solved value of the (PI operator)
% variable or function 'dvar'.
%
% INPUT
% - prog:       'struct' specifying a solved LPI optimization program
%               structure (see also 'lpiprogram'). The format should match
%               that used by SOSTOOLS.
% - dvar:       Object of type 'dpvar', 'dopvar', or 'dopvar2d', specifying 
%               an array of decision variable, a function of decision 
%               variables, or a PI operator decision variable in the LPI
%               optimization program of which to return the solution.
%
% OUTPUT
% - dvar_sol:   Solved value of the input LPI decision variable as per the
%               solved LPI program structure. 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpigetsol
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

% Check that the program structure is properly specified.
if ~isa(prog,'struct')
    error('LPI program structure should be specified as object of type ''struct''.')
end
if isempty(prog.solinfo.info)
    error('The specified LPI optimization program has not been solved; either ''lpisolve'' has not been run, or no solution was produced.')
end

% Extract the solution.
if isempty(dvar) || isa(dvar,'double') || isa(dvar,'polynomial') || isa(dvar,'opvar')
    % The variable is not really a variable at all;
    dvar_sol = dvar;
elseif isa(dvar,'dpvar') || ischar(dvar) || isa(dvar,'cell')
    % The variable is some decision variable, or linear combination of
    % decision variables.
    dvar_sol = sosgetsol(prog,dvar);
elseif isa(dvar,'dopvar') || isa(dvar,'dopvar2d')
    % The variable is a PI operator decision variable.
    dvar_sol = getsol_lpivar(prog,dvar);
else
    error('Variable of which to extract the solution should be specified as object of type ''dpvar'' or ''dopvar''.')
end

end
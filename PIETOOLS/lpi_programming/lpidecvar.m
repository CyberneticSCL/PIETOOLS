function [prog,dvar_out] = lpidecvar(prog,dvar)
% [PROG,DVAR_OUT] = LPIDECVAR(PROG,DVAR) takes an LPI program structure 
% "prog" and adds to it a decision variable "dvar", or a matrix-valued 
% decision variable of size "dvar". It returns the updated program 
% structure, as well as (if applicable) the newly generated matrix-valued
% decision variable "dvar_out", either equal to or of size "dvar". 
% The function is a shell for 'sosdecvar' and 'sospolymatrixvar'.
%
% INPUT
% - prog:       'struct' specifying an LPI optimization program structure
%               (see also 'lpiprogram'). The format should match that used
%               by SOSTOOLS.
% - dvar:       OPTION 1: nx1 array of type 'dpvar', specifying n decision 
%               variables to add to the optimization program.
%               OPTION 2: nx1 cell array specifying names of n decision 
%               variables to be generated and added to the program.
%               OPTION 3: 1x2 array of type 'double', specifying the
%               size of a new matrix-valued decision variable to be
%               generated and added to the optimization program. Defaults
%               to [1,1].
%
% OUTPUT
% - prog:       'struct' specifying the same optimization program as passed
%               by the user, but now with the decision variables "dvar_out"
%               added to the structure, specifically in the field
%               'prog.decvartable'.
% - dvar_out:   OPTION 1: nx1 array of type 'dpvar', equal to 'dvar'.
%               OPTION 2: nx1 array of type 'dpvar', with each element a
%               decision variable with name given by the corresponding
%               element of "dvar".
%               OPTION 3: Array of type 'dpvar' and of size given by
%               "dvar", specifying decision variables that have been newly
%               generated and added to the optimization program structure.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpidecvar
%
% Copyright (C)2024  PIETOOLS Team
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
% DJ, 10/19/2024: Initial coding;
% DJ, 12/14/2024: Allow function to be used to generate decision variables;

% Check that the inputs make sense.
if ~isa(prog,'struct')
    error('LPI program structure should be specified as object of type ''struct''.')
elseif ~isfield(prog,'decvartable')
    error('The LPI program structure is not well-defined; use ''lpiprogram'' to initialize your program.')
end
if nargin==1
    % Assume the user wants to generate a scalar decision variable.
    dvar = [1,1];
elseif isnumeric(dvar)                                                      % DJ, 12/14/2024
    % dvar specifies the size of a desired matrix-valued decision variable.
    if isscalar(dvar)
        % Assume a column vector of decision variables is desired.
        dvar = [dvar,1];
    elseif any(size(dvar)~=[1,2])
        error('Number of decision variables should be specified as scalar or 1x2 array.')
    end
    if any(dvar<=0) || any(round(dvar)~=dvar)
        error('Number of decision variables should be specified by strictly positive integers.')
    end
elseif ischar(dvar) || iscellstr(dvar)                                      % DJ, 12/14/2024
    % Convert variable names to actual variables.
    dvar = dpvar(dvar);
elseif ~isa(dvar,'dpvar')
    error('Decision variables to add to the optimization program should be specified as nx1 cellstr or array of type ''dpvar''.')
end

% Declare the actual decision variable(s).
if isa(dvar,'dpvar')
    % dvar corresponds to a variable to add to the LPI program
    % --> just use sosdecvar...
    prog = sosdecvar(prog,dvar);
    dvar_out = dvar;
else                                                                        % DJ, 12/14/2024
    % dvar specifies size of a desired matrix-valued decision variable
    % --> add polynomial matrix variable with just monomial 1
    [prog,dvar_out] = sospolymatrixvar(prog,1,dvar);
end

end
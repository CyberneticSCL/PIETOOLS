function prog = lpidecvar(prog,dvar)
% LPIDECVAR takes an LPI program structure 'prog' and adds to it a decision
% variable 'dvar'. The function is just a shell for 'sosdecvar'.
%
% INPUT
% - prog:       'struct' specifying an LPI optimization program structure
%               (see also 'lpiprogram'). The format should match that used
%               by SOSTOOLS.
% - dvar:       nx1 array of type 'dpvar', specifying n decision variables
%               to add to the optimization program.
%
% OUTPUT
% - prog:       'struct' specifying the same optimization program as passed
%               by the user, but now with the decision variables 'dvar'
%               added to the structure, specifically in the field
%               'prog.decvartable'.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpidecvar
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
if ~isa(dvar,'dpvar')
    error('Decision variables to add to the optimization program should be specified as nx1 array of type ''dpvar''.')
end

% Just use sosdecvar...
prog = sosdecvar(prog,dvar);

end
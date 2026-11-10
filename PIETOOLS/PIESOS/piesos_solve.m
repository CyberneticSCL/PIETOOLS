function prog = piesos_solve(prog,opts)
% PROG = PIESOS_SOLVE(PROG,OPTS) solves the PIESOS program PROG, by passing
% it to `lpisolve` (a shell for `sossolve`) with options OPTS, and returns
% the solved PIESOS program structure.
%
% INPUTS
% - prog:       'struct' specifying a PIESOS program (based on the
%               SOSTOOLS program structure) defining a distributed SOS
%               program to solve (see also `piesos_program`);
% - opts:       (optional) 'struct' of options to be passed to SOSTOOLS
%               for solving the optimization program specified by PROG;
%
% OUTPUTS
% - prog:       'struct' specifying the solved PIESOS program;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_solve
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

% The function just acts as a shell for lpisolve, which acts as shell for
% sossolve...
if nargin==1
    prog = lpisolve(prog);
else
    prog = lpisolve(prog,opts);
end

end
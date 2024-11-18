function prog_sol = lpisolve(prog,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG_SOL = LPISOLVE(PROG,OPTS) takes an LPI program structure 'prog', and
% passes this to sossolve to solve with options 'opts', returning the
% solved LPI program structure.
% 
% INPUT
% - prog:       An LPI program structure (based on SOSTOOLS structure)
%               specifying an LPI to solve;
% - opts:       'struct' of options to be passed to SOSTOOLS for solving
%               the optimization program specified by lpi;
%
% OUTPUT
% - prog:       solved LPI program() structure for desired LPI;
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% DJ, 10/18/2024 - Initial coding

% % Solve using 'sossolve' with desired settings, if specified.
if nargin>=2
    prog_sol = sossolve(prog,opts);
else
    prog_sol = sossolve(prog);
end

end
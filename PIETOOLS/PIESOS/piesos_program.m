function prog = piesos_program(vartab,decvartab)
% PROG = PIESOS_PROGRAM(VARTAB,DECVARTAB) declares a PIESOS program
% structure PROG in the distributed variables VARTAB and the decision
% variables DECVARTAB.
%
% INPUTS
% - vartab:     n x 1 'polyopvar' object, specifying the n distributed
%               polynomial variables (state variables) of the
%               optimization program;
% - decvartab:  (optional) q x 1 'dpvar' array, specifying the decision
%               variables to be used in the optimization program;
%
% OUTPUTS
% - prog:       'struct' specifying an optimization program in the format
%               of SOSTOOLS, as declared by `sosprogram`, in the spatial
%               variables of VARTAB (and the decision variables
%               DECVARTAB), with an additional field 'pievartable',
%               a cellstr with the names VARTAB.varname of the
%               distributed variables;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_program
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

% Check the input
if ~isa(vartab,'polyopvar')
    error("Distributed polynomial variables must be specified as 'polyopvar' object'.")
end
% Initialize a program structure in the spatial and decision variables
pvartab = polynomial(unique([vartab.vars.in,vartab.vars.out]));
if nargin==1
    prog = sosprogram(pvartab);
else
    if ~isa(decvartab,'dpvar')
        error("Decision variables should be specified as nx1 'dpvar' array.")
    end
    prog = sosprogram(pvartab,decvartab);
end
% Add the distributed state variables
prog.pievartable = vartab.varname;


end
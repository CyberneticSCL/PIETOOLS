function [prog,Pop,Qcell,alpha_list] = possopvar(prog,dim,vars,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,Qcell,alpha_list] = POSSOPVAR(prog,dim,vars,dom,deg,options)
% declares a positive semidefinite, self-adjoint 'sdopvar' decision operator
%
%       Pop: L_2^m[S3] -> L_2^m[S3],       Pop = Pop* >= 0
%
% parameterized as
%
%       Pop = sum_{i,j} (Z_{alpha_i})* Q_{ij} (Z_{alpha_j}),     Q >= 0.
%
% This is a wrapper. The construction lives in 'sopquadvar', which assembles
% that quadratic form from a basis list and a coefficient matrix and takes
% the variable type as an argument; 'possopvar' fixes the type to 'pos'.
% Every input and output, and every option other than 'type', is exactly as
% documented in 'sopquadvar' -- see there for dim, vars, dom, deg,
% options.psatz, options.sep, options.include and the returned Qcell and
% alpha_list.
%
% Passing options.type here is an error rather than a silent override, since
% a 'sym' variable would not be a positive operator and the name would then
% be wrong.
%
% See also SOPQUADVAR, POSCOPVAR, EQ_OPTS_SOPVAR, DEGBALANCE,
% SETTINGS2POSSOPVAR, LPI_EQ_SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% MP, 08/22/2026:  Initial coding
% MMP, 09/21/2026: Body moved to 'sopquadvar' and this reduced to a wrapper.
%                  The body never imposed positivity itself -- it assembles
%                  a quadratic form, and the whole of the positivity
%                  question was the literal 'pos' passed to 'sosquadvar' --
%                  so it belonged in a constructor that Sec. 8.4's
%                  'poscopvar' can share, rather than being duplicated. The
%                  earlier header entries for the sep option, deg.subset,
%                  the sorted-S3 cell order and the block accumulation now
%                  live with the code in 'sopquadvar'. Behaviour here is
%                  unchanged, checked against the Gram dimensions measured
%                  beforehand and the 'test_possopvar' suite.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  posmopvar -> poscopvar.

if nargin<6 || isempty(options)
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if isfield(options,'type') && ~isempty(options.type) ...
        && ~strcmp(char(options.type),'pos')
    error("'possopvar' declares a positive operator; for type '%s' call "...
          +"'sopquadvar' directly.",char(options.type))
end
options.type = 'pos';

if nargin<5
    error("Not enough input arguments.")
end
[prog,Pop,Qcell,alpha_list] = sopquadvar(prog,dim,vars,dom,deg,options);

end

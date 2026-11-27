function loc = dvar_rows(dvars_old,dmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOC = DVAR_ROWS(DVARS_OLD,DMAP) returns the row of the global decision
% variable basis that each name of DVARS_OLD occupies, LOC(i) =
% DMAP(DVARS_OLD{i}), by one lookup of DMAP per name.
%
% INPUT
% - dvars_old:  'cellstr', n x 1, decision variable names, in the order of
%               the coefficient rows of the block;
% - dmap:       'containers.Map' from name (char) to global row index
%               (double scalar), with every name of DVARS_OLD as a key.
%
% OUTPUT
% - loc:        'double', n x 1, the global row indices.
%
% NOTES
% A name of DVARS_OLD that is not a key of DMAP is an error.
%
% See also SOPQUADVAR, OPQUADVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dvar_rows
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
% Initial coding MMP, 09/22/2026

% The global row of every name, by one vectorized lookup
try
    v = values(dmap,dvars_old);
catch err
    % Missing name: report the first. Anything else is rethrown.
    tf = isKey(dmap,dvars_old);
    if all(tf),  rethrow(err);  end
    error("Internal error: unrecognized decision variable '"...
          +string(dvars_old{find(~tf,1)})+"'.")
end
loc = reshape([v{:}],[],1);     % values are scalar doubles (num2cell(1:n))

end

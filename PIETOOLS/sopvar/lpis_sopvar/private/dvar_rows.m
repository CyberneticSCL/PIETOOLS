function loc = dvar_rows(dvars_old,dmap)
% LOC = DVAR_ROWS(DVARS_OLD,DMAP) returns the row of the global decision
% variable basis that each of DVARS_OLD occupies, where DMAP is the
% containers.Map from name to global index.
%
% This is the lookup half of 'remap_dvars'. It exists so that a caller can
% run the coefficient elimination on a block's OWN decision rows and scatter
% to the global basis once at the end, rather than widening the block to the
% global row count first: the decision-variable axis is a pure batch axis, so
% widening it early pads every intermediate with zeros. Measured in
% 'sopquadvar' at two spatial variables and degree 2, the early widening made
% the sheet axis 134 to 2119 times wider than the block needed.
%
% The lookup is a hash per name rather than a search over the whole global
% list, so the cost is in the number of names the block actually carries, not
% in the size of the basis.
%
% INPUTS
% - dvars_old:  n x 1 'cellstr' of decision variable names, in the order the
%               block's coefficient rows are in;
% - dmap:       'containers.Map' from name to global row index;
%
% OUTPUTS
% - loc:        n x 1 array of global row indices.
%
% See also REMAP_DVARS, SOPQUADVAR, MOPQUADVAR.
%
% MMP, 09/22/2026: Initial coding, split out of 'remap_dvars'.

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

loc = zeros(numel(dvars_old),1);
for i = 1:numel(dvars_old)
    if ~isKey(dmap,dvars_old{i})
        error("Internal error: unrecognized decision variable '"...
              +string(dvars_old{i})+"'.")
    end
    loc(i) = dmap(dvars_old{i});
end

end

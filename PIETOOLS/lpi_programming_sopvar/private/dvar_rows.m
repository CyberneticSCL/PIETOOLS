function loc = dvar_rows(dvars_old,dmap)
% LOC = DVAR_ROWS(DVARS_OLD,DMAP) returns the row of the global decision
% variable basis that each of DVARS_OLD occupies, where DMAP is the
% containers.Map from name to global index.
%
% (was) This is the lookup half of 'remap_dvars'. It exists so that a caller can % MMP, 09/30/2026 (was)
% This lookup exists so that a caller can                                   % MMP, 09/30/2026
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
% (was) See also REMAP_DVARS, SOPQUADVAR, COPQUADVAR.                       % MMP, 09/30/2026 (was)
% See also SOPQUADVAR, COPQUADVAR.                                          % MMP, 09/30/2026
%
% MMP, 09/22/2026: Initial coding, split out of 'remap_dvars'.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change.
% MMP, 09/26/2026: One vectorized 'values' call instead of a per-name
%                  'isKey' and subsref loop; 'isKey' now runs only to name
%                  the culprit when 'values' fails. Same output. Measured
%                  in situ, 2-D Hinf build (125 calls, 7.1e5 names): 3.58 s
%                  -> 0.28 s; 2-D heavy stability (7.9e5 names, 9.1e5
%                  global): 3.40 -> 0.32 s. Standalone, 1e5 names against
%                  1e6 global: 0.58 -> 0.16 s. Still a hash per name, so
%                  linear in the names looked up and independent of the
%                  global count; 'ismember' against the sorted global list
%                  was measured and rejected: it pays for the global list on
%                  every call, 0.64 s against 0.013 s for 1e4 names of 1e6.
% MMP, 09/30/2026: 'remap_dvars' deleted: no caller since 09/22/2026 (its
%                  scatter half is built from triplets in the callers). Help
%                  no longer introduces this routine as its lookup half; the
%                  hash-per-name rationale it held is the paragraph above.
%                  Comments only.

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

% loc = zeros(numel(dvars_old),1);                                          % MMP, 09/26/2026 (was)
% for i = 1:numel(dvars_old)                                                % MMP, 09/26/2026 (was)
%     if ~isKey(dmap,dvars_old{i})                                          % MMP, 09/26/2026 (was)
%         error("Internal error: unrecognized decision variable '"...
%               +string(dvars_old{i})+"'.")                                 % MMP, 09/26/2026 (was)
%     end                                                                   % MMP, 09/26/2026 (was)
%     loc(i) = dmap(dvars_old{i});                                          % MMP, 09/26/2026 (was)
% end                                                                       % MMP, 09/26/2026 (was)
try                                                                         % MMP, 09/26/2026
    v = values(dmap,dvars_old);                                             % MMP, 09/26/2026
catch err                                                                   % MMP, 09/26/2026
    % Missing name: report the first, as before. Anything else is rethrown.
    tf = isKey(dmap,dvars_old);                                             % MMP, 09/26/2026
    if all(tf),  rethrow(err);  end                                         % MMP, 09/26/2026
    error("Internal error: unrecognized decision variable '"...
          +string(dvars_old{find(~tf,1)})+"'.")                             % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
loc = reshape([v{:}],[],1);     % values are scalar doubles (num2cell(1:n)) % MMP, 09/26/2026

end

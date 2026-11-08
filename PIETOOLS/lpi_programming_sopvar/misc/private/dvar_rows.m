function loc = dvar_rows(dvars_old,dmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOC = DVAR_ROWS(DVARS_OLD,DMAP) returns the row of the global decision
% variable basis that each name of DVARS_OLD occupies, LOC(i) =
% DMAP(DVARS_OLD{i}), by one lookup of DMAP per name.
%
% INPUT
% - dvars_old:  'cellstr', n x 1, decision variable names, in the order of
%               the block's coefficient rows;
% - dmap:       'containers.Map' from name (char) to global row index
%               (double scalar), with every name of DVARS_OLD as a key.
%
% OUTPUT
% - loc:        'double', n x 1, the global row indices.
%
% NOTES
% A name of DVARS_OLD that is not a key of DMAP is an error.
%
% See also SOPQUADVAR, COPQUADVAR.
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
% MMP, 09/22/2026: Initial coding, split out of 'remap_dvars'.
% MMP, 09/25/2026: Container classes renamed mopvar -> copvar, mdopvar ->
%                  cdopvar (mechanical, no functional change).
% MMP, 09/26/2026: One vectorized 'values' call instead of a per-name
%                  'isKey' and subsref loop; 'isKey' runs only to name the
%                  culprit when 'values' fails. Same output, 10x faster in
%                  the 2-D builds. Still a hash per name: linear in the
%                  names looked up and independent of the global count
%                  ('ismember' against the global list pays for that list
%                  on every call).
% MMP, 09/30/2026: 'remap_dvars' deleted: no caller since 09/22/2026 (its
%                  scatter half is built from triplets in the callers).
%                  Comments only.
% MMP, 10/28/2026: Header to the header protocol, step comments in the body;
%                  code unchanged. Deleted: the rationale for the lookup
%                  with a measurement, the cost of a lookup per name and two
%                  (was) lines.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed;
%                  code unchanged.

% The global row of every name, by one vectorized lookup
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

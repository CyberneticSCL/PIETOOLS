function info = verify(P)                                                   % MMP, 09/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = VERIFY(P) checks that the blocks of a 'cdopvar' agree with its     % MMP, 09/30/2026
% container metadata and with each other.
%
% OUTPUTS
% - info.true:   1 if P is consistent, 0 otherwise;                         % MMP, 09/30/2026
% - info.flags:  cell array of messages, one per inconsistency, empty when
%                P is consistent;                                           % MMP, 09/30/2026
%
% CHECKED: metadata shapes; no all-empty row or column; per block, that its
% output side matches its row and its input side matches its column, that it
% puts every variable on the registry's domain, and that an 'sdopvar' block
% carries the container's Zd.
%
% The space checks are SET equality of the names, not 'isequal' on the
% lists: the block constructors enforce vars.out = [S2,S3] and vars.in =
% [S3,S1], so blocks in one row order the same output set differently
% whenever their input spaces differ, and blocks in one column likewise.
% Comparing with 'isequal' would reject legitimate containers.
%
% The constructor already rejects anything reported here, so a 'cdopvar'
% built through 'cdopvar(C)' always verifies. VERIFY is for objects reached
% another way: properties assigned directly, the two-argument trusted
% constructor, or a block replaced in place via P.C{i,j} = ....
%
% Cost: O(M*N) block visits, spatial set operations on the registry, and
% O(q) string comparisons per 'sdopvar' block guarded by an O(1) length test.
%
% The checks are in 'verify_copvar_meta', shared with 'copvar'.             % MMP, 09/30/2026
%
% See also CDOPVAR, COPVAR, SIZE.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - verify(cdopvar)
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
% Initial coding MMP, 09/17/2026. Split out of @copvar/verify, which until
%                this date checked both cases; that file keeps the history
%                of MP's 01/19/2026 draft. This copy differs from it only by
%                admitting 'sdopvar' blocks and by checking the decision
%                variable list.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Moved from
%                  @mdopvar/ with the class.
% MMP, 09/30/2026: The body (initial coding 09/17/2026: the checks and the
%                  local 'check_dom') is deleted and moved to
%                  'verify_copvar_meta', which @copvar/verify, of which it
%                  was a copy, now calls too. One algorithm had two copies;
%                  derive_copvar_meta's NOTES name exactly this validation
%                  as what must not diverge. Verdicts and messages unchanged.
% MMP, 09/30/2026: Renamed the input Mop -> P, a name from before the
%                  09/25/2026 rename ('mopvar' is now the stub class in
%                  @mopvar). Mechanical, no functional change; each marked
%                  line differs from its old text only by that name.

% % % BEGIN body replaced by MMP, 09/30/2026 - the 09/17/2026 body, deleted
% % % here, is 'verify_copvar_meta'; see the header entry above.
info = verify_copvar_meta(P,'cdopvar',{'sopvar','sdopvar'});                % MMP, 09/30/2026

end
% % % END body replaced by MMP, 09/30/2026

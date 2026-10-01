function info = verify(Mop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = VERIFY(MOP) checks that the blocks of a 'copvar' agree with its
% container metadata and with each other.
%
% OUTPUTS
% - info.true:   1 if MOP is consistent, 0 otherwise;
% - info.flags:  cell array of messages, one per inconsistency, empty when
%                MOP is consistent;
%
% CHECKED: metadata shapes; no all-empty row or column; per block, that its
% output side matches its row and its input side matches its column, that it
% puts every variable on the registry's domain.
%
% The space checks are SET equality of the names, not 'isequal' on the
% lists: the block constructors enforce vars.out = [S2,S3] and vars.in =
% [S3,S1], so blocks in one row order the same output set differently
% whenever their input spaces differ, and blocks in one column likewise.
% Comparing with 'isequal' would reject legitimate containers.
%
% The constructor already rejects anything reported here, so a 'copvar'
% built through 'copvar(C)' always verifies. VERIFY is for objects reached
% another way: properties assigned directly, the two-argument trusted
% constructor, or a block replaced in place via P.C{i,j} = ....
%
% Cost: O(M*N) block visits, with spatial set operations on the registry.
%
% The checks are in 'verify_copvar_meta', shared with 'cdopvar'.            % MMP, 09/30/2026
%
% See also COPVAR, SIZE.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - verify(copvar)
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
% MP, 01/19/2026: Initial coding
% MMP, 09/17/2026: Replaced the body. The 01/19/2026 draft could not run:
%                  its first executable line read 'C' and 'fields' before
%                  either was assigned, 'r' was never assigned, dynamic
%                  field access was attempted with the non-field names
%                  'dim(1)' and 'vars(1)', and 'row = Mop{i,:}' takes only
%                  the first element of a comma list. It also compared
%                  variable lists for equality, which rejects legitimate
%                  containers for the canonical-order reason above. All of
%                  it is deleted rather than commented out; the 01/19/2026
%                  entry no longer describes code present here. The draft's
%                  commented-out cellfun sketches noted that empty elements
%                  should be eliminated; that is now the zero-block
%                  convention, under which [] is a legal block whose space
%                  and dimension come from the container metadata.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Moved from
%                  @mopvar/ with the class.
% MMP, 09/30/2026: The body, between the 'BEGIN/END body replaced by MMP,
%                  09/17/2026' markers (the checks and the local
%                  'check_dom'), is deleted and moved to 'verify_copvar_meta',
%                  which @cdopvar/verify, a copy of it, now calls too. One
%                  algorithm had two copies; derive_copvar_meta's NOTES name
%                  exactly this validation as what must not diverge.
%                  Verdicts and messages unchanged. The 09/17/2026 entry
%                  describes code that now lives there.

% % % BEGIN body replaced by MMP, 09/30/2026 - the 09/17/2026 body, deleted
% % % here, is 'verify_copvar_meta'; see the header entry above.
info = verify_copvar_meta(Mop,'copvar',{'sopvar'});                         % MMP, 09/30/2026

end
% % % END body replaced by MMP, 09/30/2026

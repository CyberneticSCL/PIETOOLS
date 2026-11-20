function info = verify(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = VERIFY(P) checks that the blocks of a 'cdopvar' agree with its
% container metadata and with each other.
%
% INPUTS
% - P:      'cdopvar' container to be checked;
%
% OUTPUTS
% - info.true:   1 if P is consistent, 0 otherwise;
% - info.flags:  cell array of messages, one per inconsistency, empty when
%                P is consistent;
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
% Initial coding MMP, 09/17/2026

% Check the metadata, and every 'sopvar' or 'sdopvar' block against it
info = verify_copvar_meta(P,'cdopvar',{'sopvar','sdopvar'});

end

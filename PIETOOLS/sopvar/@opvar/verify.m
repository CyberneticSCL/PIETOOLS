function info = verify(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = VERIFY(P) checks that the blocks of an 'opvar' agree with its
% container metadata and with each other.
%
% INPUTS
% - P:      'opvar' container to be checked;
%
% OUTPUTS
% - info.true:   1 if P is consistent, 0 otherwise;
% - info.flags:  cell array of messages, one per inconsistency, empty when
%                P is consistent;
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
% The constructor already rejects anything reported here, so an 'opvar'
% built through 'opvar(C)' always verifies. VERIFY is for objects reached
% another way: properties assigned directly, the two-argument trusted
% constructor, or a block replaced in place via P.C{i,j} = ....
%
% Cost: O(M*N) block visits, with spatial set operations on the registry.
%
% See also OPVAR, SIZE.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - verify(opvar)
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

% Check the metadata, and every 'sopvar' block against it
info = verify_opvar_meta(P,'opvar',{'sopvar'});

end

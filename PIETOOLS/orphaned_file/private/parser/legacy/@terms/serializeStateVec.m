function out = serializeStateVec(statevec)
% OUT = SERIALIZESTATEVEC(STATEVEC) returns in row i the byte stream
% (getByteStreamFromArray) of the state STATEVEC(i).
%
% INPUT
% - statevec:   column of 'state' objects;
%
% OUTPUT
% - out:        'uint8' array, one row per state.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding SS - 02/14/2022
% MMP, 11/19/2026: Release comment pass (comment-only).

out = zeros(length(statevec):4200);
for i=1:length(statevec)
    out(i,:) = getByteStreamFromArray(statevec(i));
end
end
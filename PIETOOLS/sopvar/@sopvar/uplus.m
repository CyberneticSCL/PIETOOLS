function Pop_in = uplus(Pop_in)
% POP_IN = UPLUS(POP_IN) returns the 'sopvar' object POP_IN itself, so that
% the unary +P is defined for a PI operator P.
%
% INPUT
% - Pop_in:     m x n 'sopvar' object;
%
% OUTPUT
% - Pop_in:     the same m x n 'sopvar' object, unchanged.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - uplus
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
% Initial coding MP, SS, DJ, AT, DB - 10/24/2026
% MMP, 11/06/2026: Help text corrected: it described (-1)*P, copied from
%                  'uminus'; +P returns P. Code unchanged.

% Nothing to do: +P is P
end
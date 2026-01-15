function Pop = uminus(Pop)
% POP = UMINUS(POP) returns the 'ndopvar' object representing the scalar
% product (-1)*P for the PI operator P defined by POP.
%
% INPUTS
% - Pop_in:     m x n 'ndopvar' object;
%
% OUTPUTS
% - Pop_out:    m x n 'ndopvar' object representing the operator defined by
%               (-1)*Pop_in;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - uminus
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
% DJ, 01/15/2026: Initial coding

% The coefficients defining (-1)*P are just (-1) times the coefficients 
% defining P;
for ii=1:numel(Pop.C)
    Pop.C{ii} = -Pop.C{ii};
end

end
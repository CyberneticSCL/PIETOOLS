function P = subsasgn(P,S)
% P = SUBSASGN_1(P,S) throws an error: 'quadPoly' objects may not be
% modified.
%
% INPUTS
% - P:      'quadPoly' object;
% - S:      'struct' array of the indexing operations;
%
% OUTPUTS
% - P:      never set, an error is thrown.
%
% NOTES
% - The function line names subsasgn, but MATLAB names a function after
%   its file, so this file defines the method subsasgn_1 and does not
%   overload subsasgn: an assignment such as P.dim = [m,n] uses the
%   built-in subsasgn and sets the property without a check.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsasgn_1
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% SS, 02/11/2026: Initial coding;

error("Modification of quadPoly objects is not permitted. Use quadPoly() constructor to create new object")
end
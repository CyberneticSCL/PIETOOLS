function n = numArgumentsFromSubscript(obj,s,indexingContext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method is automatically called when one of the following statements
% are executed:
% a) obj.s — Number of elements referenced in a statement (indexingContext = Statement) 
% b) (obj.s) — Number of elements returned in an expression (indexingContext = Expression) 
% c) [obj.s] = rhs — Number of values assigned with a comma-separated list (indexingContext = Assignment)
% For a 'state' object, every such indexing returns one output.
%
% INPUT
% - obj:              'state' object;
% - s:                'struct' array, the indexing (fields type, subs);
% - indexingContext:  matlab.mixin.util.IndexingContext, Statement,
%                     Expression or Assignment.
%
% OUTPUT
% - n:                'double', the number of outputs, always 1.

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
% Initial coding SS - 07/01/2022
% MMP, 11/19/2026: Release comment pass (comment-only).

   n=1;
end 
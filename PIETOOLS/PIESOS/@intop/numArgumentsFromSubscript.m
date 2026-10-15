function n = numArgumentsFromSubscript(~,~,~)
% N = NUMARGUMENTSFROMSUBSCRIPT(OBJ,S,INDEXINGCONTEXT) returns the number
% of values N expected from an indexing operation on the 'intop' object
% OBJ, which is always 1. This method is called automatically when one of
% the following statements is executed:
% a) obj.s - number of elements referenced in a statement
%    (indexingContext = Statement);
% b) (obj.s) - number of elements returned in an expression
%    (indexingContext = Expression);
% c) [obj.s] = rhs - number of values assigned with a comma-separated list
%    (indexingContext = Assignment);
%
% INPUTS
% - obj:              'intop' object that is indexed (not used);
% - s:                struct array specifying the indexing operation (not
%                     used);
% - indexingContext:  context of the indexing operation, Statement,
%                     Expression, or Assignment (not used);
%
% OUTPUTS
% - n:      scalar equal to 1, i.e. one value is expected in all cases;
%
% NOTES
% - The function 'numel' is overloaded for 'intop' objects, so this
%   function is needed to avoid breaking indexing and indexed assignment.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - numArgumentsFromSubscript
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
% DJ, 10/15/2026: Initial coding;

n = 1;

end
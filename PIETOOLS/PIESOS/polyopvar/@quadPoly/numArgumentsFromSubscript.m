function n = numArgumentsFromSubscript(obj,s,indexingContext)
% N = NUMARGUMENTSFROMSUBSCRIPT(OBJ,S,INDEXINGCONTEXT) returns the number
% of outputs N expected from an indexing operation on the 'quadPoly' OBJ.
% This method is automatically called when one of the following
% statements is executed:
% a) obj.s - number of elements referenced in a statement
%    (indexingContext = Statement);
% b) (obj.s) - number of elements returned in an expression
%    (indexingContext = Expression);
% c) [obj.s] = rhs - number of values assigned with a comma-separated list
%    (indexingContext = Assignment).
% Since we overloaded the function 'numel', we need this function to avoid
% breaking 'subsref' and 'subsasgn'. It returns 1 in all cases.
%
% INPUTS
% - obj:              'quadPoly' object (not used);
% - s:                'struct' array of the indexing operations (not
%                     used);
% - indexingContext:  matlab.indexing.IndexingContext, Statement,
%                     Expression or Assignment (not used);
%
% OUTPUTS
% - n:                integer, 1.

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% SS, 06/14/2022: Initial coding;

   n=1;
end
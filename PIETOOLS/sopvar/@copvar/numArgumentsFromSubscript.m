function n = numArgumentsFromSubscript(P,s,ctx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = numArgumentsFromSubscript(P,s,ctx) returns the number of outputs of
% the indexing expression s on the copvar object P, which MATLAB passes to
% the overloaded 'subsref' / 'subsasgn' as nargout.
%
% INPUT
% P:    'copvar' object
% s:    'struct' array of subscripts (fields 'type' and 'subs') of the
%       indexing expression, as passed by MATLAB
% ctx:  'matlab.indexing.IndexingContext' of the expression (statement,
%       expression or assignment)
%
% OUTPUT
% n:    nonnegative integer, the number of outputs MATLAB should expect
%
% NOTES:
% With 'subsref' overloaded, the default count of MATLAB is 1, which keeps
% only the first element of a comma list such as {P.C{:}} or [P.C{:}]. The
% builtin count returns every element, as for a class without the
% overload. A read starting with '()' is the slice P(I,J), one object,
% so it counts 1. So does a read of a field of the 1-D or 2-D layout, P.Q1,
% P.R22{2,3}, which 'subsref' computes and serves as one value; the
% builtin count errors on it, as it is not a property.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - numArgumentsFromSubscript (copvar)
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
%
% Initial coding MMP, 09/29/2026

% The slice P(I,J) is one object; assignment keeps the builtin count.
% So is a field of the 1-D or 2-D layout, which 'subsref' serves as one value.
opfields = {'P','Q1','Q2','R','I','var1','var2','dim','R00','R0x','R0y',...
    'R02','Rx0','Rxx','Rxy','Rx2','Ry0','Ryx','Ryy','Ry2','R20','R2x','R2y','R22'};
if ctx~=matlab.indexing.IndexingContext.Assignment && (strcmp(s(1).type,'()') ...
        || (strcmp(s(1).type,'.') && any(strcmp(s(1).subs,opfields))))
    n = 1;
    return
end
% Any other expression: the builtin count
n = builtin('numArgumentsFromSubscript',P,s,ctx);

end

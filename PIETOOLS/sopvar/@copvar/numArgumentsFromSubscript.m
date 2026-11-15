function n = numArgumentsFromSubscript(P,s,ctx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = numArgumentsFromSubscript(P,s,ctx) returns the number of outputs of
% the indexing expression s on the copvar object P, which MATLAB passes to
% the overloaded 'subsref' / 'subsasgn' as nargout.
%
% INPUT
% P:    copvar class object
% s:    struct array of subscripts (fields 'type' and 'subs') of the
%       indexing expression, as passed by MATLAB
% ctx:  matlab.indexing.IndexingContext of the expression (statement,
%       expression or assignment)
%
% OUTPUT
% n:    nonnegative integer, the number of outputs MATLAB should expect
%
% NOTES:
% Once 'subsref' is overloaded, MATLAB's default count is 1, so a
% comma list such as {P.C{:}} or [P.C{:}] silently kept only its first
% element. The builtin count restores the behaviour of the class without
% the overload. A read starting with '()' is the slice P(I,J), one object,
% so it counts 1. So does a read of an 'opvar' or 'opvar2d' field, P.Q1,
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
% Initial coding MMP, 09/29/2026, with 'subsref' and 'subsasgn'.
% DJ, 11/14/2026: Return 1 for ()-type subsref;
% DJ, 11/15/2026: Return 1 for a read of an 'opvar' or 'opvar2d' field,
%                which 'subsref' now serves;

% The slice P(I,J) is one object; assignment keeps the builtin count.
% So is an 'opvar'/'opvar2d' field, which 'subsref' serves as one value.   % DJ, 11/15/2026
opfields = {'P','Q1','Q2','R','I','var1','var2','dim','R00','R0x','R0y',...
    'R02','Rx0','Rxx','Rxy','Rx2','Ry0','Ryx','Ryy','Ry2','R20','R2x','R2y','R22'}; % DJ, 11/15/2026
% if strcmp(s(1).type,'()') && ctx~=matlab.indexing.IndexingContext.Assignment % DJ, 11/15/2026 (was)
if ctx~=matlab.indexing.IndexingContext.Assignment && (strcmp(s(1).type,'()') ...
        || (strcmp(s(1).type,'.') && any(strcmp(s(1).subs,opfields))))  % DJ, 11/15/2026
    n = 1;
    return
end
n = builtin('numArgumentsFromSubscript',P,s,ctx);

end

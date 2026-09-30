function n = numArgumentsFromSubscript(P,s,ctx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = NUMARGUMENTSFROMSUBSCRIPT(P,S,CTX) returns the builtin number of
% outputs of the indexing expression S on the 'cdopvar' P, which MATLAB
% passes to the overloaded 'subsref' / 'subsasgn' as nargout.
%
% Why: once 'subsref' is overloaded, MATLAB's default count is 1, so a
% comma list such as {P.C{:}} or [P.C{:}] silently kept only its first
% element (measured, R2025b). The builtin count restores the behaviour of
% the class without the overload.
%
% See also SUBSREF, SUBSASGN, COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - numArgumentsFromSubscript(cdopvar)
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

n = builtin('numArgumentsFromSubscript',P,s,ctx);

end

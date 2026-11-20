function n = numArgumentsFromSubscript(obj,s,ctx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = NUMARGUMENTSFROMSUBSCRIPT(OBJ,S,CTX) returns the number of values an
% indexing expression S on the 'sopvar' OBJ takes or yields.
%
% ASSIGNMENT, and a READ that starts with '.': the builtin count, so that a
% comma list returns or takes every element: {P.params{:}}, [P.params{:}],
% f(P.ZL{:}), [P.params{1:2}] = deal(a,b). With overloaded indexing the
% default count is 1, which keeps only the first element on reads and
% fails on assignment (measured on R2025b).
%
% A READ that starts with '()' (the component slice P(I,J)) or '{}' (an
% error): 1, since the slice is one object. A list after a slice,
% {P(I,J).ZL{:}}, is therefore not served: the subscripts after '()' are
% dropped and the slice is returned.
%
% INPUTS
% - obj:    'sopvar' object being indexed;
% - s:      struct array of subscripts, with fields 'type' ('.', '()' or
%           '{}') and 'subs';
% - ctx:    'matlab.indexing.IndexingContext' of the call, e.g. Assignment;
%
% OUTPUTS
% - n:      scalar integer, the number of values the expression takes or
%           yields;
%
% See also SUBSASGN, SUBSREF, SOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - numArgumentsFromSubscript(sopvar)
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
% Initial coding MMP, 09/29/2026
%

% Property names of the class, read once
persistent props
if strcmp(s(1).type,'.')
    % The builtin count costs about 2.3 us a call (measured), a quarter of
    % a read, so the common chains are counted here: a property; then at
    % most one field of a scalar struct property; then at most one brace
    % with a single numeric scalar, and at most one paren, last. Such a
    % chain names one value, n = 1, the builtin count. Any other chain (a
    % ':', vector or logical index, a struct array, a longer chain) is
    % counted by the builtin. A single '.' on a method name also gets 1.
    n = 1;
    L = numel(s);
    for k = 2:L
        t = s(k).type;
        if t(1)=='{'
            x = s(k).subs;
            if ~(isscalar(x) && isnumeric(x{1}) && isscalar(x{1}))
                n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
            end
        elseif t(1)=='('
            if k<L
                n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
            end
        elseif k==2
            % a field of a property: one value iff the property is a scalar
            % struct (read inside the class, so protected names too)
            if isempty(props)
                mc = ?sopvar;  props = {mc.PropertyList.Name};
            end
            if ~any(strcmp(s(1).subs,props))
                n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
            end
            v = obj.(s(1).subs);
            if ~(isstruct(v) && isscalar(v))
                n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
            end
        else
            % a second field: the builtin counts the chain
            n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
        end
    end
elseif ctx==matlab.indexing.IndexingContext.Assignment
    % Assignment: the builtin count
    n = builtin('numArgumentsFromSubscript',obj,s,ctx);
else
    % A read that starts with '()' or '{}': one value
    n = 1;
end

end

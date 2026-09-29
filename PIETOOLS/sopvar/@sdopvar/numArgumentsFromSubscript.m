function n = numArgumentsFromSubscript(obj,s,ctx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = NUMARGUMENTSFROMSUBSCRIPT(OBJ,S,CTX) returns the number of values an
% indexing expression S on the 'sdopvar' OBJ takes or yields, which MATLAB
% passes to the overloaded 'subsasgn' / 'subsref'.
%
% ASSIGNMENT, and a READ that starts with '.': the builtin count. Once
% 'subsref' / 'subsasgn' are overloaded, MATLAB's default is 1, so a comma
% list silently kept its first element on reads ({P.Zd{:}},
% [P.params.A{:}], f(P.ZL{:})) and failed on assignment ([P.Zd{1:2}] =
% deal(a,b): deal:narginNargoutMismatch), both measured on R2025b. The '.'
% case of @sdopvar/subsref serves the read by the builtin with this count.
%
% A READ that starts with '()' (the component slice P(I,J)) or '{}' (an
% error): 1, as before, since the slice is one object. A list after a
% slice, {P(I,J).Zd{:}}, therefore still holds its first element only
% (pre-existing).
%
% See also SUBSASGN, SUBSREF, SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - numArgumentsFromSubscript(sdopvar)
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
% Initial coding MMP, 09/29/2026, with 'subsasgn'. Reads that start with
%                '.' return the builtin count too (same date), with the
%                comma-list case of 'subsref'.

persistent props
if strcmp(s(1).type,'.')
    % The builtin count costs about 2.3 us a call (measured), a quarter of
    % a read, so the common chains are counted here: a property; then at
    % most one field of a scalar struct property; then at most one brace
    % with a single numeric scalar, and at most one paren, last. Such a
    % chain names one value, n = 1: the builtin count, and the count before
    % this method existed. Any other chain (a ':', vector or logical index,
    % a struct array, a longer chain) is counted by the builtin. A single
    % '.' on a method name also gets 1, as before.
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
                mc = ?sdopvar;  props = {mc.PropertyList.Name};
            end
            if ~any(strcmp(s(1).subs,props))
                n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
            end
            v = obj.(s(1).subs);
            if ~(isstruct(v) && isscalar(v))
                n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
            end
        else
            n = builtin('numArgumentsFromSubscript',obj,s,ctx);  return
        end
    end
elseif ctx==matlab.indexing.IndexingContext.Assignment
    n = builtin('numArgumentsFromSubscript',obj,s,ctx);
else
    n = 1;
end

end

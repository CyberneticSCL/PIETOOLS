function out = subsref(obj,s)
% OUT = SUBSREF(OBJ,S) extracts the elements of 'tensopvar' object OBJ
% specified by S. This function is automatically called 
% 
% INPUT
% - obj:    'tensopvar' object of which to extract certain elements;
% - s:      struct indicating which element/fields of obj to extract;
% 
% OUTPUT
% - out:    object representing the elements of 'obj' specified by 's'. If
%           s.type='.', this will be out=obj.(s.subs). If s.type.'()', the
%           output will be a tensopvar object including the elements in
%           the indices specified by s.subs;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEOOLS - subsref
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
% DJ, 01/16/2026: Initial coding

switch s(1).type
    case '.'
        % Just use the built-in function
        out = builtin('subsref',obj,s);
        
    case '()'
        % When calling C(i,j), just extract C.ops(i,j);
        ops = builtin('subsref',obj.ops,s);
        out = obj;
        out.ops = ops;
        
    case '{}'
        error("'{}' type indexing of 'tensopvar' objects is not supported.")
end
end

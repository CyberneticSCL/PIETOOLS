function out = subsasgn(obj,ref,RHS)
% OUT = SUBSASGN(OBJ,REF,RHS) assigns values RHS to field/slice REF of a
% tensopvar object OBJ
% 
% INPUT
% - obj:    'tensopvar' object of which to assign elements;
% - ref:    'struct' specifying which elements of obj to adjust;
% - RHS:    'tensopvar' object specifying the new values of the specified
%           elements;
%
% OUTPUT
% out:      'tensopvar' object with desired component value
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 01/16/2026: Initial coding


switch ref(1).type
    case '.'
        out = builtin('subsasgn',obj,ref,RHS);
    case '()'
        error("'()'-type subsassign is not supported for opvar objects.");
    case '{}'
        if ~isa(RHS,'tensopvar')
            error("For setting elements of a 'tensopvar', the new value must be specified as 'tensopvar' object as well.")
        end
        % When setting C{i,j}, just set C.ops(i,j);
        ref(1).type = '()';
        ops = builtin('subsasgn',obj.ops,ref,RHS.ops);
        out = obj;
        out.ops = ops;
end

end
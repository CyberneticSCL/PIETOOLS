function out = subsasgn(obj,ref,RHS)
% OUT = SUBSASGN(OBJ,REF,RHS) assigns values RHS to the field REF of a
% 'polyopvar' object OBJ. This function is automatically called when
% calling e.g. f.C = RHS or f.dims.in(2) = RHS for a 'polyopvar' object f,
% outside of the methods of the class.
%
% INPUTS
% - obj:    'polyopvar' object of which to assign a field;
% - ref:    'struct' specifying which field (and element of the field) of
%           obj to adjust;
% - RHS:    object specifying the new value of the field specified by ref;
%
% OUTPUTS
% - out:    'polyopvar' object with the desired field value, and with
%           out.chkval=false, indicating that the object has not been
%           checked (see initialize) since the field was changed.
%           Assigning f.chkval itself just sets the value;
%
% NOTES
% - Assignment to elements of a 'polyopvar' object, f(i,j)=RHS or
%   f{i,j}=RHS, is not supported.
% - Within the methods of the class, MATLAB uses the built-in field
%   assignment, so that the methods have to set chkval themselves.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsasgn
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 10/15/2026: Initial coding;

switch ref(1).type
    case '.'
        % Use the built-in assignment of the field (or element of the field)
        out = builtin('subsasgn',obj,ref,RHS);
        % The consistency of the fields has not been checked for the new value
        if ~strcmp(ref(1).subs,'chkval')
            out.chkval = false;
        end
    case '()'
        error("'()'-type subsasgn is not supported for polyopvar objects.")
    case '{}'
        error("'{}'-type subsasgn is not supported for polyopvar objects.")
end

end

function PIE = subsasgn(PIE,prop,val)
% SUBSASGN allows properties and elements of the "npie_struct" object PIE
% to be adjusted or set. The function is automatically called when calling
% e.g. PIE.obj = val for an "npie_struct" PIE, and sets PIE.is_initialized
% to false, since the PIE may no longer be properly specified.
%
% INPUT
% - PIE:    An npie_struct object of which to adjust a property/element.
% - prop:   A 1xq struct specifying which property/element of the input PIE
%           structure to adjust.
% - val:    Value to assign to the specified property/element of the PIE.
%
% OUTPUT
% - PIE:    An npie_struct object almost identical to the input system, but
%           with the element/prop specified by "prop" now set equal to
%           "val".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsasgn
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 10/29/2026: Initial coding;

PIE = builtin('subsasgn',PIE,prop,val);
% Any adjustment of the PIE may mean it is no longer properly specified.
if PIE.is_initialized && ~(strcmp(prop(1).type,'.') && strcmp(prop(1).subs,'is_initialized'))
    PIE.is_initialized = false;
end

end

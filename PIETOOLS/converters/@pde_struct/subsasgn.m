function PDE = subsasgn(PDE,prop,val)
% SUBSASGN allows properties and elements of the "pde_struct" object PDE to
% be adjusted or set. The function is automatically called when calling
% e.g. PDE(i,j) = val or PDE.obj = val for a "pde_struct" PDE.
%
% INPUTS:
% - PDE:    A pde_struct object of which to adjust a property/element.
% - prop:   A 1xq struct specifying which property/element of the input PDE
%           structure to adjust.
% - val:    Value to assign to the specified property/element of the PDE.
%
% OUTPUTS:
% - PDE:    A pde_struct object almost identical to the input system, but
%           with the element/prop specified by "prop" now set equal to
%           "val".
%
% EXAMPLE: calling PDE.x{1}.dom(2) = b, we have
%   PDE_in = PDE;
%   prop(1).type = '.',     prop(1).subs = 'x';
%   prop(2).type = '{}',    prop(2).subs = {[1]};
%   prop(3).type = '.',     prop(3).subs = 'dom';
%   prop(4).type = '()',    prop(4).subs = {[2]};
%   val = b;
% The returned object PDE will be the same as the input object, only with 
% PDE.x{1}.dom(2) = b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 10/17/2022
%

% At this point, we just use the built-in subsasgn function. Additional
% features such as '()' type subsasgn may be introduced in later updates.
if strcmp(prop(1).type,'()') || strcmp(prop(1).type,'{}')
    error(['''()'' and ''{}'' type subsasgn are currently not supported for "pde_struct" class objects.'])
end
PDE = builtin('subsasgn',PDE,prop,val);

% If any property of the PDE has been adjusted, the PDE is no longer
% initialized (in case it was before).
if PDE.is_initialized && ~strcmp(prop(1).subs,'is_initialized')
   PDE.is_initialized = false;
end

end
function val = subsref(PIE,prop)
% SUBSREF allows properties and elements of the "pie_struct" object PIE to
% be retrieved. The function is automatically called when calling
% e.g. PIE(i,j) or PIE.obj for a "pie_struct" PIE.
%
% INPUT
% - PIE:    A pie_struct object of which to retrieve a property/element.
% - prop:   A 1xq struct specifying which property/element of the input PIE
%           structure to retrieve.
%
% OUTPUT
% - val:    The value of the field of "PIE" specified by "prop".
%
% EXAMPLE: calling PIE.T.dim_in(2), we have
%   prop(1).type = '.',     prop(1).subs = 'T';
%   prop(2).type = '.',     prop(2).subs = 'dim_in';
%   prop(3).type = '()',    prop(3).subs = {[2]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsref
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
% Initial coding DJ - 10/20/2022
% DJ, 10/29/2026: Remove renaming of Bw, ..., Dyu to B1, ..., D22 (both are
%                   now fields);
%

% At this point, we just use the built-in subsref function. Additional
% features such as '()' type subsasgn may be introduced in later updates.
if strcmp(prop(1).type,'()') || strcmp(prop(1).type,'{}')
    error(['''()'' and ''{}'' type subsref are currently not supported for "pie_struct" class objects.'])
else
    switch prop(1).subs
        case 'T1'
            prop(1).subs = 'Tw';
        case 'T2'
            prop(1).subs = 'Tu';
    end
    val = builtin('subsref',PIE,prop);
end

end
function val = subsref(PIE,prop)
% SUBSREF allows properties and elements of the "pie_struct" object PIE to
% be retrieved. The function is automatically called when calling
% e.g. PIE(i,j) or PIE.obj for a "pie_struct" PIE.
%
% INPUTS:
% - PIE:    A pie_struct object of which to adjust a property/element.
% - prop:   A 1xq struct specifying which property/element of the input PIE
%           structure to retrieve.
%
% OUTPUTS:
% - val:    The value of the field of "PIE" specified by "prop".
%
% EXAMPLE: calling PIE.T{1}.dom(2), we have
%   PIE_in = PIE;
%   prop(1).type = '.',     prop(1).subs = 'T';
%   prop(2).type = '{}',    prop(2).subs = {[1]};
%   prop(3).type = '.',     prop(3).subs = 'dom';
%   prop(4).type = '()',    prop(4).subs = {[2]};

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
% Initial coding DJ - 10/20/2022
%

% At this point, we just use the built-in subsasgn function. Additional
% features such as '()' type subsasgn may be introduced in later updates.
if strcmp(prop(1).type,'()') || strcmp(prop(1).type,'{}')
    error(['''()'' and ''{}'' type subsref are currently not supported for "pie_struct" class objects.'])
else
    switch prop(1).subs
        case 'T1'
            prop(1).subs = 'Tw';
        case 'T2'
            prop(1).subs = 'Tu';
        case 'Bw'
            prop(1).subs = 'B1';
        case 'Bu'
            prop(1).subs = 'B2';
        case 'Cz'
            prop(1).subs = 'C1';
        case 'Cy'
            prop(1).subs = 'C2';
        case 'Dzw'
            prop(1).subs = 'D11';
        case 'Dzu'
            prop(1).subs = 'D12';
        case 'Dyw'
            prop(1).subs = 'D21';
        case 'Dyu'
            prop(1).subs = 'D22';
    end
%     if strcmp(prop(1).subs,'dom')
%         disp('WARNING: Adjusting the field "dom" of the PIE will not update the domain of the involved PI operators.')
%     end
    val = builtin('subsref',PIE,prop);
end

end
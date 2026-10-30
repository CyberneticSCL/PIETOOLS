function logval=eq(objA, objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the fields of a "pie_struct" object 'PIE', with name 'name'.
%
% INPUT
% - objA, objB:    pie_struct class objects defining PIE.
%
% OUTPUT
% - logval: 0 is objA~=objB, 1 otherwise
%
% NOTES
% This is a "disp" function, not a "display" function. It only shows which
% fields the input PIE pertains, ordered in a particular manner, without
% showing e.g. the structure of the associated system or the specific
% variables, domain or parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS
%
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


PI_list = {'T', 'A', 'Tw', 'Tu', 'B1', 'B2', 'C1','C2', 'D11', 'D12','D21', 'D22', 'dom'};

logval = 0;

if isa(objA,'pie_struct')&&isa(objB,'pie_struct')
    for i=PI_list
        if ~(objA.(i{:})==objB.(i{:}))
            logval=1;
        end
    end
else
    logval = 1;
end

logval = ~(logval);
end
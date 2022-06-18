function out = set(obj,prop,val)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = set(obj,prop,val) takes in one opvar2d object obj, and returns a
% opvar2d object out with out.prop = val
%
% Version 1.0
% Date: 07/06/21
% 
% INPUT
% obj:  opvar2d class object
% prop: a field name for the object
% val:  a value to be assigned to this field
%
% OUTPUT
% out:  opvar2d object with desired field value
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 07_06_2021
% Adjusted to convert to dopvar if dpvar is added, DJ - 12/30/2021

out = obj;
out.(prop) = val;

if isa(val,'cell')
    for j=1:length(val(:))
        if isa(val{j},'dpvar')
            out = dopvar2d(out);
        end
    end
elseif isa(val,'dpvar')
    out = dopvar2d(out);
end

% switch prop
%     case 'matdim'
%         obj.matdim = val;
%     case 'degmat'
%         obj.degmat = val;
%     case 'C'
%         obj.C = val;
%     case 'varname'
%         obj.varname = val;
%     case 'dvarname'
%         obj.dvarname = val;
%     case 'chkval'
%         error("chkval property is a dependent parameter and cannot be modified");
% end





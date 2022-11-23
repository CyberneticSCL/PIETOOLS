function obj = set(obj,prop,val)
% This function sets the property 'prop' of a dopvar object to 'val'.
% Inputs:
% prop: property to be changes, 'P', 'Q1', 'Q2','R0','R1','R2','var1','var2','I','dim'
% val: value of the property
% 
% Outputs:
% obj: dopvar object with updated property

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - set
%
% Copyright (C)2021  M. Peet, S. Shivakumar
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

if nargin~=2
    error('Incorrect number of inputs to set properties.');
end
switch prop
    case 'P'
        obj.P = val;
    case 'Q1'
        obj.Q1 = val;
    case 'Q2'
        obj.Q2 = val;
    case 'R0'
        obj.R.R0 = val;
    case 'R1'
        obj.R.R1 = val;
    case 'R2'
        obj.R.R2 = val;
    case 'var1'
        obj.var1 = val;
    case 'var2'
        obj.var2 = val;
    case 'I'
        obj.I = val;
    case 'dim'
        obj.dim = val;
    case 'dimdependent'
        error('dimdependent is a dependent property and cannot be modified');
end
end
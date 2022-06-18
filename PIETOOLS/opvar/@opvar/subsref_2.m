function Psop=subsref(Pbop,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psop=op_slice(Pbop,indr,indc) transposes an operator P: R^p x L2^q to R^m x L2^n
% Date: 7/9/2021
% Version: 2.0
% 
% INPUT
% Pbop: opvar class object to slice
% ref: subreferencing cell object
% 
% OUTPUT
% Psop: slice of the opvar object
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - ctranspose
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019
%

switch ref(1).type
    case '.' % accessing a property
        Psop = getprop(Pbop,ref);
    case '()' % accessing a slice of the opvar
        if length(ref(1).subs)==2
            indr = ref(1).subs{1};
            indc = ref(1).subs{2};
        elseif length(ref(1).subs)==1 % linear index
            sz = getprop(Pbop,struct('type','.','subs','dim'));           
            [indr,indc] = ind2sub(sz,ref(1).subs{1});
        else
            error("Incorrect subsreference. Opvar objects can be sliced only using 2 indices or 1 linear index.");
        end
        Psop = op_slice(Pbop,indr,indc);
end

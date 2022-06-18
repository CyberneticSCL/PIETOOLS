function Psop=subsref(Pbop,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psop = subsref(Pbop,ref) extracts components or certain rows or columns
% of a dopvar2d object
%
% Version 1.0
% Date: 07/12/21
% 
% INPUT
% Pbop: opvar class object to slice
% ref:  specification of the components to extract
% 
% OUTPUT
% Psop: slice of the dopvar object
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - subsref
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07_12_2021

switch ref(1).type
    case '.'
        Psop = getprop(Pbop,ref);
    case '()'
        dim = sum(Pbop.dim);
        if length(ref(1).subs)==1  % use linear indexing
            if strcmp(ref(1).subs{1},':')
                indr = 1:dim(1);
                indc = 1:dim(2);
            else
                [nr,nc] = size(Pbop);
                [indr,indc] = ind2sub([nr,nc],ref(1).subs{1});
            end
        else
            if strcmp(ref(1).subs{1},':')
                indr = 1:dim(1);
            else
                indr = ref(1).subs{1};
            end
            if strcmp(ref(1).subs{2},':')
                indc = 1:dim(2);
            else
                indc = ref(1).subs{2};
            end
        end        
        Psop = op_slice(Pbop,indr,indc);
end

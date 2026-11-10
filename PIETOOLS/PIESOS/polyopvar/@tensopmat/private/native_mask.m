function is_nat = native_mask(ops)
% IS_NAT = NATIVE_MASK(OPS) returns a logical array of the size of the cell
% OPS, with element k true if OPS{k} is a matrix of type 'double',
% 'polynomial' or 'dpvar', see `is_native_op`. This applies the test of
% `is_native_op` to all blocks of a 'tensopmat' at once, using the built-in
% forms of cellfun.
%
% INPUTS
% - ops:    cell array of blocks of a 'tensopmat' object;
%
% OUTPUTS
% - is_nat: logical array of the same size as ops, with element k true if
%           ops{k} is a nonempty matrix of type 'double', or a matrix of
%           type 'polynomial' or 'dpvar';
%
% NOTES
% - Emptiness is only tested for blocks of type 'double'.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - native_mask
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 10/15/2026: Initial coding;

is_nat = cellfun('isclass',ops,'polynomial') | cellfun('isclass',ops,'dpvar') | ...
            (cellfun('isclass',ops,'double') & ~cellfun('isempty',ops));

end

function logval = is_native_op(op)
% LOGVAL = IS_NATIVE_OP(OP) checks whether the block OP of a 'tensopmat'
% object is a nonempty matrix of type 'double', 'polynomial', or 'dpvar'.
%
% Such a block is a multiplier from R, mapping a vector x in R^n to the
% function M*x in L2^m[S], for S the spatial variables appearing in M, see
% `tensopmat`.
%
% INPUTS
% - op:     block of a 'tensopmat' object, e.g. an m x n matrix of type
%           'double', 'polynomial' or 'dpvar', a 'tensopvar' or 'intop'
%           object, or [];
%
% OUTPUTS
% - logval: boolean, true if op is a nonempty matrix of type 'double',
%           'polynomial' or 'dpvar', and false otherwise;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - is_native_op
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

logval = (isa(op,'double') || isa(op,'polynomial') || isa(op,'dpvar')) && ~isempty(op);

end

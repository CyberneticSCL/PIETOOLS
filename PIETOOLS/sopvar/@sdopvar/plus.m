function C = plus(A,B)
% C = PLUS(A,B) returns the sum A + B of two 'sdopvar' objects (or of an
% 'sdopvar' and a fixed 'sopvar' object), affine in the decision variables.
%
% INPUTS
% - A, B:   'sdopvar' objects of the same size, mapping between the same
%           spaces, or one of them a fixed 'sopvar' (promoted onto the
%           decision variables of the other), or a 'dpvar' or numeric
%           scalar (scalar*I) or 'dpvar' matrix (multiplier). Text is
%           refused;
%
% OUTPUTS
% - C:      'sdopvar' object representing A + B, on the merged decision-
%           variable and monomial bases, computed as a batch sum of two
%           operands;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding DB - 08/10/2026

% dpvar summand: scalar*I or a matrix multiplier, before
% the call to 'plus_batch' below.
% A numeric scalar is the same as a constant dpvar scalar: scalar*I.
if isa(A,'dpvar') || isa(B,'dpvar') || ...
        (isnumeric(A) && isscalar(A)) || (isnumeric(B) && isscalar(B))
    C = dpvar_op_opvar('plus',A,B);
    return
end
% The sum is the batch sum of two operands, which checks the operands,
% synchronizes their bases, and promotes a fixed 'sopvar' operand onto the
% list of its first 'sdopvar' operand; Zd_p below only sets the 'shared_Zd'
% flag. Text is an option of 'plus_batch', never a summand: refuse it here.
if ischar(A) || isstring(A) || ischar(B) || isstring(B)
    error("A summand of an 'sdopvar' must be an operator, not text.")
end
% A nonempty Zd_p is stored unchanged by 'sopvar2sdopvar', so both lists
% are then the same array. An empty one is replaced by cell(0,1), and the
% comparison is O(1) there anyway, so it is left to sync_basis.
shared_Zd = false;
if isa(A,'sopvar') || isa(B,'sopvar')
    if      isa(B,'sdopvar'),   Zd_p = B.Zd;
    elseif  isa(A,'sdopvar'),   Zd_p = A.Zd;
    else,                       Zd_p = cell(0,1);
    end
    shared_Zd = ~isempty(Zd_p);
end
if shared_Zd
    C = plus_batch(A,B,'shared_Zd');
else
    C = plus_batch(A,B);
end
end

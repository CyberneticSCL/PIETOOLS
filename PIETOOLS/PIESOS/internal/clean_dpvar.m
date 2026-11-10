function [op,is_zero] = clean_dpvar(op,ztol)
% [OP,IS_ZERO] = CLEAN_DPVAR(OP,ZTOL) removes the coefficients of absolute
% value at most ZTOL from the matrix OP of type 'double', 'polynomial', or
% 'dpvar', and returns whether all coefficients are zero.
%
% INPUTS
% - op:     m x n matrix of type 'double', 'polynomial', or 'dpvar';
% - ztol:   scalar value specifying the tolerance: coefficients of
%           absolute value at most ztol are removed (for a 'dpvar', the
%           coefficients are the entries of op.C);
%
% OUTPUTS
% - op:     m x n matrix of the same type as the input, with the
%           coefficients of absolute value at most ztol removed;
% - is_zero: logical, true if all coefficients of the cleaned op are zero;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - clean_dpvar
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

if isa(op,'double')
    op(abs(op)<=ztol) = 0;
    is_zero = ~any(op(:));
    return
end
if isa(op,'polynomial')
    op = cleanpoly(op,ztol);
    is_zero = ~any(op.coefficient(:));
    return
end
C = op.C;
[ii,jj,vv] = find(C);
is_kept = abs(vv)>ztol;
is_zero = ~any(is_kept);
if ~all(is_kept)
    op.C = sparse(ii(is_kept),jj(is_kept),vv(is_kept),size(C,1),size(C,2));
end

end

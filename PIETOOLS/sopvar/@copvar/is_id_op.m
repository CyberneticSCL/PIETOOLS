function logval = is_id_op(Pop,ztol)
% LOGVAL = IS_ID_OP(POP,ZTOL) checks whether the block matrix POP of PI
% operators is an identity operator.
%
% INPUTS
% - Pop:    p x d 'copvar' object;
% - ztol:   (optional) scalar tolerance below which coefficients are taken
%           to be zero, passed to the 'is_id_op' of each diagonal block.
%           Defaults to 1e-12;
%
% OUTPUTS
% - logval: boolean specifying whether (true) or not (false) Pop is an
%           identity operator on the function space on which it is defined;
%
% NOTES
% - Empty blocks Pop.C{i,j} represent zero operators. Pop is identity iff
%   it is square, maps from and to the same space (Pop.dim_in=Pop.dim_out
%   and Pop.space_in=Pop.space_out), every diagonal block Pop.C{i,i} is an
%   identity operator, and every off-diagonal block is empty or zero up to
%   ztol.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - is_id_op
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
% DJ, 11/06/2026: Initial coding;

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Assume the operator is not identity until we have verified
logval = false;

% Only a square block structure on the same space can be identity
[p,d] = size(Pop.C);
if p==0 || p~=d || ~isequal(Pop.dim_in,Pop.dim_out) || ~isequal(Pop.space_in,Pop.space_out)
    return
end

% Check that the diagonal blocks are identity, and the others are zero
for i=1:p
    for j=1:d
        blk = Pop.C{i,j};
        if i==j
            % An empty diagonal block is a zero operator
            if isempty(blk) || ~is_id_op(blk,ztol)
                return
            end
        elseif ~isempty(blk) && any(cellfun(@(c) ~isempty(c) && any(abs(nonzeros(c))>ztol),blk.params))
            return
        end
    end
end

% At this point, we can conclude that Pop is an identity operator
logval = true;

end

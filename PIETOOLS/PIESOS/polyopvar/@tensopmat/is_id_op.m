function logval = is_id_op(Cop,ztol)
% LOGVAL = IS_ID_OP(COP,ZTOL) checks whether the block matrix COP of
% tensor-PI operators is an identity operator.
%
% INPUTS
% - Cop:    p x d 'tensopmat' object;
% - ztol:   (optional) scalar tolerance below which coefficients are taken
%           to be zero, passed to 'clean_op' and to the 'is_id_op' of each
%           block. Defaults to 1e-12;
%
% OUTPUTS
% - logval: boolean specifying whether (true) or not (false) Cop is an
%           identity operator on the function space on which it is defined;
%
% NOTES
% - Empty blocks Cop.ops{i,j} represent zero operators. Cop is identity iff
%   it is square (dims.in = dims.out), every diagonal block Cop.ops{i,i} is
%   an identity operator, and every off-diagonal block is empty or zero up
%   to ztol.
% - Blocks are checked individually, not through Cop.vars/Cop.depmat, since
%   a zero off-diagonal block may still depend on variables that a diagonal
%   block does not.
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
% DJ, 10/15/2026: Initial coding;

% Check that the operator is properly specified, as in clean_op
if ~Cop.chkval && ~isvalid(Cop)
    error("Tensopmat object is not properly declared.")
end

% Set default tolerance, matching clean_op and the blocks' is_id_op
if nargin==1
    ztol = 1e-12;
end

% Assume the operator is not identity until we have verified
logval = false;

% Only a square block structure with matching block sizes can be identity
ops = Cop.ops;
[p,d] = size(ops);
if p==0 || p~=d || ~isequal(Cop.dims.in,Cop.dims.out)
    return
end

% Check that the diagonal blocks are identity, and the others are zero
for i=1:p
    for j=1:d
        op_ij = ops{i,j};
        if i==j
            % An empty diagonal block is a zero operator
            if isempty(op_ij)
                return
            elseif is_native_op(op_ij)
                % Matrix acting on R: identity if it is the identity matrix, and the row maps to R
                [~,is_zero_ij] = clean_dpvar(op_ij-eye(size(op_ij,1)),ztol);
                % (a row of such blocks only has the variables appearing in them, so is in R if they are not)
                is_R_row = ~any(Cop.depmat.out(i,:)) || all(cellfun(@(x) isempty(x) || is_native_op(x),ops(i,:)));
                if ~is_zero_ij || ~is_R_row
                    return
                end
            elseif ~is_id_op(op_ij,ztol)
                return
            end
        elseif ~isempty(op_ij)
            % Off-diagonal block must be zero up to ztol
            if is_native_op(op_ij)
                [~,is_zero_ij] = clean_dpvar(op_ij,ztol);
            else
                [~,is_zero_ij] = clean_op(op_ij,ztol);
            end
            if ~is_zero_ij
                return
            end
        end
    end
end

% At this point, we can conclude that Cop is an identity operator
logval = true;

end

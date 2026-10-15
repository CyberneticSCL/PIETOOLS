function Y = apply_op(Cop,X)    
% Y = APPLY_OP(COP,X) applies the 'tensopmat' object COP to the 1 x d cell
% X of inputs, one for each column-block, and returns the vertical
% concatenation Y of the outputs of the p row-blocks.
%
% INPUTS
% - Cop:    'tensopmat' object with p x d blocks;
% - X:      1 x d cell, with X{j} the input to which the blocks in column
%           j are applied. For a block of type 'double', 'polynomial' or
%           'dpvar', X{j} is a vector or matrix in R (or a 1 x 1 cell
%           containing it), which the block multiplies. For an operator
%           block, X{j} is the input passed to its 'apply_op' (to
%           'apply_sopvar' for a 'sopvar');
%
% OUTPUTS
% - Y:      polynomial array representing the vertical concatenation of
%           the outputs of the p row-blocks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_op (tensopmat)
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

if ~isa(Cop,'tensopmat')
    error("First argument must be of type 'tensopmat'.")
end
if ~Cop.chkval && ~isvalid(Cop)
    error("Tensopmat object is not properly declared.")
end
[p,d] = size(Cop.ops);
if ~isa(X,'cell') || numel(X)~=d
    error("Input must be specified as 1 x d cell, matching the number of column-blocks of Cop.")
end

Yblocks = cell(p,1);
for i=1:p
    Yi = 0;
    for j=1:d
        if isempty(Cop.ops{i,j})
            continue
        end
        % 'sopvar' has no apply_op; the others dispatch by class    
        if is_native_op(Cop.ops{i,j}) 
            % Matrix acting on R: the input X{j} is a vector (or matrix) in
            % R, or a 1x1 cell with this, the form of the inputs of the 
            % factors of an operator (see apply_poly) 
            xj = X{j}; 
            if iscell(xj) && isscalar(xj)
                xj = xj{1};
            end 
            Yij = Cop.ops{i,j}*xj;
        elseif isa(Cop.ops{i,j},'sopvar') 
            Yij = apply_sopvar(Cop.ops{i,j},X{j});
        else
            Yij = apply_op(Cop.ops{i,j},X{j});
        end
        Yi = Yi + Yij;
    end
    Yblocks{i} = Yi;
end
Y = vertcat(Yblocks{:});

end

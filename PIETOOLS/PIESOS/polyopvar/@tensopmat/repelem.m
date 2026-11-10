function B = repelem(A,r1,r2)
% B = REPELEM(A,R1,R2) returns a modified version of the operator A with
% each row repeated R1 times and each column repeated R2 times.
%
% INPUTS
% - A:      m x n 'tensopmat' object;
% - r1:     scalar integer specifying the number of times each row of A is
%           to be repeated. If r2 is not specified, r1 is the 1 x 2 array
%           [R1,R2];
% - r2:     (optional) scalar integer specifying the number of times each
%           column of A is to be repeated;
%
% OUTPUTS
% - B:      m*r1 x n*r2 'tensopmat' object representing
%           kron(A,ones(r1,r2));
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - repelem
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 10/15/2026: Initial coding;

% Check that the tensor-PI operator is appropriate
if ~isa(A,'tensopmat')
    error("First input must be of type 'tensopmat'.")
end
if ~A.chkval && ~isvalid(A)
    error("The 'tensopmat' object is not properly declared.")
end

% Check that the number of repetitions is properly specified
if nargin<2
    error("Not enough input arguments.")
elseif nargin==2
    % Specify repetition factors as row vector
    repdim = r1;
elseif nargin==3
    repdim = [r1,r2];
else
    error("Repetition is only supported along row and column directions.")
end
if any(size(repdim)~=[1,2])
    error("Repetition factors must be specified as 1x2 vector of integers, or as 2 integer scalars.")
elseif any(~isnumeric(repdim)) || any(repdim~=round(abs(repdim)))
    error("Repetition factors must be specified as nonnegative integers.")
end

% Repeat each row and column of the matrix defined by A, so B = kron(A,ones(r1,r2)):
% repeat the rows and columns of each block, which keeps its variables, so only the
% dimensions of the blocks change
B = A;
if any(repdim==0)
    % Rows or columns are removed, so no block is left to repeat: remove the blocks
    B.ops = A.ops(repelem(1:size(A.ops,1),repdim(1)),repelem(1:size(A.ops,2),repdim(2)));
    B.dims.in = repelem(A.dims.in,repdim(2),1);
    B.dims.out = repelem(A.dims.out,repdim(1),1);
    B.depmat.in = repelem(A.depmat.in,repdim(2),1);
    B.depmat.out = repelem(A.depmat.out,repdim(1),1);
else
    for k=1:numel(A.ops)
        op = A.ops{k};
        if isempty(op)
            continue        % zero operator
        elseif isa(op,'tensopvar')
            B.ops{k} = repelem(op,repdim(1),repdim(2));
        else
            % 'intop', or matrix of type 'double', 'polynomial' or 'dpvar' (builtin repelem fails for
            % a 'polynomial'): index with repeated indices
            B.ops{k} = op(repelem(1:size(op,1),repdim(1)),repelem(1:size(op,2),repdim(2)));
        end
    end
    B.dims.in = repdim(2)*A.dims.in;
    B.dims.out = repdim(1)*A.dims.out;
end
B.chkval = true;

end

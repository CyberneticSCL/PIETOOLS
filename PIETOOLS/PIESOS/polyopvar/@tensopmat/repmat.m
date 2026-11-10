function B = repmat(A,r1,r2)
% B = REPMAT(A,R1,R2) returns the concatenation of R1 x R2 copies of A
% along the row and column directions.
%
% INPUTS
% - A:      m x n 'tensopmat' object;
% - r1:     scalar integer specifying the number of times A is to be
%           repeated along the row direction. If r2 is not specified, r1
%           is the 1 x 2 array [R1,R2];
% - r2:     (optional) scalar integer specifying the number of times A is
%           to be repeated along the column direction;
%
% OUTPUTS
% - B:      m*r1 x n*r2 'tensopmat' object representing
%           kron(ones(r1,r2),A);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - repmat
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

% Simply repeat the blocks of operators defining A
B = A;
B.ops = repmat(A.ops,repdim);
B.dims.in = repmat(A.dims.in,repdim(2),1);
B.dims.out = repmat(A.dims.out,repdim(1),1);
B.depmat.in = repmat(A.depmat.in,repdim(2),1);
B.depmat.out = repmat(A.depmat.out,repdim(1),1);
B.chkval = true;

end
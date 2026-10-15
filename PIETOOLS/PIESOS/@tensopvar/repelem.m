function B = repelem(A,r1,r2)
% B = REPELEM(A,R1,R2) returns a modified version of the operator A with
% each row repeated R1 times and each column repeated R2 times.
%
% INPUTS
% - A:      m x n 'tensopvar' object;
% - r1:     scalar integer specifying the number of times each row of A is
%           to be repeated. If r2 is not specified, r1 is the 1 x 2 array
%           [R1,R2];
% - r2:     (optional) scalar integer specifying the number of times each
%           column of A is to be repeated;
%
% OUTPUTS
% - B:      m*r1 x n*r2 'tensopvar' object representing
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
if ~isa(A,'tensopvar')
    error("First input must be of type 'tensopvar'.")
end
if ~A.chkval
    A = initialize(A);
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

% For tensor product, repmat must be applied to last factor in each term
%   (R1 o ... o Rn) o Im = R1 o ... o Rn-1 o (Rn o Im)
B = A;
[ntrms,nfctrs] = size(A.ops);
repdim_arr = ones(ntrms*nfctrs,2);
repdim_arr(end-ntrms+1:end,:) = repmat(repdim,ntrms,1);
% For pointwise product, repelem must be applied to each factor in each term
%   (R1.* ... .*Rn) o Im = (R1 o Im).* ... .*(Rn o Im)
if A.otype(1)==1
    repdim_arr(1:end-ntrms,1) = repdim(1);
end
if A.otype(2)==1
    repdim_arr(1:end-ntrms,2) = repdim(2);
end
param_nums = find(any(repdim_arr>1,2));
for k=reshape(param_nums,1,[])
    B.ops{k} = repelem(A.ops{k},repdim_arr(k,:));
end

% Update the dimensions of the different factors
if A.otype(1)==1
    B.dims(:,1) = repdim(1)*A.dims(:,1);
else
    B.dims(end,1) = repdim(1)*A.dims(end,1);
end
if A.otype(2)==1
    B.dims(:,2) = repdim(2)*A.dims(:,2);
else
    B.dims(end,2) = repdim(2)*A.dims(end,2);
end
B.chkval = true;

end
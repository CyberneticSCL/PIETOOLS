function B = repelem(A,r1,r2)
% B = REPMAT(A,R1,R2) returns returns a modified version of the operator A
% with each row repeated R1 times and each column repeated R2 times
%
% INPUTS
% - A:  m x n 'sopvar' object;
% - r1: scalar integer specifying the number of times each row of A is to 
%       be repeated;
% - r2: scalar integer specifying the number of times each column of A is 
%       to be repeated;
%
% OUTPUTS
% - B:  m*r1 x n*r2 'sopvar' object representing kron(A,ones(r1,r2));
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 06/12/2026: Initial coding;

% Check the inputs
if nargin<2
    error("Not enough input arguments.")
elseif nargin==2
    % Specify repetition factors as row vector
    repdim = r1;
elseif nargin==3
    repdim = [r1,r2];
else
    error("Replication is only supported along row and column directions.")
end
if any(size(repdim)~=[1,2])
    error("Replication factors must be specified as 1x2 vector of integers, or as 2 integer scalars.")
elseif any(~isnumeric(repdim)) || any(repdim~=round(abs(repdim)))
    error("Replication factors must be specified as nonnegative integers.")
end

% By definition of the sopvar class, repmat of the operator reduces to just
% repmat of the parameters.
params = A.params;
nZL_arr = cellfun(@(a) numel(a),A.ZL);      nZL = prod([nZL_arr,1]);
nZR_arr = cellfun(@(a) numel(a),A.ZR);      nZR = prod([nZR_arr,1]);
for k=1:numel(params)
    % Extract row indices, column indices, and values of coefficients
    [idcs1,idcs2,vals] = find(params{k});
    % Split row indices into monomial index and operator row index
    ridcs = ceil(idcs1(:)./nZL);       zL_idcs = idcs1(:) - (ridcs-1)*nZL;
    % Perform repetition along the row direction
    ridcs = (ridcs-1)*r1 + (1:r1);
    zL_idcs = repmat(zL_idcs(:),1,r1);
    idcs1 = (ridcs(:)-1)*nZL + zL_idcs(:);
    idcs2 = repmat(idcs2(:),r1,1);
    vals = repmat(vals(:),r1,1);
    % Split column indices into monomial index and operator column index
    cidcs = ceil(idcs2(:)./nZR);       zR_idcs = idcs2(:) - (cidcs-1)*nZR;
    % Perform repetition along the column direction
    cidcs = (cidcs-1)*r2 + (1:r2);
    zR_idcs = repmat(zR_idcs(:),1,r2);
    idcs2 = (cidcs(:)-1)*nZR + zR_idcs(:);
    idcs1 = repmat(idcs1(:),r2,1);
    vals = repmat(vals(:),r2,1);
    % Build the updated coefficient matrix
    params{k} = sparse(idcs1,idcs2,vals,size(params{k},1)*r1,size(params{k},2)*r2);
end
% Construct the updated operator
matdim = repdim.*size(A);
B = sopvar(params,A.vars,A.ZR,A.ZL,A.dom,matdim);

end
function C = opower(A,d,type)
% C = OPOWER(A,D,TYPE) returns the D-th order tensor power of the
% tensor-PI operator A, i.e. the product A o ... o A of D copies of A, so
% that
%   (C*(Zd(x1) o ... o Zd(xD)))(s) = (A*Zd(x1))(s) .* ... .* (A*Zd(xD))(s)
% for a pointwise product, and similarly for a tensor product.
%
% INPUTS
% - A:      m x n 'tensopvar' object (or 'sopvar');
% - d:      nonnegative integer, the order of the product;
% - type:   (optional) 1 x 2 array of values 1 or 2 indicating what type of
%           product is taken along the row and column dimensions, as in
%           'otimes'. Defaults to [1,2]; a scalar specifies type(1), with
%           type(2)=2;
%
% OUTPUTS
% - C:      'tensopvar' object representing the D-th tensor power of A,
%           the same as D-1 successive calls to otimes(A,.,type);
%
% NOTES
% - If d=0, C=1 is returned, and if d=1, C=A. If A has size 0 x 0,
%   C = zeros(0,0) is returned.
% - If A has a single factor, or factors of the same type as the product,
%   the factors of d copies of each term are concatenated, so that C has
%   size(A.ops,1)^d terms. Otherwise, C has a single term, with A nested
%   as each of its d factors.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opower
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

if nargin<2
    error("Not enough inputs")
elseif ~isscalar(d) || d<0 || d~=round(d)
    error("Order of the product must be a nonnegative integer.")
end
if nargin<=2
    type = [1,2];
elseif isscalar(type)
    type = [type,2];
end
if d==0
    C = 1;
    return
end
if ~isa(A,'tensopvar')
    if isa(A,'sopvar')
        A = sopvar2tensopvar(A);
    else
        error("Input must be of type 'tensopvar' or 'sopvar'.")
    end
end
if all(size(A)==0)
    C = zeros(0,0);
    return
end
if ~A.chkval
    A = initialize(A);
end
if d==1
    C = A;
    return
end

% A single factor is of any type, as in otimes
if size(A.ops,2)==1
    A.otype = type;
end

C = A;
if all(type==A.otype)
    % Same type as the factors of A: concatenate d copies of each term, the
    % first copy varying fastest, as in otimes
    K = size(A.ops,1);
    idcs = cell(1,d);
    [idcs{:}] = ndgrid(1:K);
    C.ops = cell(K^d,size(A.ops,2)*d);
    nc = size(A.ops,2);
    for j=1:d
        C.ops(:,(j-1)*nc+(1:nc)) = A.ops(idcs{j}(:),:);
    end
    C.depmat.out = repmat(A.depmat.out,d,1);
    C.depmat.in = repmat(A.depmat.in,d,1);
    C.dims = repmat(A.dims,d,1);
else
    % Different type: nest A as each of d factors, with the dependencies of
    % A as a whole in each row of depmat
    Adeps = A.deps;
    C.ops = repmat({A},1,d);
    C.depmat.out = repmat(Adeps.out,d,1);
    C.depmat.in = repmat(Adeps.in,d,1);
    C.dims = repmat(A.matdim,d,1);
    C.otype = type;
    C.chkval = true;
end

end

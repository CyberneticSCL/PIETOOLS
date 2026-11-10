function C = opower(A,d,type)
% C = OPOWER(A,D,TYPE) returns the D-th order tensor power of the block
% matrix A of tensor-PI operators, i.e. the product A o ... o A of D copies
% of A, so that
%   (C*(Zd(x1) o ... o Zd(xD)))(s) = (A*Zd(x1))(s) .* ... .* (A*Zd(xD))(s)
% for a pointwise product, and similarly for a tensor product.
%
% INPUTS
% - A:      m x n 'tensopmat' object (or 'tensopvar', 'intop', 'sopvar');
% - d:      nonnegative integer, the order of the product;
% - type:   (optional) 1 x 2 array of values 1 or 2 indicating what type of
%           product is taken along the row and column dimensions, as in
%           'otimes'. Defaults to [1,2]; a scalar specifies type(1), with
%           type(2)=2;
%
% OUTPUTS
% - C:      'tensopmat' object representing the D-th tensor power of A,
%           with the same blocks as D-1 successive calls to
%           otimes(A,.,type) (up to the nesting of the factors in each
%           block);
%
% NOTES
% - If d=0, C=1 is returned, and if d=1, C=A.
% - If A has a single block that is a 'tensopvar', the power of this block
%   is taken using 'opower' for 'tensopvar' objects. Otherwise, the blocks
%   of different copies of A differ, and the power is computed by repeated
%   squaring, using ceil(log2(d)) products of the powers of A, instead of
%   d-1 products with A.
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
if ~isa(A,'tensopmat')
    if isa(A,'sopvar') || isa(A,'tensopvar') || isa(A,'intop')
        A = tensopmat(A);
    else
        error("Input must be of type 'tensopmat' or 'tensopvar'.")
    end
end
if ~A.chkval && ~isvalid(A)
    error("Operation not supported: The 'tensopmat' object is not properly declared.")
else
    A.chkval = true;
end
if d==1
    C = A;
    return
end

if isscalar(A.ops) && isa(A.ops{1},'tensopvar')
    % A single block: the factors of the product are all copies of it
    C = tensopmat(opower(A.ops{1},d,type));
    return
end

% Repeated squaring: P = A^(2^k), and C the product of those P for which
% the k-th bit of d is set
hasC = false;
P = A;
while true
    if mod(d,2)==1
        if hasC
            C = otimes(C,P,type);
        else
            C = P;
            hasC = true;
        end
    end
    d = floor(d/2);
    if d==0
        break
    end
    P = otimes(P,P,type);
end

end

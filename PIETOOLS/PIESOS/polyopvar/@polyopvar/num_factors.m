function [nfac,nops,is_empty] = num_factors(f)
% [NFAC,NOPS,IS_EMPTY] = NUM_FACTORS(F) returns the number of factors,
% i.e. of columns of F.Z, that each of the distributed monomials of F
% takes up.
%
% INPUTS
% - f:      'polyopvar' object defined by nZ distributed monomials
%           f.Z(i,:) = [k1,...,kd], acted on by the operators f.C.ops(:,i);
%
% OUTPUTS
% - nfac:   nZ x 1 array of integers, the number of columns of Z that
%           monomial i takes up: the number of factors in the tensor-PI
%           operators acting on it, but at least the index of the last
%           nonzero element of f.Z(i,:);
% - nops:   nZ x 1 array of integers, the number of factors in the
%           'tensopvar' objects acting on monomial i, counting also the
%           factors of nested 'tensopvar' objects. Zero if no 'tensopvar'
%           acts on monomial i, i.e. if it is acted on by 'intop' objects
%           only (a functional), or by blocks of type 'double',
%           'polynomial' or 'dpvar' (a known function, see 'tensopmat',
%           for which nfac(i)>=1);
% - is_empty: nZ x 1 logical array, true if all operators acting on
%           monomial i are empty (zero), so that it has no factors to
%           count, see factor_key;
%
% NOTES
% - Factor j of the operator acts on state f.Z(i,j), or, if f.Z(i,j)=0
%   (also for j>size(f.Z,2)), is a multiplier from R to some function
%   space, representing a known factor of the monomial. This is only so
%   if nops(i)>nnz(f.Z(i,:)), i.e. if the operator has more factors than
%   the degree of the monomial; otherwise, zeros in f.Z(i,:) are just
%   padding, and initialize moves them to the right, so that
%   nfac(i)=nnz(f.Z(i,:)) unless a factor acts on several states at once.
%   See 'polyopvar'.
% - The number of factors is read from the first term of the operator
%   in the first row of f.C.ops acting on monomial i that is a
%   'tensopvar'; initialize checks that all operators acting on the
%   monomial agree.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - num_factors
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

% Index of the last nonzero element of each row of Z
nZ = size(f.Z,1);
[ridcs,cidcs] = find(f.Z);
nfac = accumarray(ridcs(:),cidcs(:),[nZ,1],@max);

% Number of factors in the operators acting on each monomial
ops = f.C.ops;
is_tens = cellfun('isclass',ops,'tensopvar');
nops = zeros(nZ,1);
for i=find(any(is_tens,1))
    nops(i) = num_factors(ops{find(is_tens(:,i),1),i});     % method of 'tensopvar'
end
is_empty = reshape(all(cellfun('isempty',ops),1),[],1);
nfac = max(nfac,nops);
% A block of type double/polynomial/dpvar acts on R as a single multiplier, the known factor
% of a degree-0 monomial (see `tensopmat`), when no tensopvar acts on the monomial: it takes up one
% factor, so that the factors of a state in a product with it follow its own, see times. nops is not
% changed, as such a monomial is not merged with one acted on by a 'tensopvar' (see factor_key)
idcs_nat = find(~any(is_tens,1) & ~is_empty');
if ~isempty(idcs_nat)
    sub = ops(:,idcs_nat);
    is_nat = cellfun('isclass',sub,'polynomial') | cellfun('isclass',sub,'dpvar') | (cellfun('isclass',sub,'double') & ~cellfun('isempty',sub));
    nfac(idcs_nat(any(is_nat,1))) = max(nfac(idcs_nat(any(is_nat,1))),1);
end

end

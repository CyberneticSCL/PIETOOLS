function C = opower(A,d,otype)
% C = OPOWER(A,D,OTYPE) returns the 'polyopvar' object C representing the
% D-th order power of the distributed polynomial A, i.e. the product
% A o ... o A of D copies of A.
%
% INPUTS
% - A:     m x n 'polyopvar' object;
% - d:     nonnegative integer, the order of the product;
% - otype: (optional) scalar 1 or 2, specifying the type of product, as in
%          'otimes'. If otype=2 (default), the tensor product is taken, so
%          that C is a function of a copy of the variables for each factor,
%              C(s1,...,sD) = A(s1) o ... o A(sD).
%          If otype=1, the elementwise product is taken, so that C is a
%          function of a single copy of the variables,
%              C(s) = A(s).*...*A(s);
%
% OUTPUTS
% - C:     'polyopvar' object representing the D-th power of A;
%
% NOTES
% - If d=0, C=1 is returned, and if d=1, C=A.
% - For a polynomial with output variables and without constant terms, the
%   coefficients of C are the D-th power of those of A, see
%   @tensopmat/opower.m, acting on the monomials x_{i1} o ... o x_{iD},
%   for all combinations of D monomials of A, followed by combine_terms.
%   C has nZ^D monomials, for A with nZ monomials: the products of
%   monomials in a different order, e.g. x1*x2 and x2*x1, are not combined.
% - For a functional (no output variables), monomials of the same degree in
%   each state are combined, so that the number of monomials of the power
%   is much smaller than nZ^D, and the power is computed by repeated
%   squaring, combining terms after each product with otimes. The same
%   holds if A has constant terms (monomials of degree 0), as otimes
%   merges their factors with those of other monomials.
% - A.^D is the elementwise power opower(A,D,1), see @polyopvar/power.m.
% - An adjoint (see 'ctranspose') is not supported.
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
% DJ, 11/07/2026: Reject adjoints, see 'ctranspose';

if nargin<2
    error("Not enough inputs")
elseif ~isscalar(d) || ~isnumeric(d) || d<0 || d~=round(d)
    error("The order of the product, d, must be a nonnegative integer.")
end
if nargin<3
    otype = 2;
elseif ~isscalar(otype) || ~any(otype==[1,2])
    error("The type of product, otype, must be 1 (elementwise product) or 2 (tensor product).")
end
if ~isa(A,'polyopvar')
    error("Input must be of type 'polyopvar'.")
end
if A.is_adjoint
    error("polyopvar:adjoint_not_supported","'opower' is not supported for adjoint 'polyopvar' objects (u'), see 'ctranspose'.")
end
if d==0
    C = 1;
    return
end
if ~A.chkval
    A = initialize(A);
end

% Number of factors of each monomial, and which are constants (degree 0)
[nfac,~] = num_factors(A);
nfac = full(nfac(:));
is_const = nfac==1 & ~any(A.Z,2);

if ~any(is_const) && any(A.depmat.out(:)) && d>1
    % The coefficients are the d-th power of those of A, acting on the
    % monomials x_{i1} o ... o x_{id}, ordered as the columns of the power,
    % the first index running slowest. No column is a constant, so no
    % factors are merged and those of monomial k follow those of 1,...,k-1
    n = size(A.Z,1);
    % Take the tensor power of the coefficients
    C = A;
    C.C = opower(A.C,d,[otype,2]);      % sets the output space of C, see set.C
    C.chkval = true;
    % Take the tensor power of the monomials
    idx0 = (0:n^d-1)';
    shift = zeros(n^d,1);
    [ridcs,cidcs,vals] = deal(cell(d,1));
    for k=1:d
        Ik = mod(floor(idx0/n^(d-k)),n)+1;
        [r,c,v] = find(A.Z(Ik,:));
        ridcs{k} = r(:);    cidcs{k} = c(:) + shift(r(:));    vals{k} = v(:);
        shift = shift + nfac(Ik);
    end
    C.Z = sparse(vertcat(ridcs{:}),vertcat(cidcs{:}),vertcat(vals{:}),n^d,max(shift));
    C = combine_terms(C);
    return
end

% Repeated squaring: P = A^(2^k), and C the product of those P for which
% the k-th bit of d is set
hasC = false;
P = A;
while true
    if mod(d,2)==1
        if hasC
            C = otimes(C,P,otype);
        else
            C = P;
            hasC = true;
        end
    end
    d = floor(d/2);
    if d==0
        break
    end
    P = otimes(P,P,otype);
end

end

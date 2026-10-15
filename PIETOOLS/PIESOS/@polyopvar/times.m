function C = times(A,B)
% C = TIMES(A,B) returns the 'polyopvar' object C representing the
% elementwise product A.*B of two distributed polynomials, or of a
% constant with a distributed polynomial.
%
% INPUTS
% - A:     m x n 'polyopvar' object, or a constant such that A.*B.C is
%          defined, see @tensopmat/times.m;
% - B:     m x n 'polyopvar' object, or a constant such that A.C.*B is
%          defined, see @tensopmat/times.m;
%
% OUTPUTS
% - C:     'polyopvar' object representing A.*B;
%
% NOTES
% - The monomials of C are the products of each monomial of A with each
%   monomial of B, where the factors of B follow those of A, see
%   'polyopvar'.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - times
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

% Make sure the sizes are appropriate for elementwise multiplication
sz_A = size(A);     is_dim1_A = sz_A==1;
sz_B = size(B);     is_dim1_B = sz_B==1;
if any(~is_dim1_A & ~is_dim1_B & sz_A~=sz_B)
    error("Matrix dimensions must match for elementwise multiplication.")
end

% Account for case where one object is not a distributed polynomial
if ~isa(A,'polyopvar')
    % Multiplication of constant with polynomial
    C = B;
    C.C = A.*B.C;
    return
elseif ~isa(B,'polyopvar')
    % Multiplication of polynomial with constant
    C = A;
    C.C = A.C.*B;
    return
end

% Initialize
if ~A.chkval
    A = initialize(A);
end
if ~B.chkval
    B = initialize(B);
end

% Express A and B in terms of the same variables
[A,B] = common_vars(A,B);
C = A;

% Declare the full list of monomials appearing in the product of A and B
[ridcsA,cidcsA,ZvalsA] = find(A.Z);     [nrA,ncA] = size(A.Z);
[ridcsB,cidcsB,ZvalsB] = find(B.Z);     [nrB,ncB] = size(B.Z);
% The factors of B follow all factors of A, including the multipliers from R (zeros in Z), see `polyopvar`
[nfacA,nopsA] = num_factors(A);
dA = repelem(nfacA,nrB,1);
ridcsA = (ridcsA(:)-1)*nrB + (1:nrB);       % repelem(A.Z,nrB,1)
cidcsA = repmat(cidcsA(:),1,nrB);
ZvalsA = repmat(ZvalsA(:),1,nrB);
ridcsB = ridcsB(:) + (0:nrA-1)*nrB;                    % repmat(B.Z,nrA,1);
cidcsB = repmat(cidcsB(:),1,nrA) + reshape(dA(ridcsB),size(ridcsB));
ZvalsB = repmat(ZvalsB(:),1,nrA);
Znew = sparse([ridcsA(:);ridcsB(:)],[cidcsA(:);cidcsB(:)],[ZvalsA(:);ZvalsB(:)],nrA*nrB,max([ncA;nfacA(:)+ncB]));

% Declare the coefficients acting on the full vector of monomials
C.C = otimes(A.C,B.C,[1,2]);
C.dims.out = C.C.dims.out;
C.vars.out = C.C.vars.out;
C.dom.out = C.C.dom.out;
C.depmat.out = C.C.depmat.out;
C.chkval = true;

% otimes merges a constant of A (a single factor, in a monomial of degree 0) with a scalar
% factor of B without variables, e.g. a state not depending on any spatial variable, into
% one factor. B then has no factor of its own to shift: its monomial follows the (zero)
% monomial of the constant directly, as the zeros in Z are padding in that case.
is_const = full(nfacA==1 & ~any(A.Z,2));    % also a function given by a block of type double, see num_factors
if any(is_const)
    [~,nopsB] = num_factors(B);
    for i=find(is_const(:))'
        rws = (i-1)*nrB + (1:nrB);
        Csub = C;       Csub.Z = Znew(rws,:);       Csub.C.ops = C.C.ops(:,rws);
        [~,nopsC] = num_factors(Csub);
        is_merged = nopsC>0 & nopsB>0 & nopsC<1+nopsB;
        Znew(rws(is_merged),:) = 0;
        Znew(rws(is_merged),1:ncB) = B.Z(is_merged,:);
    end
end
C.Z = Znew;

% Finally, combine terms involving the same monomial
C = combine_terms(C);

end
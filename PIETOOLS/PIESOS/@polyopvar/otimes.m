function C = otimes(A,B,otype)
% C = OTIMES(A,B,OTYPE) returns the 'polyopvar' object C representing the
% tensor product of two distributed polynomials, or their elementwise
% product, or the product of a constant with a distributed polynomial.
%
% INPUTS
% - A:     m1 x n1 'polyopvar' object, or a constant such that A.*B.C is
%          defined, see @tensopmat/times.m;
% - B:     m2 x n2 'polyopvar' object, or a constant such that A.C.*B is
%          defined, see @tensopmat/times.m;
% - otype: (optional) scalar 1 or 2, specifying the type of product.
%          If otype=2 (default), the tensor product is taken, so that C is
%          the Kronecker product of A and B, with entries Ai(s1)*Bj(s2) in
%          a copy of the variables for each factor,
%              C(s1,s2) = A(s1) o B(s2).
%          If otype=1, the elementwise product A.*B is taken, requiring
%          that A and B have the same size, or that one of them is scalar,
%          so that C is a function of a single copy of the variables,
%              C(s) = A(s).*B(s).
%          For a constant, the constant must be scalar if otype=2;
%
% OUTPUTS
% - C:     'polyopvar' object representing the product of A and B;
%
% NOTES
% - The monomials of C are the products of each monomial of A with each
%   monomial of B, where the factors of B follow those of A, see
%   'polyopvar'. For otype=2, the output of such a monomial is a function
%   of a copy of the variables for each factor, for otype=1 of a single
%   copy, see 'polyopvar'.
% - A.*B calls otimes(A,B,1), see @polyopvar/times.m.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - otimes
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
% DJ, 10/22/2026: Initial coding;
% DJ, 10/28/2026: Do not merge two states without spatial variables, e.g.
%                   x.*x, into one factor;

% Check the type of product
if nargin<3
    otype = 2;
elseif ~isscalar(otype) || ~any(otype==[1,2])
    error("The type of product, otype, must be 1 (elementwise product) or 2 (tensor product).")
end

% Make sure the sizes are appropriate for elementwise multiplication
sz_A = size(A);     is_dim1_A = sz_A==1;
sz_B = size(B);     is_dim1_B = sz_B==1;
if otype==1 && any(~is_dim1_A & ~is_dim1_B & sz_A~=sz_B)
    error("Matrix dimensions must match for elementwise multiplication.")
end

% Account for case where one object is not a distributed polynomial
if ~isa(A,'polyopvar')
    % Multiplication of constant with polynomial
    if otype==2 && ~all(is_dim1_A)
        error("The tensor product of a constant with a distributed polynomial is only supported for a scalar constant.")
    end
    C = B;
    C.C = A.*B.C;
    return
elseif ~isa(B,'polyopvar')
    % Multiplication of polynomial with constant
    if otype==2 && ~all(is_dim1_B)
        error("The tensor product of a distributed polynomial with a constant is only supported for a scalar constant.")
    end
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
[nfacA,~] = num_factors(A);
dA = repelem(nfacA,nrB,1);
ridcsA = (ridcsA(:)-1)*nrB + (1:nrB);       % repelem(A.Z,nrB,1)
cidcsA = repmat(cidcsA(:),1,nrB);
ZvalsA = repmat(ZvalsA(:),1,nrB);
ridcsB = ridcsB(:) + (0:nrA-1)*nrB;                    % repmat(B.Z,nrA,1);
cidcsB = repmat(cidcsB(:),1,nrA) + reshape(dA(ridcsB),size(ridcsB));
ZvalsB = repmat(ZvalsB(:),1,nrA);
Znew = sparse([ridcsA(:);ridcsB(:)],[cidcsA(:);cidcsB(:)],[ZvalsA(:);ZvalsB(:)],nrA*nrB,max([ncA;nfacA(:)+ncB]));

% Declare the coefficients acting on the full vector of monomials: the
% product of the outputs of type otype (pointwise or Kronecker), and the
% tensor product of the inputs, the monomials
C.C = otimes(A.C,B.C,[otype,2]);
C.dims.out = C.C.dims.out;
C.vars.out = C.C.vars.out;
C.dom.out = C.C.dom.out;
C.depmat.out = C.C.depmat.out;
C.chkval = true;

% The operators of two monomials that each are a state without spatial variables (x.*x) are
% merged into one factor by otimes, as if they were constants: recompute them as two factors
if otype==1
    idcsA = find(ode_monomials(A,nfacA));
    idcsB = find(ode_monomials(B,num_factors(B)));
    if ~isempty(idcsA) && ~isempty(idcsB)
        ops = C.C.ops;
        for i=idcsA(:)'
            for j=idcsB(:)'
                for r=1:size(ops,1)
                    if isa(A.C.ops{r,i},'tensopvar') && isa(B.C.ops{r,j},'tensopvar')
                        % Columns of C are ordered A outer, B inner
                        ops{r,(i-1)*nrB+j} = otimes(A.C.ops{r,i},B.C.ops{r,j},[1,2]);
                    end
                end
            end
        end
        C.C.ops = ops;
    end
end

% otimes merges a constant of A (a single factor, in a monomial of degree 
% 0) with a scalar factor of B without variables, e.g. a state not 
% depending on any spatial variable, into one factor.
is_const = full(nfacA==1 & ~any(A.Z,2));
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function is_ode = ode_monomials(P,nfac)
% Monomials of P that consist of a single scalar state without spatial 
% variables, with one factor
is_state = P.dims.in(:)==1 & ~any(P.depmat.in,2);
[ridcs,~,vals] = find(P.Z);
is_ode = false(size(P.Z,1),1);
is_ode(ridcs(is_state(vals))) = true;
is_ode = is_ode & full(nfac(:)==1) & full(sum(P.Z~=0,2)==1);
end

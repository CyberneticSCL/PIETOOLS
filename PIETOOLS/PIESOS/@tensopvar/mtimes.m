function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'tensopvar' object C representing the product
% of the multiplier operator or tensor-PI operator defined by A with the
% tensor-PI operator defined by B.
%
% INPUTS
% - A:      m x p 'double', 'polynomial', or 'dpvar' object representing
%           some multiplier operator to multiply with B, or a 'tensopvar'
%           object (or 'sopvar', converted to 'tensopvar');
% - B:      p x n 'tensopvar' object acting on some degree-d distributed
%           monomial (or 'sopvar', converted to 'tensopvar'). A scalar B
%           that is not a 'tensopvar' is multiplied from the left, B*A;
%
% OUTPUTS
% - C:      m x n 'tensopvar' object representing the product of A with
%           the tensor-PI operator defined by B;
%
% NOTES
% - If B involves more than one factor (size(B.ops,2)>1), a multiplier A
%   must be a scalar. A general (non-scalar) matrix A rescales only the
%   first factor of each term (ops{i,1}), which represents the product A*B
%   correctly when there is only one factor, and for a scalar A
%   regardless of the number of factors (scalar multiplication commutes
%   with the tensor-product structure). For a general matrix A and a
%   multi-factor B, A would need a Kronecker structure, e.g.
%   kron(A,I,I,...), which this function does not construct.
% - A 'polynomial' multiplier is only supported if it is the identity
%   (then B is returned); other multipliers must be 'double' or 'dpvar'.
% - If A and B are both 'tensopvar' objects with the same number of factors
%   (and proper tensor products, A.otype(2)==B.otype(1)==2), the factors
%   are composed one by one, by the mixed product property
%       (A1 o ... o Ad)(B1 o ... o Bd) = A1B1 o ... o AdBd.
% - If A is nested (a factor of A is a 'tensopvar' acting on several of the
%   factors of the monomial, see 'plus') and B has a factor for each leaf
%   factor of A (num_factors(A)), the factors of B are distributed over the
%   nested factors of A, term by term, and composed recursively.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes
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

% Check that multiplication is supported
if ~isa(B,'tensopvar') && all(size(B)==1)
    % Allow right-multiplication with a scalar
    C = B*A;
    return
end
if isa(A,'sopvar')
    A = sopvar2tensopvar(A);
elseif ~isa(A,'tensopvar') && ~isa(A,'double') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("Multiplication of tensopvar objects with non-tensopvar objects is currently not supported.")
end
if isa(B,'sopvar')
    B = sopvar2tensopvar(B);
elseif ~isa(B,'tensopvar')
    error("Multiplication of tensopvar objects with non-tensopvar objects is currently not supported.")
end

% Make sure the objects are properly declared, so depmat/dims/vars below
if isa(A,'tensopvar') && ~A.chkval
    A = initialize(A);
end
if ~B.chkval
    B = initialize(B);
end

% Check that the inner dimensions of the objects match
Adim = size(A);     Bdim = size(B);
if Adim(2)~=Bdim(1) && ~all(Adim==1)
    error("Inner dimensions of the objects to multiply must match.")
end

% Nested A, with a factor for each leaf factor in B: distribute the factors of B over those of A
if isa(A,'tensopvar') && size(A.ops,2)~=num_factors(A) && size(B.ops,2)==num_factors(A)
    C = mtimes_nested(A,B);
    return
end

% Factorwise composition needs matching numbers of factors
if isa(A,'tensopvar') && size(A.ops,2)~=size(B.ops,2)
    error("Composition of tensor-PI operators is only supported if both are single-factor, or have the same number of factors (composed factor by factor); here, an operator with %d factor(s) is composed with one with %d.",size(A.ops,2),size(B.ops,2))
end
% A nested (multi-factor) operator in B cannot be composed with a
% non-nested operator in A
if isa(A,'tensopvar') && any(cellfun(@(c) isa(c,'tensopvar'),B.ops(:))) && ~any(cellfun(@(c) isa(c,'tensopvar'),A.ops(:)))
    error("Composition of a tensor-PI operator with an operator that contains nested tensor products (e.g. a product of operators acting on the same factor) is not supported.")
end

% Check that the spatial domain of the output of B matches the input of A
if isa(A,'tensopvar')
    if ~isequal(A.vars.in,B.vars.out) || ~isequal(A.depmat.in,B.depmat.out) || ~isequal(A.dims(:,2),B.dims(:,1))
        error("Function space of input of A must match that of output of B.")
    end
end

% Multiplication with an identity operator leaves the other factor
% unchanged
if isa(A,'tensopvar')
    % Identity A (or B) composes with any B (A) whose factors it matches,
    % so no mixed product property is needed.
    if size(B.ops,2)==1 || A.otype(2)==B.otype(1)
        if isequal(A.vars.out,B.vars.out) && isequal(A.dom.out,B.dom.out) &&...
              isequal(A.dom.in,B.dom.out) && isequal(A.depmat.out,B.depmat.out) &&...
                isequal(A.dims(:,1),B.dims(:,1)) && is_id_op(A)
            C = B;
            return
        elseif isequal(B.vars.in,A.vars.in) && isequal(B.dom.in,A.dom.in) &&...
                 isequal(B.dom.out,A.dom.in) && isequal(B.depmat.in,A.depmat.in) &&...
                   isequal(B.dims(:,2),A.dims(:,2)) && is_id_op(B)
            C = A;
            return
        end
    end
elseif size(A,1)==size(A,2)
    % Identity matrix (scalar 1, or identity of size(B,1))
    if isa(A,'double')
        is_id_A = isequal(A,eye(size(A,1)));
    else
        % polynomial/dpvar multiplier: A can depend on spatial variables,
        % so do not wrap A in an 'intop' (its kernels only involve its dummy
        % variables); A is the identity if A-I has no coefficients above 1e-12
        D_A = A - eye(size(A,1));
        Cval = nonzeros(D_A.C);
        is_id_A = isempty(Cval) || max(abs(Cval))<=1e-12;
    end
    if is_id_A
        C = B;
        return
    end
end

% Perform the multiplication
if isa(A,'tensopvar')
    if (size(A.ops,2)==1 && size(B.ops,2)==1) || (A.otype(2)==2 && B.otype(1)==2)
        % Support 'tensopvar' multiplication only if the operators acts on
        % monomials of degree 1, or if we have proper tensor products, so
        % that we can use the mixed product property
        %   (A o B)(C o D) = AC o BD
        C = B;
        C.ops = cell(size(A.ops,1).*size(B.ops,1),size(A.ops,2));
        for k=1:numel(C.ops)
            % Take the product of each pair of operators in A and B
            [idx1,idx2,idx3] = ind2sub([size(A.ops,1),size(B.ops,1),size(C.ops,2)],k);
            C.ops{k} = A.ops{idx1,idx3}*B.ops{idx2,idx3};
        end
        C.dims(:,1) = A.dims(:,1);
        C.vars.out = A.vars.out;
        C.dom.out = A.dom.out;
        C.depmat.out = A.depmat.out;
        C.otype(1) = A.otype(1);
    else
        error("Composition of multi-factor tensor-PI operators requires the inputs of the first and the outputs of the second to be Kronecker products (A.otype(2)==2, B.otype(1)==2); pointwise products are not supported.")
    end
else
    % @sopvar/mtimes supports numeric and 'dpvar' multipliers, but not 'polynomial' (identity returned above)
    if isa(A,'polynomial')
        error("Multiplication of a tensopvar with a 'polynomial' object is only supported if the object is the identity; other multipliers must be of type 'double' or 'dpvar'.")
    end

    % Multiply each of the operators stored in B by A
    if ~all(size(A)==1) && size(B.ops,2)>1
        error("Left-multiplication of a multi-factor tensor-PI operator by a general (non-scalar) matrix is not supported; only a scalar multiplier is supported in that case.")
    end
    C = B;
    for i=1:size(C.ops,1)
        % Multiply the first factor of each term with A
        C.ops{i,1} = A*B.ops{i,1};
    end
    if ~all(size(A)==1)
        C.dims(:,1) = size(A,1);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function C = mtimes_nested(A,B)
% C = MTIMES_NESTED(A,B) returns the composition of the nested 'tensopvar' A, with factors that are
% 'tensopvar' objects acting on several of the factors of the monomial, with the 'tensopvar' B that has one
% factor for each leaf factor of A (num_factors(A)=size(B.ops,2)), each term of B acting on the leaf factors
% of A in order. Each factor Fj of a term of A with n leaf factors is composed with the n consecutive factors
% of the term of B, which are taken together as an operator of their own (the Kronecker product of their
% outputs, as for a mixed product (F1 o ... o Fd)(B1 o ... o Bd) = F1B1 o ... o FdBd): recursively if Fj is
% itself nested. Composition is bilinear, so each pair of terms of A and B gives a term of C, with the
% outputs of the factors of A combined as in A (A.otype(1)), and the inputs of the factors of B as in B.
[KA,dA] = size(A.ops);      [KB,dB] = size(B.ops);
% B can be split between several factors of A only if it has a Kronecker product of its outputs, as for
% factors that are not nested (see the multiplication above); a single factor of A is composed with all of B
if dA>1 && ~(A.otype(2)==2 && B.otype(1)==2)
    error("Composition of multi-factor tensor-PI operators requires the inputs of the first and the outputs of the second to be Kronecker products (A.otype(2)==2, B.otype(1)==2); pointwise products are not supported.")
end
C = tensopvar();
C.otype = [A.otype(1),B.otype(2)];
C.ops = cell(KA*KB,dA);
for k=1:KA*KB
    [iA,iB] = ind2sub([KA,KB],k);
    pos = 0;        % leaf factors of B used so far
    for j=1:dA
        Fj = A.ops{iA,j};
        if isa(Fj,'tensopvar')
            n = num_factors(Fj);
        else
            n = 1;
        end
        if pos+n>dB
            error("The number of factors of B must match the number of factors in each term of A, including those of nested operators.")
        end
        if isa(Fj,'tensopvar')
            % Operator of its own acting on factors pos+1,...,pos+n of the term of B
            Bj = B;
            Bj.ops = B.ops(iB,pos+1:pos+n);
            Bj.chkval = false;
            Bj = initialize(Bj);
            C.ops{k,j} = Fj*Bj;
        else
            C.ops{k,j} = Fj*B.ops{iB,pos+1};
        end
        pos = pos+n;
    end
    if pos~=dB
        error("The number of factors of B must match the number of factors in each term of A, including those of nested operators.")
    end
end
C = initialize(C);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
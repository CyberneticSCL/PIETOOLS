function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'tensopvar' object C representing the product
% of the multiplier operator or tensor-PI operator defined by A with the
% tensor-PI operator defined by B.
%
% INPUTS
% - A:      m x p 'double' or 'polynomial' object representing
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
% - A multiplier A multiplies the first factor of each term of B. If B
%   has more than one factor (size(B.ops,2)>1), A must be a scalar.
% - A 'polynomial' multiplier is only supported if it is the identity
%   (then B is returned), and a 'dpvar' multiplier is not supported (the
%   factors cannot depend on decision variables); other multipliers must
%   be 'double'.
% - If A and B are both 'tensopvar' objects with the same number of factors
%   (and proper tensor products, A.otype(2)==B.otype(1)==2), the factors
%   are composed one by one, by the mixed product property
%       (A1 o ... o Ad)(B1 o ... o Bd) = A1B1 o ... o AdBd.
%   For different numbers of factors, the factors acting on R (without
%   spatial variables) are first merged into their neighbors; if the
%   numbers still differ, an error is thrown.
% - If A is nested (a factor of A is a 'tensopvar' acting on several of the
%   factors of the monomial) and B has a factor for each leaf factor of A
%   (num_factors(A)), the factors of B are distributed over the nested
%   factors of A, term by term, and composed recursively. If B has nested
%   factors and A has none, an error is thrown.

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

% Check that multiplication is supported
if (isnumeric(B) || isa(B,'polynomial') || isa(B,'dpvar')) && all(size(B)==1)
    % Allow right-multiplication with a scalar; a 1x1 operator ('sopvar')
    % is composed below, as B*A would compose in the wrong order
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
% --> try merging the factors that are operate on R into their neighbors
if isa(A,'tensopvar') && size(A.ops,2)~=size(B.ops,2)
    A = combine_factors(A);
    B = combine_factors(B);
end
% Throw an error if the numbers of factors still differ
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
    % Compare the variables as columns, so that empty cells of different shape (0x0, 1x0) are equal
    if ~isequal(A.vars.in(:),B.vars.out(:)) || ~isequal(A.depmat.in,B.depmat.out) || ~isequal(A.dims(:,2),B.dims(:,1))
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
        error("Multiplication of a tensopvar with a 'polynomial' object is only supported if the object is the identity; other multipliers must be of type 'double'.")
    elseif isa(A,'dpvar')
        % The factors of a 'tensopvar' cannot depend on decision variables
        error("Multiplication of a 'tensopvar' with a 'dpvar' is not supported: the factors of a 'tensopvar' cannot depend on decision variables.")
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
% C = MTIMES_NESTED(A,B) returns the composition A*B of the nested
% 'tensopvar' A, with factors that are 'tensopvar' objects acting on
% several of the factors of the monomial, with the 'tensopvar' B that has
% a factor for each leaf factor of A (num_factors(A)=size(B.ops,2)).
% Each factor Fj of a term of A with n leaf factors is composed
% (recursively, if Fj is nested) with the next n factors of a term of B,
% taken together as one 'tensopvar' with the product types B.otype, as
% for a mixed product (F1 o ... o Fd)(B1 o ... o Bd) = F1B1 o ... o FdBd.
% Each pair of terms of A and B gives a term of C, with C.otype =
% [A.otype(1),B.otype(2)]. If A has more than one factor, A.otype(2)=2
% and B.otype(1)=2 are required.
%
% INPUTS
% - A:        nested 'tensopvar' object;
% - B:        'tensopvar' object with a factor per leaf factor of A;
%
% OUTPUTS
% - C:        'tensopvar' object, the composition A*B.
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
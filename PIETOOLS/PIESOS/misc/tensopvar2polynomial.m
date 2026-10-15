function [P,is_ok] = tensopvar2polynomial(Top)
% [P,IS_OK] = TENSOPVAR2POLYNOMIAL(TOP) converts the tensor-PI operator
% TOP, if it acts on a constant monomial (R), to the 'polynomial' matrix P
% that represents it: the operator x -> P*x from R^n to the functions in
% L2^m[S] that it maps to.
%
% INPUTS
% - Top:    m x n 'tensopvar' object (or 'sopvar') acting on R, i.e.
%           without input variables, all of whose factors are multipliers
%           from R;
%
% OUTPUTS
% - P:      m x n 'polynomial' object, in the output variables of TOP,
%           such that for the factors x1,...,xd in R^n1 x ... x R^nd of a
%           distributed monomial,
%               Top(x1,...,xd) = P*(x1 o ... o xd),
%           where o is the Kronecker product (otype(2)=2, the only type
%           of product of the inputs that apply_op supports) or the
%           pointwise product of the inputs (otype(2)=1);
% - is_ok:  logical, false if TOP cannot be converted, and P=[]. If this
%           output is not requested, an error is thrown in that case;
%
% NOTES
% - TOP can be converted if no factor has input variables, and each
%   factor is a 'sopvar', or a 'tensopvar' that can be converted. Then
%   each factor is a multiplier from R^nj, the matrix of which is its
%   action on the unit vectors, apply_sopvar(Fj,eye(nj)). The product of
%   the factors in a term is taken as defined for 'tensopvar' objects:
%   pointwise (otype(1)=1) or Kronecker (otype(1)=2) along the rows,
%   Kronecker (otype(2)=2) or pointwise along the columns, and the
%   polynomials of the terms are added.
% - A Kronecker product along the rows (otype(1)=2) of factors that
%   depend on the same spatial variable is a function of several copies
%   of that variable, which a 'polynomial' cannot represent: TOP is not
%   converted.
% - P can be used as a block of a 'tensopmat', see 'tensopmat'.
%
% EXAMPLE
%   pvar s
%   vars = struct('in',{cell(1,0)},'out',{{'s'}});
%   dom = struct('in',zeros(0,2),'out',[0,1]);
%   Top = mat2tensopvar([1+s;2],vars,dom);
%   P = tensopvar2polynomial(Top);              % P = [1+s;2]
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tensopvar2polynomial
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

[P,msg] = convert_op(Top);
is_ok = isempty(msg);
if ~is_ok
    P = [];
    if nargout<2
        error(msg)
    end
end

end


function [P,msg] = convert_op(Top)
% Recursive conversion: msg is empty if P is the polynomial matrix of Top, else the reason it is not
P = [];     msg = '';
if isa(Top,'sopvar')
    Top = sopvar2tensopvar(Top);
elseif ~isa(Top,'tensopvar')
    msg = "Input must be of type 'tensopvar'.";
    return
end
if ~Top.chkval
    Top = initialize(Top);
end
vars = Top.vars;
if ~isempty(vars.in)
    msg = "The operator has input variables, so acts on a function space, not on a constant monomial.";
    return
end

ops = Top.ops;      otype = Top.otype;
[nterms,nfctrs] = size(ops);
Psum = 0;
try
    for k=1:nterms
        Pk = [];
        for j=1:nfctrs
            Fj = ops{k,j};
            if isa(Fj,'sopvar')
                Pj = apply_sopvar(Fj,eye(size(Fj,2)));
            elseif isa(Fj,'tensopvar')
                [Pj,msg] = convert_op(Fj);
                if ~isempty(msg)
                    return
                end
            else
                msg = "A factor of type '" + string(class(Fj)) + "' cannot be converted: only 'sopvar' and 'tensopvar' factors can.";
                return
            end
            if isempty(Pk)
                Pk = Pj;
            else
                Pk = combine_matrices(Pk,Pj,otype);
            end
        end
        Psum = Psum + Pk;
    end
catch ME
    msg = "The factors of the operator cannot be converted: " + string(ME.message);
    return
end
Psum = polynomial(Psum);
if ~isequal(size(Psum),size(Top))
    msg = "The polynomial matrix does not have the dimensions of the operator.";
    return
end
P = Psum;

end


function C = combine_matrices(A,B,type)
% Product of two polynomial matrices A, B of the factors of a tensor-PI operator: pointwise (type(1)=1) or
% Kronecker (type(1)=2) along the rows, and Kronecker (type(2)=2) or pointwise along the columns. The
% matrix of the products A(r1,c1)*B(r2,c2), with r1=r2 (c1=c2) for a pointwise product. Same as
% otimes_native, a local function of @tensopmat/otimes.m, which is not accessible from here.
if type(1)==2 && ~isempty(intersect(used_vars(A),used_vars(B)))
    error("A Kronecker product along the rows of factors that share a spatial variable is not a polynomial.")
end
[mA,nA] = size(A);      [mB,nB] = size(B);
if type(1)==1
    if mA~=mB
        error("Row dimensions of factors must match for a pointwise product.")
    end
    ridcsA = 1:mA;      ridcsB = 1:mB;
else
    ridcsA = repelem(1:mA,mB);      ridcsB = repmat(1:mB,1,mA);
end
if type(2)==1
    if nA~=nB
        error("Column dimensions of factors must match for a pointwise product.")
    end
    cidcsA = 1:nA;      cidcsB = 1:nB;
else
    cidcsA = repelem(1:nA,nB);      cidcsB = repmat(1:nB,1,nA);
end
C = A(ridcsA,cidcsA).*B(ridcsB,cidcsB);

end


function vars = used_vars(P)
% Names of the variables that appear in a term of the polynomial P with a nonzero coefficient
mons = find(any(P.coefficient,2));
vars = reshape(P.varname(any(full(P.degmat(mons,:))>0,1)),1,[]);

end

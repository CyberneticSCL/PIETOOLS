function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'tensopmat' object C representing the product
% of A, a multiplier operator or a matrix of tensor-PI operators, with the
% matrix of tensor-PI operators B.
%
% INPUTS
% - A:      m x p 'tensopmat' object (or 'tensopvar', 'intop' or
%           'sopvar'), or an m x p 'double' matrix, or a scalar
%           'polynomial' or 'dpvar', representing a multiplier operator;
% - B:      p x n 'tensopmat' object (or 'tensopvar', 'intop' or
%           'sopvar') representing a matrix of tensor-PI operators. If A
%           is a 'tensopmat', B may also be a scalar 'double',
%           'polynomial' or 'dpvar', in which case B*A is returned;
%
% OUTPUTS
% - C:      m x n 'tensopmat' object representing the product of A with B;
%
% NOTES
% - If A is a 'tensopmat', its input function space must match the output
%   function space of B. Multiplication with an identity operator returns
%   the other factor unchanged.
% - If A is a non-scalar 'double' matrix, C consists of a single row-block,
%   and all row-blocks of B must depend on the same output variables.
% - Multiplying a multiplier with an operator whose outputs are tensor
%   products of PI operators (B.depmat.out>1) is only supported if all
%   row-blocks of B have dimension 1. Use '.*' for elementwise
%   multiplication instead.
%

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
% DJ, 10/21/2026: Compare empty variable lists regardless of shape;

% Check that multiplication is supported
if isa(A,'sopvar') || isa(A,'tensopvar') || isa(A,'intop')
    A = tensopmat(A);
end
if isa(B,'sopvar') || isa(B,'tensopvar') || isa(B,'intop')
    B = tensopmat(B);
end
if isa(A,'double') && ~all(size(A)==1) && isa(B,'tensopmat') 
    % Build the result directly, column-by-column: for each of B's column- 
    % blocks j, sum over B's row-blocks k of (the A-columns feeding row k)  
    % times B.ops{k,j}.
    if size(B.ops,1)==1 && size(A,1)==size(A,2) && size(A,2)==size(B,1) &&...
            isequal(A,eye(size(A,1)))
        % An identity matrix leaves a single row of blocks unchanged.
        C = B;
        return
    end
    if ~all(B.dims.out==1) && any(any(B.depmat.out>1))
        error("Matrix-multiplication with tensor-products of PI operators is currently not supported. Use '.*' for elementwise multiplication instead.")
    end
    if size(A,2)~=size(B,1)
        error("Inner dimensions of the objects to multiply must match.")
    end
    [p,d] = size(B.ops);
    nn_out = [0;cumsum(B.dims.out)];
    ops_new = cell(1,d);
    % Take product of A with each block of columns of B
    for j=1:d
        Cj = [];
        for k=1:p
            if isempty(B.ops{k,j})
                % Empty operator is considered to be 0
                continue 
            end
            % Take product of columns of A acting on block k of rows in B
            Ak = A(:,nn_out(k)+1:nn_out(k+1));
            term = Ak*B.ops{k,j};
            if isempty(Cj)
                Cj = term;
            elseif is_native_op(Cj) || is_native_op(term)
                % Matrices of type double/polynomial/dpvar are added to each other, or to a tensopvar
                % acting on a constant monomial, converted to a polynomial, as in plus
                if ~is_native_op(Cj)
                    [Cj,is_ok] = tensopvar2polynomial(Cj);
                elseif ~is_native_op(term)
                    [term,is_ok] = tensopvar2polynomial(term);
                else
                    is_ok = true;
                end
                if ~is_ok
                    error("A block of type 'double', 'polynomial' or 'dpvar' can only be added to an operator that acts on R, a constant monomial, and can be converted to a 'polynomial'.")
                end
                Cj = Cj+term;
            elseif ~isequal(size(Cj),size(term)) || ~isequal(Cj.vars,term.vars) || ~isequal(Cj.dom,term.dom) || ~isequal(Cj.deps,term.deps)
                error("tensopmat:mtimes_heterogeneous_rows", ...          
                    ['Matrix-multiplication of a general matrix with a ' ...
                     'tensor-PI operator whose row-blocks depend on ' ...
                     'different (sub)sets of the output variables is ' ...
                     'currently not supported (column %d, row-block %d ' ...
                     'of the operator does not share the same variable ' ...
                     'dependency as an earlier row-block contributing to ' ...
                     'the same output).'],j,k)  
            else
                Cj = Cj+term;
            end
        end
        ops_new{1,j} = Cj;
    end
    % The d columns above populate a SINGLE row-block, so (per 
    % @tensopmat/initialize.m's row check) all columns must share the same
    % output vars/dom/deps
    Cref = [];
    for j=1:d
        if isempty(ops_new{1,j}) || is_native_op(ops_new{1,j})
            continue        % no variables to compare for a matrix, see initialize
        elseif isempty(Cref)
            Cref = ops_new{1,j};
        elseif ~isequal(Cref.vars.out,ops_new{1,j}.vars.out) || ~isequal(Cref.dom.out,ops_new{1,j}.dom.out) || ~isequal(Cref.deps.out,ops_new{1,j}.deps.out)
            error("tensopmat:mtimes_heterogeneous_columns", ...
                ['Matrix-multiplication of a general matrix with a ' ...
                 'tensor-PI operator whose row-blocks depend on ' ...
                 'different (sub)sets of the output variables is ' ...
                 'currently not supported (column %d of the result does ' ...
                 'not share the output-variable dependency established ' ...
                 'by an earlier column.'],j)
        end
    end
    if any(native_mask(ops_new))
        C = tensopmat(ops_new,B.vars,B.dom);    % domains of the variables in the matrices
    else
        C = tensopmat(ops_new);
    end
    return
end

% Check that the operators are properly declared
if isa(A,'tensopmat') 
    if ~A.chkval && ~isvalid(A)
        error("The first operator is not properly declared.")
    else
        A = clean_op(A);
    end
end
if isa(B,'tensopmat')
    if ~B.chkval && ~isvalid(B)
        error("The second operator is not properly declared.")
    else
        B = clean_op(B);
    end
end

% Check for (currently) unsupported cases
if isa(A,'tensopmat') && isa(B,'tensopmat')
    % Compare the variables as columns, so that empty cells of different shape (0x0, 1x0) are equal
    if ~isequal(A.vars.in(:),B.vars.out(:)) || ~isequal(A.depmat.in,B.depmat.out)
        error("Input function space of first factor must match output space of second factor.")
    elseif ~isequal(A.dims.in,B.dims.out)
        error("Inner dimensions of operators must match.")
    end
elseif isa(A,'tensopmat') && all(size(B)==1)
    % Allow right-multiplication with a scalar
    C = B*A;
    return
elseif isa(A,'double') || isa(A,'polynomial') || isa(A,'dpvar')
    if ~all(size(A)==1)
        error("Matrix multiplication with non-tensopmat objects is not supported.")
    end
else
    error("'tensopvar' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end


% Check that the inner dimensions of the objects match
Adim = size(A);     Bdim = size(B);
if Adim(2)~=Bdim(1) && ~all(Adim==1)
    error("Inner dimensions of the objects to multiply must match.")
end

% Multiplication with an identity operator leaves the other factor 
% unchanged.
if isa(A,'tensopmat')
    if isequal(A.vars.out,B.vars.out) && isequal(A.dom.out,B.dom.out) &&...
         isequal(A.dom.in,B.dom.out) && isequal(A.depmat.out,B.depmat.out) &&...
           isequal(A.dims.out,B.dims.out) && is_id_op(A)
        C = B;
        return
    elseif isequal(B.vars.in,A.vars.in) && isequal(B.dom.in,A.dom.in) &&...
             isequal(B.dom.out,A.dom.in) && isequal(B.depmat.in,A.depmat.in) &&...
               isequal(B.dims.in,A.dims.in) && is_id_op(B) 
        C = A;
        return
    end
else
    % A is a scalar here 
    if isa(A,'double')
        is_id_A = isequal(A,1);
    else
        D_A = A - 1;                    % polynomial/dpvar scalar
        Cval = nonzeros(D_A.C);
        is_id_A = isempty(Cval) || max(abs(Cval))<=1e-12;
    end
    if is_id_A
        C = B;
        return
    end
end

% Prohibit matrix-multiplication with actual tensor products of PI
% operators for now, since this is not quite trivial
%   (note that  A*(B1*x)(s)*(B2*x)(s) ~= (A*B1*x)(s)*(B2*x)(s))
if ~isa(A,'tensopmat') && ~all(B.dims.out==1) && any(any(B.depmat.out>1))
    error("Matrix-multiplication with tensor-products of PI operators is currently not supported. Use '.*' for elementwise multiplication instead.")
end

if isa(A,'tensopmat')
    % Take the matrix product of the tensor-PI operators
    C = tensopmat();
    C.ops = cell(size(A.ops,1),size(B.ops,2));
    for idx = 1:numel(C.ops)
        [ridx,cidx] = ind2sub(size(C.ops),idx);
        for k=1:size(A.ops,2)
            if ~isempty(A.ops{ridx,k}) && ~isempty(B.ops{k,cidx})
                % The product of the blocks, of a block of type double/polynomial/dpvar with an operator
                % if needed, see block_product
                Pk = block_product(A,ridx,k,B.ops{k,cidx});
                if isempty(C.ops{idx})
                    C.ops{idx} = Pk;
                else
                    C.ops{idx} = add_blocks(C.ops{idx},Pk);
                end
            end
        end
    end
    % Set the variables and domain of the output
    C.dims.out = A.dims.out;        C.dims.in = B.dims.in;
    C.vars.out = A.vars.out;        C.vars.in = B.vars.in;
    C.dom.out = A.dom.out;          C.dom.in = B.dom.in;
    C.depmat.out = A.depmat.out;    C.depmat.in = B.depmat.in;
    C.chkval = true;
else
    % Multiply each of the operators stored in B by A
    C = B;
    for i=1:numel(C.ops)
        if ~isempty(B.ops{i})
            C.ops{i} = A*B.ops{i};
        end
    end
    C.chkval = true;
    if ~isa(A,'double') && any(native_mask(C.ops(:)))
        C = initialize(C);      % A may introduce new variables
    end
end

end


function P = block_product(A,i,k,Bkj)
% P = BLOCK_PRODUCT(A,I,K,BKJ) returns the product of block (I,K) of the 'tensopmat' A with the block
% BKJ of the second factor. Two blocks of the same kind, both matrices of type 'double', 'polynomial'
% or 'dpvar' acting on R, or both operators, are multiplied directly. For a matrix with an operator:
% - matrix M (acting on R, mapping to its row of A) times an operator that maps to R: the matrix is
%   converted to the multiplier operator from R (native2op), and the operators are composed;
% - operator times a matrix, which maps R to functions: the operator is applied to the functions in
%   the columns of the matrix (apply_to_native), the result is a matrix.
Aik = A.ops{i,k};
is_natA = is_native_op(Aik);        is_natB = is_native_op(Bkj);
if is_natA==is_natB
    P = Aik*Bkj;
elseif is_natA
    if isa(Bkj,'intop')
        error("The product of a block of type 'double' or 'polynomial' with a functional ('intop') is not supported.")
    end
    P = native2op(Aik,A,i)*Bkj;
else
    P = apply_to_native(Aik,Bkj);
end
end

function P = add_blocks(P1,P2)
% P = ADD_BLOCKS(P1,P2) returns the sum of two blocks. A matrix of type 'double', 'polynomial' or 'dpvar'
% is added to an operator after converting the operator to a 'polynomial', if it acts on R, as in plus.
is_nat1 = is_native_op(P1);         is_nat2 = is_native_op(P2);
is_ok = true;
if is_nat1 && ~is_nat2
    [P2,is_ok] = tensopvar2polynomial(P2);
elseif is_nat2 && ~is_nat1
    [P1,is_ok] = tensopvar2polynomial(P1);
end
if ~is_ok
    error("A block of type 'double', 'polynomial' or 'dpvar' can only be added to an operator that acts on R, a constant monomial, and can be converted to a 'polynomial'.")
end
P = P1+P2;
end

function Y = apply_to_native(Aop,M)
% Y = APPLY_TO_NATIVE(AOP,M) applies the operator AOP to the functions in the columns of
% the matrix M of type 'double', 'polynomial' or 'dpvar', which maps R^k to the functions
% in L2^n[S] (see `tensopmat`), and returns the matrix Y of the same type that represents
% the composition x -> AOP*(M*x), mapping R^k to the output space of AOP.
%
% INPUTS
% - Aop:    m x n 'tensopvar' object with a single factor, or 'intop' object;
% - M:      n x k matrix of type 'double', 'polynomial' or 'dpvar';
%
% OUTPUTS
% - Y:      m x k matrix of type 'polynomial', or of type 'dpvar' if M is;
%
% NOTES
% - For a 'dpvar' M(d;s), affine in the decision variables d, the columns are
%   Y_c(d;s) = sum_{r,a} c_{r,a}(d) * AOP*(e_r s^a), for the coefficients c_{r,a}(d) of the
%   monomials s^a in the rth component of the cth column of M. So AOP is only applied to the
%   monomials appearing in M, n*na functions, and the coefficients are combined by one
%   sparse matrix product.
% - The output variables of AOP with a Kronecker product of the (single) factor, otype(1)=2,
%   would be renamed s_1 by apply_op: the product of a single factor is the same pointwise,
%   so that is used here.
%

if isa(Aop,'tensopvar')
    if size(Aop.ops,2)~=1
        error("The product of an operator with several factors with a block of type 'double', 'polynomial' or 'dpvar' is not supported.")
    end
    Aop.otype = [1,Aop.otype(2)];
end
if isa(M,'double')
    M = polynomial(M);
end
if isa(M,'polynomial')
    Y = apply_op(Aop,M);            % the columns of M are applied together
    return
end

% A matrix of type 'dpvar': apply the operator to the monomials in M, and combine the results
[n,kcols] = size(M);
nd = numel(M.dvarname);
deg = full(M.degmat);       vn = M.varname;
[na,nv] = size(deg);
if na==0 || isempty(nonzeros(M.C))
    Y = apply_op(Aop,polynomial(zeros(n,kcols)));       % M is zero
    return
end
Pcols = cell(1,n*na);
Id = eye(n);
for al=1:na
    if nv==0
        mono = polynomial(1);
    else
        mono = polynomial(1,deg(al,:),vn,[1,1]);
    end
    for r=1:n
        Pc = apply_op(Aop,Id(:,r)*mono);
        if isa(Pc,'dpvar')
            error("The product of an operator that depends on decision variables with a block of type 'dpvar' is not supported.")
        end
        Pcols{(r-1)*na+al} = polynomial(Pc);
    end
end
Pall = [Pcols{:}];                  % m x (n*na) polynomial: A applied to e_r*s^a, in column (r-1)*na+a
mout = size(Pall,1);
nbeta = size(Pall.degmat,1);        % monomials in the output

% Coefficients of the monomials of the output functions: row (i,beta), column (r,a)
[bb,ee,pv] = find(Pall.coefficient);
bb = bb(:);     ee = ee(:);     pv = pv(:);     % (:), as find returns rows for a row vector
iout = mod(ee-1,mout)+1;            rho = ceil(ee/mout);
Pm = sparse((iout-1)*nbeta+bb,rho,pv,mout*nbeta,n*na);

% Coefficients of M, rows (r outer, decision variable inner, constant first), columns (c outer, a inner):
% regroup as rows (r,a), columns (c,k) for the decision variable k, only for the (c,k) that occur
[ii,jj,cv] = find(M.C);
ii = ii(:);     jj = jj(:);     cv = cv(:);
rr = ceil(ii/(nd+1));               kk = mod(ii-1,nd+1);
cc = ceil(jj/na);                   aa = mod(jj-1,na)+1;
[ukey,~,kcol] = unique((cc-1)*(nd+1)+kk);
ukey = ukey(:);     kcol = kcol(:);
Cp = sparse((rr-1)*na+aa,kcol,cv,n*na,numel(ukey));

% Combine, and set the coefficients of the output in the layout of a 'dpvar'
Z = Pm*Cp;
[zi,zj,zv] = find(Z);
zi = zi(:);     zj = zj(:);     zv = zv(:);
i2 = ceil(zi/nbeta);                be = mod(zi-1,nbeta)+1;
kz = ukey(zj);                      c2 = floor(kz/(nd+1))+1;        k2 = mod(kz,nd+1);
CY = sparse((i2-1)*(nd+1)+k2+1,(c2-1)*nbeta+be,zv,mout*(nd+1),kcols*nbeta);
Y = dpvar(CY,full(Pall.degmat),Pall.varname,M.dvarname,[mout,kcols]);

end

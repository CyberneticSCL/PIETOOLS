function C = otimes_params(A,B,otype)
% C = OTIMES_PARAMS(A,B,OTYPE) computes the otimes product of type OTYPE
% of the parameters A and B.
%
% INPUTS
% - A:      m1 x n1 array of type 'double' (full or sparse), 'polynomial',
%           or 'dpvar'. Also supported are: a 'struct' with fields
%           'ridcs', 'cidcs', 'vals' and 'size', representing a sparse
%           matrix (B must then be of the same form); a 'cell' array of
%           parameters (see NOTES); a 'quadPoly' object; and a 'struct'
%           with field 'Zt', representing a 'quadPoly' object with
%           coefficients specified as cell;
% - B:      m2 x n2 array of type 'double', 'polynomial', or 'dpvar', or
%           another of the types listed for A;
% - otype:  (optional) 1 x 2 array of integers 1 or 2, specifying the
%           type of product to take along the row (otype(1)) and column
%           (otype(2)) dimensions, where 1 means an elementwise product
%           and 2 means a Kronecker product. A scalar value specifies
%           otype(1), and otype(2) is set to 2. Defaults to [1,2].
%           Note that otype(1)=1 requires m1=m2 (or m1=1 or m2=1), and
%           otype(2)=1 requires n1=n2 (or n1=1 or n2=1);
%
% OUTPUTS
% - C:      m x n array of the same type as the inputs, representing the
%           Kronecker or elementwise product of the input parameters along
%           the row and column dimensions, where m=m1=m2 if otype(1)=1 and
%           m=m1*m2 if otype(1)=2, and n=n1=n2 if otype(2)=1 and n=n1*n2
%           if otype(2)=2;
%
% NOTES
% - If A or B is a cell, the other is treated as a cell with a single
%   element if it is not a cell. C is then a cell of size size(A).*size(B),
%   with each element the product of an element of A with an element of
%   B (each element is a term in an operator), using the same otype.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - otimes_params
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

if nargin<=2
    % Use elementwise product along row dimension by default
    otype = [1,2];
elseif isscalar(otype)
    % Use Kronecker procuct along column dimension by default
    otype = [otype,2];
end

% Check that multiplication is supported
if otype(1)==1 && size(A,1)~=size(B,1) && ~isa(A,'cell') && ~isa(A,'struct')
    if size(A,1)==1
        A = repmat(A,size(B,1),1);
    elseif size(B,1)==1
        B = repmat(B,size(A,1),1);
    else
        error("Row dimensions must match for elementwise product.")
    end
end
if otype(2)==1 && size(A,2)~=size(B,2) && ~isa(A,'cell') && ~isa(A,'struct')
    if size(A,2)==1
        A = repmat(A,1,size(B,2));
    elseif size(B,2)==1
        B = repmat(B,1,size(A,2));
    else
        error("Column dimensions must match for elementwise product.")
    end
end

% Call the appropriate subroutine
if isa(A,'double') && isa(B,'double')
    if issparse(A) || issparse(B)
        C = otimes_sparse(A,B,otype);
    else
        C = otimes_double(A,B,otype);
    end
elseif isa(A,'struct') && isfield(A,'ridcs')
    C = otimes_sparse_struct(A,B,otype);
elseif isa(A,'cell') || isa(B,'cell')
    C = otimes_cell(A,B,otype);
elseif (isa(A,'dpvar') || isa(A,'polynomial') || isa(A,'double')) &&...
        (isa(B,'dpvar') || isa(B,'polynomial') || isa(B,'double'))
    C = otimes_dpvar(A,B,otype);
elseif isa(A,'quadPoly') || isa(B,'quadPoly')
    C = otimes_quadPoly(A,B,otype);
elseif (isa(A,'struct') && isfield(A,'Zt')) || (isa(B,'struct') && isfield(B,'Zt'))
    C = otimes_quadPolyCell(A,B,otype);
else
    error("Input type is not supported.")
end

end



%%
function C = otimes_double(A,B,type)
% % Compute the otimes product for full arrays

if type(1)==1 && type(2)==1
    % Elementwise product along both dimensions
    C = A.*B;
elseif type(1)==1
    % Kronecker product along column dimension
    Anew = repelem(A,1,size(B,2));
    Bnew = repmat(B,1,size(A,2));
    C = Anew.*Bnew;
elseif type(2)==1
    % Kronecker product along row dimension
    Anew = repelem(A,size(B,1),1);
    Bnew = repmat(B,size(A,1),1);
    C = Anew.*Bnew;
else
    % Kronecker product along both dimensions
    C = kron(A,B);
end

end


%% Sparse Matrix Implementation (explicitly using row and column indices)
function C = otimes_sparse_struct(A,B,type)
% % Compute the otimes product for sparse matrices specified as structs
ridcsA = A.ridcs(:);    ridcsB = B.ridcs(:);
cidcsA = A.cidcs(:);    cidcsB = B.cidcs(:);
valsA = A.vals(:);      valsB = B.vals(:);
mB = B.size(1);         nB = B.size(2);

% Distinguish types of products
if isempty(ridcsA) || isempty(ridcsB)
    % % One of the objects is zero
    % --> any product will be zero
    ridcs = zeros(1,0);
    cidcs = zeros(1,0);
    vals = zeros(1,0);
elseif type(1)==1 && type(2)==1
    % % Elementwise product
    % Take only product of values with the same row and column index.
    [rtnB,rtnA] = find(ridcsB==ridcsA' & cidcsB==cidcsA');
    ridcs = ridcsA(rtnA);
    cidcs = cidcsA(rtnA);
    vals = valsA(rtnA).*valsB(rtnB);

elseif type(1)==1
    % % Elementwise product along row dimension only
    % Take every combination of row and column indices
    ridcs1 = repelem(ridcsA,numel(ridcsB),1);
    ridcs2 = repmat(ridcsB,numel(ridcsA),1);
    cidcs1 = repelem(cidcsA,numel(cidcsB),1);
    cidcs2 = repmat(cidcsB,numel(cidcsA),1);
    vals1 = repelem(valsA,numel(valsB),1);
    vals2 = repmat(valsB,numel(valsA),1);


    % Take product only if row indices match
    rtn_idcs = ridcs1==ridcs2;
    ridcs = ridcs1(rtn_idcs);
    cidcs = (cidcs1(rtn_idcs)-1)*nB + cidcs2(rtn_idcs);
    vals = vals1(rtn_idcs).*vals2(rtn_idcs); 

elseif type(2)==1
    % % Elementwise product along column dimension only
    % Avoid duplicate code: just run type=[2,1] routine for transpose
    At = A;     At.ridcs = cidcsA;      At.cidcs = ridcsA; 
    Bt = B;     Bt.ridcs = cidcsB;      Bt.cidcs = ridcsB;
    At.size = A.size([2,1]);    Bt.size = B.size([2,1]);
    Ct = otimes_sparse_struct(At,Bt,[1,2]);

    ridcs = Ct.cidcs;
    cidcs = Ct.ridcs;
    vals = Ct.vals;

else
    % % Proper Kronecker product
    % Take every combination of row indices
    ridcs1 = repelem(ridcsA,numel(ridcsB),1);
    ridcs2 = repmat(ridcsB,numel(ridcsA),1);
    ridcs = (ridcs1-1)*mB + ridcs2;

    % Take every combination of column indices
    cidcs1 = repelem(cidcsA,numel(cidcsB),1);
    cidcs2 = repmat(cidcsB,numel(cidcsA),1);
    cidcs = (cidcs1-1)*nB + cidcs2;

    % Multiply values
    vals1 = repelem(valsA,numel(valsB),1);
    vals2 = repmat(valsB,numel(valsA),1);
    vals = vals1.*vals2;
end

% Return output as struct
C = struct();
C.ridcs = ridcs;
C.cidcs = cidcs;
C.vals = vals;
C.size = A.size.*(B.size.^(type-1));

end



%% Sparse Matrix Implementation
function C = otimes_sparse(A,B,type)
% % Compute the otimes product for sparse matrices

% Distinguish different types
if type(1)==1 && type(2)==1
    % % Elementwise product
    % Use the MATLAB built-in function
    C = A.*B;

elseif type(1)==2 && type(2)==2
    % % Proper Kronecker product
    % Use the MATLAB built-in function
    C = kron(A,B);

else
    % % Elementwise product along either row or column dimension
    % Use 'struct' based implementation
    [ridcsA,cidcsA,valsA] = find(A);
    [ridcsB,cidcsB,valsB] = find(B);

    Atmp = struct();        Btmp = struct();
    Atmp.ridcs = ridcsA;    Btmp.ridcs = ridcsB;
    Atmp.cidcs = cidcsA;    Btmp.cidcs = cidcsB;
    Atmp.vals = valsA;      Btmp.vals = valsB;
    Atmp.size = size(A);    Btmp.size = size(B);

    Ctmp = otimes_sparse_struct(Atmp,Btmp,type);

    % Convert to sparse matrix
    Csz = size(A).*(size(B).^(type-1));
    C = sparse(Ctmp.ridcs,Ctmp.cidcs,Ctmp.vals,Csz(1),Csz(2));
end

end



%% dpvar Implementation
function C = otimes_dpvar(A,B,otype)
% % Compute the otimes product for 'dpvar' with 'polynomial' objects

% Only one of the object can be of type 'dpvar'
if isa(A,'dpvar') && isa(B,'dpvar')
    error("Multiplication of 'dpvar' objects is not supported.")
end
if isa(A,'double')
    A = polynomial(A);
elseif isa(B,'double')
    B = polynomial(B);
end

% Decompose the objects
[mA,nA] = size(A);          [mB,nB] = size(B);
C_A = A.C;                  C_B = B.C;
Adegmat = A.degmat;         Bdegmat = B.degmat;
Apvarname = A.varname;      Bpvarname = B.varname;
nZA = size(Adegmat,1);      nZB = size(Bdegmat,1);

% Extract the decision variables, if appropriate
use_dpvar = false;
if isa(A,'dpvar')
    use_dpvar = true;
    dvarname = A.dvarname;
    ndvars = numel(dvarname);
elseif isa(B,'dpvar')
    use_dpvar = true;
    dvarname = B.dvarname;
    ndvars = numel(dvarname);
end

% Synchronize the pvarnames
[pvarnames,IA_p,~] = union(Apvarname,Bpvarname,'stable');   % returns indices in A
[~,~,IB_p] = intersect(Bpvarname,pvarnames,'stable');       % returns indices in B
np_all = length(pvarnames);

% Update the columns of the degmats accordingly
[mon_idcsA,pidcsA,degvalsA] = find(Adegmat);
[mon_idcsB,pidcsB,degvalsB] = find(Bdegmat);
pidcsA = IA_p(pidcsA);
pidcsB = IB_p(pidcsB);
Adegmat = sparse(mon_idcsA,pidcsA,degvalsA,nZA,np_all);
Bdegmat = sparse(mon_idcsB,pidcsB,degvalsB,nZB,np_all);

% Establish the basis of monomials defining the product
degmat_full = repelem(Adegmat,nZB,1) + repmat(Bdegmat,nZA,1);
[Pmat,degmat] = uniquerows_integerTable(degmat_full);    % Pmat*degmat = degmat_full
nZ = size(degmat,1);

% Deal with case where one object is zero
if nnz(C_A)==0 || nnz(C_B)==0
    % Product with 0 will always be 0
    szC = [mA,nA].*([mB,nB].^(otype-1));
    if use_dpvar
        C_C = spalloc(szC(1)*(ndvars+1),szC(2)*nZ,0);
        C = dpvar(C_C,degmat,pvarnames,dvarname,szC);
    else
        C_C = spalloc(nZ,szC(1)*szC(2),0);
        C = polynomial(C_C,degmat,pvarnames,szC);
    end
    return
end

% Establish row index, column index, dvarname and monomial associated with
% each of the coefficients defining each object
if isa(A,'dpvar')
    % A is dpvar
    [rridcsA,ccidcsA,valsA] = find(C_A);        nnzA = numel(valsA);
    ridcsA = ceil(rridcsA/(ndvars+1));
    didcsA = rridcsA - (ridcsA-1)*(ndvars+1);
    cidcsA = ceil(ccidcsA/nZA);
    ZidcsA = ccidcsA - (cidcsA-1)*nZA;
else
    % A is polynomial
    [ZidcsA,ccidcsA,valsA] = find(C_A);         nnzA = numel(valsA);
    cidcsA = ceil(ccidcsA/mA);
    ridcsA = ccidcsA - (cidcsA-1)*mA;
end
if isa(B,'dpvar')
    % B is dpvar
    [rridcsB,ccidcsB,valsB] = find(C_B);        nnzB = numel(valsB);
    ridcsB = ceil(rridcsB/(ndvars+1));
    didcsB = rridcsB - (ridcsB-1)*(ndvars+1);
    cidcsB = ceil(ccidcsB/nZB);
    ZidcsB = ccidcsB - (cidcsB-1)*nZB;
else
    % B is polynomial
    [ZidcsB,ccidcsB,valsB] = find(C_B);         nnzB = numel(valsB);
    cidcsB = ceil(ccidcsB/mB);
    ridcsB = ccidcsB - (cidcsB-1)*mB;
end

% Take all combinations of row and column indices
ridcs1 = repelem(ridcsA(:),nnzB,1);
ridcs2 = repmat(ridcsB(:),nnzA,1);
ridcs = (ridcs1-1)*mB + ridcs2;
cidcs1 = repelem(cidcsA(:),nnzB,1);
cidcs2 = repmat(cidcsB(:),nnzA,1);
cidcs = (cidcs1-1)*nB + cidcs2;

% Determine associated monomials, decision variables, and coefficients
% in each of the arrays
Zidcs1 = repelem(ZidcsA(:),nnzB,1);
Zidcs2 = repmat(ZidcsB(:),nnzA,1);
if isa(A,'dpvar')
    didcs = repelem(didcsA(:),nnzB,1);
elseif isa(B,'dpvar')
    didcs = repmat(didcsB(:),nnzA,1);
end
vals1 = repelem(valsA(:),nnzB,1);
vals2 = repmat(valsB(:),nnzA,1);

% Account for potential elementwise products
if otype(1)==1
    % Elementwise product along row dimension
    % --> only retain elements with shared rows
    rtn_idcs = ridcs1==ridcs2;
    ridcs = ridcs1(rtn_idcs);       cidcs = cidcs(rtn_idcs);
    cidcs1 = cidcs1(rtn_idcs);      cidcs2 = cidcs2(rtn_idcs);
    Zidcs1 = Zidcs1(rtn_idcs);      Zidcs2 = Zidcs2(rtn_idcs);
    vals1 = vals1(rtn_idcs);        vals2 = vals2(rtn_idcs);
    if use_dpvar
        didcs = didcs(rtn_idcs);
    end
end
if otype(2)==1
    % Elementwise product along column dimension
    % --> only retain elements with shared columns
    rtn_idcs = cidcs1==cidcs2;
    ridcs = ridcs(rtn_idcs);        cidcs = cidcs1(rtn_idcs);
    Zidcs1 = Zidcs1(rtn_idcs);      Zidcs2 = Zidcs2(rtn_idcs);
    vals1 = vals1(rtn_idcs);        vals2 = vals2(rtn_idcs);
    if use_dpvar
        didcs = didcs(rtn_idcs);
    end
end

% Determine the new monomials associated with each element
Zidcs_full = (Zidcs1-1)*nZB + Zidcs2;
cidcsP = Pmat*(1:nZ)';
Zidcs = cidcsP(Zidcs_full);

% Determine the new coefficients
vals = vals1.*vals2;

% Build the new coefficient matrix and declare the product
szC = [mA,nA].*([mB,nB].^(otype-1));
if use_dpvar
    % Output dpvar object
    rridcs = (ridcs-1)*(ndvars+1) + didcs;
    ccidcs = (cidcs-1)*nZ + Zidcs;
    C_C = sparse(rridcs,ccidcs,vals,szC(1)*(ndvars+1),szC(2)*nZ);
    C = dpvar(C_C,degmat,pvarnames,dvarname,szC);
else
    % Output polynomial object
    ccidcs = (cidcs-1)*szC(1) + ridcs;
    C_C = sparse(Zidcs,ccidcs,vals,nZ,szC(1)*szC(2));
    C = polynomial(C_C,degmat,pvarnames,szC);
end

end


%% quadPoly Implementation
function C = otimes_quadPoly(A,B,otype)
% % Compute the otimes product for 'quadPoly' objects
% A = (Im o ZL)' C (In o ZR)
% B = (Ip o ZL)' D (Iq o ZR)

if isa(A,'double')
    A = quadPoly(A);
elseif ~isa(A,'quadPoly')
    error("Both parameters must be specified as objects of type 'quadPoly'.")
end
if isa(B,'double')
    B = quadPoly(B);
elseif ~isa(B,'quadPoly')
    error("Both parameters must be specified as objects of type 'quadPoly'.")
end

% Express polynomials in terms of the same variables
if ~isequal(A.ns,B.ns) || ~isequal(A.nt,B.nt)
    [A,B] = common_vars(A,B);
end

% Augment monomial basis with constant monomial for vars not appearing 
% in the object
ZL_A = A.Zs;        ZL_B = B.Zs;    M = numel(ZL_A);
ZR_A = A.Zt;        ZR_B = B.Zt;    N = numel(ZR_A);

% Combine bases to represent all products of monomials
ZL = cell(1,M);        ZL_maps = cell(1,M);
nZL_A = zeros(1,M);    nZL_B = zeros(1,M);    nZL = ones(1,M);
for i=1:M
    % Determine the number of monomials
    nZA_i = size(ZL_A{i},1);            nZL_A(i) = nZA_i;
    nZB_i = size(ZL_B{i},1);            nZL_B(i) = nZB_i;
    % Take all combinations of monomials in A and B
    ZL_i = repelem(ZL_A{i}(:),nZB_i,1) + repmat(ZL_B{i}(:),nZA_i,1);
    % Retain only unique list
    [ZL_i,~,idcs_i] = unique(ZL_i);
    ZL{i} = ZL_i;
    nZL(i) = numel(ZL_i);
    ZL_maps{i} = idcs_i;
end
nnZL = cumprod([nZL,1],'reverse');
nnZL_A = cumprod([nZL_A,1],'reverse');
nnZL_B = cumprod([nZL_B,1],'reverse');
ZR = cell(1,N);        ZR_maps = cell(1,N);
nZR_A = zeros(1,N);    nZR_B = zeros(1,N);    nZR = ones(1,N);
for i=1:N
    % Determine the number of monomials
    nZA_i = size(ZR_A{i},1);            nZR_A(i) = nZA_i;
    nZB_i = size(ZR_B{i},1);            nZR_B(i) = nZB_i;
    % Take all combinations of monomials in A and B
    ZR_i = repelem(ZR_A{i}(:),nZB_i,1) + repmat(ZR_B{i}(:),nZA_i,1);
    % Retain only unique list
    [ZR_i,~,idcs_i] = unique(ZR_i);
    ZR{i} = ZR_i;
    nZR(i) = numel(ZR_i);
    ZR_maps{i} = idcs_i;
end
nnZR = cumprod([nZR,1],'reverse');
nnZR_A = cumprod([nZR_A,1],'reverse');
nnZR_B = cumprod([nZR_B,1],'reverse');


% Build the new coefficient matrix
% Establish row index, column index, left monomial and right monomial 
% associated with each of the coefficients defining each polynomial
szA = size(A);
[rridcsA,ccidcsA,valsA] = find(A.C);        nnzA = numel(valsA);
ridcsA = ceil(rridcsA/nnZL_A(1));
ZLidcsA = rridcsA - (ridcsA-1)*nnZL_A(1);
cidcsA = ceil(ccidcsA/nnZR_A(1));
ZRidcsA = ccidcsA - (cidcsA-1)*nnZR_A(1);

szB = size(B);
[rridcsB,ccidcsB,valsB] = find(B.C);        nnzB = numel(valsB);
ridcsB = ceil(rridcsB/nnZL_B(1));
ZLidcsB = rridcsB - (ridcsB-1)*nnZL_B(1);
cidcsB = ceil(ccidcsB/nnZR_B(1));
ZRidcsB = ccidcsB - (cidcsB-1)*nnZR_B(1);

% Take all combinations of row and column indices
ridcs1 = repelem(ridcsA(:),nnzB,1);
ridcs2 = repmat(ridcsB(:),nnzA,1);
ridcs = (ridcs1-1)*szB(1) + ridcs2;
cidcs1 = repelem(cidcsA(:),nnzB,1);
cidcs2 = repmat(cidcsB(:),nnzA,1);
cidcs = (cidcs1-1)*szB(2) + cidcs2;

% Determine associated monomials, decision variables, and coefficients
% in each of the arrays
ZRidcs1 = repelem(ZRidcsA(:),nnzB,1);
ZRidcs2 = repmat(ZRidcsB(:),nnzA,1);
ZLidcs1 = repelem(ZLidcsA(:),nnzB,1);
ZLidcs2 = repmat(ZLidcsB(:),nnzA,1);
vals1 = repelem(valsA(:),nnzB,1);
vals2 = repmat(valsB(:),nnzA,1);

% Account for potential elementwise products
if otype(1)==1
    % Elementwise product along row dimension
    % --> only retain elements with shared rows
    rtn_idcs = ridcs1==ridcs2;
    ridcs = ridcs1(rtn_idcs);       cidcs = cidcs(rtn_idcs);
    cidcs1 = cidcs1(rtn_idcs);      cidcs2 = cidcs2(rtn_idcs);
    ZLidcs1 = ZLidcs1(rtn_idcs);    ZLidcs2 = ZLidcs2(rtn_idcs);
    ZRidcs1 = ZRidcs1(rtn_idcs);    ZRidcs2 = ZRidcs2(rtn_idcs);
    vals1 = vals1(rtn_idcs);        vals2 = vals2(rtn_idcs);
end
if otype(2)==1
    % Elementwise product along column dimension
    % --> only retain elements with shared columns
    rtn_idcs = cidcs1==cidcs2;
    ridcs = ridcs(rtn_idcs);        cidcs = cidcs1(rtn_idcs);
    ZLidcs1 = ZLidcs1(rtn_idcs);    ZLidcs2 = ZLidcs2(rtn_idcs);
    ZRidcs1 = ZRidcs1(rtn_idcs);    ZRidcs2 = ZRidcs2(rtn_idcs);
    vals1 = vals1(rtn_idcs);        vals2 = vals2(rtn_idcs);
end

% Determine associated monomials in individual variables
ZLidcs = ones(numel(ZLidcs1),1);
for i=1:numel(nnZL)-1
    % Determine monomial indices in A and B associated with variable i
    ZLidcsA_i = ceil(ZLidcs1/nnZL_A(i+1));
    ZLidcsB_i = ceil(ZLidcs2/nnZL_B(i+1));
    % Determine associated monomial index in variable i in new ZL
    ZLidcs_i = (ZLidcsA_i-1)*nZL_B(i) + ZLidcsB_i;
    ZLidcs_i = ZL_maps{i}(ZLidcs_i);
    ZLidcs = ZLidcs + (ZLidcs_i-1)*nnZL(i+1);
    % Remove contribution of this monomial from the index
    ZLidcs1 = ZLidcs1 - (ZLidcsA_i-1)*nnZL_A(i+1);
    ZLidcs2 = ZLidcs2 - (ZLidcsB_i-1)*nnZL_B(i+1);    
end
ZRidcs = ones(numel(ZRidcs1),1);
for i=1:numel(nnZR)-1
    % Determine monomial indices in A and B associated with variable i
    ZRidcsA_i = ceil(ZRidcs1/nnZR_A(i+1));
    ZRidcsB_i = ceil(ZRidcs2/nnZR_B(i+1));
    % Determine associated monomial index in variable i in new ZL
    ZRidcs_i = (ZRidcsA_i-1)*nZR_B(i) + ZRidcsB_i;
    ZRidcs_i = ZR_maps{i}(ZRidcs_i);
    ZRidcs = ZRidcs + (ZRidcs_i-1)*nnZR(i+1);
    % Remove contribution of this monomial from the index
    ZRidcs1 = ZRidcs1 - (ZRidcsA_i-1)*nnZR_A(i+1);
    ZRidcs2 = ZRidcs2 - (ZRidcsB_i-1)*nnZR_B(i+1);    
end

% Determine the new coefficients
vals = vals1.*vals2;

% Build the new coefficient matrix and declare the product
szC = szA.*(szB.^(otype-1));
rridcs = (ridcs-1)*nnZL(1) + ZLidcs;
ccidcs = (cidcs-1)*nnZR(1) + ZRidcs;
C_C = sparse(rridcs,ccidcs,vals,szC(1)*nnZL(1),szC(2)*nnZR(1));

C = quadPoly(C_C,ZL,ZR,szC,A.ns,A.nt,true);

end



%% cell of quadPoly objects implementation
function C = otimes_quadPolyCell(A,B,otype)
% % Compute the otimes product for 'quadPoly' objects, where the
% coefficients are specified as cells
% A = (Im o ZL)' C (In o ZR)
% B = (Ip o ZL)' D (Iq o ZR)

if isa(A,'double')
    A = quadPoly(A);
elseif ~isa(A,'quadPoly') && ~isa(A,'struct')
    error("Both parameters must be specified as objects of type 'quadPoly'.")
end
if isa(B,'double')
    B = quadPoly(B);
elseif ~isa(B,'quadPoly') && ~isa(A,'struct')
    error("Both parameters must be specified as objects of type 'quadPoly'.")
end

% Express polynomials in terms of the same variables
if ~isequal(A.ns,B.ns) || ~isequal(A.nt,B.nt)
    [A,B] = common_vars_quadPolyCell(A,B);
end

% Augment monomial basis with constant monomial for vars not appearing 
% in the object
ZL_A = A.Zs;        ZL_B = B.Zs;    M = numel(ZL_A);
ZR_A = A.Zt;        ZR_B = B.Zt;    N = numel(ZR_A);

% Combine bases to represent all products of monomials
ZL = cell(1,M);        ZL_maps = cell(1,M);     ZL_map = 1;
nZL_A = zeros(1,M);    nZL_B = zeros(1,M);      nZL = ones(1,M);
for i=1:M
    % Determine the number of monomials
    nZA_i = size(ZL_A{i},1);            nZL_A(i) = nZA_i;
    nZB_i = size(ZL_B{i},1);            nZL_B(i) = nZB_i;
    % Take all combinations of monomials in A and B
    ZL_i = repelem(ZL_A{i}(:),nZB_i,1) + repmat(ZL_B{i}(:),nZA_i,1);
    % Retain only unique list
    [ZL_i,~,idcs_i] = unique(ZL_i);
    ZL{i} = ZL_i;
    nZL(i) = numel(ZL_i);
    ZL_maps{i} = idcs_i;
    ZL_map = (repelem(ZL_map,numel(idcs_i),1)-1)*nZL(i) + repmat(idcs_i(:),numel(ZL_map),1);
end
nnZL = cumprod([nZL,1],'reverse');
nnZL_A = cumprod([nZL_A,1],'reverse');
nnZL_B = cumprod([nZL_B,1],'reverse');
nnZL_tot = nnZL_A.*nnZL_B;
ZR = cell(1,N);        ZR_maps = cell(1,N);     ZR_map = 1;
nZR_A = zeros(1,N);    nZR_B = zeros(1,N);      nZR = ones(1,N);
for i=1:N
    % Determine the number of monomials
    nZA_i = size(ZR_A{i},1);            nZR_A(i) = nZA_i;
    nZB_i = size(ZR_B{i},1);            nZR_B(i) = nZB_i;
    % Take all combinations of monomials in A and B
    ZR_i = repelem(ZR_A{i}(:),nZB_i,1) + repmat(ZR_B{i}(:),nZA_i,1);
    % Retain only unique list
    [ZR_i,~,idcs_i] = unique(ZR_i);
    ZR{i} = ZR_i;
    nZR(i) = numel(ZR_i);
    ZR_maps{i} = idcs_i;
    ZR_map = (repelem(ZR_map,numel(idcs_i),1)-1)*nZR(i) + repmat(idcs_i(:),numel(ZR_map),1);
end
nnZR = cumprod([nZR,1],'reverse');
nnZR_A = cumprod([nZR_A,1],'reverse');
nnZR_B = cumprod([nZR_B,1],'reverse');
nnZR_tot = nnZR_A.*nnZR_B;


% Extract the coefficient matrices
A_C = A.C;          szA = A.dim;
if isnumeric(A_C)
    A_C = {A_C};
end
B_C = B.C;          szB = B.dim;
if isnumeric(B_C)
    B_C = {B_C};
end

% Initialize an empty cell of coefficients for all products of coefficients
% in A and in B
% A_C, B_C are N3-D cells, but ndims ignores trailing singleton dimensions
% (factor not depending on the last shared variable) --> pad with ones
Ndims = max(ndims(A_C),ndims(B_C));
szA_C = size(A_C,1:Ndims);
szB_C = size(B_C,1:Ndims);
szC = szA.*(szB.^(otype-1));
C_C = repmat({spalloc(szC(1)*nnZL(1),szC(2)*nnZR(1),0)},szA_C.*szB_C);

% Establish row index, column index, left monomial and right monomial 
% associated with each of the coefficients defining both polynomials
rsubs_A = cell(szA_C);      csubs_A = cell(szA_C);
vals_A = cell(szA_C);
for k=1:numel(A_C)
    [rsubs_k,csubs_k,vals_k] = find_idcs(A_C{k},nnZL_A,nnZR_A);
    rsubs_A{k} = rsubs_k;
    csubs_A{k} = csubs_k;
    vals_A{k} = vals_k;
end
rsubs_B = cell(szB_C);      csubs_B = cell(szB_C);
vals_B = cell(szB_C);
for k=1:numel(B_C)
    [rsubs_k,csubs_k,vals_k] = find_idcs(B_C{k},nnZL_B,nnZR_B);
    rsubs_B{k} = rsubs_k;
    csubs_B{k} = csubs_k;
    vals_B{k} = vals_k;
end

% Take the product of all parameters
for i=1:numel(A_C)
    % Determine which index in A_C corresponds with parameter i
    idcsA = cell(1,Ndims);
    [idcsA{:}] = ind2sub(szA_C,i);
    idcsA = cell2mat(idcsA);

    % Extract indices and coefficients associated with parameter i
    rsubs_i = rsubs_A{i};
    csubs_i = csubs_A{i};
    vals_i = vals_A{i};
    nnzA = numel(vals_i);
    if nnzA==0
        continue
    end

    for j=1:numel(B_C)
        % Determine which index in B_C corresponds with parameter j
        idcsB = cell(1,Ndims);
        [idcsB{:}] = ind2sub(szB_C,j);
        idcsB = cell2mat(idcsB);

        % Determine which element in C_C corresponds with parameter 
        % pair (i,j)
        idcsC = (idcsA-1).*szB_C + idcsB;
        idcsC = num2cell(idcsC);

        % Extract indices and coefficients associated with parameter j
        rsubs_j = rsubs_B{j};
        csubs_j = csubs_B{j};
        vals_j = vals_B{j};
        nnzB = numel(vals_j);
        if nnzB==0
            continue
        end

        % Take all combinations of row and column indices in C_A and C_B
        if all(otype==2)
            % Proper tensor product along both row and column dimension
            % --> have to consider all combinations of rows and columns
            rsubs1 = repelem(rsubs_i,nnzB,1);
            rsubs2 = repmat(rsubs_j,nnzA,1);
            csubs1 = repelem(csubs_i,nnzB,1);
            csubs2 = repmat(csubs_j,nnzA,1);
            vals1 = repelem(vals_i(:),nnzB,1);
            vals2 = repmat(vals_j(:),nnzA,1);
    
            % Deteremine row and column indices in new coefficients
            ridcs = (rsubs1(:,1)-1)*szB(1) + rsubs2(:,1);
            cidcs = (csubs1(:,1)-1)*szB(2) + csubs2(:,1);
        elseif all(otype==[1,2])
            % Elementwise product along row dimension
            % --> only retain elements with shared rows
            [rtn_idcs2,rtn_idcs1] = find_common_idcs(rsubs_j(:,1),rsubs_i(:,1));
            rsubs1 = rsubs_i(rtn_idcs1,:);   rsubs2 = rsubs_j(rtn_idcs2,:);
            csubs1 = csubs_i(rtn_idcs1,:);   csubs2 = csubs_j(rtn_idcs2,:);
            vals1 = vals_i(rtn_idcs1);       vals2 = vals_j(rtn_idcs2,:);                       

            % Deteremine row and column indices in new coefficients
            ridcs = rsubs1(:,1);
            cidcs = (csubs1(:,1)-1)*szB(2) + csubs2(:,1);
        elseif all(otype==[2,1])
            % Elementwise product along column dimension
            % --> only retain elements with shared columns
            [rtn_idcs2,rtn_idcs1] = find_common_idcs(csubs_j(:,1),csubs_i(:,1));
            rsubs1 = rsubs_i(rtn_idcs1,:);   rsubs2 = rsubs_j(rtn_idcs2,:);
            csubs1 = csubs_i(rtn_idcs1,:);   csubs2 = csubs_j(rtn_idcs2,:);
            vals1 = vals_i(rtn_idcs1);       vals2 = vals_j(rtn_idcs2,:);      

            % Deteremine row and column indices in new coefficients
            cidcs = csubs1(:,1);
            ridcs = (rsubs1(:,1)-1)*szB(1) + rsubs2(:,1);

        else
            % Elementwise product along both row and column dimensions
            % --> only retain elements with shared row and column index
            [rtn_idcs2,rtn_idcs1] = find(rsubs_j(:,1)==rsubs_i(:,1)' & csubs_j(:,1)==csubs_i(:,1)');
            rsubs1 = rsubs_i(rtn_idcs1,:);    rsubs2 = rsubs_j(rtn_idcs2,:);
            csubs1 = csubs_i(rtn_idcs1,:);    csubs2 = csubs_j(rtn_idcs2,:);
            vals1 = vals_i(rtn_idcs1);        vals2 = vals_j(rtn_idcs2,:);      

            % Deteremine row and column indices in new coefficients
            ridcs = rsubs1(:,1);
            cidcs = csubs1(:,1);
        end

        % Determine associated monomial indices in new ZL
        ZLsubs = (rsubs1(:,2:end)-1).*nZL_B + rsubs2(:,2:end);
        ZLidcs = (ZLsubs-1)*nnZL_tot(2:end)' + 1;
        ZLidcs = ZL_map(ZLidcs);        % account for uniqueness of monomials
        
        % Determine associated monomial indices in new ZR
        ZRsubs = (csubs1(:,2:end)-1).*nZR_B + csubs2(:,2:end);
        ZRidcs = (ZRsubs-1)*nnZR_tot(2:end)' + 1;
        ZRidcs = ZR_map(ZRidcs);        % account for uniqueness of monomials

        % Determine the new coefficients
        vals = vals1.*vals2;
        
        % Build the new coefficient matrix and declare the product
        rridcs = (ridcs-1)*nnZL(1) + ZLidcs;
        ccidcs = (cidcs-1)*nnZR(1) + ZRidcs;
        
        C_C{idcsC{:}} = sparse(rridcs,ccidcs,vals,szC(1)*nnZL(1),szC(2)*nnZR(1));
    end
end

% Declare the structure representing the tensor product
C = struct();
C.C = C_C;
C.Zs = ZL;
C.Zt = ZR;
C.dim = szC;
C.ns = A.ns;
C.nt = A.nt;

end



%% cell of parameters implementation
function C = otimes_cell(A,B,otype)
% Each element of the cells is assumed to represent a term in an operator,
% so that the product of the operators (independent of otype) should
% involve the product of each pair of terms in A and B.

% Assume a single parameter if object is not a cell
if ~isa(A,'cell')
    A = {A};
elseif ~isa(B,'cell')
    B = {B};
end

% Initialize a cell of the same dimensionality as A and B
if ndims(A)~=ndims(B)
    error("Dimensionality of the cells should be the same.")
end
N = ndims(A);
szA = size(A,1:N);      szB = size(B,1:N);
szC = szA.*szB;
C = cell(szC);
for k=1:numel(C)        % can be parallelized
    % Determine which index in the ND cell corresponds with parameter k
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(szC,k);
    idcs = cell2mat(idcs);
    % Determine associated index in A and B
    idcs_A = ceil(idcs./szB);
    idcs_B = idcs-(idcs_A-1).*szB;
    idcs_A = num2cell(idcs_A);
    idcs_B = num2cell(idcs_B);
    % Take the product of the appropriate parameters
    C{k} = otimes_params(A{idcs_A{:}},B{idcs_B{:}},otype);    
end

end




%%
%%
function [rsubs,csubs,vals] = find_idcs(C,nnr,nnc)
% Find indices in monomial bases associated with each nonzero coefficient
% in the matrix C
%
% INPUTS
% - C:  (sparse) array of dimensions prod(nr) x prod(nc), representing
%       coefficient matrix acting on monomials 
%           ZL = ZL{1} o ... o ZL{m}    and     ZR = ZR{1} o ... o ZR{n}
%       for o the Kronekcer product
% - nr: 1 x m array specifying the cumulative number of elements in each 
%       monomial basis ZL{i};
% - nc: 1 x n array specifying the cumulative number of elements in each 
%       monomial basis ZR{i};
%
% OUTPUTS
% - rsubs:  k x m array specifying the monomial index in ZL associated with
%           each nonzero element of C;
% - csubs:  k x m array specifying the monomial index in ZL associated with
%           each nonzero element of C;
% - vals:   k x 1 array specifying the nonzero values appearing in C
%

% Find nonzero elements of C
[rridcs,ccidcs,vals] = find(C);
vals = vals(:);

% Split row indices into indices for different monomials
M = numel(nnr);
rsubs = zeros(numel(rridcs),M);
rsubs(:,1) = ceil(rridcs(:)/nnr(1));
for i=1:M-1
    % Remove contribution of monomials in variable i
    rridcs = rridcs(:) - (rsubs(:,i)-1)*nnr(i);
    % Set indices associated with variable i+1
    rsubs(:,i+1) = ceil(rridcs/nnr(i+1));
end
% Split column indices into indices for different monomials
N = numel(nnc);
csubs = zeros(numel(ccidcs),N);
csubs(:,1) = ceil(ccidcs(:)/nnc(1));
for i=1:N-1
    % Remove contribution of monomials in variable i
    ccidcs = ccidcs(:) - (csubs(:,i)-1)*nnc(i);
    % Set indices associated with variable i+1
    csubs(:,i+1) = ceil(ccidcs/nnc(i+1));
end

end

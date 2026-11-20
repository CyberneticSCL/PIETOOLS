function [A,B] = unpack_sheets(C,mrow,ncol,nsheet,L,R,rloc,nrow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A,B] = UNPACK_SHEETS(C,MROW,NCOL,NSHEET,L,R,RLOC,NROW) returns the
% coefficients of L*C(d)*R for the matrix C(d) = C_0 + sum_s d_s C_s whose
% 1+NSHEET coefficient sheets C_0, C_1, ..., C_NSHEET, each MROW x NCOL,
% are stored side by side in C:
%
%       A = vec(L*C_0*R),       B(RLOC(s),:) = vec(L*C_s*R).',
%
% s = 1, ..., NSHEET, so that vec(L*C(d)*R) = A + B'*dd for the NROW
% decision variables dd with dd(RLOC(s)) = d_s. Rows of B that no sheet
% maps to are zero; the sheets mapped to one row add.
%
% INPUT
% - C:      'struct' with fields i, j, v (double columns) and m, n
%           (double): the triplets of sparse(i,j,v,m,n) = [C_0, C_1, ...,
%           C_nsheet], m = mrow and n = (1+nsheet)*ncol, sheet s in the
%           columns s*ncol+(1:ncol), with no repeated (i,j) and no zero v;
% - mrow:   'double', the number of rows of a sheet;
% - ncol:   'double', the number of columns of a sheet;
% - nsheet: 'double', the number of sheets after C_0;
% - L:      'double', p x mrow, sparse or full;
% - R:      'double', ncol x r, sparse or full;
% - rloc:   'double', nsheet entries, the row of B of each sheet;
% - nrow:   'double', the number of rows of B.
%
% OUTPUT
% - A:      'double', sparse, (p*r) x 1;
% - B:      'double', sparse, nrow x (p*r).
%
% NOTES
% C other than a 'struct', fewer than eight inputs, and C.m or C.n other
% than mrow or (1+nsheet)*ncol are errors.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - unpack_sheets
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
% Initial coding MMP, 08/22/2026

% Inputs: the triplets of C and all eight arguments; C of the stated size
if ~isstruct(C) || nargin<8
    error("unpack_sheets takes the triplets of int_semisep(...,'triplets') "...
          +"with L, R, RLOC and NROW.")
end
    if C.m~=mrow || C.n~=(1+nsheet)*ncol
        error("Internal error: unexpected coefficient dimensions.")
    end
    irow = C.i(:);  icol = C.j(:);  val = C.v(:);
% Sheet 0 (columns 1..ncol): the constant term, vectorized
is_A = icol<=ncol;

A = sparse((icol(is_A)-1)*mrow+irow(is_A),1,val(is_A),mrow*ncol,1);

% The other nonzeros: their sheet sg and column cc within the sheet
jcol = icol(~is_A)-ncol;
sg = floor((jcol-1)/ncol)+1;
cc = jcol - (sg-1)*ncol;

% The decision part: the product L*C_s*R and the scatter of the sheets to
% the rows RLOC, fused into one remap of the triplets where possible.
% Everything forced to columns: 'find' returns ROWS for a single-row input,
% and a row indexed into a column expands silently to a matrix.
irow = irow(:);     sg = sg(:);     cc = cc(:);     rloc = rloc(:);
vB = val(~is_A);    vB = vB(:);
% Nonzeros of L by column and of R by row; the length of vec(L*C_s*R)
[rL,cL,vL] = find(L);           % one per column: rL(cL) is its row
[cR,rR,vR] = find(R.');         % one per row of R: cR(rR) is its column
rL = rL(:);         cR = cR(:);
nout = size(L,1)*size(R,2);
% Fused route when L is 0/1 with one nonzero per column and R 0/1 with one
% per row: each coefficient of a sheet then moves to one output column kj
fused = size(L,2)==mrow && size(R,1)==ncol ...
        && isequal(cL(:),(1:mrow).') && isequal(rR(:),(1:ncol).') ...
        && all(vL==1) && all(vR==1);
if fused
    kj = rL(irow(~is_A)) + size(L,1)*(cR(cc)-1);
    % Exactness guard. B*K sums the coefficients landing on one entry in
    % ascending vec index; 'sparse' sums duplicates in an order of its own.
    % With at most TWO per entry the order is immaterial and unit weights
    % make the products exact; the count is taken on the data, since the
    % structural fan-in of the map reaches 1296.
    cnt = sparse(sg,kj,1,nsheet,nout);
    fused = all(nonzeros(cnt)<=2);
end
% A: vec(L*C*R) of the constant sheet, the mr x nc matrix multiplied and
% vectorized; B from the moved triplets (fused), or by the product
% B*kron(R.',L).' and a scatter of its rows to RLOC
X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;
A = Y(:);
if fused
    B = sparse(rloc(sg),kj,vB,nrow,nout);
else
    % Unfused: the product with kron(R.',L).', then the row scatter
    B = sparse(sg,(cc-1)*mrow+irow(~is_A),vB,nsheet,mrow*ncol);
    B = B*kron(R.',L).';
    [bi,bj,bv] = find(B);
    B = sparse(rloc(bi(:)),bj(:),bv(:),nrow,size(B,2));
end

end



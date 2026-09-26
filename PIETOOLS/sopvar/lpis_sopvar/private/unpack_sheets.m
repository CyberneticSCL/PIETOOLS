function [A,B] = unpack_sheets(C,mrow,ncol,nsheet,L,R,rloc,nrow)
% function [A,B] = unpack_sheets(C,mrow,ncol,nsheet)                        % MMP, 09/26/2026 (was)
% Inverse of 'pack_sheets' applied to the output of 'int_semisep': the
% column index of the result is (sheet,ncol,NR) with the semiseparable index
% fastest, so each sheet still occupies a contiguous block of ncol columns.
% Sheet 1 is returned as the vectorized A, the remaining nsheet sheets as
% the rows of B.
%
% The whole array is scanned once and then partitioned by column, rather
% than sliced first. Slicing a sparse matrix whose column count carries the
% decision variable dimension, as C(:,ncol+1:end) does, copies that block
% before anything useful happens, and the copy dominates once there are many
% decision variables.
%
% Optional L, R, RLOC, NROW: return instead the coefficients of L*C(d)*R    % MMP, 09/26/2026
% with sheet s moved to row RLOC(s) of NROW, i.e. A and B of                % MMP, 09/26/2026
% [A,B] = lr_multiply(L,A,B,R) followed by that row scatter. L and R from   % MMP, 09/26/2026
% 'merge_monomial_product' are 0/1 with one nonzero per column of L and per % MMP, 09/26/2026
% row of R, so B*kron(R.',L).' just moves coefficient (irow,cc) of a sheet  % MMP, 09/26/2026
% to column (colR(cc)-1)*size(L,1)+rowL(irow): one 'sparse' from the        % MMP, 09/26/2026
% triplets of C replaces B, the Kronecker product, B*K, 'find' and a second % MMP, 09/26/2026
% 'sparse'. Anything else takes the unfused route.                          % MMP, 09/26/2026
%
% C may instead be the triplet struct of int_semisep(...,'triplets'),       % MMP, 09/26/2026
% fields i, j, v, m, n, standing for sparse(i,j,v,m,n) with no repeated     % MMP, 09/26/2026
% (i,j) and no zero v. Nothing is then built or scanned over the n          % MMP, 09/26/2026
% columns, which carry the decision variables.                              % MMP, 09/26/2026
%
% MMP, 09/26/2026: Triplet input for C (above), from 'copquadvar' via
%                  int_semisep(...,'triplets'). The 'find' below was 3.6 s
%                  of a 19.7 s 2-D Hinf build (io2, light, 1146 calls),
%                  O(n) per call. Exact: every 'sparse' here depends only
%                  on the SET of triplets, not their order -- positions are
%                  distinct except in the fused B, which the guard limits
%                  to two per entry, and x1+x2 == x2+x1 -- and the triplets
%                  are the set 'find' returned. Matrix input is unchanged.
% MMP, 09/26/2026: Optional inputs L, R, RLOC, NROW, as above, for
%                  'copquadvar'. Profiled A/B, 2-D Hinf build (io2, light,
%                  1146 calls): unpack + lr_multiply + kron 5.7 + 1.3 +
%                  1.6 s -> 4.6 + 0 + 0.1 s, of which 3.5 s is the 'find'
%                  below: its cost is the column count of C (1.35e7 per
%                  call for ~3900 nonzeros), set by 'int_semisep'. Memory
%                  per call drops from O(mrow*ncol) (B's column pointers,
%                  the Kronecker product) to O(nnz(C)): peak 209 -> 5 MB
%                  here and 616 -> 5 MB in kron on 2-D heavy stability.
%                  Bit-identical: see the guard below. Without the new
%                  inputs the behaviour is unchanged ('sopquadvar').

if isstruct(C)                  % triplets: the 'find' is already done      % MMP, 09/26/2026
    if C.m~=mrow || C.n~=(1+nsheet)*ncol                                    % MMP, 09/26/2026
        error("Internal error: unexpected coefficient dimensions.")         % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
    irow = C.i(:);  icol = C.j(:);  val = C.v(:);                           % MMP, 09/26/2026
else                                                                        % MMP, 09/26/2026
if size(C,1)~=mrow || size(C,2)~=(1+nsheet)*ncol
    error("Internal error: unexpected coefficient dimensions.")
end

[irow,icol,val] = find(C);
end                                                                         % MMP, 09/26/2026
is_A = icol<=ncol;

A = sparse((icol(is_A)-1)*mrow+irow(is_A),1,val(is_A),mrow*ncol,1);

jcol = icol(~is_A)-ncol;
sg = floor((jcol-1)/ncol)+1;
cc = jcol - (sg-1)*ncol;
% B = sparse(sg,(cc-1)*mrow+irow(~is_A),val(~is_A),nsheet,mrow*ncol);       % MMP, 09/26/2026 (was)
if nargin<5                                                                 % MMP, 09/26/2026
    B = sparse(sg,(cc-1)*mrow+irow(~is_A),val(~is_A),nsheet,mrow*ncol);     % MMP, 09/26/2026
    return                                                                  % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026

% BEGIN MMP, 09/26/2026: fused lr_multiply and row scatter.
% Everything forced to columns: 'find' returns ROWS for a single-row input,
% and a row indexed into a column expands silently to a matrix.
irow = irow(:);     sg = sg(:);     cc = cc(:);     rloc = rloc(:);         % MMP, 09/26/2026
vB = val(~is_A);    vB = vB(:);                                             % MMP, 09/26/2026
[rL,cL,vL] = find(L);           % one per column: rL(cL) is its row         % MMP, 09/26/2026
[cR,rR,vR] = find(R.');         % one per row of R: cR(rR) is its column    % MMP, 09/26/2026
rL = rL(:);         cR = cR(:);                                             % MMP, 09/26/2026
nout = size(L,1)*size(R,2);                                                 % MMP, 09/26/2026
fused = size(L,2)==mrow && size(R,1)==ncol ...
        && isequal(cL(:),(1:mrow).') && isequal(rR(:),(1:ncol).') ...
        && all(vL==1) && all(vR==1);                                        % MMP, 09/26/2026
if fused                                                                    % MMP, 09/26/2026
    kj = rL(irow(~is_A)) + size(L,1)*(cR(cc)-1);                            % MMP, 09/26/2026
    % Exactness guard. B*K sums the coefficients landing on one entry in
    % ascending vec index; 'sparse' sums duplicates in an order of its own
    % (measured to differ at 1e5 nonzeros). With at most TWO per entry the
    % order is immaterial, x1+x2 == x2+x1, and unit weights make the
    % products exact. Every cx_exec build measured, 1-D and 2-D, light and
    % heavy, has at most two, though the map's structural fan-in reaches
    % 1296, so the count is taken on the data.
    cnt = sparse(sg,kj,1,nsheet,nout);                                      % MMP, 09/26/2026
    fused = all(nonzeros(cnt)<=2);                                          % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
if fused                                                                    % MMP, 09/26/2026
    A = lr_multiply(L,A,[],R);                                              % MMP, 09/26/2026
    B = sparse(rloc(sg),kj,vB,nrow,nout);                                   % MMP, 09/26/2026
else                                                                        % MMP, 09/26/2026
    % Unfused: the sequence 'copquadvar' ran before this date.
    B = sparse(sg,(cc-1)*mrow+irow(~is_A),vB,nsheet,mrow*ncol);             % MMP, 09/26/2026
    [A,B] = lr_multiply(L,A,B,R);                                           % MMP, 09/26/2026
    [bi,bj,bv] = find(B);                                                   % MMP, 09/26/2026
    B = sparse(rloc(bi(:)),bj(:),bv(:),nrow,size(B,2));                     % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
% END MMP, 09/26/2026

end



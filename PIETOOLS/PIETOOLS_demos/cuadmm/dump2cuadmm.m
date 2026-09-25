function info = dump2cuadmm(dumpfile,outdir)
% dump2cuadmm(dumpfile,outdir) -- write one ladder rung as a cuADMM problem dir.
%
% cuADMM reads five text files from a directory:
%   blk.txt      one "<type> <size>" line per block; 's' symmetric, 'u' free,
%                'l' nonnegative
%   con_num.txt  the number of constraints, m
%   At.txt       COO, 0-based, "vecidx conidx val", rows ASCENDING
%   b.txt        COO, 0-based, "idx 0 val"
%   C.txt        COO, 0-based, "idx 0 val"
% At is stored as (vec_len x m): rows index the svec'd primal vector, columns
% index constraints.
%
% SVEC CONVENTION.  cuADMM uses SDPT3's: the upper triangle taken column by
% column, diagonal unscaled, off-diagonal scaled by sqrt(2), so that
% <A,X> = svec(A)'*svec(X) for symmetric A, X.
%
% WHY THE MAP IS BUILT HERE rather than via the read_sedumi/svecADMM chain
% cuADMM ships in examples/: that chain was measured to CORRUPT these rungs --
% on P_n01 it dropped nnz(At) 1489 -> 193, returned a non-trivial block
% permutation, and produced a system whose residual at a known-good Mosek point
% was 1.2e+01 against SeDuMi's 9.9e-08.  Its output was also completely
% insensitive to the off-diagonal factor, i.e. the off-diagonal couplings were
% gone entirely.  (read_sedumi there is SDPNAL+'s copy; it is not the scaling
% that is wrong, it mishandles this data.)
%
% The A_k that sossolve emits are NOT symmetric -- on P_n01,
% max||A_k - A_k'||_F = 6.4 against max||A_k||_F = 5.2.  That is harmless for
% SeDuMi, whose PSD inner product <A_k,X> = sum_ij A_k(ij) X(ij) sees only the
% symmetric part of A_k, but it means the conversion MUST symmetrise.  T below
% does the symmetrisation and the sqrt(2) scaling in one step: the off-diagonal
% row for (i,j) carries sqrt(2)/2 on BOTH vec entries (i,j) and (j,i), which
% maps an asymmetric A_k to sqrt(2)*sym(A_k)_ij and a symmetric X to
% sqrt(2)*X_ij, as required.

S = load(dumpfile);
At = S.At;  b = full(S.b(:));  c = S.c;  K = S.K;   % SeDuMi form, At is nvar x m
if ~isfield(K,'f') || isempty(K.f), K.f = 0; end
if ~isfield(K,'l') || isempty(K.l), K.l = 0; end
if K.l ~= 0
    error('dump2cuadmm: K.l=%d unsupported here (no ladder rung has one)',K.l);
end
Ks = double(K.s(:)');
m  = numel(b);
nvar = size(At,1);
if nvar ~= K.f + sum(Ks.^2)
    error('dump2cuadmm: nvar %d ~= K.f + sum(K.s.^2) = %d',nvar,K.f+sum(Ks.^2));
end

% ---- svec map T : (vec_len x nvar), block diagonal
vec_len = K.f + sum(Ks.*(Ks+1)/2);
[Ti,Tj,Tv] = svec_triplets(K.f,Ks);
T = sparse(Ti,Tj,Tv,vec_len,nvar);

cu_At = T*At;                       % (vec_len x m)
cu_C  = T*c(:);
cu_b  = b;

blk = cell(0,2);
if K.f > 0, blk(end+1,:) = {'u',K.f}; end
for i = 1:numel(Ks), blk(end+1,:) = {'s',Ks(i)}; end     %#ok<AGROW>

if ~exist(outdir,'dir'), mkdir(outdir); end
store_blk(blk,   fullfile(outdir,'blk.txt'));
store_coo(cu_At, fullfile(outdir,'At.txt'));
store_coo(cu_b,  fullfile(outdir,'b.txt'));
store_coo(cu_C,  fullfile(outdir,'C.txt'));
fid = fopen(fullfile(outdir,'con_num.txt'),'w');
fprintf(fid,'%d\n',m);  fclose(fid);

info = struct('dir',outdir,'vec_len',vec_len,'m',m,'blk',{blk}, ...
              'Ns',S.Ns,'Kf',S.Kf,'nvar',nvar,'nnz_At',nnz(cu_At), ...
              'normb',norm(cu_b),'normC',norm(cu_C));
save(fullfile(outdir,'meta.mat'),'-struct','info');
end


function [I,J,V] = svec_triplets(Kf,Ks)
% Triplets of the svec map, assembled once per block and concatenated.  Built
% from triplets rather than by growing a matrix: nvar reaches millions on the
% large rungs, where index-assignment into a sparse matrix is O(nnz) per write.
r2 = sqrt(2)/2;
I = {};  J = {};  V = {};
if Kf > 0
    p = (1:Kf)';
    I{end+1} = p;  J{end+1} = p;  V{end+1} = ones(Kf,1);
end
bs = Kf;  bv = Kf;                      % base offsets in svec and in vec
for k = 1:numel(Ks)
    N = Ks(k);
    idx = find(triu(true(N)));          % column-major upper triangle = svec order
    [ii,jj] = ind2sub([N N],idx);
    p  = (1:numel(idx))';
    dg = (ii==jj);
    od = ~dg;
    % upper entry (i,j) for every svec row, plus the mirrored (j,i) off-diagonals
    I{end+1} = [bs+p;        bs+p(od)];                              %#ok<AGROW>
    J{end+1} = [bv+idx;      bv+sub2ind([N N],jj(od),ii(od))];       %#ok<AGROW>
    V{end+1} = [dg + od*r2;  repmat(r2,nnz(od),1)];                  %#ok<AGROW>
    bs = bs + N*(N+1)/2;
    bv = bv + N^2;
end
I = cat(1,I{:});  J = cat(1,J{:});  V = cat(1,V{:});
end


function store_blk(blk,fn)
fid = fopen(fn,'w');
for i = 1:size(blk,1)
    t = blk{i,1};
    if ~any(strcmp(t,{'s','u','l'}))
        fclose(fid);  error('unsupported block type %s',t);
    end
    fprintf(fid,'%c %d\n',t,fix(blk{i,2}));
end
fclose(fid);
end


function store_coo(mat,fn)
% 0-based COO, sorted by row -- cuADMM's COO_to_CSC assumes ascending rows.
[r,cc,v] = find(sparse(mat));
[r,idx] = sort(r-1);
cc = cc(idx)-1;  v = v(idx);
fid = fopen(fn,'w');
fprintf(fid,'%d %d %.17g\n',[r(:) cc(:) v(:)]');
fclose(fid);
end

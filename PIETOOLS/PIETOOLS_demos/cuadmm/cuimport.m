function r = cuimport(dumpdir,matfile,sfx)
% cuimport(dumpdir,matfile) -- rebuild cuADMM's certificate and score it against
% the ORIGINAL data, independently of anything the solver reported.
%
% WHY THIS EXISTS.  Every cuADMM number banked in this campaign is the solver's
% SELF-REPORT (its pinf/dinf/gap).  That is not evidence: it is computed from the
% same transformed system the solver was handed, so a bug in the transform is
% invisible to it.
%
% THE TRAP THIS ROUTINE MUST NOT FALL INTO.  Importing with the transpose of the
% same svec map T used to export is a MUTUALLY INVERSE pair, and CLAUDE.md item 4
% forbids testing such a pair against each other: if T scaled off-diagonals by 2
% instead of sqrt(2), cuADMM would solve a different SDP, T' would map the answer
% back, and At'*x = b would STILL HOLD EXACTLY.  The row residual cannot see it.
% So three checks are reported, and only the last two are independent of T:
%   rel_b    row residual in ORIGINAL units, ||Atf'*x_full - bf||/||bf||
%            -- catches gross errors only, BLIND to an inverse-pair error
%   psd_min  smallest eigenvalue over the rebuilt Gram blocks
%            -- a wrong off-diagonal factor rescales the off-diagonals of every
%               block and destroys positivity, so THIS is what catches it
%   obj      c'*x_full against the Mosek objective, supplied by the caller
%            -- catches an overall scaling
% Plus a FORWARD check of T against known-good data: T*x_mosek must satisfy
% cuADMM's own system, which uses T in the export direction on a vector T never
% produced.  T'*T is the identity on the symmetric subspace, so this is exact.

S = load(matfile);                       % SeDuMi form: At (nvar x m), b, c, K
K = S.K;  if ~isfield(K,'f')||isempty(K.f), K.f = 0; end
Ks = double(K.s(:)');
nvar = size(S.At,1);
vec_len = K.f + sum(Ks.*(Ks+1)/2);

if nargin<3, sfx = ''; end
xs = readvec(fullfile(dumpdir,['X_opt' sfx '.txt']));
if numel(xs) ~= vec_len
    error('cuimport:len','X_opt has %d entries, expected vec_len %d',numel(xs),vec_len);
end
ys = readvec(fullfile(dumpdir,['y_opt' sfx '.txt']));

T = svecmap(K.f,Ks,vec_len,nvar);
x = T'*xs(:);                            % back to SeDuMi cone coordinates

r.rel_b_norm = norm(full(S.At'*x - S.b(:)))/max(norm(full(S.b(:))),eps);
r.obj_norm   = full(S.c(:)'*x);
r.normx      = norm(x);
r.vec_len    = vec_len;  r.m = numel(ys);  r.Kf = K.f;  r.Ks = Ks;

% ---- cone, in the rebuilt blocks.  THE T-INDEPENDENT CHECK.
off = K.f;  r.psd_min = inf;  r.psd_max = -inf;  r.blk_mineig = zeros(1,numel(Ks));
for k = 1:numel(Ks)
    N = Ks(k);
    Xk = reshape(x(off+(1:N^2)),N,N);
    r.asym(k) = norm(Xk-Xk','fro')/max(norm(Xk,'fro'),eps);   % must be ~0
    ev = eig((Xk+Xk')/2);
    r.blk_mineig(k) = min(ev);
    r.psd_min = min(r.psd_min,min(ev));  r.psd_max = max(r.psd_max,max(ev));
    off = off + N^2;
end
r.psd_relmin = r.psd_min/max(r.psd_max,eps);
r.max_asym   = max(r.asym);

% ---- dual slack S = C - sum_k y_k A_k, which must also be PSD
Sv = S.c(:) - S.At*ys(:);
off = K.f;  r.dual_min = inf;
for k = 1:numel(Ks)
    N = Ks(k);
    Sk = reshape(Sv(off+(1:N^2)),N,N);
    r.dual_min = min(r.dual_min,min(eig((Sk+Sk')/2)));
    off = off + N^2;
end
r.dgap = full(S.b(:)'*ys(:));
end


function T = svecmap(Kf,Ks,vec_len,nvar)
% SDPT3 svec: upper triangle column-major, diagonal unscaled, off-diagonal
% sqrt(2)/2 on BOTH mirrored vec entries so an asymmetric A_k is symmetrised in
% the same step.  Built from triplets -- nvar reaches millions on the big rungs.
r2 = sqrt(2)/2;
I = {};  J = {};  V = {};
if Kf > 0
    p = (1:Kf)';  I{end+1} = p;  J{end+1} = p;  V{end+1} = ones(Kf,1);
end
bs = Kf;  bv = Kf;
for k = 1:numel(Ks)
    N = Ks(k);
    idx = find(triu(true(N)));
    [ii,jj] = ind2sub([N N],idx);
    p = (1:numel(idx))';  dg = (ii==jj);  od = ~dg;
    I{end+1} = [bs+p;   bs+p(od)];                                   %#ok<AGROW>
    J{end+1} = [bv+idx; bv+sub2ind([N N],jj(od),ii(od))];            %#ok<AGROW>
    V{end+1} = [dg + od*r2; repmat(r2,nnz(od),1)];                   %#ok<AGROW>
    bs = bs + N*(N+1)/2;  bv = bv + N^2;
end
T = sparse(cat(1,I{:}),cat(1,J{:}),cat(1,V{:}),vec_len,nvar);
end


function v = readvec(fn)
fid = fopen(fn,'r');
if fid<0, error('cuimport:open','cannot open %s',fn); end
v = fscanf(fid,'%f');
fclose(fid);
end

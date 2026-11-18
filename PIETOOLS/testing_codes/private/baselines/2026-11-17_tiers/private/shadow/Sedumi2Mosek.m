function prob = Sedumi2Mosek(A,b,c,K) %#ok<STOUT>
% SEDUMI2MOSEK, scratch shadow 10/28/2026: records a fingerprint of the SDP
% that sossolve hands to MOSEK in the global HD_FP and stops with the error
% 'hdshadow:stop', so no solve runs.
% INPUT  A, b, c, K: the SeDuMi data, A m x n (constraints x variables).
% Fingerprint: sizes, cone, and coefficient sums of the symmetrized A, so
% that a split of an off-diagonal coefficient between (i,j) and (j,i) and
% a permutation of rows, of PSD blocks, or within a block do not change it.
global HD_FP
[m,n] = size(A);
f = struct('m',m,'n',n,'Kf',getk(K,'f'),'Kl',getk(K,'l'),'Ks',getk(K,'s'),'Kq',getk(K,'q'));
% Column of the transposed entry of each PSD block entry (identity on the
% free, linear and Lorentz parts)
tp = (1:n)';
o = f.Kf + f.Kl + sum(f.Kq);
for s = f.Ks(:)'
    [i,j] = ndgrid(1:s,1:s);
    tp(o+(1:s^2)) = o + (i(:)-1)*s + j(:);
    o = o + s^2;
end
As = (A + A(:,tp))/2;
[~,~,v] = find(As);
f.nnzA = nnz(As);   f.sumabsA = sum(abs(v));    f.sumsqA = sum(v.^2);
f.rownnz2 = sum(full(sum(As~=0,2)).^2);
f.normb = norm(b);  f.nnzb = nnz(b);    f.sumabsb = full(sum(abs(b)));
f.nnzc = nnz(c);    f.sumc = full(sum(c));
f.scal = f.Kf + f.Kl + sum(f.Kq) + sum(f.Ks.*(f.Ks+1)/2);
HD_FP = f;
error('hdshadow:stop','SDP recorded; solve skipped.');
end

function v = getk(K,fld)
% V = GETK(K,FLD): field FLD of the cone struct K as a row, 0 or [] if absent.
v = [];
if isfield(K,fld) && ~isempty(K.(fld)),  v = double(K.(fld)(:)');  end
if isempty(v) && any(strcmp(fld,{'f','l'})),  v = 0;  end
end

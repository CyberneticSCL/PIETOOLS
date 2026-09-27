function F = bl_face(At,b,K,maxpass)                                        % CC, 09/27/2026
% BL_FACE  Forced zeros of a SeDuMi-form SDP {At'x = b, x in K}: diagonal
% entries of PSD blocks that every feasible point must have equal to zero, and
% hence whole rows/columns of those blocks.  One-sided facial reduction by
% inspection of the rows, no solve.
%
%   F = bl_face(At,b,K)           % At is nvar x m, K.f free then K.s blocks
%
% RULE.  A row with b_i = 0 whose remaining (not yet eliminated) terms are ALL
% diagonal entries of PSD blocks, all with the same sign, forces each of those
% diagonal entries to zero (a same-sign sum of nonnegative numbers is zero only
% if each is), and a PSD block with X_rr = 0 has row and column r zero.  The
% eliminated entries are removed and the rule re-applied until nothing changes
% (maxpass, default 20).  Rows with a free-variable or an off-diagonal term are
% never used, so the rule only ever finds zeros that are implied: every
% feasible X lies on the reduced face.
%
% WHY (measured 2026-09-26/27).  The pinned PIETOOLS programs have no interior
% in some Gram blocks (hinf_rd1: two eigenvalues pinned below 1e-8 x max even at
% 2 gamma*), so an exact PSD repair of a solver's point stalls at ~-1e-8 and the
% feasible certificate needed a tolerance.  On banked cuADMM iterates the rows
% of this form are exactly those with componentwise backward error omega = 1
% (b_i = 0, terms of one sign).  Restricting the repair to the face they define
% removes those directions, so the certificate's PSD test can be run on the
% reduced blocks, where strict positivity is possible.
% CC, 09/27/2026 (b) CORRECTION, measured: it is not, on the cases tried.  The
% largest margin t (Y_k >= t s_k I, s_k the block's mean eigenvalue) stays ~0
% after removing these zeros (hinf_rd1_hv 2e-9, stab_rd1_hv -3e-8), limited by
% the block the zeros came from, so more face remains there.  A diagonal-cone
% partial facial reduction LP finds exactly these zeros in 1-D and none in 2-D
% (stab2_rd_psz) or nl_fisher_opt, though those have no interior either.  Kept
% for reference; facial reduction is closed as a solver-side step.
%
% OUTPUT F: zero{k} (logical N_k: forced-zero indices of block k), elim
% (logical nvar: every variable in a forced-zero row/column), rows (the rows
% used), npass, nzero (per block).
%
% COST: O(nnz(At)) per pass on the triplets, plus O(nvar) for the index maps.

if nargin < 4 || isempty(maxpass), maxpass = 20; end
if ~isfield(K,'f') || isempty(K.f), K.f = 0; end
Ks = double(K.s(:)');  nb = numel(Ks);  nvar = size(At,1);  m = size(At,2);
b = full(b(:));
% block, row, column of every variable (0 = free)
blk = zeros(nvar,1);  rr = zeros(nvar,1);  cc = zeros(nvar,1);  off = K.f;
for k = 1:nb
    N = Ks(k);  ix = off + (1:N^2)';  p = (1:N^2)' - 1;
    blk(ix) = k;  rr(ix) = mod(p,N) + 1;  cc(ix) = floor(p/N) + 1;  off = off + N^2;
end
[vi,ci,av] = find(At);
isdg = blk(vi) > 0 & rr(vi) == cc(vi);
zero = arrayfun(@(N) false(N,1),Ks,'UniformOutput',false);
elim = false(nvar,1);  used = false(m,1);  zb = (b == 0);
for pass = 1:maxpass
    alive = ~elim(vi);
    na = accumarray(ci(alive),1,[m 1]);
    nd = accumarray(ci(alive),double(isdg(alive)),[m 1]);
    np = accumarray(ci(alive),double(isdg(alive) & av(alive) > 0),[m 1]);
    cand = zb & ~used & na > 0 & nd == na & (np == na | np == 0);
    if ~any(cand), break; end
    used = used | cand;
    sel = alive & cand(ci) & isdg;
    kk = blk(vi(sel));  ri = rr(vi(sel));
    for k = unique(kk)'
        zero{k}(ri(kk == k)) = true;
    end
    % every variable in a zeroed row or column of its block
    z = false(nvar,1);
    for k = 1:nb
        if any(zero{k})
            ix = find(blk == k);
            z(ix) = zero{k}(rr(ix)) | zero{k}(cc(ix));
        end
    end
    elim = elim | z;
end
F = struct('zero',{zero},'elim',elim,'rows',find(used),'npass',pass, ...
           'nzero',cellfun(@nnz,zero),'Ks',Ks);
end

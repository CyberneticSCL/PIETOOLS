function prob = Sedumi2Mosek(A,b,c,K)
% Phase 3 capture shadow, 11/27/2026: logs SHA-256 digests of the SeDuMi
% data A (= At'), b, c, K that sossolve hands to MOSEK, then returns what
% the original returns. Records only.
persistent h
if isempty(h),  h = d3_orig('Sedumi2Mosek');    end
d3_log('S2M',{'A','b','c','K'},{A,b,c,K},...
       sprintf('m=%d n=%d nnzA=%d nb=%.17g nKs=%d',size(A,1),size(A,2),nnz(A),norm(full(b)),...
               numel(d3_field(K,'s'))));
prob = h(A,b,c,K);
end

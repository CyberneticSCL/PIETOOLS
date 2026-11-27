function [X,Y] = MosekSol2SedumiSol(K,res)
% Phase 3 capture shadow, 11/27/2026: calls the original, then logs
% SHA-256 digests of the primal X and dual Y it returns (sossolve takes
% y = Y(1:m)) and the MOSEK status strings. Records only.
persistent h
if isempty(h),  h = d3_orig('MosekSol2SedumiSol');  end
[X,Y] = h(K,res);
st = '';
if isfield(res,'sol') && isfield(res.sol,'itr')
    st = [res.sol.itr.prosta,'/',res.sol.itr.solsta];
end
d3_log('M2S',{'x','y'},{X,Y},sprintf('nx=%.17g ny=%.17g st=%s',norm(X),norm(Y),st));
end

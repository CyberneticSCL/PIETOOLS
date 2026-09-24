function A = t1d_ipm(B,solver)                                              % CC, 09/22/2026
% PROVENANCE.  scratchpad/reach1d/r1ipm.m verbatim, renamed (B6).  This is the
% ATTRIBUTION reference of T8: where the search does not certify, the same
% program is handed to the stock interior-point path, so a miss can be charged
% to the search or to the relaxation instead of being reported as "not
% reached" with no owner.
%                                                                 % CC, 09/22/2026
% r1ipm(B) -- CURVE A: the stock interior-point path on the LPI built by
% r1build, verified by the operator gate.  This is the RELAXATION'S reach: no
% search can certify where a full IPM solve of the same program cannot.
% Verdict comes from gate1d (relP + per-block PSD), never from the solver's
% status flag, which is recorded only as a diagnostic.
if nargin<2||isempty(solver), solver='sedumi'; end
A = struct('solver',solver,'ok',false,'cert',false,'relP',NaN,'relD',NaN, ...
           'psd',false,'mineig',[],'normQ',[],'rank',[],'rank_tol',[], ...
           't_solve',NaN,'t_gate',NaN,'err','','maxPop',NaN,'maxDop',NaN);
opts.solver = solver;   opts.simplify = false;
t0 = tic;
try
    evalc('sol = lpisolve(B.prog,opts);');
catch ME
    A.err = ME.message;   A.t_solve = toc(t0);   return
end
A.t_solve = toc(t0);
x = full(sol.solinfo.RRx(:));
if numel(x) ~= B.Ntot
    A.err = sprintf('RRx length %d ~= Ntot %d',numel(x),B.Ntot);  return
end
t0 = tic;
G = gate1d(B.prog,B.H,B.P,x/B.P.nb0);
A.t_gate = toc(t0);
A.relP=G.relP; A.relD=G.rel; A.psd=G.psd; A.cert=G.cert; A.ok=true;
A.mineig=G.mineig; A.normQ=G.normQ; A.rank=G.rank;
A.maxPop=G.maxPop; A.maxDop=G.maxDop;
% numerical rank is tolerance sensitive: report the sweep, not one number
tols = [1e-4 1e-6 1e-8 1e-9 1e-12];
A.rank_tol = zeros(numel(tols),numel(B.Ns));   A.tols = tols;
for i=1:numel(B.Ns)
    Q = reshape(x(B.Kf+sum(B.Ns(1:i-1).^2)+(1:B.Ns(i)^2)),B.Ns(i),B.Ns(i));
    e = eig((Q+Q')/2);
    for j=1:numel(tols), A.rank_tol(j,i) = nnz(e > tols(j)*max(e)); end
end
A.x = x;
end

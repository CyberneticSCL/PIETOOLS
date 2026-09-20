function R = opcheck_pos2d(prog,Ptgt,Peop,P,q,use_getsol)                   % CC, 09/19/2026
% OPERATOR-LEVEL acceptance test for a 2D BARE POSITIVITY program: the 2D
% analogue of opcheck_pos.m, with opcheck_2d.m's cell enumeration.
% Push a Gram vector q (NORMALISED-b units, bm_setup convention) back through
% the dopvar2d Peop and test the LPI itself:
%     accept iff  max|Ptgt - Peop(q)| / max|Ptgt| < 1e-6  over ALL 36 leaf
%     cells, AND every Gram block PSD.
% An equality residual or a solver flag is NOT acceptance: on this study's own
% 2D family an lsqr point passed op rel = 1.7e-11 with mineig = -4.4e-7, and a
% near-zero Gram passes the equality rows of a tiny-||b|| program while
% certifying nothing.
%
% q is multiplied by P.nb0 internally (trap 2 of the 1D pm study: omitting
% nb0 reports rel = 1 for a perfectly good certificate).
% Optional use_getsol (default false): true substitutes through PIETOOLS'
% lpigetsol instead of pielr_evalgram -- the oracle/debug path only.
%
% CC, 09/19/2026: substitute via private/pielr_evalgram instead of
%          lpigetsol, as in opcheck_2d: master's (998a68cc) sosgetsol runs
%          dpvar cells through dpvar2poly, unusably slow at this package's
%          decision-variable counts; the direct sparse substitution matches
%          the branch getsol to machine precision (measured).
if nargin<6, use_getsol = false; end       % default: direct substitution   % CC, 09/19/2026
pg = prog;
pg.solinfo.info = struct('verified_by','opcheck_pos2d');
pg.solinfo.RRx  = full(q(:))*P.nb0;
%Pe  = lpigetsol(pg,Peop);                                                  % CC, 09/19/2026 (was)
if use_getsol                              % oracle path (slow on master)   % CC, 09/19/2026
    Pe = lpigetsol(pg,Peop);                                                % CC, 09/19/2026
else                                       % see pielr_evalgram             % CC, 09/19/2026
    Pe = pielr_evalgram(pg,Peop);                                           % CC, 09/19/2026
end                                                                         % CC, 09/19/2026
Res = Ptgt - Pe;
[R.maxP,R.npar] = mxop2d(Ptgt);
R.maxRes = mxop2d(Res);
R.rel    = R.maxRes/max(R.maxP,realmin);
B = numel(P.Ns);
R.normQ = zeros(1,B); R.mineig = zeros(1,B); R.maxeig = zeros(1,B); R.rank = zeros(1,B);
for i = 1:B
    Q = reshape(q(P.rows{i}),P.Ns(i),P.Ns(i));  Q = (Q+Q')/2;
    e = sort(eig(Q),'descend');
    R.normQ(i)  = norm(Q,'fro')*P.nb0;
    R.mineig(i) = e(end)*P.nb0;
    R.maxeig(i) = e(1)*P.nb0;
    R.rank(i)   = sum(e > max(e(1),eps)*1e-9);
end
R.psd = all(R.mineig >= -1e-9*max(R.maxeig,realmin));
R.ok  = (R.rel < 1e-6) && R.psd;
end

% =========================================================================
function [m,npar] = mxop2d(X)
% max|coefficient| over EVERY opvar2d/dopvar2d parameter cell (36 leaves),
% enumeration taken from properties() and ASSERTED, copied from the validated
% opcheck_2d.m: a silently missed cell is the one failure mode that turns a
% wrong answer into a passing one.
pr   = properties(X);
skip = {'I','var1','var2','dim','dimdependent'};
pars = pr(~ismember(pr,skip));
m = 0;  npar = 0;
for i = 1:numel(pars)
    v = X.(pars{i});
    if iscell(v)
        for k = 1:numel(v)
            m = max(m,mxc(v{k}));  npar = npar+1;
        end
    else
        m = max(m,mxc(v));  npar = npar+1;
    end
end
assert(numel(pars)==16 && npar==36, ...
    'opcheck_pos2d: opvar2d enumeration changed (%d fields / %d leaves)', ...
    numel(pars),npar);
end

function m = mxc(X)
if isempty(X), m = 0; return; end
if isa(X,'double')
    C = X;
elseif isa(X,'polynomial')
    C = X.coefficient;
elseif isa(X,'dpvar')
    C = X.C;
else
    m = NaN; return;
end
if isempty(C), m = 0; else, m = full(max(abs(C(:)))); end
if isempty(m), m = 0; end
end

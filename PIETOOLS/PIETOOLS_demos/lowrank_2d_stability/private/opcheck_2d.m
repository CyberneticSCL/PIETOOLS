function R = opcheck_2d(prog,H,P,q,use_getsol)                              % CC, 09/19/2026
% opcheck_2d -- OPERATOR-LEVEL verification of a 2-D LPI solution.
%
% The 2D analogue of opcheck.m.  Pushes a Gram vector q back through the real
% opvar2d operators and tests the LPI identity that was actually imposed:
%
%   Pop  = poslpivar_2d(...) + eppos*I          (positive whenever Q_LF >= 0)
%   Qop  = (Aop'*Pop*Top)' + Aop'*Pop*Top  [+ 2*epneg*Top'*Pop*Top]
%   Deop = poslpivar_2d(...)                    (the executive's Qeop)
%   identity to certify:   Qop + Deop = 0
%
% WHY THIS EXISTS.  norm(b) ~ 5e-6 on these programs, so X = 0 satisfies the
% equality rows to ~1e-9 relative while being no certificate at all.  An
% equality residual or a solver status flag therefore cannot establish 2D
% feasibility; only the operator identity plus PSD Gram blocks can.
%
% UNITS (trap 2).  q is in NORMALISED-b units, as produced by bm_setup /
% bm_resid: this routine multiplies by P.nb0 internally, exactly as opcheck
% does.  Passing original units silently reports rel = 1.
%
% INPUT
%   prog - the program from build_stab_2d_st2 (solved or not; only solinfo is used)
%   H    - operator handles from build_stab_2d_st
%   P    - bm_setup package (needs .nb0, .Ns, .rows)
%   q    - length P.Ntot vector, normalised-b units
%   use_getsol - optional, default false: true substitutes through          % CC, 09/19/2026
%          PIETOOLS' lpigetsol/getsol_lpivar_2d instead of pielr_evalgram   % CC, 09/19/2026
%          -- the oracle/debug path, and the only reason to set it          % CC, 09/19/2026
% OUTPUT
%   R.rel        relative operator residual  max|Qop+Deop| / max|Qop|   <-- the number
%   R.maxQop .maxDeop .maxRes .maxPop    the max|coeff| that ratio is built from
%   R.mineig(i)  min eig of Gram block i, ORIGINAL units
%   R.rank(i)    numerical rank of Gram block i (rel tol 1e-9)
%   R.normQ(i)   Frobenius norm of Gram block i, original units
%   R.psd        true iff every block's min eig >= -1e-8 * its own norm
%   R.ok         true iff R.rel < 1e-6 AND R.psd  (the acceptance gate)
%   R.npar       number of opvar2d leaf cells actually scanned (36; see below)
%
% CC, 09/19/2026: substitute the Gram vector via private/pielr_evalgram
%          instead of lpigetsol.  Master's (998a68cc) sosgetsol takes each
%          dpvar cell through dpvar2poly, which did not finish in 35 min on
%          this package's ~277k-dvar negativity operator (branch: ~1 s);
%          the direct sparse substitution is equivalent to the branch
%          getsol to machine precision (measured: both operators, on the
%          certificate and on random decision vectors).  use_getsol=true
%          keeps the getsol path callable as the oracle.

if nargin<5, use_getsol = false; end       % default: direct substitution   % CC, 09/19/2026
pg = prog;
pg.solinfo.info = struct('verified_by','opcheck_2d');
pg.solinfo.RRx  = full(q(:))*P.nb0;        % back to the original b scaling

%Pop  = lpigetsol(pg,H.Pop);                                                % CC, 09/19/2026 (was)
%Deop = lpigetsol(pg,H.Deop);                                               % CC, 09/19/2026 (was)
if use_getsol                              % oracle path (slow on master)   % CC, 09/19/2026
    Pop  = lpigetsol(pg,H.Pop);                                             % CC, 09/19/2026
    Deop = lpigetsol(pg,H.Deop);                                            % CC, 09/19/2026
else                                       % see pielr_evalgram             % CC, 09/19/2026
    Pop  = pielr_evalgram(pg,H.Pop);                                        % CC, 09/19/2026
    Deop = pielr_evalgram(pg,H.Deop);                                       % CC, 09/19/2026
end                                                                         % CC, 09/19/2026

% Rebuild the negativity operator from the SUBSTITUTED Pop, by the same five
% lines build_stab_2d_st2 used.  Qop is a linear map of Pop and Pop is affine in
% the decision variables, so substituting then rebuilding equals rebuilding
% then substituting -- this tests the identity semantically instead of trusting
% the assembled equality rows.
PTop  = Pop*H.Top;
APTop = H.Aop'*PTop;
if H.epneg==0
    Qop = APTop' + APTop;
else
    Qop = APTop' + APTop + 2*H.epneg*(H.Top'*PTop);
end
Qop = clean_opvar(Qop,1e-12);              % same ztol the constraint was built at

Res = Qop + Deop;

[R.maxQop ,R.npar] = mxop2d(Qop);
R.maxDeop          = mxop2d(Deop);
R.maxRes           = mxop2d(Res);
R.maxPop           = mxop2d(Pop);
R.rel = R.maxRes/max(R.maxQop,realmin);

% ---- per-block spectrum of the Gram blocks -------------------------------
B = numel(P.Ns);
R.normQ = zeros(1,B);  R.mineig = zeros(1,B);  R.rank = zeros(1,B);
for i = 1:B
    Q = reshape(q(P.rows{i}),P.Ns(i),P.Ns(i));
    Q = (Q+Q')/2;
    e = sort(eig(Q),'descend');
    R.normQ(i)  = norm(Q,'fro')*P.nb0;
    R.mineig(i) = e(end)*P.nb0;
    R.rank(i)   = sum(e > max(e(1),eps)*1e-9);
end
R.psd = all(R.mineig >= -1e-8*max(R.normQ,realmin));
R.ok  = (R.rel < 1e-6) && R.psd;
end

% =========================================================================
function [m,npar] = mxop2d(X)
% max|coefficient| over EVERY opvar2d/dopvar2d parameter cell.
%
% opvar2d has 16 parameter fields and 36 leaf cells (R00,R0x,R0y,R02, Rx0,
% Rxx{3},Rxy,Rx2{3}, Ry0,Ryx,Ryy{3},Ry2{3}, R20,R2x{3},R2y{3},R22{3,3}) --
% against 6 leaves in 1D (P,Q1,Q2,R.R0,R.R1,R.R2).  MEASURED from the live
% class, not assumed: see op2d_enum.m.
%
% The enumeration is taken from properties() rather than hard-coded, and the
% leaf count is ASSERTED, because a silently missed cell would let a nonzero
% residual read as zero -- the one failure mode that turns a wrong answer into
% a passing one.
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
    'opcheck_2d: opvar2d enumeration changed (%d fields / %d leaves); rerun op2d_enum.m', ...
    numel(pars),npar);
end

function m = mxc(X)
if isempty(X), m = 0; return; end
if     isa(X,'double'),     C = X;
elseif isa(X,'polynomial'), C = X.coefficient;
elseif isa(X,'dpvar'),      C = X.C;
else,  m = NaN; return;                       % unknown leaf type: poison, do not hide
end
if isempty(C), m = 0; else, m = full(max(abs(C(:)))); end
if isempty(m), m = 0; end
end

function out = heat3d_stability_sop(kappa,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = HEAT3D_STABILITY_SOP(KAPPA,OPTS) certifies exponential stability at
% rate KAPPA of the 3-D heat equation
%
%   u_t = u_s1s1 + u_s2s2 + u_s3s3  on [0,1]^3,  s1 DD, s2 and s3 DN,
%
% (exact rate kappa* = lambda1 = 3 pi^2/2 = 14.804407) with the heatNd
% benchmark LPI (Jagt & Peet, arXiv:2508.14840v4, Cor. 35; PIETOOLS_demos/
% sopvar_demos, preset 'bench', d = 0, ep = 0.1), and extracts and checks
% the certificate with lpigetsol_sop. Stock PIETOOLS stops at 2 spatial
% variables; the container classes pose the LPI in 3:
%
%   (15a)  P*T - ep^2 T*T - R = 0,  R >= 0   (imposed as P*T - T*P = 0 and
%          its self-adjoint part, 'split'),
%   (15b)  P*A0 + A0*P + 2 kappa P*T + Q = 0,  Q >= 0.
%
% The SDP (m 31770, 2.7e6 decision variables) is solved by heatNd_solve
% (MOSEK) on a maximal independent subset of the rows; the verdict is on
% ALL rows (+1 needs MOSEK OPTIMAL, rel_b <= 1e-6, PSD blocks). The
% measured bisection bracket at d = 0 is [14.60, 14.70) (heatNd README
% Sec. 6.3), so the default KAPPA = 14.0 is below it.
%
% CHECKS (each asserted; error id heat3d_stability_sop:check):
% (1) heatNd_solve verdict +1 (certified on all rows);
% (2) the extracted constraint operators reproduce the SDP rows: with
%     solinfo.RRx = RR*x*bscl (heatNd_sdp's map to decvartable order),
%     the (15a) and (15b) row sets have residual norms r_a, r_b with
%     r_a^2 + r_b^2 = ||At'RRx - b||^2. Every row is a canonical coefficient
%     of the extracted constraint operators, one per adjoint pair
%     ('symmetric'), so the coefficient norm c is >= r. (15b) is self-adjoint
%     for every decision value, so c_b <= sqrt(2) r_b. In (15a) P*T - T*P is
%     anti-self-adjoint and P*T - ep^2 T*T - R = S + (P*T - T*P)/2 with S
%     self-adjoint, whence c_a <= 2 r_a;
% (3) by action on test functions (opcheck_sop, 3-D quadrature): (15a)
%     P*T - ep^2 T*T = R and P*T = (P*T)', (15b) X1 + kappa X2 = -Q, each to
%     1e-5 relative; R >= 0, Q >= 0 and -(X1 + kappa X2) >= 0 (the decrease
%     condition) on the test functions (>= -1e-8, normalized);
% (4) control: RRx perturbed by 1e-3 of its rms entry breaks (3) by more
%     than 1e3 times the real residual.
%
% INPUT
% - kappa: rate (default 14.0); the SDP depends on kappa only (heatNd_lpi);
% - opts (struct, optional):
%     keep      row set of the SDP (heatNd_lindep); default: computed here,
%               one sparse QR, 363 s measured;
%     keepfile  a heatNd_sdp dump of the same family whose info.keep is
%               reused (checked: same m, K and b; At is affine in kappa, the
%               dependency structural, heatNd_lindep header);
%     nt, nq    test functions (4) and Gauss nodes per direction (6);
%     N         number of variables (default 3; 1 or 2 is a cheap dry run of
%               the same code, s1 DD, the others DN).
% OUTPUT struct: the verdict v (heatNd_solve), the row/coefficient norms of
%   (2), the residuals and margins of (3), the control, sizes, times (s)
%   and memory (MB).
%
% COST (measured, heatNd README Sec. 8): build 74-119 s; MOSEK on the
% 22226 independent rows 371-697 s at 10.7-11.45 GB; on all rows 499-921 s
% at 20.1-20.9 GB. This example, 09/29/2026, on a shared machine: build
% 133 s, MOSEK 852 s (20 iterations), process peak 12.97 GB, extraction of
% 9 operators 10.8 s, checks 57 s, 18 min in all. Price before running.
%
% Initial coding MMP, 09/29/2026. Tier 1 end-to-end example E4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(kappa),  kappa = 14.0;       end
if nargin<2 || isempty(opts),   opts = struct();    end
if ~isfield(opts,'nt'),     opts.nt = 4;    end
if ~isfield(opts,'nq'),     opts.nq = 6;    end
if ~isfield(opts,'N'),      opts.N = 3;     end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
chk = @(c,varargin) assert(c,'heat3d_stability_sop:check',varargin{:});
ep = 0.1;   d = 0;

% % % Program and SDP.
t0 = tic;
pie = heatNd_pie(opts.N,0);
[prog,meta] = heatNd_lpi(pie,d,kappa,ep,struct('preset','bench'));
t_build = toc(t0);
B = meta.base;
t0 = tic;   D = heatNd_sdp(prog);   t_sdp = toc(t0);
fprintf('E4 %d-D heat, kappa %.4f (kappa* %.6f): ndec %d, m %d, Ks %s; build %.1f s, SDP %.1f s\n',...
    opts.N,kappa,pie.exact.lambda1,numel(prog.decvartable),D.m,mat2str(D.K.s),t_build,t_sdp);

% % % Row set, then the solve (verdict on all rows).
t0 = tic;
if isfield(opts,'keep') && ~isempty(opts.keep)
    keep = opts.keep(:);    src = 'given';
elseif isfield(opts,'keepfile') && ~isempty(opts.keepfile)
    S = load(opts.keepfile,'b','K','info');
    chk(isfield(S.info,'keep') && S.info.m==D.m && isequal(S.K,D.K) && isequal(S.b,D.b),...
        '(0) keepfile is not of this SDP family (m, K or b differ)')
    keep = S.info.keep(:);  src = opts.keepfile;
else
    keep = heatNd_lindep(D.At);     src = 'heatNd_lindep';
end
t_keep = toc(t0);
fprintf('   rows: %d independent of %d (%s, %.1f s)\n',numel(keep),D.m,src,t_keep);
v = heatNd_solve(D,struct('keep',keep,'rows','keep','keepx',true));
fprintf('   heatNd_solve: st %+d (%s); rel_b %.2e (all rows), psd_relmin %.2e, MOSEK %.0f s, %d iter, peak %.0f MB\n',...
    v.st,v.why,v.rel_b,v.psd_relmin,v.t_mosek,v.iter,v.mem_peak);
chk(v.st==1,'(1) not certified: %s',v.why)

% % % Extraction: the dump's point in decvartable order is a solved program.
t0 = tic;
prog.solinfo.x = v.x*D.bscl;
prog.solinfo.RRx = full(D.RR*v.x)*D.bscl;
prog.solinfo.info = struct('route','heatNd_solve','st',v.st,'rel_b',v.rel_b);
v.x = [];
Ps = lpigetsol_sop(prog,B.P);   Rs = lpigetsol_sop(prog,B.R);   Qs = lpigetsol_sop(prog,B.Q);
PT = B.P'*pie.T;                                % as heatNd_lpi composes it
Ea1 = PT - PT';     Ea2 = B.E1 - B.R;   Eb = B.X1 + meta.kappa*B.X2 + B.Q;
Ea1s = lpigetsol_sop(prog,Ea1);     Ea2s = lpigetsol_sop(prog,Ea2);     Ebs = lpigetsol_sop(prog,Eb);
E1s = lpigetsol_sop(prog,B.E1);     PTs = lpigetsol_sop(prog,PT);
Ls = lpigetsol_sop(prog,B.X1 + meta.kappa*B.X2);
t_extract = toc(t0);
chk(all(cellfun(@(X) isa(X,'copvar'),{Ps,Rs,Qs,Ea1s,Ea2s,Ebs,E1s,PTs,Ls})),'(1) extracted operators not fixed')

% (2) extraction against the SDP rows.
na = B.prog.expr.num;
ra = row_res(prog,1:na);    rb_ = row_res(prog,na+1:prog.expr.num);
rall = norm(D.At'*v_x(D,prog) - D.b)*D.bscl;
ca = hypot(opcheck_sop('coef',Ea1s),opcheck_sop('coef',Ea2s));     cb = opcheck_sop('coef',Ebs);
fprintf('   rows (15a) %.3e, (15b) %.3e, all %.3e (hypot %.3e); ||coef|| / ||rows||: (15a) %.5f, (15b) %.5f\n',...
    ra,rb_,rall,hypot(ra,rb_),ca/ra,cb/rb_);
chk(abs(hypot(ra,rb_)/rall-1)<=1e-6,'(2) the two row sets are not the SDP''s rows (%.6g vs %.6g)',hypot(ra,rb_),rall)
chk(ca>=ra*(1-1e-6) && ca<=2*ra*(1+1e-6),'(2) (15a): coefficient norm %.6g vs row residual %.6g',ca,ra)
chk(cb>=rb_*(1-1e-6) && cb<=sqrt(2)*rb_*(1+1e-6),'(2) (15b): coefficient norm %.6g vs row residual %.6g',cb,rb_)

% (3) by action, 3-D quadrature.
o = struct('nt',opts.nt,'deg',2,'nqi',opts.nq,'nqo',opts.nq);
t0 = tic;
r15a = opcheck_sop('res',E1s,Rs,o);             % P*T - ep^2 T*T = R
rPT = opcheck_sop('res',PTs,PTs',o);            % P*T self-adjoint
r15b = opcheck_sop('res',Ls,-Qs,o);             % X1 + kappa X2 = -Q
pR = opcheck_sop('psd',Rs,o);   pQ = opcheck_sop('psd',Qs,o);   pL = opcheck_sop('psd',-Ls,o);
t_check = toc(t0);
fprintf(['   by action: (15a) %.2e, P*T = (P*T)'' %.2e, (15b) %.2e | min <x,Rx> %.2e, <x,Qx> %.2e, '...
         '-<x,(X1+kX2)x> %.2e (normalized); %.1f s\n'],r15a,rPT,r15b,pR,pQ,pL,t_check);
chk(max([r15a,rPT,r15b])<=1e-5,'(3) the extracted operators violate (15a)/(15b) (%.3g, %.3g, %.3g)',r15a,rPT,r15b)
chk(pR>=-1e-8 && pQ>=-1e-8 && pL>=-1e-8,'(3) positivity fails on test functions (%.3g, %.3g, %.3g)',pR,pQ,pL)

% (4) control: RRx perturbed by 1e-3 of its rms entry.
bad = prog;     x = prog.solinfo.RRx;   s = rng;    rng(7);
bad.solinfo.RRx = x + 1e-3*norm(x)/sqrt(numel(x))*randn(size(x));  rng(s);
r15ax = opcheck_sop('res',lpigetsol_sop(bad,B.E1),lpigetsol_sop(bad,B.R),o);
r15bx = opcheck_sop('res',lpigetsol_sop(bad,B.X1 + meta.kappa*B.X2),-lpigetsol_sop(bad,B.Q),o);
fprintf('   control RRx + 1e-3 rms noise: (15a) %.2e, (15b) %.2e\n',r15ax,r15bx);
chk(max(r15ax,r15bx)>1e3*max([r15a,r15b,1e-12]),'(4) a perturbed RRx is not detected')

out = struct('kappa',kappa,'kstar',pie.exact.lambda1,'v',v,'nkeep',numel(keep),'m',D.m,...
    'ndec',numel(prog.decvartable),'Ks',D.K.s,'rows',[ra rb_ rall],'coef_rows',[ca/ra cb/rb_],...
    'res',[r15a rPT r15b],'psd',[pR pQ pL],'control',[r15ax r15bx],'mem',meta.mem,...
    't',struct('build',t_build,'sdp',t_sdp,'keep',t_keep,'mosek',v.t_mosek,'extract',t_extract,'check',t_check));
fprintf('E4 PASS: kappa %.4f certified, certificate extracted and checked\n',kappa);
end


% ------------------------------------------------------------------------
function x = v_x(D,prog)
% The solved cone vector, normalized as D.b (heatNd_sdp: x_orig = x*bscl).
x = prog.solinfo.x/D.bscl;
end

function r = row_res(prog,rows)
% ||At_i' x - b_i|| over the expressions ROWS, x = solinfo.RRx
% (decvartable order, the order of At rows).
x = prog.solinfo.RRx(:);    r2 = 0;
for i = rows
    At = prog.expr.At{i};   b = prog.expr.b{i};
    r2 = r2 + sum((At'*x(1:size(At,1)) - b).^2);
end
r = sqrt(full(r2));
end

function out = stability_2d_sop(frac,level)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = STABILITY_2D_SOP(FRAC,LEVEL) solves the 2-D PIE stability LPI of
% executives/2D/PIETOOLS_stability_2D.m (direct form) on the container
% path for the 2-D heat equation u_t = u_s1s1 + u_s2s2 + FRAC*2 pi^2 u,
% Dirichlet on [0,1]^2 (cx_plant('rd2',FRAC)), and tests the EXTRACTION of
% its solution with lpigetsol_sop:
%
%   P = P0 + diag(eppos) >= 0,   Q = (A'(P T))' + A'(P T) = -Qe,  Qe >= 0.
%
% NOT A CERTIFICATE. Psatz-free 2-D stability is limited by solver
% accuracy (MOSEK fails at every fraction; SeDuMi stops at rel_b ~ 7e-6 at
% FRAC = 0.1, above the 1e-6 gate: memory notes dopvar2d-vs-sopvar-2d-
% precision, psatz-generators-2d). So this solves with SeDuMi and shows
% that the extracted operators reproduce the solver's own residual level;
% it does not certify stability.
%
% The program is cx_stability_2D's (sopvar/Testfolder/sdopvar/
% claude_tests/cx_exec), which returns no operators and is frozen by the
% maintainer. It is transcribed below (BUILD_2D) with lpiprogram_sop and
% lpi_eq_sop in place of lpiprogram and lpi_eq_cdopvar, and returns P, Q,
% Qe; the transcription must build cx_stability_2D's program exactly (At,
% b, variables, decvartable), which cx_stability_2D, called unmodified,
% provides as the reference.
%
% CHECKS (each asserted; error id stability_2d_sop:check):
% (0) the transcription builds cx_stability_2D's program exactly;
% (1) SeDuMi returns a point (not X = 0; rel_b < 1); pinf, numerr,
%     feasratio, rel_b and the PSD margin are reported;
% (2) the extracted constraint operator E = Qe + Q reproduces the SDP
%     residual: every row of the program is a canonical coefficient of E,
%     one per adjoint pair ('symmetric'), so ||r_rows|| <= ||coef E(sol)||
%     <= sqrt(2) ||r_rows|| (1e-6 slack);
% (3) getsol(Q) equals Q rebuilt from getsol(P) by the operator algebra
%     to 1e-10 (by action);
% (4) by action on test functions: Qe + Q = 0 to within the solver's
%     accuracy class (relative residual <= 1e-3; the value is reported
%     beside rel_b); P and Qe >= -1e-6 on the test functions;
% (5) control: RRx perturbed by 1e-3 of its rms entry breaks (4) by more
%     than 10 times (solinfo.x equals RRx here: no free variable, no
%     inequality).
%
% INPUT  frac (default 0.1), level (lpisettings name, default 'light').
% OUTPUT struct: SeDuMi info, rel_b, psd margins, the ratio of (2), the
%   residuals of (3)-(5), the SDP shape and the times (s).
% Cost (measured 09/29/2026, FRAC 0.1, light): m 3456, 179840 decision
% variables, Ks [424 8]; build 1.8-2.7 s each (transcription and
% reference); SeDuMi 215-273 s; about 4-5 min in all.
%
% Initial coding MMP, 09/29/2026. Tier 1 end-to-end example E3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(frac),   frac = 0.1;         end
if nargin<2 || isempty(level),  level = 'light';    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
chk = @(c,varargin) assert(c,'stability_2d_sop:check',varargin{:});
st = cx_settings(level,'sedumi');
PIE = initialize(cx_plant('rd2',frac));

% % % Reference program (frozen file, unmodified) and the transcription.
t0 = tic;   pref = cx_stability_2D(PIE,st);     t_ref = toc(t0);
t0 = tic;   [prog,ops] = build_2d(PIE,st);      t_build = toc(t0);
same = isequal(pref.expr.At,prog.expr.At) && isequal(pref.expr.b,prog.expr.b) && ...
       isequal(pref.var.idx,prog.var.idx) && isequal(pref.var.type,prog.var.type) && ...
       isequal(pref.decvartable,prog.decvartable) && isequal(pref.dom,prog.dom);
S = cx_shape(prog);
fprintf('E3 rd2 %.2f/%s: m %d, ndv %d, Ks %s; build %.1f s (reference %.1f s); transcription = cx_stability_2D: %d\n',...
    frac,level,S.m,S.ndv,mat2str(S.Ks),t_build,t_ref,same);
chk(same,'(0) the transcription does not build cx_stability_2D''s program')

% % % Solve (SeDuMi).
t0 = tic;   prog = lpisolve(prog,st.sos_opts);  t_solve = toc(t0);
info = prog.solinfo.info;
[rb,pmin,prel,triv] = cx_resid(prog);
fprintf('   SeDuMi: pinf %d, dinf %d, numerr %d, feasratio %.4f; rel_b %.2e, psd_min %.2e (rel %.2e); solve %.1f s\n',...
    info.pinf,info.dinf,info.numerr,info.feasratio,rb,pmin,prel,t_solve);
chk(~triv && rb<1,'(1) the solver returned no point (trivial %d, rel_b %.3g)',triv,rb)

% % % Extraction.
t0 = tic;
Ps = lpigetsol_sop(prog,ops.P);     Qs = lpigetsol_sop(prog,ops.Q);
Qes = lpigetsol_sop(prog,ops.Qe);   Es = lpigetsol_sop(prog,ops.Qe + ops.Q);
t_extract = toc(t0);
chk(all(cellfun(@(X) isa(X,'copvar'),{Ps,Qs,Qes,Es})),'(1) extracted operators not fixed')

% (2) the extracted constraint operator against the SDP rows.
rall = row_res(prog,1:prog.expr.num);
ce = opcheck_sop('coef',Es);
fprintf('   ||coef E(sol)|| %.4e / ||r_rows|| %.4e = %.7f (in [1, sqrt(2)]); ||b|| %.3e\n',ce,rall,ce/rall,rall/rb);
chk(ce>=rall*(1-1e-6) && ce<=sqrt(2)*rall*(1+1e-6),'(2) coefficient norm %.6g vs row residual %.6g',ce,rall)

% (3) getsol commutes with the algebra.
o = struct('nt',6,'nqi',10,'nqo',10);
PT = Ps*ops.T;  APT = ops.A'*PT;    Qr = APT' + APT;
if isfield(ops,'epneg') && ops.epneg~=0,    Qr = Qr + 2*ops.epneg*(ops.T'*PT);   end
rQ = opcheck_sop('res',Qs,Qr,o);

% (4) the relation and positivity, by action.
rE = opcheck_sop('res',Qes,-Qs,o);
pP = opcheck_sop('psd',Ps,o);   pQe = opcheck_sop('psd',Qes,o);
fprintf('   by action: getsol(Q) vs Q(getsol(P)) %.2e | Qe + Q %.2e (rel_b %.2e) | min <x,Px> %.2e, <x,Qe x> %.2e\n',...
    rQ,rE,rb,pP,pQe);
chk(rQ<=1e-10,'(3) getsol(Q) differs from Q rebuilt from getsol(P) (%.3g)',rQ)
chk(rE<=1e-3,'(4) Qe + Q = 0 fails beyond the solver''s accuracy (%.3g)',rE)
chk(pP>=-1e-6 && pQe>=-1e-6,'(4) P or Qe not positive on test functions (%.3g, %.3g)',pP,pQe)

% (5) control: RRx perturbed by 1e-3 of its rms entry. (solinfo.x is no
% control here: with no free variable and no inequality it equals RRx,
% measured ||x - RRx|| = 0.)
bad = prog;     nd = numel(prog.decvartable);
x = prog.solinfo.RRx(1:nd);     s = rng;    rng(7);
bad.solinfo.RRx = x + 1e-3*norm(x)/sqrt(nd)*randn(nd,1);    rng(s);
rEx = opcheck_sop('res',lpigetsol_sop(bad,ops.Qe),-lpigetsol_sop(bad,ops.Q),o);
nxr = min(nd,numel(prog.solinfo.x));
fprintf('   control RRx + 1e-3 rms noise: Qe + Q %.2e (||x - RRx|| = %.2e)\n',rEx,...
    norm(prog.solinfo.x(1:nxr)-x(1:nxr)));
chk(rEx>10*rE,'(5) a perturbed RRx is not detected (%.3g vs %.3g)',rEx,rE)

out = struct('info',info,'rel_b',rb,'psd_min',pmin,'psd_relmin',prel,'coef_rows',ce/rall,...
    'rows_res',rall,'res_Qgetsol',rQ,'res_E',rE,'psd',[pP pQe],'control',rEx,'shape',S,...
    't',struct('ref',t_ref,'build',t_build,'solve',t_solve,'extract',t_extract));
fprintf('E3 PASS (extraction; not a certificate)\n');
end


% ------------------------------------------------------------------------
function [prog,ops] = build_2d(PIE,st)
% cx_stability_2D.m (MMP 09/25/2026), body transcribed line for line with
% two substitutions, lpiprogram -> lpiprogram_sop and lpi_eq_cdopvar ->
% lpi_eq_sop, and the operators returned. Keep in step with it by hand;
% check (0) catches a divergence.
if ~isfield(st,'is2D') || ~st.is2D
    sos_opts = st.sos_opts;     st = st.settings_2d;    st.sos_opts = sos_opts;
end
if ~isfield(st,'eppos'),    eppos = [1e-4;1e-6;1e-6;1e-6];
else,                       eppos = st.eppos;   if numel(eppos)==1,  eppos = eppos*ones(4,1);  end
end
epneg = 0;  if isfield(st,'epneg'),    epneg = st.epneg;   end
LF_use_psatz = st.LF_use_psatz;
LF_deg_psatz = psatz_cell(st.LF_deg_psatz,LF_use_psatz);
LF_opts_psatz = psatz_cell(st.LF_opts_psatz,LF_use_psatz);
if st.use_sosineq
    error('stability_2d_sop:gap','use_sosineq=1 needs lpi_ineq_2d, which has no container counterpart.')
end
eq_opts = st.eq_opts;   eq_deg = st.eq_deg;     eq_use_psatz = st.eq_use_psatz;
eq_deg_psatz = psatz_cell(st.eq_deg_psatz,eq_use_psatz);
eq_opts_psatz = psatz_cell(st.eq_opts_psatz,eq_use_psatz);

Tm = opvar2d2copvar(PIE.T);     Am = opvar2d2copvar(PIE.A);
prog = lpiprogram_sop(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);

[prog,Pm] = cx_stability_pos2d(prog,Tm,'out',PIE.dom,st.LF_deg,st.LF_opts);
for j = 1:length(LF_use_psatz)
    if LF_use_psatz(j)~=0
        [prog,P2m] = cx_stability_pos2d(prog,Tm,'out',PIE.dom,LF_deg_psatz{j},LF_opts_psatz{j});
        Pm = Pm + P2m;
    end
end
if ~all(eppos==0)
    np = PIE.T.dim(:,1);
    Ip = blkdiag(eppos(1)*eye(np(1)),eppos(2)*eye(np(2)),eppos(3)*eye(np(3)),eppos(4)*eye(np(4)));
    Pm = Pm + opvar2d2copvar(opvar2d(Ip,[np,np],PIE.dom,PIE.vars));
end

PTm = Pm*Tm;    APTm = Am'*PTm;
if epneg==0
    Qm = APTm' + APTm;
else
    Qm = APTm' + APTm + 2*epneg*(Tm'*PTm);
end

ztol = 1e-12;
eq_opts = cx_stability_eq_opts_2D(Qm,eq_opts,ztol);
[prog,Qem] = cx_stability_pos2d(prog,Qm,'out',PIE.dom,eq_deg,eq_opts);
for j = 1:length(eq_use_psatz)
    if eq_use_psatz(j)~=0
        eq_opts_psatz{j}.exclude = eq_opts_psatz{j}.exclude | eq_opts.exclude;
        eq_opts_psatz{j}.sep = eq_opts_psatz{j}.sep | eq_opts.sep;
        [prog,Qe2m] = cx_stability_pos2d(prog,Qm,'out',PIE.dom,eq_deg_psatz{j},eq_opts_psatz{j});
        Qem = Qem + Qe2m;
    end
end
prog = lpi_eq_sop(prog,Qem+Qm,'symmetric');
ops = struct('T',Tm,'A',Am,'P',Pm,'Q',Qm,'Qe',Qem,'epneg',epneg);
end

function outcell = psatz_cell(incell,use_psatz)
% cx_stability_2D's local function of the same name, unchanged.
if all(use_psatz==0),                   outcell = {};   return, end
if isa(incell,'struct'),                outcell = repmat({incell},1,length(use_psatz));
elseif numel(incell)==1,                outcell = repmat(incell(1),1,length(use_psatz));
elseif numel(incell)==length(use_psatz),outcell = incell;
elseif numel(incell)>=max(use_psatz),   outcell = incell(use_psatz);
else,   error('For each element of ''use_psatz'', a field should be defined.')
end
end

function r = row_res(prog,rows)
% ||At_i' x - b_i|| over the expressions ROWS, x = solinfo.RRx.
x = prog.solinfo.RRx(:);    r2 = 0;
for i = rows
    At = prog.expr.At{i};   b = prog.expr.b{i};
    r2 = r2 + sum((At'*x(1:size(At,1)) - b).^2);
end
r = sqrt(full(r2));
end

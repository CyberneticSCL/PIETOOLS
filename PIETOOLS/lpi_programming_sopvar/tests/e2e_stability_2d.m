function out = e2e_stability_2d(frac,level)                                 % MMP, 10/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = E2E_STABILITY_2D(FRAC,LEVEL) builds, on the container path, the
% stability LPI of the PIE  T x_t = A x  of the 2-D heat equation
%
%   u_t = u_s1s1 + u_s2s2 + FRAC*2 pi^2 u  on [0,1]^2,  u = 0 on the
%   boundary,
%
%   P = P0 (+ psatz terms) + diag(eppos) >= 0,
%   Q = (A'(P T))' + A'(P T) (+ 2 epneg T'(P T)) = -Qe,  Qe >= 0,
%
% with the 2-D settings of lpisettings(LEVEL) on an lpiprogram_sop program
% (BUILD_2D), solves it with SeDuMi, extracts P, Q, Qe and Qe + Q with
% lpigetsol_sop, and checks the extraction, erroring at the first failed
% check. The solve is held to no residual bound, so the test does not
% certify stability:
%   (0) the program equals, in expr.At, expr.b, var.idx, var.type,
%       decvartable and dom, the program of cx_stability_2D;
%   (1) SeDuMi returns a point: ||solinfo.x|| > 1e-12 and
%       rel_b < 1 - 1e-6, rel_b = ||At'x - b||/||b|| at x = solinfo.RRx;
%       the extracted operators are 'copvar' objects;
%   (2) the coefficient norm of the extracted Qe + Q lies in [r, sqrt(2) r]
%       (1e-6 slack), r the residual norm of all rows at solinfo.RRx;
%   (3) by the action on 6 random polynomial test functions (opcheck_sop,
%       10 Gauss nodes per direction), the extracted Q equals
%       (A'(P T))' + A'(P T) (+ 2 epneg T'(P T)) at the extracted P to
%       1e-10 relative;
%   (4) on the same test functions, Qe + Q = 0 to 1e-3 relative, and
%       <x, P x> and <x, Qe x> at least -1e-6 (normalized);
%   (5) with solinfo.RRx perturbed by 1e-3 times its rms entry (Gaussian,
%       seed 7), the residual of Qe + Q = 0 exceeds 10 times its
%       unperturbed value.
% The pinf, dinf, numerr and feasratio of SeDuMi, rel_b and the PSD margin
% are printed.
%
% INPUT
% - frac:   (optional) 'double', the reaction coefficient over 2 pi^2 (the
%           PDE is stable for FRAC < 1); default 0.1;
% - level:  (optional) 'char', the name of the LPI settings (lpisettings);
%           default 'light'.
%
% OUTPUT
% - out:    'struct' with fields info ('struct', the SeDuMi report),
%           rel_b, psd_min and psd_relmin ('double'), coef_rows ('double',
%           the coefficient norm of (2) over r), rows_res ('double', r),
%           res_Qgetsol and res_E ('double', the residuals of (3) and (4)),
%           psd (1x2 'double', the least <x,Px> and <x,Qe x>), control
%           ('double', the residual of (5)), shape ('struct' of cx_shape:
%           ndv, Kf, Ks, m, nvar), t ('struct' of seconds: ref, build,
%           solve, extract); printed lines with the values.
%
% NOTES
% A failed check errors with identifier 'stability_2d_sop:check'; 2-D
% settings with use_sosineq = 1 are an error ('stability_2d_sop:gap').
% Requires SeDuMi.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - e2e_stability_2d
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/29/2026. Tier 1 end-to-end example E3.
% MMP, 10/21/2026: moved from examples/stability_2d_sop.m to
%                tests/e2e_stability_2d.m: an end-to-end check built on test
%                helpers (cx_exec, opcheck_sop), not an example of standard
%                usage; examples/stability_2d_sop.m is now that. Code
%                unchanged; the error identifiers keep the former name.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was the measured solver
%                accuracy and cost, references to memory notes, the
%                provenance of the program, the rationale for the bounds of
%                check (2) and notes on keeping BUILD_2D in step with its
%                reference.
%                The (was) line of the former signature moved here from
%                under the function line, so that help starts with the
%                description:
% function out = stability_2d_sop(frac,level)                               % MMP, 10/21/2026 (was)

% Defaults; no warnings for non-canonical multipliers; a check asserts
% with one identifier
if nargin<1 || isempty(frac),   frac = 0.1;         end
if nargin<2 || isempty(level),  level = 'light';    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
chk = @(c,varargin) assert(c,'stability_2d_sop:check',varargin{:});
% Settings with SeDuMi; the PIE of the 2-D heat equation
st = cx_settings(level,'sedumi');
PIE = initialize(cx_plant('rd2',frac));

% % % Reference program (frozen file, unmodified) and the transcription.
t0 = tic;   pref = cx_stability_2D(PIE,st);     t_ref = toc(t0);
t0 = tic;   [prog,ops] = build_2d(PIE,st);      t_build = toc(t0);
% (0) the two programs must be equal; the shape of the SDP
same = isequal(pref.expr.At,prog.expr.At) && isequal(pref.expr.b,prog.expr.b) && ...
       isequal(pref.var.idx,prog.var.idx) && isequal(pref.var.type,prog.var.type) && ...
       isequal(pref.decvartable,prog.decvartable) && isequal(pref.dom,prog.dom);
S = cx_shape(prog);
fprintf('E3 rd2 %.2f/%s: m %d, ndv %d, Ks %s; build %.1f s (reference %.1f s); transcription = cx_stability_2D: %d\n',...
    frac,level,S.m,S.ndv,mat2str(S.Ks),t_build,t_ref,same);
chk(same,'(0) the transcription does not build cx_stability_2D''s program')

% % % Solve (SeDuMi).
t0 = tic;   prog = lpisolve(prog,st.sos_opts);  t_solve = toc(t0);
% (1) solver report; rel_b, PSD margin and whether x = 0
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

% Output
out = struct('info',info,'rel_b',rb,'psd_min',pmin,'psd_relmin',prel,'coef_rows',ce/rall,...
    'rows_res',rall,'res_Qgetsol',rQ,'res_E',rE,'psd',[pP pQe],'control',rEx,'shape',S,...
    't',struct('ref',t_ref,'build',t_build,'solve',t_solve,'extract',t_extract));
fprintf('E3 PASS (extraction; not a certificate)\n');
end


% ------------------------------------------------------------------------
function [prog,ops] = build_2d(PIE,st)
% [PROG,OPS] = BUILD_2D(PIE,ST) builds, unsolved, the LPI
% P = P0 (+ psatz terms) + diag(eppos),  Q = (A'(P T))' + A'(P T)
% (+ 2 epneg T'(P T)),  Qe + Q = 0  with P0, Qe and the psatz terms
% positive operators (cx_stability_pos2d) at the degrees and options of
% the 2-D settings, Qe pruned from the zero structure of Q.
% INPUT  PIE: 'pie_struct' object of a 2-D PIE whose state is in
%             L2[s1,s2];
%        st:  'struct' of LPI settings with fields sos_opts and
%             settings_2d, or 2-D settings with is2D true; use_sosineq = 1
%             is an error.
% OUTPUT prog: 'struct', the LPI program; ops: 'struct' with fields T, A
%        ('copvar'), P, Q, Qe ('cdopvar') and epneg ('double').
% 2-D settings, with the solver options of ST
if ~isfield(st,'is2D') || ~st.is2D
    sos_opts = st.sos_opts;     st = st.settings_2d;    st.sos_opts = sos_opts;
end
% Margins (one per space of T) and decay rate
if ~isfield(st,'eppos'),    eppos = [1e-4;1e-6;1e-6;1e-6];
else,                       eppos = st.eppos;   if numel(eppos)==1,  eppos = eppos*ones(4,1);  end
end
epneg = 0;  if isfield(st,'epneg'),    epneg = st.epneg;   end
% Degrees and options of the psatz terms of P
LF_use_psatz = st.LF_use_psatz;
LF_deg_psatz = psatz_cell(st.LF_deg_psatz,LF_use_psatz);
LF_opts_psatz = psatz_cell(st.LF_opts_psatz,LF_use_psatz);
if st.use_sosineq
    error('stability_2d_sop:gap','use_sosineq=1 needs lpi_ineq_2d, which has no container counterpart.')
end
% Degrees and options of Qe and of its psatz terms
eq_opts = st.eq_opts;   eq_deg = st.eq_deg;     eq_use_psatz = st.eq_use_psatz;
eq_deg_psatz = psatz_cell(st.eq_deg_psatz,eq_use_psatz);
eq_opts_psatz = psatz_cell(st.eq_opts_psatz,eq_use_psatz);

% The PIE operators as containers; the program
Tm = opvar2d2copvar(PIE.T);     Am = opvar2d2copvar(PIE.A);
prog = lpiprogram_sop(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);

% P = P0 (+ psatz terms) + diag(eppos) on the output space of T
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

% Q = (A'(P T))' + A'(P T) (+ 2 epneg T'(P T))
PTm = Pm*Tm;    APTm = Am'*PTm;
if epneg==0
    Qm = APTm' + APTm;
else
    Qm = APTm' + APTm + 2*epneg*(Tm'*PTm);
end

% Qe >= 0, pruned from the zero structure of Q (tolerance 1e-12), and its
% psatz terms, which exclude what Qe excludes
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
% Qe + Q = 0, one parameter per adjoint pair; the operators
prog = lpi_eq_sop(prog,Qem+Qm,'symmetric');
ops = struct('T',Tm,'A',Am,'P',Pm,'Q',Qm,'Qe',Qem,'epneg',epneg);
end

function outcell = psatz_cell(incell,use_psatz)
% OUTCELL = PSATZ_CELL(INCELL,USE_PSATZ) returns one entry of INCELL per
% psatz term: {} when every entry of USE_PSATZ is 0; INCELL repeated when
% it is a 'struct' or has one entry; INCELL itself when it has one entry
% per term; INCELL(USE_PSATZ) when it has at least max(USE_PSATZ)
% entries; any other INCELL is an error.
% INPUT  incell: 'struct' or 'cell' of degrees or options; use_psatz:
%        'double' vector, the psatz terms.
% OUTPUT outcell: 'cell', one entry per element of USE_PSATZ.
if all(use_psatz==0),                   outcell = {};   return, end
if isa(incell,'struct'),                outcell = repmat({incell},1,length(use_psatz));
elseif numel(incell)==1,                outcell = repmat(incell(1),1,length(use_psatz));
elseif numel(incell)==length(use_psatz),outcell = incell;
elseif numel(incell)>=max(use_psatz),   outcell = incell(use_psatz);
else,   error('For each element of ''use_psatz'', a field should be defined.')
end
end

function r = row_res(prog,rows)
% R = ROW_RES(PROG,ROWS) returns sqrt(sum_i ||At_i'*x - b_i||^2) over the
% expressions i in ROWS, x = prog.solinfo.RRx (decvartable order, the
% order of the rows of At_i).
% INPUT  prog: 'struct', a solved LPI program; rows: 'double' vector of
%        indices into prog.expr.At and prog.expr.b.
% OUTPUT r: 'double'.
x = prog.solinfo.RRx(:);    r2 = 0;
for i = rows
    At = prog.expr.At{i};   b = prog.expr.b{i};
    r2 = r2 + sum((At'*x(1:size(At,1)) - b).^2);
end
r = sqrt(full(r2));
end

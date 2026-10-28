function out = e2e_heat3d_stability(kappa,opts)                             % MMP, 10/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = E2E_HEAT3D_STABILITY(KAPPA,OPTS) certifies exponential stability
% at rate KAPPA of the heat equation in N = OPTS.N spatial variables,
%
%   u_t = u_s1s1 + ... + u_sNsN  on [0,1]^N,  u = 0 at s1 = 0 and s1 = 1,
%   u = 0 at si = 0 and u_si = 0 at si = 1 for i >= 2,
%
% whose exact decay rate is pi^2 + (N-1) pi^2/4 (3 pi^2/2 for N = 3), by
% the LPI of Jagt & Peet, arXiv:2508.14840v4, Cor. 35 (heatNd_lpi, preset
% 'bench', d = 0, ep = 0.1): for the PIE  T v_t = A0 v  (heatNd_pie), find
% a decision operator P of degree d, R >= 0 and Q >= 0 with
%
%   (15a)  P'*T - (P'*T)' = 0  and  P'*T - ep^2 T'*T - R = 0,
%   (15b)  X1 + KAPPA*X2 + Q = 0,  X1 = P'*A0 + (P'*A0)',
%          X2 = P'*T + (P'*T)',
%
% each equality on one parameter per adjoint pair ('symmetric'). The SDP
% (heatNd_sdp) is solved with MOSEK on a maximal linearly independent
% subset of its equality rows (heatNd_solve), P, R, Q and the constraint
% operators are extracted with lpigetsol_sop, and they are checked,
% erroring at the first failed check:
%   (0) a row set from OPTS.keepfile belongs to an SDP with the same m, K
%       and b;
%   (1) heatNd_solve returns the verdict +1: MOSEK OPTIMAL, rel_b =
%       ||At'x - b||/||b|| <= 1e-6 on all rows, the least eigenvalue of the
%       PSD blocks at least -1e-8 times the largest, x ~= 0; the extracted
%       operators are 'copvar' objects;
%   (2) at solinfo.RRx = RR*x*bscl (x the MOSEK point, RR and bscl from
%       heatNd_sdp), the residual norms r_a of the (15a) rows and r_b of
%       the (15b) rows satisfy sqrt(r_a^2 + r_b^2) = r_all to 1e-6
%       relative, r_all = ||At'x - b||*bscl; the coefficient norm of the
%       extracted (15a) operators lies in [r_a, 2 r_a] and that of (15b)
%       in [r_b, sqrt(2) r_b] (1e-6 slack);
%   (3) by the action on OPTS.nt random polynomial test functions of
%       degree at most 2 (opcheck_sop, OPTS.nq Gauss nodes per direction):
%       P'*T - ep^2 T'*T = R, P'*T = (P'*T)' and X1 + KAPPA*X2 = -Q to 1e-5
%       relative; R, Q and -(X1 + KAPPA*X2) at least -1e-8 (normalized);
%   (4) with solinfo.RRx perturbed by 1e-3 times its rms entry (Gaussian,
%       seed 7), the larger of the residuals of P'*T - ep^2 T'*T = R and
%       X1 + KAPPA*X2 = -Q exceeds 1e3 times the larger unperturbed one.
%
% INPUT
% - kappa:  (optional) 'double', the rate; default 14.0;
% - opts:   (optional) 'struct' with optional fields
%             keep      'double' vector, the equality rows to solve on;
%             keepfile  'char', a .mat file saved by heatNd_sdp whose
%                       info.keep is the row set;
%             nt        'double', the number of test functions (default 4);
%             nq        'double', Gauss nodes per direction (default 6);
%             N         'double', the number of spatial variables
%                       (default 3);
%           without keep and keepfile, the row set is computed by one
%           sparse QR (heatNd_lindep).
%
% OUTPUT
% - out:    'struct' with fields kappa and kstar ('double', the rate and
%           the exact rate), v ('struct', the verdict of heatNd_solve,
%           field x emptied), nkeep, m, ndec ('double', the rows solved on,
%           all rows, the decision variables), Ks ('double', the PSD block
%           sizes), rows ([r_a r_b r_all]), coef_rows (1x2 'double', the
%           coefficient norms over r_a and r_b), res (1x3 'double', the
%           residuals of (3)), psd (1x3 'double', the least values for R, Q
%           and -(X1 + KAPPA*X2)), control (1x2 'double', the residuals of
%           (4)), mem ('struct', the memory record of heatNd_lpi), t
%           ('struct' of seconds: build, sdp, keep, mosek, extract, check);
%           printed lines with the values.
%
% NOTES
% A failed check errors with identifier 'heat3d_stability_sop:check'. By
% Cor. 35, a certificate exists only for KAPPA at most the exact rate.
% Requires MOSEK.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - e2e_heat3d_stability
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
% Initial coding MMP, 09/29/2026. Tier 1 end-to-end example E4.
% MMP, 10/21/2026: moved from examples/heat3d_stability_sop.m to
%                tests/e2e_heat3d_stability.m: an end-to-end check built on
%                the heatNd benchmark functions and opcheck_sop, not an
%                example of standard usage. No standard-usage 3-D example
%                replaces it: the PDE front end converts at most two
%                spatial variables. Code unchanged; the error identifiers
%                keep the former name.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was the measured size,
%                cost and bisection bracket of the SDP, a comparison with
%                stock PIETOOLS, the derivation of the bounds of check (2)
%                and references to a README.
%                The (was) line of the former signature moved here from
%                under the function line, so that help starts with the
%                description:
% function out = heat3d_stability_sop(kappa,opts)                           % MMP, 10/21/2026 (was)

% Defaults; no warnings for non-canonical multipliers; a check asserts
% with one identifier
if nargin<1 || isempty(kappa),  kappa = 14.0;       end
if nargin<2 || isempty(opts),   opts = struct();    end
if ~isfield(opts,'nt'),     opts.nt = 4;    end
if ~isfield(opts,'nq'),     opts.nq = 6;    end
if ~isfield(opts,'N'),      opts.N = 3;     end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
chk = @(c,varargin) assert(c,'heat3d_stability_sop:check',varargin{:});
% epsilon of (15a), degree of P
ep = 0.1;   d = 0;

% % % Program and SDP.
% The PIE (r = 0), the LPI at KAPPA and its operators (meta.base), then
% the SDP in SeDuMi form
t0 = tic;
pie = heatNd_pie(opts.N,0);
[prog,meta] = heatNd_lpi(pie,d,kappa,ep,struct('preset','bench'));
t_build = toc(t0);
B = meta.base;
t0 = tic;   D = heatNd_sdp(prog);   t_sdp = toc(t0);
fprintf('E4 %d-D heat, kappa %.4f (kappa* %.6f): ndec %d, m %d, Ks %s; build %.1f s, SDP %.1f s\n',...
    opts.N,kappa,pie.exact.lambda1,numel(prog.decvartable),D.m,mat2str(D.K.s),t_build,t_sdp);

% % % Row set, then the solve (verdict on all rows).
% Rows given, read from a dump of the same SDP family (check (0)), or
% found by one sparse QR
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
% (1) MOSEK on the kept rows, the verdict on all rows
v = heatNd_solve(D,struct('keep',keep,'rows','keep','keepx',true));
fprintf('   heatNd_solve: st %+d (%s); rel_b %.2e (all rows), psd_relmin %.2e, MOSEK %.0f s, %d iter, peak %.0f MB\n',...
    v.st,v.why,v.rel_b,v.psd_relmin,v.t_mosek,v.iter,v.mem_peak);
chk(v.st==1,'(1) not certified: %s',v.why)

% % % Extraction: the dump's point in decvartable order is a solved program.
t0 = tic;
% Cone vector and decvartable-order point, both scaled back by bscl
prog.solinfo.x = v.x*D.bscl;
prog.solinfo.RRx = full(D.RR*v.x)*D.bscl;
prog.solinfo.info = struct('route','heatNd_solve','st',v.st,'rel_b',v.rel_b);
v.x = [];
% P, R, Q, the constraint operators of (15a), (15b) and the parts that
% check (3) uses
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

% Output
out = struct('kappa',kappa,'kstar',pie.exact.lambda1,'v',v,'nkeep',numel(keep),'m',D.m,...
    'ndec',numel(prog.decvartable),'Ks',D.K.s,'rows',[ra rb_ rall],'coef_rows',[ca/ra cb/rb_],...
    'res',[r15a rPT r15b],'psd',[pR pQ pL],'control',[r15ax r15bx],'mem',meta.mem,...
    't',struct('build',t_build,'sdp',t_sdp,'keep',t_keep,'mosek',v.t_mosek,'extract',t_extract,'check',t_check));
fprintf('E4 PASS: kappa %.4f certified, certificate extracted and checked\n',kappa);
end


% ------------------------------------------------------------------------
function x = v_x(D,prog)
% X = V_X(D,PROG) returns the solved cone vector in the scaling of D.b,
% prog.solinfo.x/D.bscl.
% INPUT  D: 'struct' of heatNd_sdp (field bscl); prog: 'struct', a solved
%        LPI program (field solinfo.x).
% OUTPUT x: 'double' column.
x = prog.solinfo.x/D.bscl;
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

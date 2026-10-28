function out = e2e_hinf_gain_1d(plant,level)                                % MMP, 10/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = E2E_HINF_GAIN_1D(PLANT,LEVEL) computes, on the container path, an
% upper bound gam on the H-infinity gain from w to z of the 1-D PIE
% T x_t = A x + B_w w,  z = C_z x + D_zw w  of PLANT, with gam a decision
% variable minimized in one SDP,
%
%   min gam  s.t.  gam >= 0,  R >= 0,  T'Q - R = 0,
%       K = [-gam*I_w, D_zw',    B_w'Q    ]
%           [ D_zw,    -gam*I_z, C_z      ]  <= 0,  i.e. K + N = 0, N >= 0,
%           [ Q'B_w,   C_z',     A'Q + Q'A]
%
% on an lpiprogram_sop program: gam from lpidecvar, gam >= 0 by lpi_ineq;
% R and N the positive operators of the 'lf' and 'slack' settings of
% lpisettings(LEVEL) (poslpivar_settings_sop); Q free (lpivar_cdopvar at
% the degrees get_lpivar_degs_sop reads off R); I_w, I_z from
% eye_copvar_sop; the two equalities by lpi_eq_sop, the second
% 'symmetric'. It solves the LPI with MOSEK, extracts gam, Q, R, N, K and
% the two constraint operators with lpigetsol_sop, and checks them,
% erroring at the first failed check:
%   (1) MOSEK numerr 0 and pinf 0; rel_b = ||At'x - b||/||b|| <= 1e-6 at
%       x = solinfo.RRx; the least eigenvalue of the Gram blocks at least
%       -1e-8 times the largest; the extracted operators are 'copvar'
%       objects;
%   (2) gam equals the bound of stock PIETOOLS_Hinf_gain on the same PIE
%       and settings to 1e-4 relative; for PLANT {'io1'} at 'light', gam
%       lies in [0.182503, 0.182626], widened by 1e-5 relative;
%   (3) with r1, r2 the residual norms at solinfo.RRx of the rows of
%       T'Q - R = 0 and K + N = 0, the coefficient norm of the extracted
%       T'Q - R equals r1 to 1e-6 relative, that of K + N lies in
%       [r2, sqrt(2) r2] (1e-6 slack);
%   (4) by the action on 6 random polynomial test functions x, y
%       (opcheck_sop, quadrature of the kernels): <T y, Q x> = <y, R x>
%       and N + K = 0 to 1e-5 relative; the extracted K equals K rebuilt
%       from the extracted gam and Q to 1e-10 relative; <x, R x> and
%       <x, N x> at least -1e-8, and <x, K x> at most the larger of 1e-8
%       and the residual of N + K (normalized);
%   (5) with the first entries of solinfo.x in place of solinfo.RRx, the
%       larger of the residuals of <T y, Q x> = <y, R x> and N + K = 0
%       exceeds 1e3 times the larger one at solinfo.RRx.
%
% INPUT
% - plant:  (optional) 'cell', the arguments of cx_plant naming a 1-D
%           plant with Tw = 0; default {'io1'}, the PIE of
%           u_t = u_ss + 0.5 pi^2 u + w(s) on [0,1], u(0) = u(1) = 0,
%           z = int_0^1 u ds;
% - level:  (optional) 'char', the name of the LPI settings (lpisettings);
%           default 'light'.
%
% OUTPUT
% - out:    'struct' with fields gam and gam_stock ('double', the container
%           and the stock bound), rel_b, psd_min and psd_relmin ('double',
%           of the container solve), coef_rows (1x2 'double', the
%           coefficient norms over r1 and r2), rows_res ([r1 r2]),
%           weak_TQR, res_KN and res_Kgetsol ('double', the residuals of
%           (4)), psd (1x3 'double', the least <x,Rx> and <x,Nx> and the
%           largest <x,Kx>), control (1x2 'double', the residuals of (5)),
%           shape_stock and shape ('struct' of cx_shape: ndv, Kf, Ks, m,
%           nvar), t ('struct' of seconds: stock, build, solve, extract);
%           printed lines with the values.
%
% NOTES
% A 2-D plant or Tw ~= 0 is an error ('hinf_gain_1d_sop:plant'); a failed
% check errors with identifier 'hinf_gain_1d_sop:check'. Requires MOSEK.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - e2e_hinf_gain_1d
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
% Initial coding MMP, 09/29/2026. Tier 1 end-to-end example E1.
% MMP, 10/06/2026: cx_hinf_lf -> poslpivar_settings_sop ('lf'), cx_hinf_qdeg
%                -> get_lpivar_degs_sop, cx_space_list -> copvar_space_list,
%                and the re-typed cx_hinf_slack block -> poslpivar_settings_sop
%                ('slack'), so the example runs on the library translators
%                (lpi_programming_sopvar); programs unchanged (126-program
%                bit-identity check). 'cx_hinf_lf', 'cx_hinf_qdeg' and
%                'cx_hinf_slack' above now refer to those library routines.
%                Check (0) ('checked to build cx_hinf_slack's program
%                exactly', above) is superseded and removed: the example and
%                the transcriptions now share one routine.
% MMP, 10/21/2026: moved from examples/hinf_gain_1d_sop.m to
%                tests/e2e_hinf_gain_1d.m: an end-to-end check built on test
%                helpers (cx_exec, opcheck_sop), not an example of standard
%                usage; examples/hinf_gain_1d_sop.m is now that. Code
%                unchanged; the error identifiers keep the former name.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and the header of the local function; code
%                unchanged. The deleted header text was the stock executive
%                the LPI reproduces, the provenance of the program and of
%                the bracket of check (2), the rationale for the
%                tolerances, and the path of the helpers. The references of
%                the 10/06/2026 entry to the header text above no longer
%                apply.
%                The (was) line of the former signature moved here from
%                under the function line, so that help starts with the
%                description:
% function out = hinf_gain_1d_sop(plant,level)                              % MMP, 10/21/2026 (was)

% Defaults; no warnings for non-canonical multipliers; settings with
% MOSEK; the PIE, 1-D with Tw = 0; a check asserts with one identifier
if nargin<1 || isempty(plant),  plant = {'io1'};    end
if nargin<2 || isempty(level),  level = 'light';    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
st = cx_settings(level,'mosek');
PIE = initialize(cx_plant(plant{:}));
if PIE.dim~=1 || ~(PIE.Tw==0)
    error('hinf_gain_1d_sop:plant','A 1-D PIE with Tw = 0 is required (the executive''s non-coercive branch).')
end
chk = @(c,varargin) assert(c,'hinf_gain_1d_sop:check',varargin{:});

% % % Stock executive, unmodified.
t0 = tic;
[prog_s,~,gam_s] = PIETOOLS_Hinf_gain(PIE,st);
t_stock = toc(t0);

% % % Container program, gamma a decision variable.
t0 = tic;
% The PIE operators as containers; identity operators on the spaces of w
% and z
Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.Bw);  Cz = opvar2copvar(PIE.Cz);  Dzw = opvar2copvar(PIE.Dzw);
% [spw,dmw] = cx_space_list(Bw,'in');     [spz,dmz] = cx_space_list(Cz,'out'); % MMP, 10/06/2026 (was)
[spw,dmw] = copvar_space_list(Bw,'in');     [spz,dmz] = copvar_space_list(Cz,'out'); % MMP, 10/06/2026
Iw = eye_copvar_sop(dmw,spw,PIE.dom);   Iz = eye_copvar_sop(dmz,spz,PIE.dom);
% Program; gam >= 0 a decision variable, the objective
prog = lpiprogram_sop(Tm);                          % registry and dom of T
[prog,gam] = lpidecvar(prog,'gam');
prog = lpi_ineq(prog,gam);                          % gam >= 0: legacy sosineq
prog = lpisetobj(prog,gam);
% R >= 0 ('lf' settings) and Q free on the output spaces of T; T'Q - R = 0,
% its rows recorded
% [prog,Rm] = cx_hinf_lf(prog,Tm,PIE,st);             % R >= 0              % MMP, 10/06/2026 (was)
% [sp,dm] = cx_space_list(Tm,'out');                                        % MMP, 10/06/2026 (was)
% [prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,cx_hinf_qdeg(Rm));          % MMP, 10/06/2026 (was)
[prog,Rm] = poslpivar_settings_sop(prog,Tm,st,'lf','out',PIE.dom); % R >= 0 % MMP, 10/06/2026
[sp,dm] = copvar_space_list(Tm,'out');                                      % MMP, 10/06/2026
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,get_lpivar_degs_sop(Rm));     % MMP, 10/06/2026
E1 = Tm'*Qm - Rm;
e0 = prog.expr.num;
prog = lpi_eq_sop(prog,E1);                         % NOT symmetric (stock line 157)
rows1 = e0+1:prog.expr.num;
% KYP operator K, affine in gam and Q
Km = [-(gam*Iw),   Dzw',        Bw'*Qm;
       Dzw,        -(gam*Iz),   Cz;
       Qm'*Bw,     Cz',         Am'*Qm + Qm'*Am];
% BEGIN MMP, 10/06/2026: N by the library poslpivar_settings_sop ('slack'),
% which returns N, in place of the re-typed cx_hinf_slack block (kept below,
% commented out). prog0 served check (0) only, which is removed.
% % N >= 0 on K's spaces, as cx_hinf_slack (which does not return N).       % MMP, 10/06/2026 (was)
% prog0 = prog;                                                             % MMP, 10/06/2026 (was)
% [spk,dmk] = cx_space_list(Km,'out');                                      % MMP, 10/06/2026 (was)
% [deg,co] = cx_hinf_posdeg(st.dd2,st.options2,spk);                        % MMP, 10/06/2026 (was)
% [prog,Nm] = poscopvar(prog,dmk,spk,PIE.dom,deg,co);                       % MMP, 10/06/2026 (was)
% if st.override2~=1                                                        % MMP, 10/06/2026 (was)
%     [deg,co] = cx_hinf_posdeg(st.dd3,st.options3,spk);                    % MMP, 10/06/2026 (was)
%     [prog,N2] = poscopvar(prog,dmk,spk,PIE.dom,deg,co);                   % MMP, 10/06/2026 (was)
%     Nm = Nm + N2;                                                         % MMP, 10/06/2026 (was)
% end                                                                       % MMP, 10/06/2026 (was)
% N >= 0 on K's spaces (dd2/options2, + dd3/options3 when override2 ~= 1).  % MMP, 10/06/2026
[prog,Nm] = poslpivar_settings_sop(prog,Km,st,'slack','out',PIE.dom);       % MMP, 10/06/2026
% END MMP, 10/06/2026
% K + N = 0, one parameter per adjoint pair, its rows recorded
E2 = Nm + Km;
e0 = prog.expr.num;
prog = lpi_eq_sop(prog,E2,'symmetric');             % Deop + Dop = 0
rows2 = e0+1:prog.expr.num;
t_build = toc(t0);
% Check (0) removed: the transcriptions now build N by the same routine.    % MMP, 10/06/2026
% pchk = cx_hinf_slack(prog0,Km,PIE,st);                                    % MMP, 10/06/2026 (was)
% chk(isequal(pchk.expr,prog.expr) && isequal(pchk.var,prog.var) && isequal(pchk.decvartable,prog.decvartable),...
%     '(0) the slack transcription does not build cx_hinf_slack''s program') % MMP, 10/06/2026 (was)

% % % One solve.
t0 = tic;
prog = lpisolve(prog,st.sos_opts);
t_solve = toc(t0);
% (1) solver flags, rel_b and the PSD margin; the shapes of both SDPs
info = prog.solinfo.info;
[rb,pmin,prel] = cx_resid(prog);
Ss = cx_shape(prog_s);  Sc = cx_shape(prog);
fprintf('E1 %s/%s: container gamma solve: numerr %d, pinf %d, rel_b %.2e, psd_relmin %.2e; build %.1f s, solve %.1f s\n',...
    strjoin(cellfun(@num2str,plant,'uni',0),','),level,info.numerr,info.pinf,rb,prel,t_build,t_solve);
fprintf('   shape stock: ndv %d, m %d, Ks %s | container: ndv %d, m %d, Ks %s\n',...
    Ss.ndv,Ss.m,mat2str(Ss.Ks),Sc.ndv,Sc.m,mat2str(Sc.Ks));
chk(info.numerr==0 && info.pinf==0 && rb<=1e-6 && prel>=-1e-8,...
    '(1) solve not certified: numerr %d, pinf %d, rel_b %.3g, psd_relmin %.3g',info.numerr,info.pinf,rb,prel)

% % % Extraction.
t0 = tic;
g = double(lpigetsol_sop(prog,gam));
Qs = lpigetsol_sop(prog,Qm);    Rs = lpigetsol_sop(prog,Rm);
Ns = lpigetsol_sop(prog,Nm);    Ks = lpigetsol_sop(prog,Km);
E1s = lpigetsol_sop(prog,E1);   E2s = lpigetsol_sop(prog,E2);
t_extract = toc(t0);
chk(all(cellfun(@(X) isa(X,'copvar'),{Qs,Rs,Ns,Ks,E1s,E2s})),'(1) extracted operators not fixed')

% (2) gamma against stock and the bisection bracket.
fprintf('   gamma: container %.8f, stock %.8f (rel diff %.1e)\n',g,gam_s,abs(g-gam_s)/gam_s);
chk(abs(g-gam_s)<=1e-4*gam_s,'(2) gamma %.8g differs from stock %.8g',g,gam_s)
if isequal(plant,{'io1'}) && strcmp(level,'light')
    chk(g>=0.182503*(1-1e-5) && g<=0.182626*(1+1e-5),'(2) gamma %.8g outside the bisection bracket',g)
end

% (3) extraction against the SDP rows.
r1 = row_res(prog,rows1);   r2 = row_res(prog,rows2);
c1 = opcheck_sop('coef',E1s);   c2 = opcheck_sop('coef',E2s);
fprintf('   ||coef E(sol)|| / ||r_rows||: T''Q-R %.6f (rows %.2e), K+N %.4f (rows %.2e)\n',c1/r1,r1,c2/r2,r2);
chk(abs(c1/r1-1)<=1e-6,'(3) T''Q-R: coefficient norm %.6g vs row residual %.6g',c1,r1)
chk(c2>=r2*(1-1e-6) && c2<=sqrt(2)*r2*(1+1e-6),'(3) K+N: coefficient norm %.6g vs row residual %.6g',c2,r2)

% (4) semantic certificate.
o = struct('nt',6,'nqi',12,'nqo',14);
w1 = opcheck_sop('weak',{1,Tm,Qs; -1,[],Rs},Rs,o);     % <Ty,Qx> - <y,Rx>
rN = opcheck_sop('res',Ns,-Ks,o);                       % N = -K
Kf = [(-g)*Iw,     Dzw',       Bw'*Qs;                  % K(gamma, Q) rebuilt
      Dzw,         (-g)*Iz,    Cz;
      Qs'*Bw,      Cz',        Am'*Qs + Qs'*Am];
rK = opcheck_sop('res',Ks,Kf,o);
pR = opcheck_sop('psd',Rs,o);   pN = opcheck_sop('psd',Ns,o);   pK = -opcheck_sop('psd',-Ks,o);
fprintf(['   by action: T''Q=R weak %.2e | K+N %.2e | getsol(K) vs K(gamma,getsol(Q)) %.2e | '...
         'min <x,Rx> %.2e, min <x,Nx> %.2e, max <x,Kx> %.2e (normalized)\n'],w1,rN,rK,pR,pN,pK);
chk(w1<=1e-5 && rN<=1e-5,'(4) the extracted operators violate the equalities (%.3g, %.3g)',w1,rN)
chk(rK<=1e-10,'(4) getsol(K) differs from K rebuilt at the extracted gamma and Q (%.3g)',rK)
chk(pR>=-1e-8 && pN>=-1e-8 && pK<=max(1e-8,rN),'(4) positivity fails on test functions (%.3g, %.3g, %.3g)',pR,pN,pK)

% (5) control: the solver's cone vector in RRx's place.
bad = prog;     nd = numel(prog.decvartable);
bad.solinfo.RRx = prog.solinfo.x(1:nd);
w1x = opcheck_sop('weak',{1,Tm,lpigetsol_sop(bad,Qm); -1,[],lpigetsol_sop(bad,Rm)},Rs,o);
rNx = opcheck_sop('res',lpigetsol_sop(bad,Nm),-lpigetsol_sop(bad,Km),o);
fprintf('   control solinfo.x: T''Q=R weak %.2e, K+N %.2e\n',w1x,rNx);
chk(max(w1x,rNx)>1e3*max([w1,rN,1e-12]),'(5) reading solinfo.x is not detected')

% Output
out = struct('gam',g,'gam_stock',gam_s,'rel_b',rb,'psd_min',pmin,'psd_relmin',prel,...
    'coef_rows',[c1/r1, c2/r2],'rows_res',[r1 r2],'weak_TQR',w1,'res_KN',rN,'res_Kgetsol',rK,...
    'psd',[pR pN pK],'control',[w1x rNx],'shape_stock',Ss,'shape',Sc,...
    't',struct('stock',t_stock,'build',t_build,'solve',t_solve,'extract',t_extract));
fprintf('E1 PASS: gamma %.8f (stock %.8f)\n',g,gam_s);
end


% ------------------------------------------------------------------------
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

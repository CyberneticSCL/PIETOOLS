function R = test_dpvar_hinf_chain_sop(plant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_DPVAR_HINF_CHAIN_SOP(PLANT) poses, on the container path, the
% H-infinity gain LPI of the 1-D PIE  T x_t = A x + B_w w,
% z = C_z x + D_zw w  of PLANT, with the gain gamma a 'dpvar' decision
% variable minimized in one SDP,
%
%   min gam  s.t.  gam >= 0,  R >= 0,  T'Q - R = 0,  K + N = 0,  N >= 0,
%   K = [-gam*I_w, D_zw',    B_w'Q;
%         D_zw,    -gam*I_z, C_z;
%         Q'B_w,   C_z',     A'Q + Q'A],
%
% on an lpiprogram program: gam from lpidecvar, gam >= 0 by lpi_ineq; R
% and N the positive operators of the 'lf' and 'slack' settings
% (poslpivar_settings_sop); Q free (lpivar_cdopvar at the degrees
% get_lpivar_degs_sop reads off R); I_w, I_z from eye_copvar_sop; the two
% equalities by lpi_eq_cdopvar, the second 'symmetric'. It solves the LPI
% and stock PIETOOLS_Hinf_gain on the same PIE, both at
% lpisettings('light') with MOSEK, and checks:
%   (1) K is a 'cdopvar' with gam among its decision variables; when z is
%       in R^1, K equals (eq) the operator with the 'dpvar' -gam in place
%       of -gam*I_z;
%   (2) the two SDPs have the same number of decision variables and the
%       same PSD block sizes;
%   (3) both solves end with MOSEK numerr 0, the container solve also with
%       pinf 0;
%   (4) the two values of gamma agree to 1e-4 relative.
% The rel_b and PSD margins of both solves and the shapes of both SDPs
% (decision variables, free variables, rows, PSD blocks) are printed.
%
% INPUT
% - plant:  (optional) 'cell', the arguments of cx_plant naming a 1-D plant
%           with disturbance w and output z; default {'io1'}, the PIE of
%           u_t = u_ss + 0.5 pi^2 u + w(s) on [0,1], u(0) = u(1) = 0,
%           z = int_0^1 u ds.
%
% OUTPUT
% - R:      'struct' with fields gam_s and gam_c ('double', the stock and
%           container gamma), rel_b and psd_min ('double', of the container
%           solve), Ss and Sc ('struct' of cx_shape for the stock and the
%           container SDP: fields ndv, Kf, Ks, m, nvar); one printed line
%           per failed check, then the counts.
%
% NOTES
% After the last check, a failed check is an error
% ('test_dpvar_hinf_chain_sop:failed'). Requires MOSEK.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_dpvar_hinf_chain_sop
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
% Initial coding MMP, 09/29/2026
% MMP, 10/06/2026: cx_hinf_lf, cx_hinf_qdeg, cx_space_list and
%                  cx_hinf_slack -> library poslpivar_settings_sop ('lf';
%                  'slack' + the equality cx_hinf_slack imposed),
%                  get_lpivar_degs_sop, copvar_space_list, so the test runs
%                  on the library translators (lpi_programming_sopvar);
%                  programs unchanged (126-program bit-identity check).
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), the licence block added, step comments
%                  in the body and the header of the nested function; code
%                  unchanged. The deleted header text was the provenance of
%                  the program, design notes on the identity operators and
%                  a measured difference of the row counts.

% Default plant; no warnings for non-canonical multipliers; check counters
if nargin<1 || isempty(plant),  plant = {'io1'};    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
npass = 0;  nfail = 0;
    function ck(c,msg)
        % CK(C,MSG) counts a passed check when C is true, else a failed
        % one, and prints MSG.
        % INPUT  c: 'logical'; msg: 'char'.
        % OUTPUT none; npass or nfail of the enclosing function grows by 1.
        if c,   npass = npass+1;
        else,   nfail = nfail+1;    fprintf('  FAIL: %s\n',msg);
        end
    end

% Settings ('light', MOSEK, no simplification) and the PIE
st = cx_settings('light','mosek');
PIE = initialize(cx_plant(plant{:}));

% % % Stock executive, unmodified.
t = tic;
[prog_s,~,gam_s] = PIETOOLS_Hinf_gain(PIE,st);
ts = toc(t);
% Its solver report, rel_b and PSD margin
info_s = prog_s.solinfo.info;
[rb_s,pm_s] = cx_resid(prog_s);

% % % Container program, gamma a decision variable.
t = tic;
% The PIE operators as containers; identity operators on the spaces of w
% and z
Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.Bw);  Cz = opvar2copvar(PIE.Cz);
Dzw = opvar2copvar(PIE.Dzw);
% [spw,dmw] = cx_space_list(Bw,'in');     [spz,dmz] = cx_space_list(Cz,'out'); % MMP, 10/06/2026 (was)
[spw,dmw] = copvar_space_list(Bw,'in');     [spz,dmz] = copvar_space_list(Cz,'out'); % MMP, 10/06/2026
Iw = eye_copvar_sop(dmw,spw,PIE.dom);   Iz = eye_copvar_sop(dmz,spz,PIE.dom);
% Program on the variables of the PIE; gam >= 0 a decision variable, the
% objective
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
gam = dpvar('gam');
prog = lpidecvar(prog,gam);
prog = lpi_ineq(prog,gam);
prog = lpisetobj(prog,gam);
% R >= 0 ('lf' settings) and Q free on the output spaces of T; T'Q - R = 0
% [prog,Rm] = cx_hinf_lf(prog,Tm,PIE,st);                                   % MMP, 10/06/2026 (was)
% Qdeg = cx_hinf_qdeg(Rm);                                                  % MMP, 10/06/2026 (was)
% [sp,dm] = cx_space_list(Tm,'out');                                        % MMP, 10/06/2026 (was)
[prog,Rm] = poslpivar_settings_sop(prog,Tm,st,'lf','out',PIE.dom);          % MMP, 10/06/2026
Qdeg = get_lpivar_degs_sop(Rm);                                             % MMP, 10/06/2026
[sp,dm] = copvar_space_list(Tm,'out');                                      % MMP, 10/06/2026
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);
% (1) KYP operator K; with z in R^1, the form with the raw 'dpvar' entry
% -gam must equal it
Km = [-gam*Iw,     Dzw',       Bw'*Qm;
       Dzw,        -gam*Iz,    Cz;
       Qm'*Bw,     Cz',        Am'*Qm + Qm'*Am];
if all(cellfun(@isempty,spz)) && sum(dmz)==1
    Km2 = [-gam*Iw,     Dzw',   Bw'*Qm;
            Dzw,        -gam,   Cz;
            Qm'*Bw,     Cz',    Am'*Qm + Qm'*Am];
    ck(eq(Km,Km2),'KYP operator: -gam*Iz equals the raw dpvar entry -gam');
end
ck(isa(Km,'cdopvar') && any(strcmp(Km.Zd,'gam')),'KYP operator carries gam');
% N >= 0 ('slack' settings) on the spaces of K; N + K = 0
% prog = cx_hinf_slack(prog,Km,PIE,st);                                     % MMP, 10/06/2026 (was)
[prog,Nm] = poslpivar_settings_sop(prog,Km,st,'slack','out',PIE.dom);       % MMP, 10/06/2026
prog = lpi_eq_cdopvar(prog,Nm + Km,'symmetric');   % Deop + Dop = 0         % MMP, 10/06/2026
ta = toc(t);
% Solve the container program
t = tic;
prog = lpisolve(prog,st.sos_opts);
tc = toc(t);
% Both shapes post-solve: sossolve adds the slack block of an 'ineq' (gam >= 0).
Ss = cx_shape(prog_s);  Sc = cx_shape(prog);
% Solver report, rel_b, PSD margins and gamma of the container solve
info_c = prog.solinfo.info;
[rb_c,pm_c,pr_c] = cx_resid(prog);
gam_c = double(lpigetsol(prog,gam));

% (2) SDP shapes, (3) solver flags, (4) gamma
ck(Sc.ndv==Ss.ndv,sprintf('decision variables: container %d, stock %d',Sc.ndv,Ss.ndv));
ck(isequal(Sc.Ks,Ss.Ks),sprintf('PSD blocks: container %s, stock %s',mat2str(Sc.Ks),mat2str(Ss.Ks)));
ck(info_s.numerr==0,sprintf('stock numerr %d',info_s.numerr));
ck(info_c.numerr==0 && info_c.pinf==0,sprintf('container numerr %d pinf %d',info_c.numerr,info_c.pinf));
ck(abs(gam_c-gam_s)<=1e-4*abs(gam_s),sprintf('gamma: container %.8g, stock %.8g',gam_c,gam_s));
% Print the values, the shapes and the counts; the output; an error when a
% check failed
fprintf(['%s: gamma stock %.8g (numerr %d, rel_b %.2e, psd_min %.2e, %.1f s) | container %.8g ' ...
         '(numerr %d, rel_b %.2e, psd_min %.2e, psd_relmin %.2e; assembly %.2f s, solve %.1f s)\n'], ...
        strjoin(cellfun(@num2str,plant,'UniformOutput',false),','),gam_s,info_s.numerr,rb_s,pm_s,ts, ...
        gam_c,info_c.numerr,rb_c,pm_c,pr_c,ta,tc);
fprintf('  shape stock: ndv %d, Kf %d, m %d, Ks %s\n',Ss.ndv,Ss.Kf,Ss.m,mat2str(Ss.Ks));
fprintf('  shape cont.: ndv %d, Kf %d, m %d, Ks %s\n',Sc.ndv,Sc.Kf,Sc.m,mat2str(Sc.Ks));
fprintf('test_dpvar_hinf_chain_sop: %d passed, %d failed\n',npass,nfail);
R = struct('gam_s',gam_s,'gam_c',gam_c,'rel_b',rb_c,'psd_min',pm_c,'Ss',Ss,'Sc',Sc);
if nfail>0
    error('test_dpvar_hinf_chain_sop:failed','%d checks failed.',nfail);
end
end

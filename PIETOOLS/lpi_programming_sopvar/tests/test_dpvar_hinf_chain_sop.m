function R = test_dpvar_hinf_chain_sop(plant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_DPVAR_HINF_CHAIN_SOP(PLANT) chain test of Tier 1b/1c through
% lpi_programming and an executive: the 1-D primal H-infinity gain LPI of
% PIETOOLS_Hinf_gain, posed on the container path with gamma a DECISION
% variable and minimized in ONE SDP, against the stock executive.
%
% The container program is cx_Hinf_gain (cx_exec) with the fixed double
% gamma replaced exactly as the stock executive declares it:
%   dpvar gam; lpidecvar; lpi_ineq(prog,gam); lpisetobj(prog,gam)
%   Km = [-gam*Iw, Dzw', Bw'*Q; Dzw, -gam*Iz, Cz; Q'*Bw, Cz', A'*Q+Q'*A]
% where -gam*Iw is dpvar*copvar and Iw, Iz come from eye_copvar_sop, not
% opvar2copvar(mat2opvar(...)). A second KYP operator writes the z entry as
% a raw dpvar in [ ] ([Dzw, -gam, Cz]) and must equal the first (container
% eq, affine in the decision variables, matched by name).
%
% Checks: the two KYP forms are equal; the container program has the stock
% shape (decision variables, PSD block sizes; the row counts are reported:
% the container has 5 fewer, as in cx_run_1d); both solves
% end with MOSEK numerr = 0, and the container's rel_b on solinfo.RRx rows
% and PSD margin are reported (cx_resid); the two gammas agree to 1e-4
% relative.
%
% PLANT: cx_plant argument cell, default {'io1'}.
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(plant),  plant = {'io1'};    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
npass = 0;  nfail = 0;
    function ck(c,msg)
        if c,   npass = npass+1;
        else,   nfail = nfail+1;    fprintf('  FAIL: %s\n',msg);
        end
    end

st = cx_settings('light','mosek');
PIE = initialize(cx_plant(plant{:}));

% % % Stock executive, unmodified.
t = tic;
[prog_s,~,gam_s] = PIETOOLS_Hinf_gain(PIE,st);
ts = toc(t);
info_s = prog_s.solinfo.info;
[rb_s,pm_s] = cx_resid(prog_s);

% % % Container program, gamma a decision variable.
t = tic;
Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.Bw);  Cz = opvar2copvar(PIE.Cz);
Dzw = opvar2copvar(PIE.Dzw);
[spw,dmw] = cx_space_list(Bw,'in');     [spz,dmz] = cx_space_list(Cz,'out');
Iw = eye_copvar_sop(dmw,spw,PIE.dom);   Iz = eye_copvar_sop(dmz,spz,PIE.dom);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
gam = dpvar('gam');
prog = lpidecvar(prog,gam);
prog = lpi_ineq(prog,gam);
prog = lpisetobj(prog,gam);
[prog,Rm] = cx_hinf_lf(prog,Tm,PIE,st);
Qdeg = cx_hinf_qdeg(Rm);
[sp,dm] = cx_space_list(Tm,'out');
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);
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
prog = cx_hinf_slack(prog,Km,PIE,st);
ta = toc(t);
t = tic;
prog = lpisolve(prog,st.sos_opts);
tc = toc(t);
% Both shapes post-solve: sossolve adds the slack block of an 'ineq' (gam >= 0).
Ss = cx_shape(prog_s);  Sc = cx_shape(prog);
info_c = prog.solinfo.info;
[rb_c,pm_c,pr_c] = cx_resid(prog);
gam_c = double(lpigetsol(prog,gam));

ck(Sc.ndv==Ss.ndv,sprintf('decision variables: container %d, stock %d',Sc.ndv,Ss.ndv));
ck(isequal(Sc.Ks,Ss.Ks),sprintf('PSD blocks: container %s, stock %s',mat2str(Sc.Ks),mat2str(Ss.Ks)));
ck(info_s.numerr==0,sprintf('stock numerr %d',info_s.numerr));
ck(info_c.numerr==0 && info_c.pinf==0,sprintf('container numerr %d pinf %d',info_c.numerr,info_c.pinf));
ck(abs(gam_c-gam_s)<=1e-4*abs(gam_s),sprintf('gamma: container %.8g, stock %.8g',gam_c,gam_s));
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

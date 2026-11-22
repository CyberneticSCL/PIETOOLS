function pc_run(tag)
% PC_RUN(TAG) runs the jobs of TAG for the library examples that no tested
% setting certifies (problem cases, maintainer 11/20/2026: "create and run a
% plan for trying different settings variations and trades to certify these
% problem cases"; plan cleaning/private/problem_cases_1120/PLAN.md), in a
% worktree copy of the 11/18 study harness.
%
% INPUT
%   tag:    char:
%     a tag of pc_jobs.csv  stability LPI of a library row per job (diag_amt:
%                           one build and solve; variant 'A' = the row as is);
%                           results_<TAG>.csv, jobs already there skipped
%     'EXS_<name>'          ex_search over the starts pc_starts_<name>.csv
%                           (H2 dual rows; descents from a certified start)
%     'TRADE3'              ex_search trades on the one-variable rows at
%                           which TRADE1 stopped (statement split, fixed in
%                           apply_mods 11/20/2026)
%
% NOTES
%   - lb_setup first: restoredefaultpath, then this checkout alone on the
%     path, single resolution asserted. diag_amt and ex_search set no path,
%     and the saved pathdef holds the main checkout.
%   - pc_jobs.csv columns: tag, key, form, tier, mods, variant, param,
%     ndvcap; mods are MATLAB statements on st without commas; param empty
%     for the library default; ndvcap the variable cap (default 2e6).
%   - mods naming the solver 'sdpt3' or 'sedumi' put SDPT3 4.0 or SeDuMi
%     1.3 on the path first (neither survives lb_setup).
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
fprintf('## PC_RUN %s %s\n## lpisettings_sop: %s\n',tag,char(datetime('now')),which('lpisettings_sop'));
if strcmp(tag,'TRADE3')
    sf = fullfile(here,'..','..','..','..','..','cleaning','private','example_settings_1119','phaseB3_starts.csv');
    T = readtable(fullfile(here,'exs_TRADE1_final.csv'),'Delimiter',',','TextType','char');
    keys = T.key(contains(T.reason,'incomplete'))';
    ex_search(sf,'TRADE3',struct('keys',{keys},'statedir','_trade3','maxbuild',120));
    return
end
if startsWith(tag,'EXS_')
    sf = fullfile(here,['pc_starts_',tag(5:end),'.csv']);
    T = readtable(sf,'Delimiter',',','TextType','char');
    nb = 40;   if contains(tag,'50d'),  nb = 16;  end     % H2 dual builds: up to 412 s each
    ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',45000,'ndvcap',2e6, ...
              'rowmin',40,'maxbuild',nb,'statedir',['_',lower(tag)]));
    return
end
if strcmp(tag,'DIAG41')             % PDE 41: which part of the system defeats the LPI
    pc_diag41();
    return
end
if strcmp(tag,'DIAG41B')            % PDE 41 with an explicit delay state: settings
    pc_diag41b();
    return
end
if strcmp(tag,'DIAG44')             % PDE 44: which part of the formulation defeats the LPI
    pc_diag44();
    return
end
if strcmp(tag,'DIAG69')             % 69, 70 in one variable: position or energy states
    pc_diag69();
    return
end
if strcmp(tag,'DIAG44D')            % the PDE 44 variants in the dual forms
    pc_diag44({'Pd',''; 'Qd',''; 'Pd','st.epneg = 1e-3'},'DIAG44D');
    return
end
if strcmp(tag,'DIAG70')             % 70 with the velocity of order (2,2), as PDE 38
    pc_diag70();
    return
end
if strcmp(tag,'DIAG38')             % the order of the velocity alone, on the wave of PDE 38
    pc_diag38();
    return
end
if strcmp(tag,'DIAG44E')            % the minimal form of 44 in 'Pd': margins, tiers, X3 layout
    pc_diag44e();
    return
end
if strcmp(tag,'DIAG44F')            % X3 and 44 minimal forms: exponential rows; a damped wave in 'Q'
    pc_diag44f();
    return
end
if strcmp(tag,'DIAG2C')             % which coupling defeats the 2-D 'Q' test of 69 and 70
    pc_diag2c();
    return
end
if strcmp(tag,'WITNESS')            % the numerical witness of the position forms against N
    pc_witness();
    return
end
if strcmp(tag,'DIAG2M')             % 69, 70 in a mixed position-strain form
    pc_diag2m();
    return
end
if strcmp(tag,'DIAG56')             % 43 and 56 with the delay state of order 0 in space
    pc_diag56();
    return
end
if strcmp(tag,'DIAG69M')            % 69 in the mixed form with the strain of order (1,2)
    pc_diag69m();
    return
end
if strcmp(tag,'DIAG69N')            % settings around the 69mix2 certificate (measure 0.722)
    pc_diag69m({'st.lpi.ineq.dD = 0; st.eppos = 1e-3'; 'st.lpi.ineq.dD = 0; st.eppos = 1e-2'; 'st.lpi.ineq.dD = 0'; 'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto'''; 'st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1'}, ...
               'DIAG69N',{'light','light','stripped','light','light'});
    return
end
if strcmp(tag,'DIAG70N')            % setting cuts on the 70mix certificate (187,116 at light dD 0)
    pc_diag2m({'Q','st.lpi.ineq.dD = 0','stripped'; 'Q','st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''','light'; 'Q','st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = ''none''','light'; 'Q','st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''','stripped'; 'Q','','light'}, ...
              'DIAG70N',{'70mix'});
    return
end
if strcmp(tag,'DIAG69P')            % other programs for 69mix2 (light dD 0: measure 0.722)
    pc_diag69m({''; 'st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = ''faces1'''; 'st.lpi.ineq.dD = 0'},'DIAG69P',{'light','light','heavy'});
    return
end
if strcmp(tag,'DIAG56T')            % how far below tau* the default degree certifies 43 and 56
    pc_diag56t();
    return
end
if strcmp(tag,'DIAGPD')             % the dual form against the primal on self-adjoint PIEs
    pc_diagpd();
    return
end
if strcmp(tag,'DIAGRD')             % the 2-D reaction-diffusion side result of DIAGPD
    pc_diagrd();
    return
end
if strcmp(tag,'DIAG70T')            % storage-degree cuts on the 70mix certificate (187,116)
    pc_diag2m({'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{1} = [2 2 2]','light'; 'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{1} = [2 1 1]','light'; 'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0','light'; 'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{1} = [2 2 2]; st.sop.deg_Q{2} = [0 1 1]','light'}, ...
              'DIAG70T',{'70mix'});
    return
end
if strcmp(tag,'DIAG69T')            % storage-degree cuts on the 69mix2 certificate (267,060)
    pc_diag69m({'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0'; 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [1 0 0]'},'DIAG69T');
    return
end
if strcmp(tag,'DIAGX3')             % X3 in the energy (Riemann) states
    pc_diagx3();
    return
end
if startsWith(tag,'EXV')            % st_run verification of a plan (EXVT3: the TRADE3 cuts)
    st_run(tag);
    return
end
if strcmp(tag,'BLK_T3')             % the DDE and NDS/DDF blocks among the verified TRADE3 cuts
    blockcheck('DDE',{'DDE02','DDE08','DDE13','DDE16','DDE21','DDE22','DDE23'},fullfile(here,'pc_blk_T3_results.csv'));
    blockcheck('NDSDDF',{'NDS12'},fullfile(here,'pc_blk_T3_results.csv'));
    return
end
if strcmp(tag,'DIAG44G')            % the mu sweep of 44 across mu = 1/2 in 'Pd' (soundness)
    pc_diag44g();
    return
end
if strcmp(tag,'DIAG44H')            % 44 at tau = 1 past mu = 0.6: where the witness turns positive
    pc_diag44h();
    return
end
if strcmp(tag,'DIAG70S')            % 70mix with the edge damping reversed must fail (soundness)
    pc_diag70s();
    return
end
if strcmp(tag,'DIAG69S')            % 69mix2 with an unstable heat part must fail (soundness)
    pc_diag69s();
    return
end
if strcmp(tag,'DIAG44I')            % heavier settings for 44 past mu = 0.49
    pc_diag44i();
    return
end
if strcmp(tag,'DIAG70P')            % 70mix across the edge oscillator parameters
    pc_diag70p();
    return
end
if strcmp(tag,'SND70')              % the 70mix storage at epneg > 0 must fail (not uniformly stable)
    pc_diag2m({'Q','st.lpi.ineq.dD = 0; st.epneg = 1e-3','light'; 'Q','st.lpi.ineq.dD = 0; st.epneg = 1e-2','light'},'SND70',{'70mix'});
    return
end
if strcmp(tag,'SND69')              % the 69mix2 storage at epneg > 0 must fail (not exponentially stable)
    pc_diag69m({'st.lpi.ineq.dD = 0; st.epneg = 1e-3'; 'st.lpi.ineq.dD = 0; st.epneg = 1e-2'},'SND69');
    return
end
if strcmp(tag,'DIAGSED')            % the new certificates with SeDuMi in place of MOSEK
    pc_diagsed();
    return
end
if strcmp(tag,'SED69')              % 69mix2 with SeDuMi in place of MOSEK
    addpath(genpath('C:/YALMIP/yalmip_install/tbxmanager/toolboxes/sedumi/1.3/pcwin64/sedumi_1_3_pcwin64'));
    pc_diag69m({'st.lpi.ineq.dD = 0; st.sos_opts.solver = ''sedumi'''},'SED69');
    return
end
if startsWith(tag,'LIB_')           % the library path: example file -> its st -> auto executor
    libcheck(fullfile(here,['pc_lib_',tag(5:end),'.csv']),fullfile(here,['pc_lib_',tag(5:end),'_results.csv']));
    return
end
J = readtable(fullfile(here,'pc_jobs.csv'),'Delimiter',',','TextType','char');
J = J(strcmp(J.tag,tag),:);
if isempty(J),  error('pc_run:tag','no jobs with tag %s',tag);  end
if any(contains(J.mods,'sdpt3'))     % SDPT3 4.0, never on the saved path
    D = 'C:/Users/mpeet/solvers/SDPT3-4.0';
    addpath(D,fullfile(D,'Solver'),fullfile(D,'Solver','Mexfun'),fullfile(D,'HSDSolver'));
end
if any(contains(J.mods,'sedumi'))    % SeDuMi 1.3: on the saved path only, which lb_setup resets
    addpath(genpath('C:/YALMIP/yalmip_install/tbxmanager/toolboxes/sedumi/1.3/pcwin64/sedumi_1_3_pcwin64'));
end
jobs = cell(1,height(J));
for i = 1:height(J)
    v = J.param(i);   if isnan(v),  v = [];  end
    jobs{i} = {J.key{i},J.form{i},J.tier{i},J.mods{i},J.variant{i},v};
end
cap = J.ndvcap(~isnan(J.ndvcap));
ndvcap = 2e6;   if ~isempty(cap),  ndvcap = max(cap);  end
diag_amt(jobs,tag,ndvcap);
end

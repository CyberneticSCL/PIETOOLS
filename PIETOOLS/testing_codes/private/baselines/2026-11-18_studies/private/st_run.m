function st_run(tag)
% ST_RUN(TAG) runs one configuration of cleaning/private/settings_1116/
% PLAN_1118_studies.md through this study copy of the 11/17 library
% harness (lb_run, path 'sop'), into results_<TAG>.csv in this folder.
%
% INPUT
%   tag:    char, the configuration:
%     'L1ctl'   light, 1-D rows (PDE 1-32, DDE, NDSDDF): control of B2/D,
%               and check of this copy against results_ddeQ4.csv (11/17)
%     'B2dD0'   light, 1-D rows, lpi.ineq.dD = 0
%     'H1ctl'   heavy, 1-D rows: control of B3
%     'B3q522'  heavy, 1-D rows, sop.deg_Q{1} = [5 2 2]
%     'B3q622'  heavy, 1-D rows, sop.deg_Q{1} = [6 2 2]
%     'D1P'     light, DDE and NDSDDF rows, stability forms 'Q'->'P',
%               'Qd'->'Pd'
%     'D2st'    light, the five DDE rows with an example st (eec4f07d,
%               EXAMPLE_SETTINGS.md Section 7) in their form 'Q'
%     'E1str'   stripped, PDE 61-72;  'E2l', 'E2h': light, heavy, PDE 53-56
%     'C1ext2', 'C2ext2'  extreme, two-variable rows: as is; stripped 1-D entries
%     'H2ctl', 'H2T3'     heavy, two-variable rows: as is; with T3 (B1 follow-up)
%     'D1Ppde', 'D1Ppdeh' light, heavy: 1-D PDE stability rows in form 'P'
%     'D2Pl'    light, two-variable stability rows in form 'P'
%     'GallQh', 'GallPh'  heavy, every H-infinity gain row in form 'Q', 'P'
%     'M1L', 'M2L'        light, after the fast-forward to origin 7764389f
%               (opvar size, 7d0acf2c): the 1-D rows (control L1ctl) and the
%               two-variable rows of B1c0 under 40 s (control B1c0)
%     Phase B of the per-example settings (maintainer 11/19: "all examples
%     should have optimized settings files"; cleaning/private/
%     example_settings_1119/PHASE_A.md), MMP, 11/19/2026:
%     'XSP', 'XSPd0', 'XLPd0', 'XSQ', 'XSQd0'  1-D rows at stripped (S) or
%               light (L), with dD 0 ('d0'), in form 'P' (stability and gain
%               executives mapped 'Q'->'P', 'Qd'->'Pd') or as is ('Q')
%     'X2SP', 'X2SPd0', 'X2LPd0', 'X2SQ', 'X2SQd0'  the same on the
%               two-variable rows with a Phase A choice ('P': those light 'P'
%               certifies; 'X2SQ': those without a stripped run on this tree)
%     'X2HPd0'  PDE 45-47 (H-infinity gains), heavy, dD 0, form 'P'
%     'XVP8'    DDE08 at veryheavy in form 'P' (light, heavy 'P' fall short)
%     'XVH'     veryheavy on this tree: Phase A rows measured on the 11/17
%               tree only (DDE13, DDE18, PDE 31, 32, 49/H2_norm_dual)
%
% NOTES
%   - Size gate 2e6 variables, m 45,000 (gate/lpisolve.m); rows started
%     after 15 minutes are left for a rerun (lb_run resumes), so a MATLAB
%     session stays near 20 minutes.
%   - Hooks (lb_exec of this copy): LB_ST_HOOK, @(S,PIE,exec) -> S;
%     LB_FORM_MAP, struct of stability forms; LB_ROW_KEY set per row.
% MMP, 11/18/2026.
% MMP, 11/23/2026: 'VDA' takes the root from the checkout that holds this
%   file (was the live repo: a run from an export or a worktree put the live
%   tree on the path; stress test of 11/22/2026). 'CTLCHK', variant 'Ai':
%   Adi{1} moved to Ai{1} only where the block sets Adi (blocks/DDE19.m sets
%   Ai{1} since 99f7dff5, and rmfield of the absent Adi errored). 'DIAGA1',
%   'DIAGA3': no PDE 16 rows; the example left the library on 11/19/2026
%   and each such row could only record an error.
here = fileparts(mfilename('fullpath'));    cd(here);
if strcmp(tag,'A2size'),  st_d2size();  return,  end     % study A2: D2 sizes, no solve
if strcmp(tag,'VDA')      % D-a verified: PIETOOLS_DDE with Example B ('P') and Example E ('Q')
%   R = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';   cd(R);  pietools_path_update;  % MMP, 11/23/2026 (was)
    % R: the checkout holding this file, the first ancestor of here with    % MMP, 11/23/2026
    % pietools_path_update.m (the loop ends at the drive root if none)      % MMP, 11/23/2026
    R = here;                                                               % MMP, 11/23/2026
    while ~isfile([R '/pietools_path_update.m']) && ~strcmp(R,fileparts(R)) % MMP, 11/23/2026
        R = fileparts(R);                                                   % MMP, 11/23/2026
    end                                                                     % MMP, 11/23/2026
    cd(R);  pietools_path_update;                                           % MMP, 11/23/2026
    run_script(fullfile(here,'vda','vda_B.m'));   % each in its own workspace:
    run_script(fullfile(here,'vda','vda_E.m'));   % the scripts begin with clear
    return
end
% study A3: D2 solved at the cap-1 rung: 'A3cap1' at 0.5 tau*, 'A3cap1f09' at
% 0.9 tau*, 'A3cap1f11' at 1.1 tau* (control)
% 'A3unc': no cap with the triangular cap (2,134,428 variables, m 90,777),
% equal in certificates to the uncapped program; 'A3cap1ci': cap 1 with it
switch tag
    case 'A3cap1',     st_d2solve(1,0.5);  return
    case 'A3cap1f09',  st_d2solve(1,0.9);  return
    case 'A3cap1f11',  st_d2solve(1,1.1);  return
    case 'A3unc',      st_d2solve(NaN,0.5,'auto');  return
    case 'A3cap1ci',   st_d2solve(1,0.5,'auto');  return
    % the small rungs of the A2 sweep (dD 0 and/or deg_Q{3} 0), at 0.5 tau*
    case 'A3s1',       st_d2solve(NaN,0.5,'auto',0,0);  return   % 102,582 / m 10,569
    case 'A3s2',       st_d2solve(1,0.5,[],0,0);  return         %  83,574 / m  8,733
    case 'A3s3',       st_d2solve(1,0.5,'auto',0,1);  return     % 387,420 / m 29,348
    case 'A3s4',       st_d2solve(1,0.5,[],0,1);  return         % 655,708 / m 35,881
    case 'A3s5',       st_d2solve(1,0.5,[],1,0);  return         % 764,582 / m 42,665
    case 'A3s6',       st_d2solve(0,0.5,[],0,1);  return         % 155,004 / m 23,286
    % the storage form 'P' at light, uncapped with the triangular cap (9,108 /
    % m 12,826, sized by the Fable session, d2/d2_size3.log): a certificate
    % settles D2 at 0.5 tau* in the stronger form
    case 'A3P1',       st_d2solve(NaN,0.5,'auto',1,1,'P');  return
    case 'A3P1f09',    st_d2solve(NaN,0.9,'auto',1,1,'P');  return
    case 'A3P1f11',    st_d2solve(NaN,1.1,'auto',1,1,'P');  return
end
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
% Session length: minutes in the file ST_STOPMIN of this folder (15 while
% the maintainer is online, 60 while the maintainer is offline), default 15
stopmin = 15;
if exist(fullfile(here,'ST_STOPMIN'),'file'),  stopmin = str2double(strtrim(fileread(fullfile(here,'ST_STOPMIN'))));  end
if ~isfinite(stopmin) || stopmin<=0,  stopmin = 15;  end
opts = struct('path','sop','results',fullfile(here,['results_',tag,'.csv']),'max2d',600, ...
    'stopat',char(datetime('now')+minutes(stopmin),'yyyy-MM-dd HH:mm'));
lib1 = {'PDE','DDE','NDSDDF'};      % with tiers2d {}: the 1-D rows only
switch tag
    case 'L1ctl',   ids = lib1;  tiers = {'light'};  opts.tiers2d = {};
    case 'M1L',     ids = lib1;  tiers = {'light'};  opts.tiers2d = {};     % L1ctl after 7764389f
    case 'M2L'      % B1c0 rows under 40 s, after 7764389f
        ids = {'PDE:33','PDE:34','PDE:35','PDE:37','PDE:40','PDE:42','PDE:43','PDE:46', ...
               'PDE:57','PDE:59','PDE:60','PDE:62','PDE:63','PDE:66'};  tiers = {'light'};
    case 'B2dD0',   ids = lib1;  tiers = {'light'};  opts.tiers2d = {};
        setappdata(0,'LB_ST_HOOK',@hook_dD0);
    case 'H1ctl',   ids = lib1;  tiers = {'heavy'};  opts.tiers2d = {};
    case 'B3q522',  ids = lib1;  tiers = {'heavy'};  opts.tiers2d = {};
        setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_degQ1(S,[5 2 2]));
    case 'B3q622',  ids = lib1;  tiers = {'heavy'};  opts.tiers2d = {};
        setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_degQ1(S,[6 2 2]));
    case 'D1P',     ids = {'DDE','NDSDDF'};  tiers = {'light'};  opts.tiers2d = {};
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd'));
    case 'D1Ph',    ids = {'DDE','NDSDDF'};  tiers = {'heavy'};  opts.tiers2d = {};   % 'P' at heavy
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd'));
    case 'D2st',    ids = {'DDE:DDE03','DDE:DDE04','DDE:DDE05','DDE:DDE07','DDE:DDE08'};
        tiers = {'light'};  opts.tiers2d = {};
        setappdata(0,'LB_ST_HOOK',@hook_blockst);
    case {'G45P','G46P','G47P'}   % 2-D gains, exact 0.11101 / 0.17110 / 0.20668, in form 'P'
        ids = {['PDE:',tag(2:3)]};  tiers = {'light'};
        setappdata(0,'LB_FORM_MAP',struct('Q','P','only','PIETOOLS_Hinf_gain_sop'));
    case {'D1Ppde','D1Ppdeh'}   % D on the 1-D PDE stability rows (D100: 'P' wherever it certifies):
                                % light (heavy: 'h'), 'Q'->'P', 'Qd'->'Pd'; controls L1ctl, H1ctl
        S = load(fullfile(here,'catalog.mat'));
        sel = strcmp({S.C.lib},'PDE') & [S.C.dim]==1 & contains({S.C.exec},'stability');
        ids = {S.C(sel).key};  tiers = {'light'};  opts.tiers2d = {};
        if tag(end)=='h',  tiers = {'heavy'};  end
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd'));
    case 'D2Pl'     % D on the two-variable stability rows at light, 'Q'->'P', 'Qd'->'Pd'; control B1c0
        S = load(fullfile(here,'catalog.mat'));
        sel = [S.C.dim]==2 & contains({S.C.exec},'stability');
        ids = {S.C(sel).key};  tiers = {'light'};
        opts.skip = {'PDE:65@light'};    % as in B1: 2.85M variables at light in 'Q' (size gate)
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd'));
    case {'GallQ','GallP','GallQh','GallPh'}   % every H-infinity gain row at light (heavy: 'h'), form
                                               % 'Q' (control) or 'P' ('Pd' dual)
        S = load(fullfile(here,'catalog.mat'));
        ids = {S.C(contains({S.C.exec},'Hinf_gain')).key};  tiers = {'light'};
        if tag(end)=='h',  tiers = {'heavy'};  end
        if startsWith(tag,'GallP')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only','PIETOOLS_Hinf_gain_sop'));
        end
    case {'B1c0','B1c1','B1c2','B1c3'}   % B1 candidates on every catalog row in two variables
        S = load(fullfile(here,'catalog.mat'));     ids = {S.C([S.C.dim]==2).key};  tiers = {'light'};
        % PDE 65 (Turing pair) not run: under dD 0 with the degree raise of Q it
        % held 46 GB at 09:30 with 5.7 GB free on the machine (B1c2), stopped
        opts.skip = {'PDE:65@light'};
        switch tag(4)
            case '1',  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,0,[]));
            case '2',  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,0,0));
            case '3',  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,1,0));
        end
    case 'B1c1g'    % the gain rows of B1c1 again, on the committed gain executive: the B1c1
                    % session of 07:47 ran them on the Fable session's uncommitted candidate
                    % (tree edit 07:04-07:54, D121); their lines are deleted, then rerun
        opts.results = fullfile(here,'results_B1c1.csv');
        ids = {'PDE:45','PDE:46','PDE:47','PDE:71'};  tiers = {'light'};
        drop_rows(opts.results,{'PDE,45,','PDE,46,','PDE,47,','PDE,71,'});
        setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,0,[]));
    case {'C1ext','C2ext','C3ext'}       % C: extreme as is; with the stripped 1-D entries; and dD 0
        ids = lib1;  tiers = {'extreme'};  opts.tiers2d = {};
        if ~strcmp(tag,'C1ext'),  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_ext(S,PIE,strcmp(tag,'C3ext')));  end
    case {'C1ext2','C2ext2'}             % C, step 2: extreme as is, and with the stripped 1-D entries
                                         % (C2ext, the best 1-D candidate), on the two-variable rows
        S = load(fullfile(here,'catalog.mat'));     ids = {S.C([S.C.dim]==2).key};  tiers = {'extreme'};
        if strcmp(tag,'C2ext2'),  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_ext(S,PIE,false));  end
    case {'H2ctl','H2T3'}                % B1 follow-up: heavy on the two-variable rows, as is and with
                                         % the B1 winner T3 (dD 1, deg_Q{2} 0 raised by degree twice);
                                         % both files seeded with the light rows of B1c0 (seed_light2),
                                         % so the 2-D time projection skips the same rows in both
        S = load(fullfile(here,'catalog.mat'));     ids = {S.C([S.C.dim]==2).key};  tiers = {'heavy'};
        opts.skip = {'PDE:65@heavy'};    % 8.08M variables at heavy (size gate, 11/17); 46 GB under raises (B1c2)
        seed_light2(opts.results,fullfile(here,'results_B1c0.csv'));
        if strcmp(tag,'H2T3'),  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,1,0));  end
    case {'B4ctl','B4q866'}              % B4: veryheavy on the ten 1-D bisection rows; Q at [8 6 6]
        ids = {'DDE:DDE01','DDE:DDE02','DDE:DDE03','DDE:DDE04','DDE:DDE05','DDE:DDE06','DDE:DDE07', ...
               'DDE:DDE08','NDSDDF:NDS05','NDSDDF:NDS08'};  tiers = {'veryheavy'};  opts.tiers2d = {};
        if strcmp(tag,'B4q866'),  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_degQ1(S,[8 6 6]));  end
    case 'E3vh',    ids = arrayfun(@(k) sprintf('PDE:%d',k),56:72,'UniformOutput',false);  tiers = {'veryheavy'};
    case 'E1str',   ids = arrayfun(@(k) sprintf('PDE:%d',k),61:72,'UniformOutput',false);  tiers = {'stripped'};
    case 'E2l',     ids = arrayfun(@(k) sprintf('PDE:%d',k),53:56,'UniformOutput',false);  tiers = {'light'};
    case 'E2h',     ids = arrayfun(@(k) sprintf('PDE:%d',k),53:56,'UniformOutput',false);  tiers = {'heavy'};
    % BEGIN MMP, 11/19/2026: Phase B of the per-example settings (header)
    case {'XSP','XSPd0','XLPd0','XSQ','XSQd0', 'X2SP','X2SPd0','X2LPd0','X2SQ','X2SQd0'}
        if tag(2)=='2'      % two-variable rows with a Phase A choice
            if contains(tag,'P')    % rows that light 'P' certifies
                ids = {'PDE:33','PDE:34','PDE:35','PDE:37','PDE:42','PDE:59','PDE:60','PDE:61', ...
                       'PDE:62','PDE:63','PDE:64','PDE:67','PDE:68','PDE:71'};
            elseif strcmp(tag,'X2SQ')   % no stripped run on this tree (E1str has 61-72)
                ids = {'PDE:33','PDE:34','PDE:35','PDE:37','PDE:38','PDE:40','PDE:42','PDE:43', ...
                       'PDE:57','PDE:58','PDE:59','PDE:60'};
            else
                ids = {'PDE:33','PDE:34','PDE:35','PDE:37','PDE:38','PDE:40','PDE:42','PDE:43', ...
                       'PDE:57','PDE:58','PDE:59','PDE:60','PDE:61','PDE:62','PDE:63','PDE:64', ...
                       'PDE:66','PDE:67','PDE:68','PDE:71'};
            end
        else
            ids = lib1;  opts.tiers2d = {};      % the 1-D rows
        end
        tiers = {'stripped'};  if contains(tag,'LP'),  tiers = {'light'};  end
        if endsWith(tag,'d0'),  setappdata(0,'LB_ST_HOOK',@hook_dD0);  end
        if contains(tag,'P')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd', ...
                'only',{{'PIETOOLS_stability_sop','PIETOOLS_Hinf_gain_sop'}}));
        end
    case 'STOCK2D'      % the New Posopvar session (11/19): PDE 65 certified on 11/17 in 'P'
                        % with the stock 2-D storage (2,102,408 variables, 740 s); PDE 39
                        % untried in 'P' with it (about 0.9M). light + settings_2d_stock
        m = 'st.settings_2d = settings_2d_stock(''light'')';
        diag_amt({{'PDE:39','P','light',m,'A',[]},{'PDE:65','P','light',m,'A',[]}},tag,2.2e6);
        return
    case 'EXS8'         % PDE 8 (Das, flux sign reversed): the knob search from 'Q' stripped,
                        % the primal at its certified threshold 9.7924981, the dual
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB_das_starts.csv');
        ex_search(sf,tag,struct('keys',{{'PDE:8','PDE:8/stability_dual'}},'maxbuild',30, ...
                  'rowmin',15,'statedir','_das'));
        return
    case {'D8P','D8Q'}  % PDE 8 (Das) with the flux sign reversed (11/19/2026): tier passes,
                        % bisection on lam toward the documented pi^2 (the threshold is
                        % 9.7946), form 'P' (map) or as is ('Q')
        ids = {'PDE:8','PDE:8/stability_dual'};  tiers = {'stripped','light','heavy'};
        opts.tiers2d = {};
        if strcmp(tag,'D8P')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
        end
    case {'TQ2ctl','TQ2011'}  % tier pattern (PHASE_C.md section 6): the two-variable stability
                        % rows at light in 'Q', as is and with deg_Q{2} = [0 1 1], the most
                        % frequent two-variable change of the tuned settings besides dD
        ids = {'PDE:33','PDE:34','PDE:35','PDE:37','PDE:38','PDE:40','PDE:42','PDE:43', ...
               'PDE:57','PDE:58','PDE:59','PDE:60','PDE:61','PDE:62','PDE:63','PDE:64', ...
               'PDE:66','PDE:67','PDE:68'};
        tiers = {'light'};
        if strcmp(tag,'TQ2011'),  setappdata(0,'LB_ST_HOOK',@hook_degQ2011);  end
    case {'TQ2hctl','TQ2hcap','TG2ctl','TG2cap'}  % the triangular cap at heavy on the same
                        % stability rows ('Q'), and on the two-variable gain rows at light ('P')
        if startsWith(tag,'TQ')
            ids = {'PDE:33','PDE:34','PDE:35','PDE:37','PDE:38','PDE:40','PDE:42','PDE:43', ...
                   'PDE:57','PDE:58','PDE:59','PDE:60','PDE:61','PDE:62','PDE:63','PDE:64', ...
                   'PDE:66','PDE:67','PDE:68'};
            tiers = {'heavy'};
        else
            ids = {'PDE:45','PDE:46','PDE:47','PDE:71'};  tiers = {'light'};
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only','PIETOOLS_Hinf_gain_sop'));
        end
        if endsWith(tag,'cap')
            setappdata(0,'LB_ST_HOOK',@(S,~,~) setfield(S,'lpi',setfield(S.lpi,'ineq',setfield(S.lpi.ineq,'cap_inc','auto'))));
        end
    case {'TQ2cap','TQ2psn'}  % the same rows with cap_inc 'auto' (30% of the tuned rows), and
                        % with psatz{2} 'none' (25%)
        ids = {'PDE:33','PDE:34','PDE:35','PDE:37','PDE:38','PDE:40','PDE:42','PDE:43', ...
               'PDE:57','PDE:58','PDE:59','PDE:60','PDE:61','PDE:62','PDE:63','PDE:64', ...
               'PDE:66','PDE:67','PDE:68'};
        tiers = {'light'};
        if strcmp(tag,'TQ2cap'),  setappdata(0,'LB_ST_HOOK',@(S,~,~) setfield(S,'lpi',setfield(S.lpi,'ineq',setfield(S.lpi.ineq,'cap_inc','auto'))));
        else,  setappdata(0,'LB_ST_HOOK',@hook_psatz2none);  end
    case {'T1Pc0','T1Pi0'}  % the one-variable rows at stripped 'P' (control: XSP) with
                        % len_com{1} = 0 (49% of the tuned stability rows), len_inc{1} = 0 (32%)
        ids = lib1;  opts.tiers2d = {};  tiers = {'stripped'};
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd', ...
            'only',{{'PIETOOLS_stability_sop','PIETOOLS_Hinf_gain_sop'}}));
        if strcmp(tag,'T1Pc0'),  setappdata(0,'LB_ST_HOOK',@(S,~,~) hook_len1(S,'len_com'));
        else,  setappdata(0,'LB_ST_HOOK',@(S,~,~) hook_len1(S,'len_inc'));  end
    case 'X2HPd0'       % H-infinity gains 45-47: light 'P' 0.4-0.9% above exact, heavy 'P' exact
        ids = {'PDE:45','PDE:46','PDE:47'};  tiers = {'heavy'};
        setappdata(0,'LB_ST_HOOK',@hook_dD0);
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only','PIETOOLS_Hinf_gain_sop'));
    case 'XVP8'         % DDE08: 'P' reaches 3.5471558 at heavy, 'Q' [8 6 6] 3.9144165
        ids = {'DDE:DDE08'};  tiers = {'veryheavy'};  opts.tiers2d = {};
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd'));
    case 'XVH'          % Phase A choices measured on the 11/17 tree only
        ids = {'DDE:DDE13','DDE:DDE18','PDE:31','PDE:32','PDE:49/H2_norm_dual'};
        tiers = {'veryheavy'};  opts.tiers2d = {};
    case {'EXT','EXS1','EXS2'}  % the knob search of each row (ex_search): a test on
                                % four rows; the 1-D rows; the 2-D rows
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB_choices_starts.csv');
        S = load(fullfile(here,'catalog.mat'));
        T = readtable(sf,'Delimiter',',','TextType','char');
        dims = arrayfun(@(k) S.C(find(strcmp({S.C.key},T.key{k}),1)).dim, 1:height(T));
        switch tag
            case 'EXT',   keys = {'PDE:1','DDE:DDE03','DDE:DDE10','PDE:33'};
            case 'EXS1',  keys = T.key(dims==1)';
            case 'EXS2',  keys = T.key(dims==2)';
        end
        ex_search(sf,tag,struct('keys',{keys}));
        return
    case 'NEWCAT'       % catalog rows of the example files outside the library (lb_meta X1-X3)
        % lb_catalog_add({'X1','X2','X3'});   (11/19/2026: X2 moved to cleaning/private/dead_code)
        lb_catalog_add({'X1','X3'});
        return
    case {'NEWT','NEWTP','NEWTd0','NEWTPd0'}   % tier passes on those rows: as is, form 'P',
                                               % and both with dD 0
        S = load(fullfile(here,'catalog.mat'));
        ids = {S.C(startsWith({S.C.key},'PDE:X') & [S.C.run] & ~strcmp({S.C.key},'PDE:X2')).key};  % X2 moved 11/19/2026
        tiers = {'stripped','light','heavy'};
        if endsWith(tag,'d0'),  tiers = {'stripped','light'};  setappdata(0,'LB_ST_HOOK',@hook_dD0);  end
        if contains(tag,'P')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd', ...
                'only',{{'PIETOOLS_stability_sop','PIETOOLS_Hinf_gain_sop'}}));
        end
    case {'TRADE1','TRADE2'}    % trades between storage and slack (ex_search), from
                                % the final settings of the knob search; 1-D, 2-D rows
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB3_starts.csv');
        S = load(fullfile(here,'catalog.mat'));
        T = readtable(sf,'Delimiter',',','TextType','char');
        dims = arrayfun(@(k) S.C(find(strcmp({S.C.key},T.key{k}),1)).dim, 1:height(T));
        nd = 1 + strcmp(tag,'TRADE2');
        ex_search(sf,tag,struct('keys',{T.key(dims==nd)'},'statedir','_trade', ...
                  'maxbuild',120 - 80*(nd==2)));
        return
    case {'EXV1','EXV2','EXV1B','EXVX'}  % verification of the planned st and forms of each example
                            % (EXVX: the rows of the example files outside the library)
                            % (EXV1B: the 1-D plan again, after the trade pass and the
                            % verified overrides)
                            % (phaseC_plan.csv) on every row of it: bisections for
                            % thresholds, values for gains; 1-D rows, then 2-D rows
        pf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseC_plan.csv');
        P = readtable(pf,'Delimiter',',','TextType','char');
        P = P(strcmp(P.status,'ready') | strcmp(P.status,'partial'),:);
        S = load(fullfile(here,'catalog.mat'));
        ids = {};  plan = containers.Map();
        for q = 1:height(P)
            rk = strsplit(P.rows{q},' ');
            for r = rk
                plan(r{1}) = struct('tier',P.tier{q},'stm',P.statements{q}, ...
                    'sf',P.stability_form{q},'gf',P.Hinf_gain_form{q});
                ids{end+1} = r{1}; %#ok<AGROW>
            end
        end
        k = cellfun(@(x) find(strcmp({S.C.key},x),1), ids);
        dims = [S.C(k).dim];  ids0 = ids;
        if any(strcmp(tag,{'EXV1','EXV1B'})),  ids = ids(dims==1);  opts.tiers2d = {};  else,  ids = ids(dims==2);  end
        if strcmp(tag,'EXVX'),  ids = ids0(startsWith(ids0,'PDE:X'));  opts.tiers2d = {};  end  % the rows outside the library
        setappdata(0,'EXV_PLAN',plan);
        setappdata(0,'LB_ST_HOOK',@hook_plan);
        tiers = {'light'};      % a label: hook_plan builds the settings of each row
        opts.results = fullfile(here,['results_',tag,'.csv']);
    case 'EXU'          % the rows no measured setting certifies: alternative starts
                        % (cross terms, the other form, heavier tiers), m cap 30,000
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        keys = unique(T.key,'stable')';
        ex_search(sf,tag,struct('keys',{keys},'mcap',30000,'rowmin',25,'maxbuild',30));
        return
    case 'EXX'          % the knob search on the rows of the example files outside the
                        % library that a tier certifies (X1; X2 and X3 certify at none)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseBX_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'}));
        return
    case 'EXU3D'        % the 3-D rows whose light start the m cap of EXU refused: that
                        % start alone, m cap 45,000, run exclusive (the Schur complement
                        % at m = 32,846 is 8.6 GB)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_3d_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',45000,'rowmin',40, ...
                  'maxbuild',8,'statedir','_3d'));
        return
    case {'EXS2X','EXS2Y'}  % the two-variable rows whose search ended on the build budget
                        % (40), resumed from the saved state, 30 more builds: stability
                        % rows (X); gain rows (Y; PDE 48 left out, no cut in 40 builds)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB_choices_starts.csv');
        keys = {'PDE:64','PDE:63','PDE:61','PDE:62','PDE:58','PDE:57'};
        if strcmp(tag,'EXS2Y'),  keys = {'PDE:46','PDE:45','PDE:71','PDE:47'};  end
        for q = 1:numel(keys)
            f = fullfile(here,'exs_state',[strrep(strrep(keys{q},':','_'),'/','_'),'.mat']);
            L = load(f);  st = L.st;
            if st.done && strcmp(st.reason,'build budget'),  st.done = false;  st.reason = '';  save(f,'st');  end
        end
        ex_search(sf,tag,struct('keys',{keys},'maxbuild',70,'rowmin',60));
        return
    case 'EXU2DU'       % the two-variable starts of EXU that the caps refused (PDE 39: m
                        % 30,721; PDE 65: 2.04M variables), uncapped to m 45,000, 2.5M
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_2du_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',45000,'ndvcap',2.5e6, ...
                  'rowmin',40,'maxbuild',8,'statedir','_2du'));
        return
    case 'EXSED'        % the open rows 15, 16, 36, 41, 44 with SeDuMi in place of MOSEK
                        % (41 and 36 end inconclusive in MOSEK at every start)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_sed_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',30000,'rowmin',30, ...
                  'maxbuild',4,'statedir','_sed'));
        return
    case 'EXS54'        % PDE 54 in 'P' stripped (the setting that certifies PDE 53 at
                        % 580,612; PDE 54 light 'Q' ended at feasibility ratio 1.69)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_54_starts.csv');
        ex_search(sf,tag,struct('keys',{{'PDE:54'}},'mcap',45000,'ndvcap',3e6, ...
                  'rowmin',55,'maxbuild',12,'statedir','_54'));
        return
    case 'EXS53'        % PDE 53 from 'P' stripped (580,612 variables, certified; the written
                        % light 'Q' setting: 668,443): the knob search, 12 builds
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_53_starts.csv');
        ex_search(sf,tag,struct('keys',{{'PDE:53'}},'mcap',45000,'ndvcap',3e6, ...
                  'rowmin',55,'maxbuild',12,'statedir','_53'));
        return
    case 'EXP3C'        % PDE 55 in the direct form with the triangular cap (D151: the cap
                        % cuts the 2-D direct form 1.28-1.54x and keeps its certificates)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_3dc_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',45000,'ndvcap',3e6, ...
                  'rowmin',50,'maxbuild',1,'statedir','_3dc'));
        return
    case 'EXP3D'        % probes before a Sol job for PDE 55: its 'Q' stripped with dD 0 and
                        % deg_Q{3} 0; 'P' stripped with dD 1 on the one-field PDE 53
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_3dp_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',45000,'ndvcap',3e6, ...
                  'rowmin',50,'maxbuild',1,'statedir','_3dp'));
        return
    case 'EXU3D3'       % PDE 55, 'P' stripped: the lift of dD 1 under the triangular cap
                        % (dD 1 alone: 7.34M variables; dD 0: infeasible at 114,704)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_3d3_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',45000,'ndvcap',2.5e6, ...
                  'rowmin',50,'maxbuild',8,'statedir','_3d3'));
        return
    case 'EXU3D2'       % PDE 55 in the form 'P' with the slack lift degree of the tiers
                        % (its light 'Q' start: 13 million variables)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB2_3d2_starts.csv');
        T = readtable(sf,'Delimiter',',','TextType','char');
        ex_search(sf,tag,struct('keys',{unique(T.key,'stable')'},'mcap',45000,'rowmin',40, ...
                  'maxbuild',8,'statedir','_3d2'));
        return
    case 'BLOCKCHECK'   % the DDE and NDS/DDF blocks with planned settings, each run through
                        % PIETOOLS_DDE.m or PIETOOLS_DDF.m uncommented (blockcheck.m)
        pf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseC_plan.csv');
        P = readtable(pf,'Delimiter',',','TextType','char');
        P = P(strcmp(P.status,'ready') & (startsWith(P.group,'DDE:') | startsWith(P.group,'NDSDDF:')),:);
        of = fullfile(here,'blockcheck_results.csv');  if exist(of,'file'),  delete(of);  end
        for lib = {'DDE','NDSDDF'}
            g = P.group(startsWith(P.group,[lib{1},':']));
            blockcheck(lib{1},regexprep(g,'^\w+:','')',of);
        end
        return
    case 'BLOCKCHECK3'  % DDE19 and DDF03 with their regulated outputs and st (11/19/2026)
        blockcheck('DDE',{'DDE19'},fullfile(here,'blockcheck_results.csv'));
        blockcheck('NDSDDF',{'DDF03'},fullfile(here,'blockcheck_results.csv'));
        return
    case 'EXSFIX'       % the search on the rows corrected 11/19/2026 (maintainer: "Apply all
                        % four"): PDE 15 with its bending term (dual test), PDE 4.5 at
                        % sig2 = 1.048, DDF03 with tau(i) = i (best: heavy, 0.42133672)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB_fix_starts.csv');
        ex_search(sf,tag,struct('keys',{{'PDE:15','PDE:4.5','NDSDDF:DDF03'}},'maxbuild',30, ...
                  'rowmin',15,'statedir','_fix'));
        return
    case 'DIAGB'        % PDE 15 with its bending term: the damping c of the source (c > 0
                        % certified in Peet et al. 2019), c = 0.1 and 1; c = 0 in the
                        % primal forms; light
        J = {};
        for fm = {'Q','Qd','P','Pd'}
            for cv = [0.1 1]
                J{end+1} = {'PDE:15',fm{1},'light','','A',cv,'c=%g;'}; %#ok<AGROW>
            end
        end
        J = [J, {{'PDE:15','Q','light','','A',[]},{'PDE:15','P','light','','A',[]}}];
        diag_amt(J,tag);
        return
    case 'EXSFIX2'      % DDE19 with the discrete delay: its st of the distributed delay
                        % (265 variables) ends unknown (ratio 0.594); search again
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB_fix2_starts.csv');
        ex_search(sf,tag,struct('keys',{{'DDE:DDE19'}},'maxbuild',30,'rowmin',15, ...
                  'statedir','_fix2'));
        return
    case 'LIBCHECK65'   % the library run path of case 65 with its recommended st (stock 2-D
                        % storage, form 'P'; 1,635,656 variables in the harness)
        cf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','libcheck_fix3_cases.csv');
        libcheck(cf,fullfile(here,'libcheck_results.csv'));
        return
    case 'LIBCHECK8'    % the library run path of case 8 (Das, flux sign reversed, its st)
        cf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','libcheck_fix2_cases.csv');
        libcheck(cf,fullfile(here,'libcheck_results.csv'));
        return
    case 'CATFIX'       % catalog rows of examples removed 11/19/2026 (PDE 16: orphaned_file;
                        % X2: cleaning/private/dead_code) marked not to run
        cf = fullfile(here,'catalog.mat');  S = load(cf);
        k = ismember({S.C.key},{'PDE:16','PDE:X2'});
        fprintf('## CATFIX: %d rows (%s)\n',nnz(k),strjoin({S.C(k).key},', '));
        [S.C(k).run] = deal(false);  save(cf,'-struct','S');
        return
    case 'LIBCHECKFIX'  % the library run path of cases 15 and 4.5 after their corrections
        cf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','libcheck_fix_cases.csv');
        libcheck(cf,fullfile(here,'libcheck_results.csv'));
        return
    case 'BLOCKCHECK5'  % DDE19 with the discrete delay (Ai), its st unchanged
        blockcheck('DDE',{'DDE19'},fullfile(here,'blockcheck_results.csv'));
        return
    case 'BLOCKCHECK4'  % DDF03 again, after blockcheck matches active blocks as they are
        blockcheck('NDSDDF',{'DDF03'},fullfile(here,'blockcheck_results.csv'));
        return
    case 'BLOCKCHECK2'  % the H2 blocks DDE25, DDE26 (BLOCKCHECK lacked the H2 flags)
        blockcheck('DDE',{'DDE25','DDE26'},fullfile(here,'blockcheck_results.csv'));
        return
    case 'LIBCHECK'     % the library run path of every example file with an output st
                        % (libcheck.m): the example st and forms reach the executive
        cf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','libcheck_cases.csv');
        libcheck(cf,fullfile(here,'libcheck_results.csv'));
        return
    case {'DIAGA1','DIAGA2','DIAGA3','DIAGA4'}  % rows no setting certifies, with A
                        % unchanged, A = 0, A = -T and A - T (diag_amt.m): storage
                        % side or lpi_ineq side
        J = {};  % rows under 50,000 variables at light, forms 'Q' and 'P'; the
                 % bisection rows (44, X3) at their tested parameter
%       R1 = {{'PDE:15',[]},{'PDE:16',[]},{'PDE:41',[]},{'PDE:44',0.374625},{'PDE:X3',0.375}}; % MMP, 11/23/2026 (was)
        % no PDE 16: the example left the library on 11/19/2026             % MMP, 11/23/2026
        R1 = {{'PDE:15',[]},{'PDE:41',[]},{'PDE:44',0.374625},{'PDE:X3',0.375}}; % MMP, 11/23/2026
        R2 = {{'PDE:69','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0; st.lpi.ineq.wcross = 1'}, ...
              {'PDE:70','st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1'}, ...
              {'PDE:36','st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1'}, ...
              {'PDE:39','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [2 1 1]; st.lpi.ineq.wcross = 1'}, ...
              {'PDE:65','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0; st.lpi.ineq.wcross = 1'}};
        vs = {'A','A0','AmT'};  if any(strcmp(tag,{'DIAGA3','DIAGA4'})),  vs = {'AsT'};  end
        if any(strcmp(tag,{'DIAGA1','DIAGA3'}))
            for k = R1
                for fm = {'Q','P'},  for vr = vs
                    J{end+1} = {k{1}{1},fm{1},'light','',vr{1},k{1}{2}}; %#ok<AGROW>
                end,  end
            end
            % DIAGA3 PDE 16 rows removed with the example (11/19/2026)      % MMP, 11/23/2026
%           if strcmp(tag,'DIAGA3')     % PDE 16 'P' fails at A = -T: a larger slack lift % MMP, 11/23/2026 (was)
%               J{end+1} = {'PDE:16','P','light','st.lpi.ineq.dD = 2','AmT',[]}; % MMP, 11/23/2026 (was)
%               J{end+1} = {'PDE:16','P','light','st.lpi.ineq.Dmax = true','AmT',[]}; % MMP, 11/23/2026 (was)
%           end                                                             % MMP, 11/23/2026 (was)
        else      % two-variable rows in 'Q' at their first EXU start, whose A is known
            if strcmp(tag,'DIAGA2'),  vs = {'A0','AmT'};  end
            % DIAGA4: not 69, 65, whose start (deg_Q{2} 0) fails at A = 0 (DIAGA2)
            if strcmp(tag,'DIAGA4'),  R2 = R2(2:4);  end
            for k = R2
                for vr = vs
                    J{end+1} = {k{1}{1},'Q','light',k{1}{2},vr{1},[]}; %#ok<AGROW>
                end
            end
        end
        diag_amt(J,tag);
        return
    case 'DIAGA5'       % PDE 69 at the default deg_Q (its cap start, A infeasible at
                        % 1,005,460 variables): A = 0 and A = -T
        m = 'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''; st.lpi.ineq.wcross = 1';
        diag_amt({{'PDE:69','Q','light',m,'A0',[]},{'PDE:69','Q','light',m,'AmT',[]}},tag);
        return
    case 'DIAGA6'       % PDE 44 (abscissa -0.29 at mu 0.4): the shift sweep A - c*T in
                        % 'Q' (A - T certifies, A does not, DIAGA3), light and heavy;
                        % X3 at shift 2 (fails at shift 1)
        J = {};
        for c = {'AsT0.5','AsT0.25','AsT0.1','AsT0.05'}
            J{end+1} = {'PDE:44','Q','light','',c{1},0.374625}; %#ok<AGROW>
        end
        for c = {'A','AsT0.25','AsT0.1'}
            J{end+1} = {'PDE:44','Q','heavy','',c{1},0.374625}; %#ok<AGROW>
        end
        J{end+1} = {'PDE:X3','Q','light','','AsT2',0.375};
        diag_amt(J,tag);
        return
    case 'DIAGA7'       % PDE 70: shifts 2, 4 (fails at 1, no cancellation); PDE 69 at its
                        % cap start: A - T; PDE 41 (unknown at A and A - T): lam 0, tau
                        % 0.5 and 2, the reaction and the delay length
        m70 = 'st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1';
        m69 = 'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''; st.lpi.ineq.wcross = 1';
        J = {{'PDE:70','Q','light',m70,'AsT2',[]},{'PDE:70','Q','light',m70,'AsT4',[]}, ...
             {'PDE:41','Q','light','','A',0,'lam=%g;'},{'PDE:41','Q','light','','A',0.5,'tau=%g;'}, ...
             {'PDE:41','Q','light','','A',2,'tau=%g;'},{'PDE:69','Q','light',m69,'AsT',[]}};
        diag_amt(J,tag);
        return
    case 'CTLCHK'       % DDE19 and DDF03 with their regulated outputs (11/19/2026):
                        % DDE19 as written (distributed delay Adi) and with the
                        % discrete delay of [5] (Ai); DDF03 as written (tau(i) = n)
                        % and with the delays of the source (tau(i) = i)
%       dis = 'DDE.Ai{1}=DDE.Adi{1}; DDE=rmfield(DDE,''Adi'');';            % MMP, 11/23/2026 (was)
        % 'Ai': the discrete delay of [5] in Ai{1}, from Adi{1} where the   % MMP, 11/23/2026
        % block sets Adi; blocks/DDE19.m sets Ai{1} since 99f7dff5 (11/19), % MMP, 11/23/2026
        % so there the variant is the block as written ('tau_i' likewise)   % MMP, 11/23/2026
        dis = ['if isfield(DDE,''Adi''), DDE.Ai{1}=DDE.Adi{1}; ', ...
               'DDE=rmfield(DDE,''Adi''); end'];                            % MMP, 11/23/2026
        J = {};
        for tr = {'stripped','light','heavy'}
            J = [J, {{'DDE:DDE19','Adi',tr{1},''},{'DDE:DDE19','Ai',tr{1},dis}, ...
                     {'NDSDDF:DDF03','tau_n',tr{1},''},{'NDSDDF:DDF03','tau_i',tr{1},'DDF.tau=1:n;'}}]; %#ok<AGROW>
        end
        ctl_check(J,tag);
        return
    case 'EXSCTL'       % the knob search on DDE19 and DDF03 once they have a regulated
                        % output (value rows: within 1e-3 of the best bound, heavy for
                        % DDF03, the flat stripped bound for DDE19)
        sf = fullfile(here,'..','..','..','..','..','cleaning','private', ...
                      'example_settings_1119','phaseB_ctl_starts.csv');
        ex_search(sf,tag,struct('keys',{{'DDE:DDE19','NDSDDF:DDF03'}},'maxbuild',30, ...
                  'rowmin',15,'statedir','_ctl'));
        return
    case 'DIAGA9'       % the delay rows 41, 44, X3 (a delay adds a second variable: mixed
                        % domains), never run with cross terms: wcross 1, and len_cross{2}
                        % 2 with it, at light and heavy, forms 'Q' and 'P', A unchanged
        R = {{'PDE:41',[]},{'PDE:44',0.374625},{'PDE:X3',0.375}};
        M = {'st.lpi.ineq.wcross = 1','st.lpi.pos.len_cross{2} = 2; st.lpi.ineq.wcross = 1'};
        J = {};
        for k = R,  for tr = {'light','heavy'},  for fm = {'Q','P'},  for mm = M
            J{end+1} = {k{1}{1},fm{1},tr{1},mm{1},'A',k{1}{2}}; %#ok<AGROW>
        end,  end,  end,  end
        diag_amt(J,tag);
        return
    case 'DIAGA8'       % PDE 56 (Sol job held since 11/18): its uncapped light 'Q' start with
                        % dD 0 (infeasible, 840,156 variables) at A = 0 and A = -T
        J = {{'PDE:56','Q','light','st.lpi.ineq.dD = 0','A0',[]}, ...
             {'PDE:56','Q','light','st.lpi.ineq.dD = 0','AmT',[]}};
        diag_amt(J,tag);
        return
    % END MMP, 11/19/2026
    otherwise
        % B1 tightness: 'T<k>_<system>', k = 0 light, 1 dD 0, 2 dD 0 and
        % deg_Q{2} 0 raised by degree (2 raises), 3 deg_Q{2} 0 raised by
        % degree; the system bisected on its parameter toward its known
        % margin (bracket [0.5, 1.1] x margin, reltol 1e-3 of the margin).
        tok = regexp(tag,'^T([0-3])_(\w+)$','tokens','once');
        if isempty(tok),  error('st_run:tag','Unknown tag %s.',tag);  end
        opts.catalog = tight_catalog(here);
        ids = {['PDE:',tight_key(tok{2}),'/reach']};  tiers = {'light'};
        switch tok{1}
            case '1',  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,0,[]));
            case '2',  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,0,0));
            case '3',  setappdata(0,'LB_ST_HOOK',@(S,PIE,ex) hook_b1(S,PIE,1,0));
        end
end
fprintf('## ST_RUN %s start %s: %s, %s\n',tag,char(datetime('now')),strjoin(tiers,' '),strjoin(ids,' '));
lb_run(ids,tiers,opts);
fprintf('## ST_RUN %s DONE %s\n',tag,char(datetime('now')));
end

function drop_rows(f,prefixes)
% deletes from the results file F the lines that start with any of PREFIXES
if ~exist(f,'file'),  return,  end
L = strsplit(fileread(f),newline);
keep = true(size(L));
for k = 2:numel(L)
    keep(k) = ~any(startsWith(L{k},prefixes));
end
fid = fopen(f,'w');  fprintf(fid,'%s',strjoin(L(keep),newline));  fclose(fid);
fprintf('st_run: dropped %d lines from %s\n',sum(~keep),f);
end

function seed_light2(f,src)
% copies the header and the light rows of SRC into the results file F when F
% does not exist yet (a resumed session keeps F as it is); lb_run skips rows
% already in F, and its 2-D projection reads the light rows as the lighter tier
if exist(f,'file'),  return,  end
L = splitlines(fileread(src));
keep = [L(1); L(contains(L,',light,'))];
fid = fopen(f,'w');  fprintf(fid,'%s\n',keep{:});  fclose(fid);
fprintf('st_run: seeded %s with %d light rows of %s\n',f,numel(keep)-1,src);
end

function run_script(f)
% runs the script F in this function's workspace
try
    run(f);
catch ME
    fprintf('@@ VDA %s: ERROR %s\n',f,strrep(ME.message,newline,' '));
end
end

function S = hook_dD0(S,~,~)
% B2: the lift raise dD = 0 (the 1-D rows only are run)
S.lpi.ineq.dD = 0;
end

function S = hook_degQ2011(S,~,~)
% TQ2011: Q of the direct form at degree [0 1 1] in two variables
S.sop.deg_Q{2} = [0 1 1];
end

function S = hook_psatz2none(S,~,~)
% TQ2psn: no psatz terms in two variables
S.lpi.ineq.psatz{2} = 'none';
end

function S = hook_len1(S,f)
% T1Pc0, T1Pi0: one storage length of the one-variable entry set to 0
S.lpi.pos.(f){1} = 0;
end

function S = hook_degQ1(S,v)
% B3: the 1-D degrees of Q in the direct stability form
S.sop.deg_Q{1} = v;
end

function S = hook_b1(S,PIE,dD,degQ2)
% B1: on a PIE in two variables only, the lift raise dD and, when degQ2 is
% given, Q of the direct form at that degree, raised by degree twice
if PIE.dim~=2,  return,  end
S.lpi.ineq.dD = dD;
if ~isempty(degQ2)
    S.sop.deg_Q{2} = degQ2;     S.sop.raise_Q = 'deg';  S.sop.max_raise_Q(2) = 2;
end
end

function S = hook_ext(S,~,dD0)
% C: the 1-D entries of stripped in extreme (storage row, Q, free
% operator), and with dD0 the lift raise 0
T = settings_sopvar_stripped();
for f = fieldnames(T.lpi.pos)'
    v = T.lpi.pos.(f{1});
    if iscell(v),  S.lpi.pos.(f{1}){1} = v{1};  end
end
S.sop.deg_Q{1} = T.sop.deg_Q{1};    S.lpi.var.deg{1} = T.lpi.var.deg{1};
if dD0,  S.lpi.ineq.dD = 0;  end
end

function k = tight_key(name)
% B1 tightness systems: library case of each
T = struct('kiss','34','para','37','s18','43','fo2s','57','tl','64','age','66','s17','72');
k = T.(name);
end

function f = tight_catalog(here)
% The catalog with one bisection row per B1 tightness system, key
% '<case>/reach', on the parameter whose margin the library documents
% (examples_PDE_library_PIETOOLS headers); written once.
f = fullfile(here,'catalog_tight.mat');
if exist(f,'file'),  return,  end
S = load(fullfile(here,'catalog.mat'));     C = S.C;
% case, parameter, setter, margin
R = {'34','lam','lam=%.10g;',2*pi^2
     '37','a','a=%.10g;',0.5*pi^2
     '43','tau','tau=%.10g;',1.034726
     '57','tau','tau=%.10g;',2.881401
     '64','lam','lam=%.10g;',13.686938
     '66','b','b=%.10g;',2.302685
     '72','r2','r2=%.10g;',pi^2};
for r = 1:size(R,1)
    k = find(strcmp({C.key},['PDE:',R{r,1}]),1);
    assert(~isempty(k),'tight_catalog: no row PDE:%s',R{r,1});
    c = C(k);
    c.key = ['PDE:',R{r,1},'/reach'];
    c.param = R{r,2};   c.setter = R{r,3};  c.docnum = R{r,4};
    c.lo = 0.5*R{r,4};  c.hi = 1.1*R{r,4};  c.bisect = true;
    if isfield(c,'run'),  c.run = true;  end
    if isfield(c,'doc'),  c.doc = sprintf('margin %.10g',R{r,4});  end
    C(end+1) = c; %#ok<AGROW>
end
save(f,'C');
end

function S = hook_plan(~,~,exec)
% EXV: the planned settings of the row (phaseC_plan.csv), from
% lpisettings_sop(tier) and the statements, and the form of this executive
% call: stability_form for the stability tests, Hinf_gain_form for the gain
key = getappdata(0,'LB_ROW_KEY');  plan = getappdata(0,'EXV_PLAN');  p = plan(key);
st = lpisettings_sop(p.tier,0,'','mosek'); %#ok<NASGU>
if ~isempty(p.stm),  eval(p.stm);  end
S = st;
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
if contains(exec,'stability') && ~isempty(p.sf)
    setappdata(0,'LB_FORM_MAP',struct('Q',p.sf,'Qd',[p.sf 'd']));
elseif contains(exec,'Hinf_gain') && ~isempty(p.gf)
    setappdata(0,'LB_FORM_MAP',struct('Q',p.gf,'Qd',[p.gf 'd'],'only','PIETOOLS_Hinf_gain_sop'));
end
end
function S = hook_blockst(S,~,~)
% D2: the st of the DDE library blocks committed in eec4f07d (form 'Q'),
% from EXAMPLE_SETTINGS.md Section 7, on top of light
switch getappdata(0,'LB_ROW_KEY')
    case 'DDE:DDE03',  S.sop.deg_Q{1} = [5 2 2];  S.sop.max_raise_Q(1) = 0;
    case 'DDE:DDE04',  S.sop.deg_Q{1} = [3 2 2];  S.sop.max_raise_Q(1) = 0;
    case 'DDE:DDE05',  S.sop.deg_Q{1} = [4 2 2];  S.sop.max_raise_Q(1) = 0;
    case 'DDE:DDE07',  S.sop.deg_Q{1} = [4 2 2];  S.sop.max_raise_Q(1) = 0;
    case 'DDE:DDE08',  S.sop.deg_Q{1} = [8 6 6];
end
end

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
%               (copvar size, 7d0acf2c): the 1-D rows (control L1ctl) and the
%               two-variable rows of B1c0 under 40 s (control B1c0)
%
% NOTES
%   - Size gate 2e6 variables, m 45,000 (gate/lpisolve.m); rows started
%     after 15 minutes are left for a rerun (lb_run resumes), so a MATLAB
%     session stays near 20 minutes.
%   - Hooks (lb_exec of this copy): LB_ST_HOOK, @(S,PIE,exec) -> S;
%     LB_FORM_MAP, struct of stability forms; LB_ROW_KEY set per row.
% MMP, 11/18/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
if strcmp(tag,'A2size'),  st_d2size();  return,  end     % study A2: D2 sizes, no solve
if strcmp(tag,'VDA')      % D-a verified: PIETOOLS_DDE with Example B ('P') and Example E ('Q')
    R = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';   cd(R);  pietools_path_update;
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

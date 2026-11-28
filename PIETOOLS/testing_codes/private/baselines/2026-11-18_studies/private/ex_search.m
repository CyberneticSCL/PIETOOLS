function ex_search(startsfile, tag, opts)
% EX_SEARCH(STARTSFILE, TAG, OPTS) searches, for each library row of
% STARTSFILE, example-specific changes of the settings that reduce the
% program while the row keeps its accuracy (Phase B of the per-example
% settings, maintainer 11/19/2026; cleaning/private/example_settings_1119).
%
% From the start of a row (tier, changes, form) a greedy descent tries one
% change of one knob at a time, among the knobs the executive of the row
% reads; a change is kept when the program has fewer decision variables
% and the row keeps its accuracy:
%   stability rows - certified under both gates at the library parameter;
%   bisection rows - certified at the best measured threshold less 1e-3 of
%                    it (one solve, not a bisection);
%   value rows     - certified, with a value within 1e-3 of the best.
% Both gates: harness status 'certified' and the executive status optimal
% or inaccurate (lb_exec). The size gate of gate/lpisolve.m is set below
% the current size, so a change that does not shrink the program is
% refused after its build, before the solver runs.
%
% INPUT
%   startsfile: CSV, columns key, kind, form, tier, mods, best (written by
%               cleaning/private/example_settings_1119/phaseA.pl); mods are
%               MATLAB statements on st, or the keywords stripped1d (the
%               1-D entries of stripped, study C2) and blockst (the DDE
%               block st of eec4f07d, study D2)
%   tag:        name of the run: log exs_<TAG>.csv, final exs_<TAG>_final.csv
%   opts:       struct, optional fields
%       .keys      cellstr, rows to search (default: every row of the file)
%       .stopmin   minutes after which no build starts (default ST_STOPMIN)
%       .maxbuild  builds per row, solved or refused (default 40)
%       .rowmin    minutes per row (default 20)
%       .ndvcap    decision-variable cap of a start (default 2e6)
%
% OUTPUT (files in this folder)
%   exs_<TAG>.csv        one line per build: time, key, step, move, outcome,
%                        status, exst, eta, fr, passes, ndv, m, gam, seconds,
%                        note
%   exs_state/<key>.mat  the state of a row, so a stopped session resumes it
%   exs_<TAG>_final.csv  one line per finished row: key, form, tier, the st
%                        statements (from lpisettings(tier)), start and
%                        final decision variables and m, builds, reason
%
% NOTES
%   - Moves by form: 'Q' stability (deg_Q, slack); 'P' stability, gains and
%     synthesis (storage row lpi.pos, slack, dw_feas or dw_P; lpi.var.deg
%     for synthesis); in two or more variables, cross terms (wcross 1,
%     len_cross 2) combined with each reduction refused or failed alone.
%   - The search is greedy: a kept change restarts the list of moves.
% MMP, 11/19/2026.
% MMP, 11/20/2026: apply_mods splits the statements of a start at a ';'
%                   outside brackets only: the split at every ';' cut a
%                   column vector (len_sos{1} = [2;1;2]), and the trade
%                   pass TRADE1 stopped at the start of 80 of 115 rows.
% MMP, 11/21/2026: opts.post true (the rerun of cleaning/private/
%                   settings_rerun_1121/PLAN.md): the gate of 11/20 (harness
%                   'certified', executive optimal, eta <= 1e-7), and a build
%                   refused for numerical reasons only (executive optimal,
%                   inaccurate or unknown; no infeasibility certificate) is
%                   checked a posteriori (post_verify: storage fixed, slacks
%                   sought) and kept when that check passes. The row ends on
%                   the smallest direct certificate when it costs at most
%                   opts.ratio (3) times the variables of the smallest one,
%                   and is not near the memory limit (opts.near, [1.2e6 27000]
%                   variables and rows, 60% of the size gate); else on the
%                   smallest (maintainer, 11/21: prefer certificates which
%                   directly certify if the cost is less than 3 times, not
%                   counting cases near the 32GB limit). The final file gains
%                   the columns path, ndv_direct, m_direct. Start keyword
%                   examplest: the st an example file returns. 'FQ' takes the
%                   moves of the 'Q' family.
% MMP, 11/21/2026: a row with no certified start keeps the starts it tried
%                   (st.altdone) and a resumed session skips them: each
%                   session retried every start from the first, so a row
%                   whose first start takes most of a session to build
%                   (PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq, light 'P': 28M
%                   variables, 13.6 min, refused by the size gate) never
%                   reached its later starts.
if nargin<3,  opts = struct();  end
post = isfield(opts,'post') && opts.post;                                   % MMP, 11/21/2026
setappdata(0,'EXS_POST',post);                                              % MMP, 11/21/2026
ratio = 3;  if isfield(opts,'ratio'),  ratio = opts.ratio;  end             % MMP, 11/21/2026
nearv = [1.2e6 27000];  if isfield(opts,'near'),  nearv = opts.near;  end   % MMP, 11/21/2026
setappdata(0,'EXS_RULE',struct('ratio',ratio,'near',nearv));                % MMP, 11/21/2026
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'gate'));     % the size gate first on the path, as in lb_run
stopmin = 15;
if exist(fullfile(here,'ST_STOPMIN'),'file')
    stopmin = str2double(strtrim(fileread(fullfile(here,'ST_STOPMIN'))));
end
if isfield(opts,'stopmin'),  stopmin = opts.stopmin;  end
maxbuild = 40;  if isfield(opts,'maxbuild'),  maxbuild = opts.maxbuild;  end
mcap = 45000;   if isfield(opts,'mcap'),      mcap = opts.mcap;          end
setappdata(0,'EXS_MCAP',mcap);      % m cap of every build (session length)
rowmin = 20;    if isfield(opts,'rowmin'),    rowmin = opts.rowmin;      end
ndvcap = 2e6;   if isfield(opts,'ndvcap'),    ndvcap = opts.ndvcap;      end
setappdata(0,'EXS_NDVCAP',ndvcap);
tstop = datetime('now') + minutes(stopmin);
logf = fullfile(here,['exs_',tag,'.csv']);
finf = fullfile(here,['exs_',tag,'_final.csv']);
sfx = '';       if isfield(opts,'statedir'),  sfx = opts.statedir;  end   % a pass of its own
sdir = fullfile(here,['exs_state',sfx]);  if ~exist(sdir,'dir'),  mkdir(sdir);  end
% A flag file EXS_YIELD_<tag> ends the run at once, without the stop-time
% line, so the queue runner takes its next tag; the rows keep their state
yieldf = fullfile(here,['EXS_YIELD_',tag]);
if exist(yieldf,'file'),  fprintf('## EX_SEARCH %s yields (flag file)\n',tag);  return,  end
if ~exist(logf,'file')
    fid = fopen(logf,'w');
    fprintf(fid,'time,key,step,move,outcome,status,exst,eta,fr,passes,ndv,m,gam,seconds,note\n');
    fclose(fid);
end
if ~exist(finf,'file')
    fid = fopen(finf,'w');
%   fprintf(fid,'key,form,tier,statements,ndv0,m0,ndv,m,builds,reason\n');  fclose(fid); % MMP, 11/21/2026 (was)
    fprintf(fid,'key,form,tier,statements,ndv0,m0,ndv,m,builds,reason,path,ndv_direct,m_direct\n');  fclose(fid); % MMP, 11/21/2026
end
T = readtable(startsfile,'Delimiter',',','TextType','char','ReadVariableNames',true);
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
keys = T.key;  if isfield(opts,'keys'),  keys = opts.keys;  end
fprintf('## EX_SEARCH %s start %s: %d rows, stop %s\n',tag,char(datetime('now')),numel(keys),char(tstop));
for i = 1:numel(keys)
    key = keys{i};  jj = find(strcmp(T.key,key));
    if isempty(jj),  fprintf('## %s: not in the starts file\n',key);  continue,  end
    j = jj(1);
    sf = fullfile(sdir,[strrep(strrep(key,':','_'),'/','_'),'.mat']);
    if exist(sf,'file'),  L = load(sf);  st = L.st;  if st.done,  continue,  end
    else,  st = [];  end
    if datetime('now') >= tstop
        fprintf('\nlb_run: stop time %s reached; rerun to resume.\n',char(tstop));  return
    end
    if exist(yieldf,'file'),  fprintf('## EX_SEARCH %s yields (flag file)\n',tag);  return,  end
    k = find(strcmp({Cat.key},key),1);
    c = Cat(k);  c.exec = ['sop:',c.exec];
    row = table2struct(T(j,:));
    alts = table2struct(T(jj,:));   % alternative starts, tried in order
    fprintf('\n##### %s  [%s %s %s]  %s\n',key,row.tier,row.form,row.mods,char(datetime('now')));
    try
        st = search_row(c,alts,st,sf,logf,tstop,maxbuild,rowmin);
        if ~isempty(st.S),  row = alts(st.alt);  end
    catch ME
        fprintf('##### %s  error: %s\n',key,ME.message);
        lg(logf,key,0,'row','error','error',strrep(ME.message,',',';'));
        st = struct('done',true,'S',[],'S0',[],'cost0',[NaN NaN],'cost',[NaN NaN], ...
                    'nbuild',0,'reason',['error: ',ME.message]);  save(sf,'st');
    end
    if ~st.done
        fprintf('\nlb_run: stop time %s reached; rerun to resume.\n',char(tstop));  return
    end
    % The final setting: the choice between the smallest direct and the     % MMP, 11/21/2026
    % smallest certificate (opts.post), else the end of the descent         % MMP, 11/21/2026
    [Sf,cf,pth,cd] = final_choice(st);                                      % MMP, 11/21/2026
    stm = '';
%   if ~isempty(st.S),  stm = strjoin(st_diff(lpisettings_sop(row.tier,0,'','mosek'),st.S,'st'),' ');  end % MMP, 11/21/2026 (was)
    if ~isempty(Sf),  stm = strjoin(st_diff(lpisettings(row.tier,0,'','mosek'),Sf,'st'),' ');  end % MMP, 11/21/2026
    fid = fopen(finf,'a');
%   fprintf(fid,'%s,%s,%s,"%s",%d,%d,%d,%d,%d,"%s"\n',key,row.form,row.tier,stm, ... % MMP, 11/21/2026 (was)
%       st.cost0(1),st.cost0(2),st.cost(1),st.cost(2),st.nbuild,st.reason); % MMP, 11/21/2026 (was)
    fprintf(fid,'%s,%s,%s,"%s",%d,%d,%d,%d,%d,"%s",%s,%d,%d\n',key,row.form,row.tier,stm, ...
        st.cost0(1),st.cost0(2),cf(1),cf(2),st.nbuild,st.reason,pth,cd(1),cd(2)); % MMP, 11/21/2026
    fclose(fid);
%   fprintf('##### %s  done: %d -> %d variables  (%s)\n',key,st.cost0(1),st.cost(1),stm); % MMP, 11/21/2026 (was)
    fprintf('##### %s  done: %d -> %d variables, %s  (%s)\n',key,st.cost0(1),cf(1),pth,stm); % MMP, 11/21/2026
end
fprintf('## EX_SEARCH %s DONE %s\n',tag,char(datetime('now')));
end

% ------------------------------------------------------------------------
function st = search_row(c,alts,st,sf,logf,tstop,maxbuild,rowmin)
% The greedy descent of one row from the first of its alternative starts
% ALTS that certifies; st is the saved state ([] for a new row).
setappdata(0,'LB_ROW_KEY',c.key);
t0 = tic;
if isempty(st)                      % the start: tier, changes, form
    st = struct('S',[],'S0',[],'cost0',[NaN NaN],'cost',[NaN NaN],'tried',{{}}, ...
                'blocked',{{}},'nbuild',0,'step',0,'done',false,'reason','','alt',0, ...
                'inert',{{}});
end
if ~isfield(st,'inert'),  st.inert = {};  end   % moves that leave the size unchanged
% The smallest direct certificate and the path of the current point         % MMP, 11/21/2026
if ~isfield(st,'Sdir'),  st.Sdir = [];  st.costdir = [Inf Inf];  st.path = '';  end % MMP, 11/21/2026
if ~isfield(st,'altdone'),  st.altdone = [];  end   % starts tried          % MMP, 11/21/2026
if st.alt == 0
    for a = 1:numel(alts)
        if ismember(a,st.altdone),  continue,  end  % tried before          % MMP, 11/21/2026
        if datetime('now') >= tstop,  save(sf,'st');  return,  end
        row = alts(a);  row.exec = c.exec;
        [PIE,~] = row_pie(c,row);
        S = lpisettings(row.tier,0,'','mosek');
%       S = apply_mods(S,row.mods,c.key);                                   % MMP, 11/21/2026 (was)
        if strcmp(strtrim(row.mods),'examplest')                            % MMP, 11/21/2026
            S = example_st(c);              % the st of the example file    % MMP, 11/21/2026
        else                                                                % MMP, 11/21/2026
            S = apply_mods(S,row.mods,c.key);                               % MMP, 11/21/2026
        end                                                                 % MMP, 11/21/2026
        setappdata(0,'LB_NDV_CAP',getappdata(0,'EXS_NDVCAP'));   setappdata(0,'LB_M_CAP',getappdata(0,'EXS_MCAP'));
        [ok,r,note] = try_solve(c,PIE,S,row);
        st.nbuild = st.nbuild + 1;
        lg(logf,c.key,0,sprintf('start%d %s %s %s',a,row.form,row.tier,strrep(row.mods,',',';')), ...
           ternary(ok,'start','start fails'),r,note);
        if ok
            st.S = S;  st.S0 = S;  st.cost0 = [r.ndec r.m];  st.cost = st.cost0;  st.alt = a;
            st.path = r.path;                                               % MMP, 11/21/2026
            if strcmp(r.path,'direct'),  st.Sdir = S;  st.costdir = st.cost0;  end % MMP, 11/21/2026
            st.ref = r.gam;         % value rows: reductions keep the start value
            save(sf,'st');  break
        end
        st.altdone(end+1) = a;                                              % MMP, 11/21/2026
        save(sf,'st');
    end
    if st.alt == 0
        st.done = true;  st.reason = 'no start certified';  save(sf,'st');  return
    end
end
row = alts(st.alt);  row.exec = c.exec;
if strcmp(row.kind,'value') && isfield(st,'ref') && isfinite(st.ref)
    row.best = min(num(row.best),st.ref);    % a loose start target: keep its value
end
PIE = row_pie(c,row);
nv = c.dim;                         % a delayed PDE has a space in nv+1 variables
try,  nv = max(nv,double(PIE.dim));  catch,  end
while true
    if st.nbuild >= maxbuild,  st.reason = 'build budget';  break,  end
    if toc(t0) > 60*rowmin,  st.reason = 'row time';  break,  end
    if datetime('now') >= tstop,  save(sf,'st');  return,  end      % resume later
    M = moves(st.S,row,nv);
    M = M(~ismember({M.name},[st.tried,st.inert]));
    if isempty(M)                   % reductions exhausted: cross terms
        M = cross_moves(st.S,row,nv,st.blocked);
        M = M(~ismember({M.name},[st.tried,st.inert]));
    end
    if isempty(M)                   % then trades between storage and slack
        M = trade_moves(st.S,row,nv,st.blocked);
        M = M(~ismember({M.name},[st.tried,st.inert]));
        if isempty(M),  st.reason = 'converged';  break,  end
    end
    mv = M(1);
    setappdata(0,'LB_NDV_CAP',st.cost(1)-1);  setappdata(0,'LB_M_CAP',getappdata(0,'EXS_MCAP'));
    [ok,r,note] = try_solve(c,PIE,mv.S,row);
    st.nbuild = st.nbuild + 1;  st.step = st.step + 1;
    if ok && strcmp(r.path,'direct') && r.ndec < st.costdir(1)              % MMP, 11/21/2026
        st.Sdir = mv.S;  st.costdir = [r.ndec r.m];                         % MMP, 11/21/2026
    end                                                                     % MMP, 11/21/2026
    if ok && ~(r.ndec < st.cost(1))     % kept only when strictly smaller
        ok = false;  note = 'certified, not smaller';
    end
    if ok
        st.S = mv.S;  st.cost = [r.ndec r.m];  st.tried = {};  st.blocked = {};
        st.path = r.path;                                                   % MMP, 11/21/2026
        lg(logf,c.key,st.step,mv.name,'kept',r,note);
    else
        st.tried{end+1} = mv.name;
        if contains(note,'LB_SIZE') && r.ndec == st.cost(1)
            st.inert{end+1} = mv.name;      % the knob does not act on this PIE
        elseif ~startsWith(mv.name,'x:')
            st.blocked{end+1} = mv.name;
        end
        lg(logf,c.key,st.step,mv.name,ternary(contains(note,'LB_SIZE'),'not smaller','fails'),r,note);
    end
    save(sf,'st');
end
st.done = true;  save(sf,'st');
end

% ------------------------------------------------------------------------
function [PIE,v] = row_pie(c,row)
% The PIE of a row at its test parameter, and the form map of its start:
% 'P'/'Pd' for stability and gains, 'oco'/'cco' for the H2 norm.
for a = {'LB_ST_HOOK','LB_FORM_MAP'},  if isappdata(0,a{1}),  rmappdata(0,a{1});  end,  end
if any(strcmp(row.form,{'P','Pd'}))
    setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd', ...
        'only',{{'PIETOOLS_stability','PIETOOLS_Hinf_gain'}}));
elseif any(strcmp(row.form,{'oco','cco'}))
    setappdata(0,'LB_FORM_MAP',struct('o','oco','c','cco','only','PIETOOLS_H2_norm'));
elseif contains(row.form,'>')       % the executive of the row in another form, e.g.
    ab = strsplit(row.form,'>');    % 'Qd>Q': a dual row tested in the primal form
    setappdata(0,'LB_FORM_MAP',struct(ab{1},ab{2}));
end
v = [];                             % the library parameter
if strcmp(row.kind,'bisect')        % the best threshold less 1e-3 of it
    b = num(row.best);  v = b - 1e-3*abs(b);
end
persistent cache
if isempty(cache),  cache = containers.Map();  end
ck = sprintf('%s@%.12g',c.key,v);  if isempty(v),  ck = c.key;  end
if cache.isKey(ck),  PIE = cache(ck);  return,  end
PIE = lb_build_system(c, v, true);
if cache.Count > 4,  remove(cache,cache.keys);  end
cache(ck) = PIE;
end

% ------------------------------------------------------------------------
function [ok,r,note] = try_solve(c,PIE,S,row)
% One executive call; ok by both gates and, for value rows, the best value.
% With EXS_POST (opts.post): the gate of 11/20, and a build refused for     % MMP, 11/21/2026
% numerical reasons only kept when post_verify passes; r.path 'direct'      % MMP, 11/21/2026
% or 'post' tells which                                                     % MMP, 11/21/2026
note = '';  ok = false;
r = struct('status','error','exst','','eta',NaN,'fr',NaN,'passes',NaN,'ndec',NaN,'m',NaN,'gam',NaN,'sec',NaN);
r.path = '';                                                                % MMP, 11/21/2026
post = isappdata(0,'EXS_POST') && isequal(getappdata(0,'EXS_POST'),true);   % MMP, 11/21/2026
try
    setappdata(0,'LB_KEEP_LAST',post);                                      % MMP, 11/21/2026
%   t = tic;  r = lb_exec(c.exec, PIE, S);  r.sec = toc(t);                 % MMP, 11/21/2026 (was)
    t = tic;  r = lb_exec(c.exec, PIE, S);  r.sec = toc(t);  r.path = '';   % MMP, 11/21/2026
    setappdata(0,'LB_KEEP_LAST',false);                                     % MMP, 11/21/2026
    ok = strcmp(r.status,'certified') && (isempty(r.exst) || any(strcmp(r.exst,{'optimal','inaccurate'})));
    if post                         % the gate of 11/20                     % MMP, 11/21/2026
        ok = strcmp(r.status,'certified') && strcmp(r.exst,'optimal') && r.eta <= 1e-7; % MMP, 11/21/2026
    end                                                                     % MMP, 11/21/2026
    if ok,  r.path = 'direct';  end                                         % MMP, 11/21/2026
    if post && ~ok && isappdata(0,'LB_LAST')                                % MMP, 11/21/2026
        [ok,r,note] = try_post(r);                                          % MMP, 11/21/2026
    end                                                                     % MMP, 11/21/2026
    if isappdata(0,'LB_LAST'),  rmappdata(0,'LB_LAST');  end                % MMP, 11/21/2026
    if ok && strcmp(row.kind,'value')
        b = num(row.best);
        ok = isfinite(r.gam) && r.gam <= b + 1e-3*abs(b);
    end
catch ME
    setappdata(0,'LB_KEEP_LAST',false);                                     % MMP, 11/21/2026
    if isappdata(0,'LB_LAST'),  rmappdata(0,'LB_LAST');  end                % MMP, 11/21/2026
    note = strrep(ME.message,',',';');
    tok = regexp(note,'LB_SIZE ndv=(\d+) m=(\d+)','tokens','once');
    if ~isempty(tok),  r.ndec = str2double(tok{1});  r.m = str2double(tok{2});  end
end
end

% ------------------------------------------------------------------------  % MMP, 11/21/2026
function [ok,r,note] = try_post(r)
% [OK,R,NOTE] = TRY_POST(R) checks the solve kept in appdata LB_LAST a
% posteriori (post_verify) when the gate refused it for numerical reasons
% only: executive optimal, inaccurate or unknown, and no infeasibility
% certificate; stability forms and the H2 form 'c' (the forms post_verify
% covers). OK true: kept, r.path 'post' (r.gam the verified bound for H2).
% MMP, 11/21/2026.
ok = false;     note = '';
L = getappdata(0,'LB_LAST');
frm = '';   if ~isempty(L.form),  frm = L.form{1};  end
numer = any(strcmp(r.exst,{'optimal','inaccurate','unknown'})) && isequal(r.pinf,0) && isequal(r.dinf,0);
can = strcmp(L.fn,'PIETOOLS_stability') || (strcmp(L.fn,'PIETOOLS_H2_norm') && strcmp(frm,'c'));
if ~numer || ~can,  return,  end
% The check is a verification, not a candidate: the size gate of the
% descent (just below the current size) does not apply to its programs
cap0 = getappdata(0,'LB_NDV_CAP');  mcap0 = getappdata(0,'LB_M_CAP');
setappdata(0,'LB_NDV_CAP',getappdata(0,'EXS_NDVCAP'));  setappdata(0,'LB_M_CAP',getappdata(0,'EXS_MCAP'));
try
    [ok,r2,pre] = post_verify(L.PIE,L.S,frm,L.prog,L.info,0);
    note = sprintf('post %s: stage 2 %s fr %.3g exst %s eta %.2g ndec %d m %d; %s', ...
        ternary(ok,'passes','fails'),r2.status,r2.fr,r2.exst,r2.eta,r2.ndec,r2.m,pre);
    if ok
        r.path = 'post';
        if isfield(r2,'gam'),  r.gam = r2.gam;  end
    end
catch ME
    note = ['post error: ',ME.message];
end
setappdata(0,'LB_NDV_CAP',cap0);    setappdata(0,'LB_M_CAP',mcap0);     % the descent's gate again
note = strrep(strrep(note,',',';'),newline,' ');
end

function [Sf,cf,pth,cd] = final_choice(st)
% [SF,CF,PTH,CD] = FINAL_CHOICE(ST) the final setting of a row: with
% opts.post, the smallest direct certificate when it costs at most
% EXS_RULE.ratio times the variables of the smallest certificate and is not
% near the memory limit (EXS_RULE.near: variables, rows); else the end of
% the descent. PTH 'direct', 'post' or '' (no certificate); CD the cost of
% the smallest direct certificate ([NaN NaN] if none).
% MMP, 11/21/2026.
Sf = st.S;  cf = st.cost;   pth = '';   cd = [NaN NaN];
if isfield(st,'path'),  pth = st.path;  end
if ~isfield(st,'Sdir') || isempty(st.Sdir),  return,  end
cd = st.costdir;
R = getappdata(0,'EXS_RULE');
near = cd(1) >= R.near(1) || cd(2) >= R.near(2);
if strcmp(pth,'post') && ~near && cd(1) <= R.ratio*cf(1)
    Sf = st.Sdir;   cf = cd;    pth = 'direct';
end
end

function S = example_st(c)
% S = EXAMPLE_ST(C) the settings the example file of the catalog row C
% returns as its last output, [PDE_t,st] = fn([index,]GUI,params).
% MMP, 11/21/2026.
if ~strcmp(c.lib,'PDE'),  error('ex_search:examplest','examplest: %s is not a PDE example.',c.key);  end
% st is the last output: [PDE_t,st], or [PDE_t,PDE_b,st] in the files that
% also return the batch format
outs = cell(1,nargout(c.fn));
if c.idxarg,  [outs{:}] = feval(c.fn,str2double(c.id),0,{});
else,         [outs{:}] = feval(c.fn,0,{});
end
S = outs{end};
if isempty(S) || ~isstruct(S) || ~isfield(S,'lpi')
    error('ex_search:examplest','%s returns no settings (its last output is not an lpisettings structure).',c.fn);
end
end

% ------------------------------------------------------------------------
function M = moves(S,row,nv)
% Reductions of the knobs the executive of the row reads, as a struct
% array with fields name, S (most promising first).
M = struct('name',{},'S',{});
fform = row.form;  if contains(fform,'>'),  fform = extractAfter(fform,'>');  end
% famQ = any(strcmp(fform,{'Q','Qd'}));   % the form the executive runs in  % MMP, 11/21/2026 (was)
famQ = any(strcmp(fform,{'Q','Qd','FQ'}));  % the form the executive runs in % MMP, 11/21/2026
famP = ~famQ;                      % 'P' forms, H2, synthesis: the storage
syn  = contains(row.exec,{'control','estimator'});
stab = contains(row.exec,'stability');
K = 1:min(nv,3);
for k = K                           % slack lift and terms
    M = add(M,'dD0',S,'lpi.ineq.dD',0,[]);
end
for k = K
    if famQ,  M = vecdown(M,S,'sop.deg_Q',k,'degQ');  end
    if syn,   M = vecdown(M,S,'lpi.var.deg',k,'vdeg');  end
    if famP
        M = vecdown(M,S,'lpi.pos.len_com',k,'d2');
        M = vecdown(M,S,'lpi.pos.len_inc',k,'d1');
        M = vecdown(M,S,'lpi.pos.len_sos',k,'d3');
        M = patdown(M,S,k);
    end
    M = psatzdown(M,S,k);
    M = rowdown(M,S,'lpi.ineq.wmin',k,'wmin');
    if famQ,  M = rowdown(M,S,'sop.max_raise_Q',k,'raiseQ');  end
    if famP && stab
        M = rowdown(M,S,'sop.dw_feas',k,'dwfeas');  M = rowdown(M,S,'sop.max_raise_feas',k,'raisefeas');
    elseif famP
        M = rowdown(M,S,'sop.dw_P',k,'dwP');        M = rowdown(M,S,'sop.max_raise',k,'raise');
    end
end
if isempty(getf(S,'lpi.ineq.cap_inc'))
    M = add(M,'capinc',S,'lpi.ineq.cap_inc','auto',[]);
end
% one entry per name (the dD0 move repeats over k)
[~,u] = unique({M.name},'stable');  M = M(u);
end

function M = cross_moves(S,row,nv,blocked)
% Cross terms (two or more variables) with each reduction that was refused
% or failed alone; names 'x:<cross>+<reduction>'.
M = struct('name',{},'S',{});
if nv < 2,  return,  end
X = struct('name',{},'S',{});
if getf(S,'lpi.ineq.wcross')==0,  X(end+1) = struct('name','wcross1','S',setf(S,'lpi.ineq.wcross',1));  end
if ~any(strcmp(row.form,{'Q','Qd'}))
    for k = 2:min(nv,3)
        lc = getf(S,'lpi.pos.len_cross');
        if iscell(lc) && numel(lc)>=k && isequal(lc{k},1)
            X(end+1) = struct('name',sprintf('lencross%d',k),'S',setc(S,'lpi.pos.len_cross',k,2)); %#ok<AGROW>
        end
    end
end
for x = 1:numel(X)
    R = moves(X(x).S,row,nv);
    R = R(ismember({R.name},blocked));
    for q = 1:numel(R)
        M(end+1) = struct('name',['x:',X(x).name,'+',R(q).name],'S',R(q).S); %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------------------
function M = trade_moves(S,row,nv,blocked)
% Trades: raise one knob of one family (storage or slack), then lower a
% knob of the other family that failed alone; names 't:<raise>+<reduction>'.
% A trade is kept, like any move, only when the program shrinks and the row
% keeps its accuracy.
M = struct('name',{},'S',{});
U = raises(S,row,nv);
for u = 1:numel(U)
    R = moves(U(u).S,row,nv);
    R = R(ismember({R.name},blocked));
    for q = 1:numel(R)
        if strcmp(famof(R(q).name),U(u).fam),  continue,  end     % the other family
        M(end+1) = struct('name',['t:',U(u).name,'+',R(q).name],'S',R(q).S); %#ok<AGROW>
    end
end
end

function U = raises(S,row,nv)
% One-step raises: storage (lengths of the storage row, deg_Q, var.deg) and
% slack (dD, psatz, wmin, raises, dw), with the family of each.
U = struct('name',{},'S',{},'fam',{});
fform = row.form;  if contains(fform,'>'),  fform = extractAfter(fform,'>');  end
famQ = any(strcmp(fform,{'Q','Qd'}));  syn = contains(row.exec,{'control','estimator'});
stab = contains(row.exec,'stability');
K = 1:min(nv,3);
for k = K
    if famQ,  U = vecup(U,S,'sop.deg_Q',k,'degQ+');
    else
        U = vecup(U,S,'lpi.pos.len_com',k,'d2+');  U = vecup(U,S,'lpi.pos.len_sos',k,'d3+');
        U = vecup(U,S,'lpi.pos.len_inc',k,'d1+');
    end
    if syn,  U = vecup(U,S,'lpi.var.deg',k,'vdeg+');  end
end
d = getf(S,'lpi.ineq.dD');
if isnumeric(d) && d < 3,  U(end+1) = struct('name','dD+','S',setf(S,'lpi.ineq.dD',d+1),'fam','slack');  end
P = getf(S,'lpi.ineq.psatz');  if ischar(P),  P = {P,P,P};  end
up = struct('none','faces1','faces1','auto');
for k = K
    p = P{min(k,end)};
    if ischar(p) && isfield(up,p)
        U(end+1) = struct('name',sprintf('psatz+%d',k),'S',setc(S,'lpi.ineq.psatz',k,up.(p)),'fam','slack'); %#ok<AGROW>
    end
    U = rowup(U,S,'lpi.ineq.wmin',k,'wmin+');
    if famQ,       U = rowup(U,S,'sop.max_raise_Q',k,'raiseQ+');
    elseif stab,   U = rowup(U,S,'sop.dw_feas',k,'dwfeas+');
    else,          U = rowup(U,S,'sop.dw_P',k,'dwP+');
    end
end
end

function U = vecup(U,S,path,k,tag)
C = getf(S,path);  if ~iscell(C) || isempty(C),  return,  end
x = C{min(k,end)};  if any(strcmp(tag,{'degQ+','vdeg+'})) && isscalar(x),  x = x*[1 1 1];  end
for e = 1:numel(x)
    y = x;  y(e) = y(e)+1;
    U(end+1) = struct('name',sprintf('%s%d.%d',tag,k,e),'S',setc(S,path,k,y),'fam','storage'); %#ok<AGROW>
end
end

function U = rowup(U,S,path,k,tag)
x = getf(S,path);  if isempty(x) || ~isnumeric(x),  return,  end
kk = min(k,numel(x));  y = x;  y(kk) = y(kk)+1;
U(end+1) = struct('name',sprintf('%s%d',tag,k),'S',setf(S,path,y),'fam','slack');
end

function f = famof(name)
% the knob family of a reduction move
if ~isempty(regexp(name,'^(d[123]|pat|degQ|vdeg)','once')),  f = 'storage';  else,  f = 'slack';  end
end

function M = add(M,name,S,path,val,k)
if isempty(k),  S1 = setf(S,path,val);  else,  S1 = setc(S,path,k,val);  end
if ~isequaln(S1,S),  M(end+1) = struct('name',name,'S',S1);  end   % NaN fields (cap_delta)
end

function M = vecdown(M,S,path,k,tag)
% Each positive entry of the cell entry k (scalars expand to 3 entries
% for deg_Q and var.deg) lowered by one.
C = getf(S,path);  if ~iscell(C) || isempty(C),  return,  end
x = C{min(k,end)};
if any(strcmp(tag,{'degQ','vdeg'})) && isscalar(x),  x = x*[1 1 1];  end
for e = 1:numel(x)
    if x(e) > 0
        y = x;  y(e) = y(e)-1;
        M(end+1) = struct('name',sprintf('%s%d.%d',tag,k,e),'S',setc(S,path,k,y)); %#ok<AGROW>
    end
end
end

function M = rowdown(M,S,path,k,tag)
x = getf(S,path);  if isempty(x) || ~isnumeric(x),  return,  end
kk = min(k,numel(x));
if x(kk) > 0
    y = x;  y(kk) = y(kk)-1;
    M(end+1) = struct('name',sprintf('%s%d',tag,k),'S',setf(S,path,y));
end
end

function M = patdown(M,S,k)
P = getf(S,'lpi.pos.patterns');  if ~iscell(P),  return,  end
p = P{min(k,end)};  nxt = struct('all','separable','separable','uniform');
if ischar(p) && isfield(nxt,p)
    M(end+1) = struct('name',sprintf('pat%d',k),'S',setc(S,'lpi.pos.patterns',k,nxt.(p)));
end
end

function M = psatzdown(M,S,k)
P = getf(S,'lpi.ineq.psatz');
if ischar(P),  P = {P,P,P};  S = setf(S,'lpi.ineq.psatz',P);  end
if ~iscell(P),  return,  end
p = P{min(k,end)};  nxt = struct('auto','faces1','faces1','none','faces','faces1');
if ischar(p) && isfield(nxt,p)
    M(end+1) = struct('name',sprintf('psatz%d',k),'S',setc(S,'lpi.ineq.psatz',k,nxt.(p)));
end
end

% ------------------------------------------------------------------------
function b = num(x)
% a number of the starts file (readtable gives double, or char when quoted)
if ischar(x),  b = str2double(x);  else,  b = double(x);  end
end

function S = apply_mods(S,mods,key)
% The changes of a start: MATLAB statements on st, or keywords.
if ~ischar(mods) || isempty(mods),  return,  end
st = S; %#ok<NASGU>
% parts = strtrim(strsplit(mods,';'));                                      % MMP, 11/20/2026 (was)
% A ';' inside [], {} or () separates entries, not statements               % MMP, 11/20/2026
parts = strtrim(regexp(mods,';(?![^\[\{\(]*[\]\}\)])','split'));            % MMP, 11/20/2026
for p = parts
    s = p{1};  if isempty(s),  continue,  end
    switch s
        case 'stripped1d'
            Tt = settings_PIETOOLS_stripped();
            for f = fieldnames(Tt.lpi.pos)'
                vv = Tt.lpi.pos.(f{1});
                if iscell(vv),  st.lpi.pos.(f{1}){1} = vv{1};  end
            end
            st.sop.deg_Q{1} = Tt.sop.deg_Q{1};    st.lpi.var.deg{1} = Tt.lpi.var.deg{1};
        case 'blockst'
            switch key
                case 'DDE:DDE03',  st.sop.deg_Q{1} = [5 2 2];  st.sop.max_raise_Q(1) = 0;
                case 'DDE:DDE04',  st.sop.deg_Q{1} = [3 2 2];  st.sop.max_raise_Q(1) = 0;
                case 'DDE:DDE05',  st.sop.deg_Q{1} = [4 2 2];  st.sop.max_raise_Q(1) = 0;
                case 'DDE:DDE07',  st.sop.deg_Q{1} = [4 2 2];  st.sop.max_raise_Q(1) = 0;
                case 'DDE:DDE08',  st.sop.deg_Q{1} = [8 6 6];
            end
        otherwise
            eval([s,';']);
    end
end
S = st;
end

function v = getf(S,path)
v = [];  f = strsplit(path,'.');
for i = 1:numel(f)
    if ~isstruct(S) || ~isfield(S,f{i}),  v = [];  return,  end
    S = S.(f{i});
end
v = S;
end

function S = setf(S,path,v)
f = strsplit(path,'.');  S = setfield(S,f{:},v);
end

function S = setc(S,path,k,v)
% Entry k of a cell field; a shorter cell is padded with its last entry.
C = getf(S,path);  if ~iscell(C),  C = {C};  end
while numel(C) < k,  C{end+1} = C{end};  end %#ok<AGROW>
C{k} = v;  S = setf(S,path,C);
end

function lg(logf,key,step,mv,outcome,r,note,varargin)
if ~isstruct(r)
    r = struct('status',r,'exst','','eta',NaN,'fr',NaN,'passes',NaN,'ndec',NaN,'m',NaN,'gam',NaN,'sec',NaN);
end
fid = fopen(logf,'a');
fprintf(fid,'%s,%s,%d,%s,%s,%s,%s,%.3g,%.3g,%g,%g,%g,%.8g,%.1f,%s\n',char(datetime('now','Format','HH:mm:ss')), ...
    key,step,mv,outcome,fv(r,'status'),fv(r,'exst'),nv_(r,'eta'),nv_(r,'fr'),nv_(r,'passes'), ...
    nv_(r,'ndec'),nv_(r,'m'),nv_(r,'gam'),nv_(r,'sec'),strrep(note,newline,' '));
fclose(fid);
fprintf('   %3d %-24s %-12s %s/%s ndv %g m %g gam %.6g (%.1f s)\n',step,mv,outcome,fv(r,'status'), ...
    fv(r,'exst'),nv_(r,'ndec'),nv_(r,'m'),nv_(r,'gam'),nv_(r,'sec'));
end

function s = fv(r,f),  s = '';  if isfield(r,f) && ischar(r.(f)),  s = r.(f);  end,  end
function x = nv_(r,f),  x = NaN;  if isfield(r,f) && isnumeric(r.(f)) && isscalar(r.(f)),  x = double(r.(f));  end,  end
function y = ternary(c,a,b),  if c,  y = a;  else,  y = b;  end,  end

% ------------------------------------------------------------------------
function L = st_diff(A,B,pre)
% MATLAB statements that turn the settings A into B (fields of B only).
L = {};
if isstruct(B) && isstruct(A) && isscalar(B) && isscalar(A)
    for f = fieldnames(B)'
        if isfield(A,f{1})
            L = [L, st_diff(A.(f{1}),B.(f{1}),[pre,'.',f{1}])]; %#ok<AGROW>
        else
            L{end+1} = sprintf('%s.%s = %s;',pre,f{1},lit(B.(f{1}))); %#ok<AGROW>
        end
    end
elseif iscell(B) && iscell(A) && numel(A)==numel(B)
    for i = 1:numel(B)
        if ~isequaln(A{i},B{i}),  L{end+1} = sprintf('%s{%d} = %s;',pre,i,lit(B{i}));  end %#ok<AGROW>
    end
elseif ~isequaln(A,B)
    L{end+1} = sprintf('%s = %s;',pre,lit(B));
end
end

function s = lit(v)
if ischar(v),  s = ['''',v,''''];
elseif isnumeric(v) || islogical(v),  s = mat2str(v);
elseif iscell(v),  s = ['{',strjoin(cellfun(@lit,v,'UniformOutput',false),','),'}'];
else,  s = '<struct>';
end
end

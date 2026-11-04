function lb_run(ids, tiers, opts)
% LB_RUN(IDS, TIERS, OPTS) runs the library baseline: for each selected
% catalog row and each settings tier, the stock executive on the stock
% settings file, one results.csv line per (example, tier).
%
% INPUT
%   ids:    'all', a library name ('PDE','DDE','NDSDDF'), or a cellstr of
%           catalog keys / library names, e.g. {'PDE:6','DDE:DDE05',
%           'PDE:20/Hinf_gain_dual'}; rows with run=no are never run
%   tiers:  'all' or a cellstr of lpisettings types, from {'extreme',
%           'stripped','light','heavy','veryheavy','custom'}
%   opts:   optional struct
%       .results    - results file (default results.csv in this folder)
%       .max2d      - skip a 2-D row whose projected time per solve
%                     exceeds this [s] (default 600)
%       .f2d        - floor on the growth of the time per solve per 2-D
%                     tier step into custom, heavy or veryheavy (default
%                     10: PDE 35 took 7.8 s on light, 731 s on heavy, two
%                     2-D tier steps, measured 10/28/2026); the floor
%                     into stripped and light is 1
%       .maxsolves  - bisection solve cap (default 14)
%       .reltol     - bisection stop: hi-lo <= reltol*|documented|, or
%                     reltol*(hi0-lo0) if the documented value is 0 or
%                     not numeric (default 1e-3)
%       .dry        - true: list the rows to run, solve nothing
%       .stopat     - date-time string, e.g. '2026-10-29 08:00'; no row is
%                     started after it (default '': no limit)
%       .tiers2d    - tiers on which 2-D rows run (default: all six);
%                     2-D rows on other tiers are not run and get no line
%       .skip       - cellstr of 'KEY@tier' not to run (no line), e.g.
%                     {'DDE:DDE16@veryheavy'} (default {})
%       .path       - 'stock' (default): the stock executives with
%                     lpisettings; 'sop': the _sop executives with
%                     lpisettings_sop (lb_exec maps each stock executive
%                     to its _sop form). Each path has its own results,
%                     started and PAUSE files.
%
% OUTPUT (file)
%   results.csv, columns: library, id, settings, status, result,
%   documented, n_solves, seconds, decision_variables, note.
%     status  - single solve: certified / inaccurate / infeasible /
%               unknown / nosolinfo (lb_exec); bisection: converged /
%               max_solves / lo_fail / hi_pass; error; skipped
%     result  - stability: stock solve_val (1, 0.5, 0); gain or norm: the
%               value the executive returns; bisection: the largest
%               certified parameter value (NaN for lo_fail)
%     seconds - wall time of the row, setup and all solves
%     decision_variables - numel(prog.decvartable) of the last solve
%     note    - executive, solver exit, SDP size (m rows, nx cone
%               entries), bisection steps value:C/I/F/U/E
%   runs.log, one line per call: time, tiers, rows, git state, MATLAB
%   and MOSEK versions.
%
% NOTES
%   - Resumable: a (library, id, settings) already in results.csv is
%     skipped; delete its line to rerun it.
%   - Settings: lpisettings(tier,0,'','mosek'), unmodified. Tiers run
%     outer, from light to heavy, so that each tier completes first and a
%     2-D row is projected from the lighter tiers measured before it.
%   - A bisection step passes iff lb_exec returns 'certified'.
%   - Crash recovery: each row start is appended to started.txt; a row
%     started twice without a results line (MATLAB died during it twice)
%     is recorded with status 'crash' and not started again.

if nargin<3, opts = struct(); end
here = fileparts(mfilename('fullpath'));
def = struct('results',fullfile(here,'results.csv'),'max2d',600,'f2d',10, ...
             'maxsolves',14,'reltol',1e-3,'dry',false,'stopat','', ...
             'tiers2d',{{'extreme','stripped','light','heavy','veryheavy','custom'}}, ...
             'skip',{{}},'path','stock');
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts,fn{k}), opts.(fn{k}) = def.(fn{k}); end
end
order1 = {'extreme','stripped','light','heavy','veryheavy','custom'};
order2 = {'extreme','stripped','light','custom','heavy','veryheavy'};
% settings_sopvar: custom is light, not a tier step between light and heavy
if strcmp(opts.path,'sop'), order2 = order1; end
if ischar(tiers) && strcmp(tiers,'all'), tiers = order1; end
tiers = cellstr(tiers);

info = lb_setup();
S = load(fullfile(here,'catalog.mat'));  C = S.C;
sel = lb_select(C, ids);
done = lb_done(opts.results);
[rp,rn] = fileparts(opts.results);
startf = fullfile(rp,'started.txt');                          % results.csv
if ~strcmp(rn,'results'), startf = fullfile(rp,[rn,'_started.txt']); end
started = lb_started(startf);
pausef = fullfile(here,'PAUSE');
if strcmp(opts.path,'sop'), pausef = fullfile(here,'PAUSE_sop'); end

if ~opts.dry                       % a dry run solves nothing: not logged
    lg = fopen(fullfile(here,'runs.log'),'a');
    fprintf(lg,'%s | path %s | tiers %s | %d rows | git %s | %s | %s | threads %d\n', ...
        char(datetime('now')), opts.path, strjoin(tiers,' '), numel(sel), info.git, ...
        info.version, info.mosek, info.threads);
    fclose(lg);
end

% Tiers outer, so that each tier completes for all rows before the next
% starts. Order1 runs extreme, stripped and light first, so a 2-D row on
% a heavier tier is projected from them; lb_project ranks 2-D tiers by
% order2, in which custom lies between light and heavy.
tt = order1(ismember(order1,tiers));
for t = 1:numel(tt)
    tier = tt{t};
    for k = sel
        c = C(k);
        if strcmp(opts.path,'sop'), c.exec = ['sop:',c.exec]; end   % lb_exec maps it
        if c.dim>=2, ord = order2; else, ord = order1; end
        rowid = c.key(numel(c.lib)+2:end);
        key = strjoin({c.lib,rowid,tier},'|');
        if done.isKey(key), continue, end
        if c.dim>=2 && ~ismember(tier,opts.tiers2d), continue, end   % not selected, no row
        if ismember([c.key,'@',tier],opts.skip), continue, end           % not selected, no row
        if opts.dry
            fprintf('would run %s %s\n', c.key, tier);   continue
        end
        if ~isempty(opts.stopat) && datetime('now') >= datetime(opts.stopat)
            fprintf('\nlb_run: stop time %s reached; rerun to resume.\n', opts.stopat);
            return
        end
        if exist(pausef,'file')
            fprintf('\nlb_run: %s present; delete it and rerun to resume.\n', pausef);
            return
        end
        fprintf('\n##### %s  [%s]  %s\n', c.key, tier, char(datetime('now')));
        doc = lb_docfor(c.doc, tier);
        % Two starts without a results line: MATLAB died twice in this row
        if started.isKey(key) && started(key)>=2
            lb_append(opts.results, {c.lib, rowid, tier, 'crash', 'NaN', doc, '0', 'NaN', ...
                'NaN', [c.exec,' | MATLAB exited during this row twice (see full.log)']});
            done(key) = true;
            fprintf('##### %s  [%s]  crash (twice): recorded, not started\n', c.key, tier);
            continue
        end
        fid = fopen(startf,'a');  fprintf(fid,'%s\n',key);  fclose(fid);
        t0 = tic;
        out = struct('status','','result',NaN,'n',0,'ndec',NaN,'note','');
        skip = '';
        if c.dim>=2
            skip = lb_project(opts.results, c, tier, ord, opts);
        end
        if ~isempty(skip)
            out.status = 'skipped';   out.note = skip;
        else
            try
                if strcmp(opts.path,'sop')
                    S0 = lpisettings_sop(tier,0,'','mosek');
                else
                    S0 = lpisettings(tier,0,'','mosek');
                end
                if c.bisect
                    out = lb_bisect(c, S0, opts);
                else
                    PIE = lb_build_system(c, [], true);
                    r = lb_exec(c.exec, PIE, S0);
                    out.status = r.status;   out.n = 1;   out.ndec = r.ndec;
                    if isnan(r.gam), out.result = r.val; else, out.result = r.gam; end
                    out.note = lb_rnote(r);
                end
            catch ME
                out.status = 'error';
                out.note = lb_errnote(ME, info.root);
            end
        end
        sec = toc(t0);
        out.note = strtrim([c.exec,' | ',out.note]);
        lb_append(opts.results, {c.lib, rowid, tier, out.status, num2str(out.result,8), ...
            doc, num2str(out.n), sprintf('%.2f',sec), num2str(out.ndec), out.note});
        done(key) = true;
        fprintf('##### %s  [%s]  %s  result %s  (%.1f s)\n', c.key, tier, out.status, ...
            num2str(out.result,8), sec);
    end
end
end

% ------------------------------------------------------------------------
function sel = lb_select(C, ids)
% Indices of catalog rows selected by ids, rows with run=no excluded.
if ischar(ids), ids = {ids}; end
keys = {C.key};  libs = {C.lib};
sel = false(1,numel(C));
for k = 1:numel(ids)
    if strcmp(ids{k},'all')
        sel(:) = true;
    elseif any(strcmp(ids{k},{'PDE','DDE','NDSDDF'}))
        sel = sel | strcmp(libs,ids{k});
    else
        hit = strcmp(keys,ids{k});
        if ~any(hit), error('lb_run:key','No catalog row %s.',ids{k}); end
        sel = sel | hit;
    end
end
sel = find(sel & [C.run] & ~strcmp({C.exec},'ERROR'));
end

function done = lb_done(file)
% Keys library|id|settings of the rows already in the results file.
done = containers.Map('KeyType','char','ValueType','logical');
if ~exist(file,'file'), return, end
L = splitlines(fileread(file));
for k = 2:numel(L)
    f = strsplit(L{k},',','CollapseDelimiters',false);
    if numel(f)>=3, done(strjoin(f(1:3),'|')) = true; end
end
end

function started = lb_started(file)
% Number of starts of each key library|id|settings recorded in FILE.
started = containers.Map('KeyType','char','ValueType','double');
if ~exist(file,'file'), return, end
L = splitlines(strtrim(fileread(file)));
for k = 1:numel(L)
    if isempty(L{k}), continue, end
    if started.isKey(L{k}), started(L{k}) = started(L{k})+1; else, started(L{k}) = 1; end
end
end

function lb_append(file, fields)
% Append one line; commas and newlines inside a field become ';' and ' '.
newf = ~exist(file,'file');
fid = fopen(file,'a');
if newf
    fprintf(fid,'library,id,settings,status,result,documented,n_solves,seconds,decision_variables,note\n');
end
for k = 1:numel(fields)
    fields{k} = strrep(strrep(char(fields{k}),',',';'),newline,' ');
end
fprintf(fid,'%s\n',strjoin(fields,','));
fclose(fid);
end

function d = lb_docfor(doc, tier)
% Documented value for this tier: the 'tier=value' entry if the library
% gives one, otherwise the whole documented string.
d = doc;
tok = regexp(doc,['(^|;)',tier,'=([^;]+)'],'tokens','once');
if ~isempty(tok), d = tok{2}; end
end

function s = lb_rnote(r)
s = sprintf('fr=%.3g pinf=%g dinf=%g numerr=%g iter=%g cpu=%.1fs m=%g nx=%g', ...
    r.fr, r.pinf, r.dinf, r.numerr, r.iter, r.cpusec, r.m, r.nx);
if isfield(r,'called') && startsWith(r.called,'PIETOOLS_') && contains(r.called,'_sop')
    s = [r.called,' ',s];
end
end

function s = lb_errnote(ME, root)
% First line of the error, its identifier, and the deepest frame inside
% the PIETOOLS root (file:line), for reproduction.
s = strtok(ME.message, newline);
if ~isempty(ME.identifier), s = [ME.identifier,': ',s]; end
for k = 1:numel(ME.stack)
    if startsWith(ME.stack(k).file, root)
        [~,f] = fileparts(ME.stack(k).file);
        s = sprintf('%s @ %s:%d', s, f, ME.stack(k).line);
        break
    end
end
end

function out = lb_bisect(c, S0, opts)
% Largest certified parameter value in [c.lo, c.hi]: lo must certify, hi
% must not; then halve until hi-lo <= tol or opts.maxsolves solves.
lo = c.lo;  hi = c.hi;
if isfinite(c.docnum) && c.docnum~=0
    tol = opts.reltol*abs(c.docnum);
else
    tol = opts.reltol*(hi-lo);
end
steps = {};  n = 0;  ndec = NaN;  last = [];
[ok, r, e] = lb_step(c, S0, lo);  n = n+1;  steps{end+1} = lb_stp(lo,r,e);
if ~isempty(r), ndec = r.ndec; last = r; end
if ~isempty(e) && isempty(r)
    rethrow(e);                      % the default bracket itself fails to run
end
if ~ok
    out = lb_bout('lo_fail', NaN, n, ndec, steps, last, c.lo, c.hi, lo, hi);  return
end
[ok, r, e] = lb_step(c, S0, hi);  n = n+1;  steps{end+1} = lb_stp(hi,r,e);
if ~isempty(r), ndec = r.ndec; last = r; end
if ok
    out = lb_bout('hi_pass', hi, n, ndec, steps, last, c.lo, c.hi, lo, hi);  return
end
while hi-lo > tol && n < opts.maxsolves
    mid = (lo+hi)/2;
    [ok, r, e] = lb_step(c, S0, mid);  n = n+1;  steps{end+1} = lb_stp(mid,r,e); %#ok<AGROW>
    if ~isempty(r), ndec = r.ndec; last = r; end
    if ok, lo = mid; else, hi = mid; end
end
if hi-lo <= tol, st = 'converged'; else, st = 'max_solves'; end
out = lb_bout(st, lo, n, ndec, steps, last, c.lo, c.hi, lo, hi);
end

function [ok, r, e] = lb_step(c, S0, v)
% One bisection solve at parameter value v; an error counts as a failure.
r = [];  e = [];
try
    PIE = lb_build_system(c, v, true);
    r = lb_exec(c.exec, PIE, S0);
    ok = strcmp(r.status,'certified');
catch e
    ok = false;
end
end

function s = lb_stp(v, r, e)
code = struct('certified','C','inaccurate','I','infeasible','F','unknown','U','nosolinfo','N');
if ~isempty(e) || isempty(r), s = sprintf('%.6g:E',v);
else, s = sprintf('%.6g:%s(%.0fs)',v,code.(r.status),r.sec); end
end

function out = lb_bout(st, res, n, ndec, steps, last, lo0, hi0, lo, hi)
out.status = st;  out.result = res;  out.n = n;  out.ndec = ndec;
out.note = sprintf('bracket [%.6g %.6g] final [%.6g %.6g] | steps %s', lo0, hi0, lo, hi, strjoin(steps,' '));
if ~isempty(last), out.note = [out.note,' | last: ',lb_rnote(last)]; end
end

function skip = lb_project(file, c, tier, ord, opts)
% Projected time per solve of a 2-D row from the lighter tiers already in
% the results file; '' to run, otherwise the reason to skip.
skip = '';
if ~exist(file,'file'), return, end
rowid = c.key(numel(c.lib)+2:end);
L = splitlines(fileread(file));
tsolve = nan(1,numel(ord));  bad = false(1,numel(ord));
for k = 2:numel(L)
    f = strsplit(L{k},',','CollapseDelimiters',false);
    if numel(f)<10 || ~strcmp(f{1},c.lib) || ~strcmp(f{2},rowid), continue, end
    j = find(strcmp(ord,f{3}));
    if isempty(j), continue, end
    % A lighter tier that was skipped or ran out of memory rules out the
    % heavier ones.
    bad(j) = strcmp(f{4},'skipped') || (strcmp(f{4},'error') && contains(lower(f{10}),'memory'));
    ns = str2double(f{7});
    if ns>0, tsolve(j) = str2double(f{8})/ns; end
end
jt = find(strcmp(ord,tier));
jb = find(bad(1:jt-1),1,'last');
if ~isempty(jb)
    skip = sprintf('skipped: lighter tier %s skipped or out of memory',ord{jb});   return
end
jm = find(~isnan(tsolve(1:jt-1)));
if isempty(jm), return, end                  % nothing lighter measured: run
j1 = jm(end);
% Growth per tier step: the larger of this row's measured growth and
% opts.f2d. The growth between the light tiers understates the jump to
% heavy (PDE 35: 94x from light to heavy), so f2d acts as a floor.
% Floor per step INTO tier ord{j}: 1 into stripped and light (PDE 35:
% 30.8 s on extreme, 7.8 s on light), opts.f2d into custom, heavy and
% veryheavy. A floor of 10 from extreme skipped every 2-D row over 60 s on
% extreme on all heavier tiers, light included (10/28/2026).
fl = opts.f2d*ones(1,numel(ord));
fl(ismember(ord,{'extreme','stripped','light'})) = 1;
gm = 1;
if numel(jm)>=2
    j0 = jm(end-1);
    gm = (tsolve(j1)/tsolve(j0))^(1/(j1-j0));
end
proj = tsolve(j1)*prod(max(fl(j1+1:jt),gm));
if proj > opts.max2d
    skip = sprintf('skipped: projected %.0f s per solve > %.0f s (from %s %.0f s, factor %.3g)', ...
        proj, opts.max2d, ord{j1}, tsolve(j1), proj/tsolve(j1));
end
end

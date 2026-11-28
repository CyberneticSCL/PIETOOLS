function d3run(wt, qfile, outdir, batch)
% D3RUN(WT,QFILE,OUTDIR,BATCH), 11/27/2026: Phase 3 identity runner of the
% opvar rename; strun1122 (C:\Users\mpeet\cuadmm_out\stress_1122\runner)
% adapted. WT holds PIETOOLS/, an export of a commit without cleaning/,
% orphaned_file/ and .claude/. Runs, one after another in this MATLAB,
% every row of the tab-separated QFILE whose first column is BATCH:
%   batch, id, folder, cmd [, extra]
% folder relative to WT/PIETOOLS; extra: folders separated by ';', each
% added with genpath, relative to WT/PIETOOLS or absolute.
% Per row: path reset (restoredefaultpath, pietools_path_update of the
% export, the extra folders, SeDuMi, MOSEK), checked: nothing of the
% checkout, no cleaning/orphaned_file code, the core routines and the
% container class of the export (copvar before the rename, opvar after)
% resolve inside it. Then base workspace, globals and functions cleared;
% figures invisible; rng(1); the capture folder (shadows of Sedumi2Mosek,
% MosekSol2SedumiSol, sedumi) first on the path, D3_ROW = id and D3_LOG =
% OUTDIR/sdp.tsv; cwd = folder; cmd evaluated in the base workspace.
% Appends "batch id status seconds message" to OUTDIR/results.tsv after
% each row; output of each row in OUTDIR/logs/<batch>_<id>.log. Restart:
% rows already in results.tsv are skipped; the row named in
% OUTDIR/running_<batch>.txt is recorded as CRASH. A file OUTDIR/STOP ends
% the batch before the next row.
SEDUMI = fullfile('C:','YALMIP','yalmip_install','tbxmanager','toolboxes','sedumi','1.3','pcwin64');
MOSEK = 'C:\Program Files\Mosek\11.0\toolbox\r2019b';
REPO = 'C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS';
CAP = fullfile(fileparts(mfilename('fullpath')),'capture');
P = fullfile(wt,'PIETOOLS');
cls = 'copvar';     if exist(fullfile(P,'sopvar','@opvar'),'dir'),    cls = 'opvar';  end
% _sop rename, 11/28/2026: the executive and the settings function checked below carry
% the suffix _sop before the rename and none after
ex1 = 'PIETOOLS_stability_sop';  ex2 = 'lpisettings_sop';
if exist(fullfile(P,'executives','PIETOOLS_stability.m'),'file'), ex1 = 'PIETOOLS_stability';  ex2 = 'lpisettings';  end
if ~exist(fullfile(outdir,'logs'),'dir'), mkdir(fullfile(outdir,'logs')); end
if ~exist(fullfile(outdir,'ws'),'dir'), mkdir(fullfile(outdir,'ws')); end
res = fullfile(outdir,'results.tsv');
mk = fullfile(outdir,['running_' batch '.txt']);
nrm = @(s) lower(strrep(s,'/','\'));
rows = {};
L = strsplit(fileread(qfile),newline);
for i = 1:numel(L)
    f = strsplit(strtrim(L{i}),sprintf('\t'),'CollapseDelimiters',false);
    if numel(f) >= 4 && strcmp(f{1},batch), rows{end+1} = f; end      %#ok<AGROW>
end
if exist(mk,'file')
    m = strsplit(strtrim(fileread(mk)),sprintf('\t'));
    append_row(res,batch,m{1},'CRASH',NaN,'MATLAB ended during the row');
    delete(mk);
end
done = {};
if exist(res,'file')
    R = strsplit(fileread(res),newline);
    for i = 1:numel(R)
        f = strsplit(R{i},sprintf('\t'));
        if numel(f) >= 2 && strcmp(f{1},batch), done{end+1} = f{2}; end    %#ok<AGROW>
    end
end
fprintf('D3 batch %s on %s (class %s): %d rows, %d already recorded\n',batch,P,cls,numel(rows),numel(done));
tb = tic;   % budget: no new row after 17 min, so that one process stays under 20 min
for r = 1:numel(rows)
    f = rows{r};  id = f{2};  folder = f{3};  cmd = f{4};
    if ismember(id,done), continue; end
    if exist(fullfile(outdir,'STOP'),'file'), fprintf('D3 STOP file: batch %s ended\n',batch); break; end
    if toc(tb) > 1020, fprintf('D3 BUDGET: batch %s paused after %.0f s\n',batch,toc(tb)); break; end
    fid = fopen(mk,'w');  fprintf(fid,'%s\t%s\n',id,char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));  fclose(fid);
    restoredefaultpath;
    st = 'OK';  msg = '';  el = 0;  other = {};
    try
        old = cd(P);  evalc('pietools_path_update');  cd(old);
        warning('off','MATLAB:mpath:privateDirectoriesNotAllowedOnPath');
        cl = fullfile(P,'cleaning');
        if exist(cl,'dir'), warning('off','MATLAB:rmpath:DirNotFound');  rmpath(genpath(cl)); end
        ex = {};
        if numel(f) >= 5 && ~isempty(strtrim(f{5}))
            ex = strtrim(strsplit(f{5},';'));  ex = ex(~cellfun(@isempty,ex));
            for k = 1:numel(ex)
                if isempty(regexp(ex{k},'^[A-Za-z]:','once')), ex{k} = fullfile(P,ex{k}); end
                if ~exist(ex{k},'dir'), error('d3run:extra','extra folder "%s" not found',ex{k}); end
                addpath(genpath(ex{k}));
            end
        end
        addpath(genpath(SEDUMI));
        if isempty(which('mosekopt')) && exist(MOSEK,'dir'), addpath(MOSEK); end
        pe = strsplit(path,pathsep);
        for k = 1:numel(pe)
            e = nrm(pe{k});
            if startsWith(e,nrm(REPO)), error('d3run:path','path entry in the checkout: "%s"',pe{k}); end
            if (contains(e,'\cleaning') || contains(e,'orphaned_file')) && hascode(pe{k})
                error('d3run:path','path entry with cleaning/orphaned_file code: "%s"',pe{k});
            end
            known = startsWith(e,nrm(P)) || startsWith(e,nrm(matlabroot)) || ...
                    startsWith(e,nrm(SEDUMI)) || startsWith(e,nrm(MOSEK)) || ...
                    any(cellfun(@(x) startsWith(e,nrm(x)),ex));
            if ~known, other{end+1} = pe{k}; end                             %#ok<AGROW>
        end
        for fn = {'lpiprogram','sopvar',cls,'polyopvar','pde_var','sossolve',ex1,'PIESIM',ex2,'Sedumi2Mosek','MosekSol2SedumiSol'}
            w = which(fn{1});
            if isempty(w) || ~startsWith(nrm(w),nrm(P))
                error('d3run:path','%s resolves outside the export: "%s"',fn{1},w);
            end
        end
    catch ME
        st = 'ERR';  msg = ['path: ' ME.message];
    end
    lg = fullfile(outdir,'logs',[batch '_' id '.log']);
    if exist(lg,'file'), delete(lg); end
    if strcmp(st,'OK')
        evalin('base','clear');  clear global;  clear functions;          %#ok<CLFUNC>
        set(0,'DefaultFigureVisible','off');  rng(1);
        addpath(CAP);
        setenv('D3_ROW',id);    setenv('D3_LOG',fullfile(outdir,'sdp.tsv'));
        setenv('D3_WS',fullfile(outdir,'ws',[id '.mat']));
        cd(fullfile(P,folder));
        diary(lg);
        fprintf('D3_START %s | %s | %s | %s\n',id,folder,cmd,char(datetime('now')));
        for k = 1:numel(other), fprintf('D3_PATH_OTHER %s\n',other{k}); end
        t0 = tic;
        try
            evalin('base',cmd);
        catch ME
            st = 'ERR';  msg = ME.message;
            fprintf('%s\n',getReport(ME,'extended','hyperlinks','off'));
        end
        el = toc(t0);
        pe = strsplit(path,pathsep);  pn = cellfun(nrm,pe,'UniformOutput',false);
        bad = pe(startsWith(pn,nrm(REPO)) | ((contains(pn,'\cleaning') | contains(pn,'orphaned_file')) & cellfun(@hascode,pe)));
        if ~isempty(bad)
            fprintf('D3_PATH_LEAK %s\n',strjoin(bad,' ; '));
            if strcmp(st,'OK'), st = 'PATHLEAK';  msg = ['path after the row holds ' bad{1}]; end
        end
        fprintf('D3_RESULT %s %s %.1f\n',id,st,el);
        diary off;
        setenv('D3_LOG','');
        close all force;
    end
    append_row(res,batch,id,st,el,msg);
    delete(mk);
end
fprintf('D3 batch %s done\n',batch);
end

function tf = hascode(d)
% True when folder D holds MATLAB code or class/package folders
tf = false;
if ~exist(d,'dir'), return; end
for pat = {'*.m','*.p','*.mex*','*.mlx','@*','+*'}
    if ~isempty(dir(fullfile(d,pat{1}))), tf = true; return; end
end
end

function append_row(res,batch,id,st,el,msg)
% One line per row: batch, id, status, seconds, message (240 characters at most)
msg = regexprep(msg,'\s+',' ');  msg = msg(1:min(end,240));
fid = fopen(res,'a');
fprintf(fid,'%s\t%s\t%s\t%.1f\t%s\n',batch,id,st,el,msg);
fclose(fid);
end

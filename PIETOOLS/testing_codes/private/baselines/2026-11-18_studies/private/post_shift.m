function post_shift(tag,sfx,startsfile,deltas,keys)
% POST_SHIFT(TAG,SFX,STARTSFILE,DELTAS) the soundness rows of the settings
% rerun (cleaning/private/settings_rerun_1121/PLAN.md): every row of
% exs_<TAG>_final.csv kept through stage 2 (path 'post') is solved again at
% its final settings (st.S of the row state exs_state<SFX>/<key>.mat) and
% checked by post_verify at each shift in DELTAS: the PIE T x_t = (A +
% delta*T) x, whose spectrum is the original one moved right by delta. A
% sound check passes at 0 and fails once delta exceeds the decay rate.
%
% INPUT
% - tag:        char, the tag of the ex_search run;
% - sfx:        char, its opts.statedir;
% - startsfile: char, its starts file (kind and best of each row: the test
%               parameter of 'bisect' rows, as ex_search sets it);
% - deltas:     (optional) row of shifts, default [0 10];
% - keys:       (optional) cellstr, only these rows (default: every row
%               with path 'post'); the results are appended.
% OUTPUT
% results_SHIFT_<TAG>.csv: key, form, delta, the stage-1 exit, ok and the
% stage-2 exit of post_verify, its note. H2 'c' has no shifted check (its
% stage 2 does not read delta) and is run at 0 only.
% MMP, 11/21/2026.
if nargin<4 || isempty(deltas),  deltas = [0 10];  end
if nargin<5,  keys = {};  end
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);   addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
F = readtable(fullfile(here,['exs_',tag,'_final.csv']),'Delimiter',',','TextType','char');
T = readtable(startsfile,'Delimiter',',','TextType','char');
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
of = fullfile(here,['results_SHIFT_',tag,'.csv']);
if ~exist(of,'file')
    fid = fopen(of,'w');
    fprintf(fid,'key,form,delta,s1_status,s1_fr,s1_exst,s1_eta,s1_ndec,ok,s2_status,s2_fr,s2_exst,s2_eta,s2_ndec,note\n');  fclose(fid);
end
fprintf('## POST_SHIFT %s %s\n## %s\n',tag,char(datetime('now')),which('post_verify'));
for i = 1:height(F)
    if ~strcmp(F.path{i},'post'),  continue,  end
    key = F.key{i};    form = F.form{i};
    if ~isempty(keys) && ~any(strcmp(keys,key)),  continue,  end
    L0 = load(fullfile(here,['exs_state',sfx],[strrep(strrep(key,':','_'),'/','_'),'.mat']));
    c = Cat(find(strcmp({Cat.key},key),1));     c.exec = ['sop:',c.exec];
    % the form map and test parameter of ex_search (row_pie)
    for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_LAST'},  if isappdata(0,a{1}),  rmappdata(0,a{1});  end,  end
    if any(strcmp(form,{'P','Pd'}))
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop','PIETOOLS_Hinf_gain_sop'}}));
    elseif any(strcmp(form,{'oco','cco'}))
        setappdata(0,'LB_FORM_MAP',struct('o','oco','c','cco','only','PIETOOLS_H2_norm_sop'));
    elseif contains(form,'>')
        ab = strsplit(form,'>');    setappdata(0,'LB_FORM_MAP',struct(ab{1},ab{2}));
    end
    j = find(strcmp(T.key,key),1);   v = [];
    if strcmp(T.kind{j},'bisect')
        b = T.best(j);  if iscell(b),  b = str2double(b{1});  end
        v = b - 1e-3*abs(b);
    end
    fprintf('\n##### SHIFT %s %s  %s\n',key,form,char(datetime('now')));
    PIE = lb_build_system(c,v,true);
    setappdata(0,'LB_KEEP_LAST',true);
    r1 = lb_exec(c.exec,PIE,L0.st.S);
    setappdata(0,'LB_KEEP_LAST',false);
    if ~isappdata(0,'LB_LAST'),  fprintf('##### SHIFT %s: no solve kept\n',key);  continue,  end
    L = getappdata(0,'LB_LAST');    frm = L.form{1};
    dd = deltas;    if strcmp(frm,'c'),  dd = 0;  end
    for d = dd
        note = '';  ok = false;  r2 = struct();
        try
            [ok,r2,note] = post_verify(L.PIE,L.S,frm,L.prog,L.info,d);
        catch ME
            note = ['error: ',ME.message];
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,%g,%s,%.4g,%s,%.3g,%g,%d,%s,%.4g,%s,%.3g,%g,"%s"\n',key,frm,d,f(r1,'status'),g(r1,'fr'), ...
            f(r1,'exst'),g(r1,'eta'),g(r1,'ndec'),ok,f(r2,'status'),g(r2,'fr'),f(r2,'exst'),g(r2,'eta'),g(r2,'ndec'), ...
            strrep(strrep(note,',',';'),'"',''''));
        fclose(fid);
        fprintf('##### SHIFT %s delta %g: ok %d, stage 2 %s fr %.3g eta %.2g\n',key,d,ok,f(r2,'status'),g(r2,'fr'),g(r2,'eta'));
    end
    rmappdata(0,'LB_LAST');
end
fprintf('\n## POST_SHIFT %s DONE %s\n',tag,char(datetime('now')));
end

function s = f(r,k),  s = '';   if isstruct(r) && isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isstruct(r) && isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

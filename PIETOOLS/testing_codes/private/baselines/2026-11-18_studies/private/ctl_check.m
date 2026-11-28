function ctl_check(jobs, tag)
% CTL_CHECK(JOBS, TAG) runs the executive of a DDE or NDS/DDF block of the
% catalog at a tier, with an optional statement run after the block code (a
% variant of the block), and writes results_<TAG>.csv in this folder. Built
% for the regulated outputs of DDE19 (wind tunnel) and DDF03 (showering
% people), 11/19/2026.
%
% INPUT
%   jobs:   cell array of rows {key, label, tier, post}: key as in
%           catalog.mat (e.g. 'DDE:DDE19'); label names the variant; post
%           MATLAB statements run after the block ('' for none)
%   tag:    results_<TAG>.csv; a row (key, label, tier) already in it is
%           skipped
% MMP, 11/19/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
of = fullfile(here,['results_',tag,'.csv']);
done = {};
if exist(of,'file')
    L = splitlines(strtrim(fileread(of)));
    for i = 2:numel(L),  x = strsplit(L{i},',');  done{end+1} = strjoin(x(1:3),','); end %#ok<AGROW>
else
    fid = fopen(of,'w');  fprintf(fid,'key,label,tier,exec,status,gam,ndec,m,fr,sec,note\n');  fclose(fid);
end
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
for i = 1:numel(jobs)
    [key,label,tier,post] = jobs{i}{:};
    if any(strcmp(done,strjoin({key,label,tier},','))),  continue,  end
    c = Cat(find(strcmp({Cat.key},key),1));
    fprintf('\n##### %s %s %s  %s\n',key,label,tier,char(datetime('now')));
    note = '';  r = struct();  v = [];
    if ~isempty(post),  c.setter = post;  v = 1;  end     % lb_build_system runs sprintf(setter,v)
    try
        PIE = lb_build_system(c, v, true);
        S = lpisettings(tier,0,'','mosek');
        r = lb_exec(['sop:',c.exec], PIE, S);
    catch ME
        note = strrep(strrep(ME.message,',',';'),newline,' ');  r.status = 'error';
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%s,%s,%s,%s,%.8g,%.0f,%.0f,%.4g,%.1f,"%s"\n',key,label,tier,c.exec,f(r,'status'), ...
        g(r,'gam'),g(r,'ndec'),g(r,'m'),g(r,'fr'),g(r,'sec'),note);
    fclose(fid);
    fprintf('##### %s %s %s: %s gam %.6g ndec %g\n',key,label,tier,f(r,'status'),g(r,'gam'),g(r,'ndec'));
end
fprintf('## CTL_CHECK %s DONE %s\n',tag,char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

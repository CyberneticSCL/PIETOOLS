function diag_amt(jobs, tag, ndvcap)
% DIAG_AMT(JOBS, TAG) separates, for library rows that no setting certifies,
% the storage side of the stability LPI from its dynamics. Each job builds
% the stability LPI of a row with the PIE generator A replaced:
%   'A'    - unchanged (the control);
%   'A0'   - A = 0: the derivative term vanishes, so the program asks only
%            for a storage that meets the positivity constraints of the
%            form; infeasible means no admissible storage exists in the
%            parameterization at all;
%   'AmT'  - A = -T: x_t = -x, which every coercive storage certifies
%            (P = I in 'P'; Q = T in 'Q'); infeasible means the slack of the
%            negativity constraint cannot carry even 2T'T;
%   'AsT'  - A - T: the dynamics of the row shifted by rate 1. A storage V of
%            the program certifies A - T exactly when dV/dt <= 2V under A.
%            Infeasible means the program cannot even bound the growth of
%            the row by rate 1; certified, with A infeasible, means the
%            program fails only between growth rate 1 and the marginal case;
%   'AsT<c>' - A - c*T, e.g. 'AsT0.25': the shift of a sweep. The least
%            certified c bounds the growth rate the program can certify.
% Maintainer 11/19/2026: "The known storage for the beam is very simple. The
% problem we have is that oftentimes the complexity lies in lpi_ineq."
%
% INPUT
%   jobs:   cell array of rows {key, form, tier, mods, variant, param}:
%           key as in catalog.mat; form 'Q' or 'P'; tier of lpisettings_sop;
%           mods MATLAB statements on st ('' for none); variant 'A', 'A0' or
%           'AmT'; param the library parameter ([] for the default)
%   tag:    results_<TAG>.csv in this folder; a job already in it is skipped
%
% NOTES
%   - Size gate 2e6 variables, m 45,000 (gate/lpisolve.m). No job starts
%     after ST_STOPMIN minutes; the queue runner reruns the tag to resume.
% MMP, 11/19/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
if nargin < 3,  ndvcap = 2e6;  end     % optional: a larger variable cap for a known size
setappdata(0,'LB_NDV_CAP',ndvcap);   setappdata(0,'LB_M_CAP',45000);
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
stopmin = 15;
if exist(fullfile(here,'ST_STOPMIN'),'file'),  stopmin = str2double(strtrim(fileread(fullfile(here,'ST_STOPMIN'))));  end
tstop = datetime('now') + minutes(stopmin);
of = fullfile(here,['results_',tag,'.csv']);
done = {};
if exist(of,'file')
    L = splitlines(strtrim(fileread(of)));
    for i = 2:numel(L),  x = strsplit(L{i},',');  done{end+1} = strjoin(x([1 2 3 5]),','); end %#ok<AGROW>
else
    fid = fopen(of,'w');
    fprintf(fid,'key,form,variant,tier,mods,status,fr,pinf,dinf,numerr,ndec,m,sec,note\n');  fclose(fid);
end
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
for i = 1:numel(jobs)
    [key,form,tier,mods,variant,v] = jobs{i}{1:6};
    stt = '';  if numel(jobs{i}) > 6,  stt = jobs{i}{7};  end  % setter override, e.g. 'lam=%g;'
    id = strjoin({key,form,[variant,sprintf(stt,v)],['"',mods,'"']},',');   % mods: no commas
    if any(strcmp(done,id)),  continue,  end
    if datetime('now') >= tstop
        fprintf('\nlb_run: stop time %s reached; rerun to resume.\n',char(tstop));  return
    end
    fprintf('\n##### %s %s %s  [%s %s]  %s\n',key,form,variant,tier,mods,char(datetime('now')));
    c = Cat(find(strcmp({Cat.key},key),1));
    if ~isempty(stt),  c.setter = stt;  end
    note = '';  r = struct();
    try
        PIE = lb_build_system(c, v, true);
        switch variant                  % the generator under test
            case 'A0',   PIE.A = 0*PIE.A;
            case 'AmT',  PIE.A = -PIE.T;
            case 'AsT',  PIE.A = PIE.A - PIE.T;
            otherwise    % 'AsT<c>': A - c*T, the shift c of a sweep
                if startsWith(variant,'AsT'),  PIE.A = PIE.A - str2double(variant(4:end))*PIE.T;  end
        end
        S = lpisettings_sop(tier,0,'','mosek');
        if ~isempty(mods),  eval(strrep(mods,'st.','S.'));  end
        if strcmp(form,'P')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
        elseif any(strcmp(form,{'Qd','Pd'}))       % the dual test in that form
            setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
        elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
        end
        r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, S);
    catch ME
        note = strrep(strrep(ME.message,',',';'),newline,' ');  r.status = 'error';
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%s,%s,%s,"%s",%s,%.4g,%g,%g,%g,%.0f,%.0f,%.1f,"%s"\n',key,form,[variant,sprintf(stt,v)],tier, ...
        mods,f(r,'status'),g(r,'fr'),g(r,'pinf'),g(r,'dinf'),g(r,'numerr'),g(r,'ndec'),g(r,'m'), ...
        g(r,'sec'),note);
    fclose(fid);
    fprintf('##### %s %s %s: %s fr %.3g ndec %g\n',key,form,variant,f(r,'status'),g(r,'fr'),g(r,'ndec'));
end
fprintf('## DIAG_AMT %s DONE %s\n',tag,char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

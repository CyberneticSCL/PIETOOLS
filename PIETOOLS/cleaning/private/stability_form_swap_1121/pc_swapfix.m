function pc_swapfix(part,mins)
% PC_SWAPFIX(PART,MINS) tests the change of 11/21/2026 to PIETOOLS_stability_sop:
% 'Qd' and 'FQ' posed directly in two or more variables (before: replaced
% by 'Pd' and 'FP' without a message), and the warning when PIE.K turns
% 'Pd' or 'Qd' into 'P' or 'Q'.
%
% INPUT
%   part:  'b0' - one-variable rows ('Qd', 'FQ') on the executive before the change;
%          'b1' - the same rows after the change (the programs must be identical);
%          'a'  - rows in two and three variables after the change, each new
%                 'Qd'/'FQ' row beside the 'Pd'/'FP' row it replaced; a 'Qd'/'FQ'
%                 row that passes the gate is rerun at an unstable parameter;
%          'k'  - periodic heat (PIE.K nonempty) in 'Pd' and 'Qd': the warning,
%                 the form posed and the verdict.
%   mins:  (optional) minutes after which no row starts (default 18).
% results_SWAP_<PART>.csv in this folder. Maintainer 11/21/2026: the silent
% substitution "sounds like a bug to fix".
% MMP, 11/21/2026.
if nargin<2,  mins = 18;  end
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
fprintf('## PC_SWAPFIX %s %s\n## executive: %s\n',part,char(datetime('now')),which('PIETOOLS_stability_sop'));
tstop = datetime('now') + minutes(mins);
of = fullfile(here,['results_SWAP_',part,'.csv']);
fid = fopen(of,'w');
fprintf(fid,'key,form,tier,mods,param,status,fr,pinf,dinf,numerr,exst,eta,ndec,m,sec,gate,note\n');  fclose(fid);
if strcmp(part,'k'),  run_k(of);  return,  end
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
ts = 1.034726;                      % delay threshold of PDE 43 and 56 (library)
e1 = 'st.epneg = 1';                % 'FQ', 'FP': homogeneous in (P, epneg)
d0 = 'st.lpi.ineq.dD = 0';
s40 = 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0';                         % st of case 40
s43 = 'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''; st.sop.deg_Q{2} = [1 2 2]';  % st of case 43
% Rows {key, form, tier, mods, setter, parameter ([] library), unstable parameter}
switch part
    case {'b0','b1'}
        J = {'PDE:5','Qd','light','','',[],[];
             'PDE:5','FQ','light',e1,'',[],[];
             'PDE:8','Qd','light','','',[],[]};
    case 'a'
        J = {'PDE:40','Q','light',s40,'',[],[];
             'PDE:40','Qd','light',s40,'',[],[];
             'PDE:40','Pd','light',s40,'',[],[];
             'PDE:40','FQ','light',[s40,'; ',e1],'',[],[];
             'PDE:40','FP','light',[s40,'; ',e1],'',[],[];
             'PDE:43','Qd','light',s43,'tau=%.10g;',[],1.1*ts;
             'PDE:43','FQ','light',[s43,'; ',e1],'tau=%.10g;',[],1.1*ts;
             'PDE:56','FQ','stripped',[d0,'; ',e1],'tau=%.10g;',[],1.1*ts;
             'PDE:55','FQ','stripped',[d0,'; ',e1],'lam=%g;',[],17};
    case 'a2'   % 'a' without case 55 (its 'FQ' build, 29.8M variables, stops at the size gate)
        J = {'PDE:40','Q','light',s40,'',[],[];
             'PDE:40','Qd','light',s40,'',[],[];
             'PDE:40','Pd','light',s40,'',[],[];
             'PDE:40','FQ','light',[s40,'; ',e1],'',[],[];
             'PDE:40','FP','light',[s40,'; ',e1],'',[],[];
             'PDE:43','Qd','light',s43,'tau=%.10g;',[],1.1*ts;
             'PDE:43','FQ','light',[s43,'; ',e1],'tau=%.10g;',[],1.1*ts;
             'PDE:56','FQ','stripped',[d0,'; ',e1],'tau=%.10g;',[],1.1*ts};
    case 'c'    % case 40 at K = -0.5 (params set K; k only sets its default): A + B*K = 0.5*I, unstable
        J = {'PDE:40','Qd','light',s40,'K=%g;',-0.5,[];
             'PDE:40','FQ','light',[s40,'; ',e1],'K=%g;',-0.5,[];
             'PDE:40','Q','light',s40,'K=%g;',-0.5,[]};
    otherwise
        error('pc_swapfix:part','Unknown part %s.',part);
end
for i = 1:size(J,1)
    [key,form,tier,mods,stt,v,vu] = J{i,:};
    ok = false;
    for pass = 1:2              % 2: the unstable parameter, after a row that passes
        if pass==2 && (~ok || isempty(vu) || ~any(strcmp(form,{'Qd','FQ'}))),  break,  end
        if datetime('now') >= tstop
            fprintf('\n## PC_SWAPFIX stop time %s reached\n',char(tstop));  return
        end
        p = v;  if pass==2,  p = vu;  end
        fprintf('\n##### SWAP %s %s %s [%s] param %s  %s\n',key,form,tier,mods,num2str(p),char(datetime('now')));
        note = '';  r = struct();  t0 = tic;
        try
            c = Cat(find(strcmp({Cat.key},key),1));
            if isempty(c),  error('pc_swapfix:key','No catalog row %s.',key);  end
            if ~isempty(stt),  c.setter = stt;  end
            PIE = lb_build_system(c, p, true);
            S = lpisettings_sop(tier, double(any(strcmp(form,{'FP','FQ'}))), '', 'mosek');
            if ~isempty(mods),  eval([strrep(mods,'st.','S.'),';']);  end
            setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
            r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, S);
        catch ME
            note = strrep(strrep(ME.message,',',';'),newline,' ');  r.status = 'error';
        end
        r.sec = toc(t0);
        % The gate: harness verdict, executive optimal, eta <= 1e-7
        ok = strcmp(f(r,'status'),'certified') && strcmp(f(r,'exst'),'optimal') && g(r,'eta') <= 1e-7;
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,%s,"%s",%s,%s,%.4g,%g,%g,%g,%s,%.3g,%.0f,%.0f,%.1f,%d,"%s"\n',key,form,tier, ...
            strrep(mods,'''',''),num2str(p),f(r,'status'),g(r,'fr'),g(r,'pinf'),g(r,'dinf'),g(r,'numerr'), ...
            f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'),ok,note);
        fclose(fid);
        fprintf('##### SWAP %s %s %s param %s: %s fr %.3g exst %s eta %.3g ndec %g m %g gate %d\n',key,form,tier, ...
            num2str(p),f(r,'status'),g(r,'fr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),ok);
    end
end
fprintf('\n## PC_SWAPFIX %s DONE %s\n',part,char(datetime('now')));
end

function run_k(of)
% RUN_K(OF) poses the periodic heat equation x_t = x_ss - x on [0,1],
% x(0) = x(1), x_s(0) = x_s(1) (stable; PIE.K nonempty) in 'Pd' and 'Qd',
% and appends the warning, the form posed and the verdict to the file OF.
s = pvar('s');  th = pvar('theta');
PDE = pde_struct();  PDE.vars = [s, th];  PDE.dom = [0 1];
PDE.x{1} = struct('vars',[s th],'dom',[0 1],'size',1,'diff',2, ...
    'term',{{struct('x',1,'D',2,'C',1), struct('x',1,'C',-1)}});
B = [1,0,-1,0; 0,1,0,-1];           % acts on [x(0); x_s(0); x(1); x_s(1)]
for i = 1:2
    terms = {};
    for block = 1:4
        if B(i,block)~=0
            terms{end+1} = struct('x',1,'C',B(i,block),'D',mod(block-1,2),'loc',PDE.dom(1+(block>2))); %#ok<AGROW>
        end
    end
    PDE.BC{i} = struct('vars',[],'size',1,'term',{terms});
end
PDE = initialize(PDE,true);
PIE = convert_PIETOOLS_PDE(PDE,'silent');
for form = {'Pd','Qd'}
    lastwarn('');  t0 = tic;  stt = 'error';
    try
        S = lpisettings_sop('light');
        [~,~,info] = PIETOOLS_stability_sop(PIE,S,form{1});
        [msg,id] = lastwarn;
        note = sprintf('posed %s; K rows %d; warning [%s] %s',info.form,size(PIE.K,1),id,msg);
        stt = info.status;
    catch ME
        note = ME.message;
    end
    fid = fopen(of,'a');
    fprintf(fid,'periodic_heat,%s,light,"",,%s,,,,,,,,,%.1f,,"%s"\n',form{1},stt,toc(t0),strrep(note,'"',''''));
    fclose(fid);
    fprintf('##### SWAP periodic heat %s: %s; %s\n',form{1},stt,note);
end
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

function pc_fe(part)
% PC_FE(PART) tests the open problem cases with the finite-energy forms of
% PIETOOLS_stability_sop (maintainer 11/20/2026: "maybe test if there is a
% certificate using finite energy stability test?"):
%   'FP'  P >= 0,  T'*P*A + A'*P*T + 2*epneg*T'*T <= 0;
%   'FQ'  (1-D only) P >= 0,  T'*Q = P,  A'*Q + Q'*A + 2*epneg*T'*T <= 0;
% and the dual H2 test of PDE 50 without the margin of W ('cco' at
% eppos 0): the Gramian of a finite-rank B1 is compact, which no
% W >= eppos I with eppos > 0 contains; with W >= 0 the bound needs the
% stability certificate of the primal test.
% epneg = 1 in 'FP', 'FQ': both are homogeneous in (P, epneg), so
% feasibility does not depend on the value of epneg > 0, while a small
% epneg leaves P = 0 within the solver tolerance. A row that passes the gate
% is rerun at an unstable parameter.
%
% INPUT
%   part:   'a' - PDE 43 (the 1-D twin of 56) in 'FQ' and 'FP', then 56 in 'FP';
%           'b' - PDE 50 dual H2 in 'cco' at eppos 0, then PDE 55 in 'FP'.
% results_FE_<PART>.csv in this folder; no job starts after 18 minutes.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
fprintf('## PC_FE %s %s\n## lpisettings_sop: %s\n',part,char(datetime('now')),which('lpisettings_sop'));
tstop = datetime('now') + minutes(18);
of = fullfile(here,['results_FE_',part,'.csv']);
fid = fopen(of,'w');
fprintf(fid,'key,form,tier,mods,param,status,fr,pinf,dinf,numerr,exst,eta,ndec,m,sec,gam,gate,note\n');  fclose(fid);
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
ts = 1.034726;                  % delay threshold of PDE 43 and 56 (library)
d0 = 'st.lpi.ineq.dD = 0';
e0 = 'st.lpi.ineq.dD = 0; st.eppos = 0';
% Rows {key, form, tier, mods, setter, parameter ([] library), unstable parameter}
if strcmp(part,'a')
    J = {'PDE:43','FQ','light','', 'tau=%.10g;',[],1.1*ts;
         'PDE:43','FQ','light',d0, 'tau=%.10g;',[],1.1*ts;
         'PDE:43','FP','light','', 'tau=%.10g;',[],1.1*ts;
         'PDE:43','FP','light',d0, 'tau=%.10g;',[],1.1*ts;
         'PDE:43','FQ','heavy','', 'tau=%.10g;',[],1.1*ts;
         'PDE:43','FP','heavy','', 'tau=%.10g;',[],1.1*ts;
         'PDE:56','FP','stripped',d0,'tau=%.10g;',[],1.1*ts;
         'PDE:56','FP','stripped',d0,'tau=%.10g;',0.5*ts,1.1*ts};
elseif strcmp(part,'c')     % the 50 rows of 'b' again (exec prefix fixed)
    J = {'PDE:50/H2_norm_dual','cco','stripped',[e0,'; st.lpi.ineq.cap_inc = ''auto'''],'',[],[];
         'PDE:50/H2_norm_dual','cco','stripped',e0,'',[],[];
         'PDE:50/H2_norm_dual','cco','light',e0,'',[],[]};
else
    J = {'PDE:50/H2_norm_dual','cco','stripped',[e0,'; st.lpi.ineq.cap_inc = ''auto'''],'',[],[];
         'PDE:50/H2_norm_dual','cco','stripped',e0,'',[],[];
         'PDE:55','FP','stripped',d0,'lam=%g;',[],17;
         'PDE:50/H2_norm_dual','cco','light',e0,'',[],[]};
end
for i = 1:size(J,1)
    [key,form,tier,mods,stt,v,vu] = J{i,:};
    ok = false;
    for pass = 1:2              % 2: the unstable parameter, after a row that passes
        if pass==2 && (~ok || isempty(vu)),  break,  end
        if datetime('now') >= tstop
            fprintf('\n## PC_FE stop time %s reached\n',char(tstop));  return
        end
        p = v;  if pass==2,  p = vu;  end
        fprintf('\n##### FE %s %s %s [%s] param %s  %s\n',key,form,tier,mods,num2str(p),char(datetime('now')));
        note = '';  r = struct();  t0 = tic;
        try
            c = Cat(find(strcmp({Cat.key},key),1));
            if isempty(c),  error('pc_fe:key','No catalog row %s.',key);  end
            if ~isempty(stt),  c.setter = stt;  end
            PIE = lb_build_system(c, p, true);
            S = lpisettings_sop(tier, double(any(strcmp(form,{'FP','FQ'}))), '', 'mosek');
            if ~isempty(mods),  eval([strrep(mods,'st.','S.'),';']);  end
            if strcmp(form,'cco')
                setappdata(0,'LB_FORM_MAP',struct('c','cco','only','PIETOOLS_H2_norm_sop'));
                r = lb_exec(['sop:',c.exec], PIE, S);   % the catalog exec without 'sop:', as in ex_search
            else
                setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
                r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, S);
            end
        catch ME
            note = strrep(strrep(ME.message,',',';'),newline,' ');  r.status = 'error';
        end
        r.sec = toc(t0);
        % The gate of the run: harness verdict, executive optimal, eta <= 1e-7
        ok = strcmp(f(r,'status'),'certified') && strcmp(f(r,'exst'),'optimal') && g(r,'eta') <= 1e-7;
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,%s,"%s",%s,%s,%.4g,%g,%g,%g,%s,%.3g,%.0f,%.0f,%.1f,%.8g,%d,"%s"\n',key,form,tier, ...
            strrep(mods,'''',''),num2str(p),f(r,'status'),g(r,'fr'),g(r,'pinf'),g(r,'dinf'),g(r,'numerr'), ...
            f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'),g(r,'gam'),ok,note);
        fclose(fid);
        fprintf('##### FE %s %s %s param %s: %s fr %.3g exst %s eta %.3g ndec %g gate %d\n',key,form,tier, ...
            num2str(p),f(r,'status'),g(r,'fr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),ok);
    end
end
fprintf('\n## PC_FE %s DONE %s\n',part,char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

function pc_diag70()
% PC_DIAG70() tests PDE 70 (2-D wave with an acoustic boundary at s2 = 0,
% Beale 1976; infeasible with a certificate at every setting of 11/19-11/20)
% with the velocity ps of order (2,2), as the velocity of PDE 38 (2-D
% undamped wave, certified in 'Q' light dD 0 at deg_Q{2} = [0 1 1]); the
% example declares ps of order (2,1). The second boundary condition of ps
% in s2 is the time derivative of ph_s2(s1,0) = -et with
% et_t = -et - de - ps(s1,0):
%   ps_s2(s1,0) = et + de + ps(s1,0),
% which every solution of the example with compatible data satisfies.
% Variant '70m', each in 'Q' light dD 0 with deg_Q{2} [0 1 1], 'Q' light
% dD 0 and 'Pd' light dD 0; results_DIAG70.csv in this folder. A
% diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG70.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]';
     'Q','st.lpi.ineq.dD = 0';
     'Pd','st.lpi.ineq.dD = 0'};
for v = {'70m'}
    try,  PIE = build70(v{1});
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v{1},strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG70 %s: build error %s\n',v{1},ME.message);  continue
    end
    for i = 1:size(J,1)
        [form,mods] = J{i,:};
        if strcmp(form,'P')         % the primal 'Q' of the executor posed in form 'P'
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
        elseif any(strcmp(form,{'Pd','Qd'}))    % the dual test in that form
            setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
        elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
        end
        st = lpisettings_sop('light',0,'','mosek');
        if ~isempty(mods),  eval([mods,';']);  end
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},form,mods,r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAG70 %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v{1},form,mods,r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG70 DONE %s\n',char(datetime('now')));
end

function PIE = build70(~)
pvar s1 s2 t
clear stateNameGenerator
ph = pde_var([s1;s2],[0,1;0,1],[2;2]);
ps = pde_var([s1;s2],[0,1;0,1],[2;2]);          % 2nd order in s1 and s2, as the velocity of PDE 38
de = pde_var(s1,[0,1],2);
et = pde_var(s1,[0,1],2);
PDE_t = [diff(ph,'t')==ps;
         diff(ps,'t')==diff(ph,s1,2)+diff(ph,s2,2);
         diff(de,'t')==et;
         diff(et,'t')==-et-de-subs(ps,s2,0);
         subs([ph;ps],s1,0)==0;     subs([ph;ps],s1,1)==0;
         subs(diff(ph,s2),s2,0)==-et;
         subs(diff(ps,s2),s2,0)==et+de+subs(ps,s2,0);   % d/dt of the line above
         subs(ph,s2,1)==0;          subs(ps,s2,1)==0;
         subs([de;et],s1,0)==0;     subs([de;et],s1,1)==0];
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

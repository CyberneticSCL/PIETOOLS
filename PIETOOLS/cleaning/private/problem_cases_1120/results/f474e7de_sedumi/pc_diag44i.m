function pc_diag44i()
% PC_DIAG44I() heavier settings for the PDE 44 certificates past mu = 0.49 (k = 1, tau = 1; exact boundary mu =
% 2/3, REPORT finding 8): '44min' in 'Pd' and 'riemann' in 'Pd' and 'P' at heavy and veryheavy, mu = 0.55, 0.6,
% 0.65. results_DIAG44I.csv. A diagnostic. MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG44I.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,mu,form,tier,status,fr,numerr,exst,eta,ndec\n');  fclose(fid);
J = {};
for mu = [0.55 0.6 0.65]
    for tier = {'heavy','veryheavy'}
        J(end+1,:) = {'44min',mu,'Pd',tier{1}};  J(end+1,:) = {'riemann',mu,'Pd',tier{1}};  J(end+1,:) = {'riemann',mu,'P',tier{1}};
    end
end
for i = 1:size(J,1)
    [v,mu,form,tier] = J{i,:};
    if strcmp(form,'P'),  setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
    else,                 setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
    end
    PIE = build44i(v,mu);
    st = lpisettings_sop(tier,0,'','mosek');
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%g,%s,%s,%s,%.4g,%g,%s,%.3g,%.0f\n',v,mu,form,tier,r.status,g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'));
    fclose(fid);
    fprintf('##### DIAG44I %s mu %g %s %s: %s fr %.3g ndec %g\n',v,mu,form,tier,r.status,g(r,'fr'),g(r,'ndec'));
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG44I DONE %s\n',char(datetime('now')));
end

function PIE = build44i(v,mu)
k = 1;  tau = 1;
pvar s t
clear stateNameGenerator
if strcmp(v,'riemann')
    p = pde_var(s,[0,1],1);  q = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
    PDE_t = [diff(p,'t')==diff(q,s);   diff(q,'t')==diff(p,s);   diff(z,'t')==-(1/tau)*diff(z,s);
             subs(q,s,0)==0;   subs(p,s,1)==-k*(1-mu)*subs(q,s,1)-k*mu*subs(z,s,1);   subs(z,s,0)==subs(q,s,1)];
else
    u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
    PDE_t = [diff(u,'t')==v1;   diff(v1,'t')==diff(u,s,2);   diff(z,'t')==-(1/tau)*diff(z,s);
             subs(u,s,0)==0;  subs(v1,s,0)==0;
             subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,s,1);   subs(z,s,0)==subs(v1,s,1)];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

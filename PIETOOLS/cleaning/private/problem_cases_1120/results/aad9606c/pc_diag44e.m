function pc_diag44e()
% PC_DIAG44E() follows DIAG44D (11/20/2026): the minimal form of PDE 44 –
% u (order 2), v = u_t (order 1) and the delayed velocity z = u_t(1,t-tau s)
% (order 1), without the ODE state and the delayed position of the example
% file – certifies in the dual form 'Pd' at epneg 0 (2,520 variables,
% measure 1.000) and misses the gates at epneg 1e-3 (measure 1.316, eta
% 2.05e-7). Rows {variant, mu, form, tier, mods}:
%   '44min'  that system on s in [0,1] (k = 1, tau = 1): 'Pd' with epneg
%            1e-4 and 1e-2 at light, at stripped with epneg 0 and 1e-3, at
%            heavy with epneg 1e-3
%   'x3min'  the same system with z on a second variable s2 (the layout of
%            X3, PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D) at mu = 0.375
%            (the catalog parameter of X3) and 0.4: 'Pd' light at epneg 0
%            and 1e-3, 'Pd' stripped
% results_DIAG44E.csv in this folder. A diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG44E.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,mu,form,tier,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'44min',0.4,'Pd','light','st.epneg = 1e-4';
     '44min',0.4,'Pd','light','st.epneg = 1e-2';
     '44min',0.4,'Pd','stripped','';
     '44min',0.4,'Pd','stripped','st.epneg = 1e-3';
     '44min',0.4,'Pd','heavy','st.epneg = 1e-3';
     'x3min',0.375,'Pd','light','';
     'x3min',0.375,'Pd','light','st.epneg = 1e-3';
     'x3min',0.375,'Pd','stripped','';
     'x3min',0.4,'Pd','light','';
     'x3min',0.4,'Pd','light','st.epneg = 1e-3'};
for i = 1:size(J,1)
    [v,mu,form,tier,mods] = J{i,:};
    try,  PIE = build44e(v,mu);
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,%g,%s,%s,"%s","build error: %s"\n',v,mu,form,tier,mods,strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG44E %s: build error %s\n',v,ME.message);  continue
    end
    setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));   % the dual test in that form
    st = lpisettings_sop(tier,0,'','mosek');
    if ~isempty(mods),  eval([mods,';']);  end
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%g,%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v,mu,form,tier,mods,r.status, ...
        g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAG44E %s mu %g %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v,mu,form,tier,mods,r.status, ...
        g(r,'fr'),f(r,'exst'),g(r,'ndec'));
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG44E DONE %s\n',char(datetime('now')));
end

function PIE = build44e(v,mu)
k = 1;  tau = 1;
pvar s s2 t
clear stateNameGenerator
u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);
if strcmp(v,'44min'),  z = pde_var(s,[0,1],1);   zs = s;     % delay on the same variable (example 44)
else,                  z = pde_var(s2,[0,1],1);  zs = s2;    % delay on a second variable (X3)
end
PDE_t = [diff(u,'t')==v1;
         diff(v1,'t')==diff(u,s,2);
         diff(z,'t')==-(1/tau)*diff(z,zs);
         subs(u,s,0)==0;  subs(v1,s,0)==0;
         subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,zs,1);
         subs(z,zs,0)==subs(v1,s,1)];
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

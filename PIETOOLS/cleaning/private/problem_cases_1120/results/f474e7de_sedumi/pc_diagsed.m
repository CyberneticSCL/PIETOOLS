function pc_diagsed()
% PC_DIAGSED() confirms the new certificates of the problem-case run with SeDuMi 1.3 in place of MOSEK: the minimal
% forms of 44 (mu = 0.4, 'Pd' light, epneg 0 and 1e-1) and X3 (mu = 0.4, 'Pd' light, epneg 0 and 1e-2), then the
% mixed form of 70 ('Q' light dD 0, 187,116 variables with MOSEK). results_DIAGSED.csv. A diagnostic. MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
addpath(genpath('C:/YALMIP/yalmip_install/tbxmanager/toolboxes/sedumi/1.3/pcwin64/sedumi_1_3_pcwin64'));  % SeDuMi 1.3 (lb_setup resets the path)
fprintf('## sedumi: %s\n',which('sedumi'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAGSED.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'44min','Pd',''; '44min','Pd','st.epneg = 1e-1'; 'x3min','Pd',''; 'x3min','Pd','st.epneg = 1e-2'; '70mix','Q','st.lpi.ineq.dD = 0'};
for i = 1:size(J,1)
    [v,form,mods] = J{i,:};
    if strcmp(form,'Pd'),  setappdata(0,'LB_FORM_MAP',struct('Q','Pd','only',{{'PIETOOLS_stability_sop'}}));
    elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
    end
    PIE = buildsed(v);
    st = lpisettings_sop('light',0,'','sedumi');
    st.sos_opts.solver = 'sedumi';
    if ~isempty(mods),  eval([mods,';']);  end
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v,form,mods,r.status,g(r,'fr'),g(r,'numerr'),f(r,'exst'), ...
        g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAGSED %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v,form,mods,r.status,g(r,'fr'),f(r,'exst'),g(r,'ndec'));
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAGSED DONE %s\n',char(datetime('now')));
end

function PIE = buildsed(v)
k = 1;  tau = 1;  mu = 0.4;
pvar s s2 s1 t
clear stateNameGenerator
switch v
    case {'44min','x3min'}
        u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);
        if strcmp(v,'44min'),  z = pde_var(s,[0,1],1);   zs = s;
        else,                  z = pde_var(s2,[0,1],1);  zs = s2;
        end
        PDE_t = [diff(u,'t')==v1;   diff(v1,'t')==diff(u,s,2);   diff(z,'t')==-(1/tau)*diff(z,zs);
                 subs(u,s,0)==0;  subs(v1,s,0)==0;
                 subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,zs,1);   subs(z,zs,0)==subs(v1,s,1)];
    case '70mix'
        ph = pde_var([s1;s2],[0,1;0,1],[2;0]);   e2 = pde_var([s1;s2],[0,1;0,1],[0;1]);
        ps = pde_var([s1;s2],[0,1;0,1],[2;1]);   de = pde_var(s1,[0,1],2);   et = pde_var(s1,[0,1],2);
        PDE_t = [diff(ph,'t')==ps;   diff(ps,'t')==diff(ph,s1,2)+diff(e2,s2);   diff(e2,'t')==diff(ps,s2);
                 diff(de,'t')==et;   diff(et,'t')==-et-de-subs(ps,s2,0);
                 subs([ph;ps],s1,0)==0;   subs([ph;ps],s1,1)==0;   subs(ps,s2,1)==0;   subs(e2,s2,0)==-et;
                 subs([de;et],s1,0)==0;   subs([de;et],s1,1)==0];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

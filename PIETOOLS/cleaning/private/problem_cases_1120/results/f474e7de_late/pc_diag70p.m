function pc_diag70p()
% PC_DIAG70P() the mixed position-strain form of PDE 70 ('70mix', certified in 'Q' light dD 0) across the edge
% oscillator parameters: et_t = -d*et - kk*de - ps(s1,0), de_t = kk*et (energy (|ps|^2 + |ph_s1|^2 + |e2|^2 +
% kk|de|^2... scaled so that dE/dt = -d*|et|^2), for (d, kk) = (0.1, 1), (10, 1), (1, 10), (1, 1) (the example).
% A robustness check of the certificate. results_DIAG70P.csv. A diagnostic. MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
of = fullfile(here,'results_DIAG70P.csv');
fid = fopen(of,'w');  fprintf(fid,'d,kk,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
for P = [0.1 1; 10 1; 1 10; 1 1]'
    d = P(1);  kk = P(2);
    pvar s1 s2 t
    clear stateNameGenerator
    ph = pde_var([s1;s2],[0,1;0,1],[2;0]);
    e2 = pde_var([s1;s2],[0,1;0,1],[0;1]);
    ps = pde_var([s1;s2],[0,1;0,1],[2;1]);
    de = pde_var(s1,[0,1],2);
    et = pde_var(s1,[0,1],2);
    PDE_t = [diff(ph,'t')==ps;
             diff(ps,'t')==diff(ph,s1,2)+diff(e2,s2);
             diff(e2,'t')==diff(ps,s2);
             diff(de,'t')==kk*et;
             diff(et,'t')==-d*et-kk*de-subs(ps,s2,0);
             subs([ph;ps],s1,0)==0;     subs([ph;ps],s1,1)==0;
             subs(ps,s2,1)==0;
             subs(e2,s2,0)==-et;
             subs([de;et],s1,0)==0;     subs([de;et],s1,1)==0];
    PIE = convert(PDE_t,'pie');
    st = lpisettings_sop('light',0,'','mosek');
    st.lpi.ineq.dD = 0;
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%g,%g,%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',d,kk,r.status,g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAG70P d %g kk %g: %s fr %.3g ndec %g\n',d,kk,r.status,g(r,'fr'),g(r,'ndec'));
end
fprintf('## DIAG70P DONE %s\n',char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

function pc_diag70s()
% PC_DIAG70S() a soundness test of the mixed position-strain form of PDE 70
% ('70mix', certified in 'Q' light dD 0 at 187,116 variables): with the
% damping of the edge oscillator reversed, et_t = +c*et - de - ps(s1,0),
% the energy grows, dE/dt = c*|et|^2 >= 0, and no non-increasing storage
% exists; the program must fail. Rows: c = 1 and c = 0.1 in 'Q' light dD 0;
% c = -1 (the example, control). results_DIAG70S.csv. A diagnostic.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
of = fullfile(here,'results_DIAG70S.csv');
fid = fopen(of,'w');  fprintf(fid,'c,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
for c = [1 0.1 -1]
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
             diff(de,'t')==et;
             diff(et,'t')==c*et-de-subs(ps,s2,0);
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
    fprintf(fid,'%g,%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',c,r.status,g(r,'fr'),g(r,'numerr'),f(r,'exst'), ...
        g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAG70S c %g: %s fr %.3g ndec %g\n',c,r.status,g(r,'fr'),g(r,'ndec'));
end
fprintf('## DIAG70S DONE %s\n',char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

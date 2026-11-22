function pc_diag69s()
% PC_DIAG69S() a soundness test of the mixed position-strain form of PDE 69
% ('69mix2', certified in 'Q' light dD 0 at 267,060 variables): with a
% reaction term in the heat equation, w_t = w_s1s1 + w_s2s2 + r*w, r = 30
% (above 2*pi^2 = 19.74, the decay rate of the heat part with Dirichlet data
% on [0,1]^2, and the heat part here has Dirichlet data on three sides
% only), the heat part grows and the coupled system has no non-increasing
% storage; the program must fail. Control: r = 0 (the example).
% results_DIAG69S.csv. A diagnostic. MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
of = fullfile(here,'results_DIAG69S.csv');
fid = fopen(of,'w');  fprintf(fid,'r,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
for r0 = [30 0]
    pvar s1 s2 t
    clear stateNameGenerator
    u = pde_var([s1;s2],[0,1;0,1],[0;2]);
    e1 = pde_var([s1;s2],[0,1;0,1],[1;2]);
    vv = pde_var([s1;s2],[0,1;0,1],[1;2]);
    w = pde_var([s1;s2],[0,1;0,1],[2;2]);
    PDE_t = [diff(u,'t')==vv;
             diff(vv,'t')==diff(e1,s1)+diff(u,s2,2);
             diff(e1,'t')==diff(vv,s1);
             diff(w,'t')==diff(w,s1,2)+diff(w,s2,2)+r0*w;
             subs(u,s2,0)==0;           subs(u,s2,1)==0;
             subs(e1,s1,0)==-subs(diff(w,s1),s1,0);
             subs(e1,s2,0)==0;          subs(e1,s2,1)==0;
             subs(vv,s1,1)==0;          subs(vv,s2,0)==0;    subs(vv,s2,1)==0;
             subs(w,s1,1)==0;           subs(w,s1,0)==subs(vv,s1,0);
             subs(w,s2,0)==0;           subs(w,s2,1)==0];
    PIE = convert(PDE_t,'pie');
    st = lpisettings_sop('light',0,'','mosek');
    st.lpi.ineq.dD = 0;
    try,  res = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  res = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%g,%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',r0,res.status,g(res,'fr'),g(res,'numerr'),f(res,'exst'), ...
        g(res,'eta'),g(res,'ndec'),g(res,'m'),g(res,'sec'));
    fclose(fid);
    fprintf('##### DIAG69S r %g: %s fr %.3g ndec %g\n',r0,res.status,g(res,'fr'),g(res,'ndec'));
end
fprintf('## DIAG69S DONE %s\n',char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

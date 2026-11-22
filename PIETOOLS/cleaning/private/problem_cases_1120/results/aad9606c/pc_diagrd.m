function pc_diagrd()
% PC_DIAGRD() follows a side result of DIAGPD (11/20/2026): the 2-D
% reaction-diffusion equation u_t = u_s1s1 + u_s2s2 + r*u on [0,1]^2, u = 0
% on every edge, is not certified at r = 15 (76% of the threshold 2*pi^2 =
% 19.74) in 'P' or 'Pd' at light dD 0 (5,745 / 37,745 variables). Rows
% {r, form, tier, mods}: r = 15 in 'P' light (dD 1), 'P' heavy, 'Q' light dD
% 0, 'Q' light; r = 10 (51%) and 5 (25%) in 'P' light dD 0.
% results_DIAGRD.csv in this folder. A diagnostic.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAGRD.csv');
fid = fopen(of,'w');  fprintf(fid,'r,form,tier,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {15,'P','light','';  15,'P','heavy','';  15,'Q','light','st.lpi.ineq.dD = 0';  15,'Q','light','';
     10,'P','light','st.lpi.ineq.dD = 0';  5,'P','light','st.lpi.ineq.dD = 0'};
for i = 1:size(J,1)
    [r,form,tier,mods] = J{i,:};
    pvar s1 s2 t
    clear stateNameGenerator
    u = pde_var([s1;s2],[0,1;0,1],[2;2]);
    PDE_t = [diff(u,'t')==diff(u,s1,2)+diff(u,s2,2)+r*u;
             subs(u,s1,0)==0;   subs(u,s1,1)==0;   subs(u,s2,0)==0;   subs(u,s2,1)==0];
    PIE = convert(PDE_t,'pie');
    if strcmp(form,'P')
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
    elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
    end
    st = lpisettings_sop(tier,0,'','mosek');
    if ~isempty(mods),  eval([mods,';']);  end
    try,  res = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  res = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%g,%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',r,form,tier,mods,res.status, ...
        g(res,'fr'),g(res,'numerr'),f(res,'exst'),g(res,'eta'),g(res,'ndec'),g(res,'m'),g(res,'sec'));
    fclose(fid);
    fprintf('##### DIAGRD r %g %s %s [%s]: %s fr %.3g ndec %g\n',r,form,tier,mods,res.status,g(res,'fr'),g(res,'ndec'));
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAGRD DONE %s\n',char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

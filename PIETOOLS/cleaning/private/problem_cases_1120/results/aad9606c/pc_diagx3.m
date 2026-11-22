function pc_diagx3()
% PC_DIAGX3() X3 (PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D: the system of
% example 44 with its delay states on a second variable s2) in the energy
% (Riemann) states w_s, w_t on s1 and the delayed velocity on s2, against
% the example as written; mu = 0.375 (the catalog parameter) and 0.4, k = 1,
% tau = 1; 'Q' and 'P' at light. PC_DIAG44 (11/20/2026): example 44 is
% infeasible in every form that keeps the position u as a state and
% certifies in the Riemann states. Results results_DIAGX3.csv in this folder.
% A diagnostic, not a setting. MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAGX3.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,mu,form,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
for mu = [0.375 0.4]
    for v = {'example','riemann'}
        try,  PIE = buildx3(v{1},mu);
        catch ME
            fid = fopen(of,'a');  fprintf(fid,'%s,%g,,"build error: %s"\n',v{1},mu,strrep(ME.message,',',';'));  fclose(fid);  continue
        end
        for form = {'Q','P'}
            if strcmp(form{1},'P')
                setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
            elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
            end
            st = lpisettings_sop('light',0,'','mosek');
            try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
            catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
            end
            fid = fopen(of,'a');
            fprintf(fid,'%s,%g,%s,%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},mu,form{1},r.status, ...
                g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
            fclose(fid);
            fprintf('##### DIAGX3 %s mu %g %s: %s fr %.3g ndec %g\n',v{1},mu,form{1},r.status,g(r,'fr'),g(r,'ndec'));
        end
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAGX3 DONE %s\n',char(datetime('now')));
end

function PIE = buildx3(v,mu)
k = 1;  tau = 1;
if strcmp(v,'example')
    PDE_t = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D(0,{'k=1;','tau=1;',sprintf('mu=%g;',mu)});
    PIE = convert(PDE_t,'pie');   return
end
pvar s1 s2 t
clear stateNameGenerator
p = pde_var(s1,[0,1],1);  q = pde_var(s1,[0,1],1);  z = pde_var(s2,[0,1],1);
PDE_t = [diff(p,'t')==diff(q,s1);
         diff(q,'t')==diff(p,s1);
         diff(z,'t')==-(1/tau)*diff(z,s2);
         subs(q,s1,0)==0;
         subs(p,s1,1)==-k*(1-mu)*subs(q,s1,1)-k*mu*subs(z,s2,1);
         subs(z,s2,0)==subs(q,s1,1)];
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

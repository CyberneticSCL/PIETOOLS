function pc_diag41()
% PC_DIAG41() tests which part of PDE 41 (heat equation with a delayed
% Dirichlet boundary input, library example 41) keeps its stability LPI
% inconclusive. Every tested setting of 11/19-11/20 ends MOSEK 'unknown' with
% the dual objective settling above 0 while the primal residual vanishes:
% a weakly infeasible program. Variants of the system, the delay written as
% the transport x2_t = x2_s2 on s2 in [1, 1+tau]:
%   'as is'    x1(1) = x2(1)          (the example, delay state explicit)
%   'flux'     x1_s(1) = x2(1)        (the delayed input enters as a flux)
%   'decoupled' x1(1) = 0             (heat and transport not coupled)
%   'example'  the example file itself (delay expanded by the parser)
% each in the forms 'Q' and 'P' at light; results_DIAG41.csv in this folder.
% A diagnostic of the system, not a setting: the example is not changed.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG41.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
lam = 1;  tau = 1;
for v = {'as is','flux','decoupled','example'}
    PIE = build41(v{1},lam,tau);
    for form = {'Q','P'}
        if strcmp(form{1},'P')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
        elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
        end
        S = lpisettings_sop('light',0,'','mosek');
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, S);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},form{1},r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAG41 %s %s: %s fr %.3g exst %s eta %.3g ndec %g\n',v{1},form{1},r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'));
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG41 DONE %s\n',char(datetime('now')));
end

function PIE = build41(v,lam,tau)
% The PIE of variant V of PDE 41 at lam, tau.
if strcmp(v,'example')
    PDE_t = PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input(0,{sprintf('lam=%g;',lam),sprintf('tau=%g;',tau)});
    PIE = convert(PDE_t,'pie');   return
end
pvar s1 s2 t
clear stateNameGenerator
pde_var x1 x2 control u
x1.vars = s1;   x1.dom = [0,1];
x2.vars = s2;   x2.dom = [1,1+tau];
switch v
    case 'as is',      bc = subs(x1,s1,1)==subs(x2,s2,1);
    case 'flux',       bc = subs(diff(x1,s1),s1,1)==subs(x2,s2,1);
    case 'decoupled',  bc = subs(x1,s1,1)==0;
end
PDE_t = [diff(x1,t)==diff(x1,s1,2)+lam*x1;
         diff(x2,t)==diff(x2,s2);
         subs(x1,s1,0)==0;  bc;
         subs(x2,s2,1+tau)==u];
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

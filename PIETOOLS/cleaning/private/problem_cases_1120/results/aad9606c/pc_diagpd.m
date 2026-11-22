function pc_diagpd()
% PC_DIAGPD() tests the dual stability form 'Pd' of PIETOOLS_stability_sop
% against the primal 'P' on self-adjoint PIEs. For u_t = c*Lap(u) + r*u with
% u = 0 on the boundary, T (the Green operator of the fundamental state) and
% A are self-adjoint, so the dual condition T*P*A' + A*P*T' <= 0 and the
% primal condition A'*P*T + T'*P*A <= 0 coincide: the two forms must agree.
% In two variables 'Pd' failed on every system tried (PDE 38, 50, 56, 69,
% 70) where 'P' or 'Q' certified some. Variants:
%   'heat1d'  u_t = u_ss on [0,1], u(0) = u(1) = 0 (control)
%   'heat2d'  u_t = u_s1s1 + u_s2s2 on [0,1]^2, u = 0 on every edge
%   'rd2d'    as 'heat2d' with + 15*u (threshold 2*pi^2 = 19.74)
% each in 'P' and 'Pd' light dD 0; results_DIAGPD.csv in this folder. A
% diagnostic.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAGPD.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
for v = {'heat1d','heat2d','rd2d'}
    PIE = buildpd(v{1});
    for form = {'P','Pd'}
        if strcmp(form{1},'P')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
        else
            setappdata(0,'LB_FORM_MAP',struct('Q','Pd','only',{{'PIETOOLS_stability_sop'}}));
        end
        st = lpisettings_sop('light',0,'','mosek');
        st.lpi.ineq.dD = 0;
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,"st.lpi.ineq.dD = 0",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},form{1},r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAGPD %s %s: %s fr %.3g exst %s ndec %g\n',v{1},form{1},r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAGPD DONE %s\n',char(datetime('now')));
end

function PIE = buildpd(v)
pvar s s1 s2 t
clear stateNameGenerator
switch v
    case 'heat1d'
        u = pde_var(s,[0,1],2);
        PDE_t = [diff(u,'t')==diff(u,s,2);  subs(u,s,0)==0;  subs(u,s,1)==0];
    case {'heat2d','rd2d'}
        r = 0;   if strcmp(v,'rd2d'),  r = 15;  end
        u = pde_var([s1;s2],[0,1;0,1],[2;2]);
        PDE_t = [diff(u,'t')==diff(u,s1,2)+diff(u,s2,2)+r*u;
                 subs(u,s1,0)==0;   subs(u,s1,1)==0;   subs(u,s2,0)==0;   subs(u,s2,1)==0];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

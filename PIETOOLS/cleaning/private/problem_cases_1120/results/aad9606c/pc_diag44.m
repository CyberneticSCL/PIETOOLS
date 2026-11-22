% function pc_diag44()                                                      % MMP, 11/20/2026 (11:33, was)
function pc_diag44(J, tag)                                                  % MMP, 11/20/2026 (11:33)
% PC_DIAG44() tests which part of the formulation of PDE 44 (wave equation,
% partly delayed boundary damping, Xu 2006; library example 44) defeats the
% stability LPI: infeasible with a certificate at every setting of
% 11/19-11/20, while the system is exponentially stable at the default
% mu = 0.4 (spectral abscissa -0.29) and its Riemann form (example 59)
% certifies. Variants, k = 1, tau = 1, mu = 0.4:
%   'example'  the example file (states x1 ODE, x2 = u, x3 = u_t of order 2,
%              x4 = u(1,t-tau s), x5 = u_t(1,t-tau s))
%   'orders'   the same states, x2 of order 2 and x3 of order 1
%   'minimal'  u (order 2), u_t (order 1), u_t(1,t-tau s): no ODE state and
%              no delayed position
%   'riemann'  w_s, w_t (order 1) and w_t(1,t-tau s): the energy states
% each in 'Q' and 'P' at light, and 'Q' light with epneg 1e-3 or eppos 1e-2;
% results_DIAG44.csv in this folder. A diagnostic, not a setting.
% Optional: J, rows {form, mods} in place of those above (form 'Pd' or 'Qd'
% the dual test); tag, the results file results_<tag>.csv (default DIAG44).
% MMP, 11/20/2026.
% MMP, 11/20/2026 (11:33): inputs J and tag; the dual forms mapped as in
%                diag_amt (DIAG69: the 1-D position forms of 69 and 70
%                certify in 'Pd' only).
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if nargin < 2 || isempty(tag),  tag = 'DIAG44';  end                       % MMP, 11/20/2026 (11:33)
% of = fullfile(here,'results_DIAG44.csv');                                 % MMP, 11/20/2026 (11:33, was)
of = fullfile(here,['results_',tag,'.csv']);                                % MMP, 11/20/2026 (11:33)
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
% J = {'Q',''; 'P',''; 'Q','st.epneg = 1e-3'; 'Q','st.eppos = 1e-2'};       % MMP, 11/20/2026 (11:33, was)
if nargin < 1 || isempty(J)                                                 % MMP, 11/20/2026 (11:33)
    J = {'Q',''; 'P',''; 'Q','st.epneg = 1e-3'; 'Q','st.eppos = 1e-2'};     % MMP, 11/20/2026 (11:33)
end                                                                         % MMP, 11/20/2026 (11:33)
for v = {'example','orders','minimal','riemann'}
    try,  PIE = build44(v{1});
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v{1},strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG44 %s: build error %s\n',v{1},ME.message);  continue
    end
    for i = 1:size(J,1)
        [form,mods] = J{i,:};
        if strcmp(form,'P')
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
        elseif any(strcmp(form,{'Pd','Qd'}))    % the dual test in that form  % MMP, 11/20/2026 (11:33)
            setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}})); % MMP, 11/20/2026 (11:33)
        elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
        end
        st = lpisettings_sop('light',0,'','mosek');
        if ~isempty(mods),  eval([mods,';']);  end
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},form,mods,r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAG44 %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v{1},form,mods,r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG44 DONE %s\n',char(datetime('now')));
end

function PIE = build44(v)
k = 1;  tau = 1;  mu = 0.4;
if strcmp(v,'example')
    PDE_t = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay(0,{'k=1;','tau=1;','mu=0.4;'});
    PIE = convert(PDE_t,'pie');   return
end
pvar s t
clear stateNameGenerator
switch v
    case 'orders'       % the example states, x2 of order 2 and x3 of order 1
        x1 = pde_var();  x2 = pde_var(s,[0,1],2);  x3 = pde_var(s,[0,1],1);
        x4 = pde_var(s,[0,1]);  x5 = pde_var(s,[0,1]);
        PDE_t = [diff(x1,'t')==subs(diff(x2,s),s,1);
                 diff(x2,'t')==x3;
                 diff(x3,'t')==diff(x2,s,2);
                 diff(x4,'t')==-(1/tau)*diff(x4,s);
                 diff(x5,'t')==-(1/tau)*diff(x5,s);
                 subs([x2;x3],s,0)==0;
                 x1==-k*(1-mu)*subs(x2,s,1)-k*mu*subs(x4,s,1);
                 subs(diff(x2,s),s,1)==-k*(1-mu)*subs(x3,s,1)-k*mu*subs(x5,s,1);
                 subs(x4,s,0)==subs(x2,s,1);
                 subs(x5,s,0)==subs(x3,s,1)];
    case 'minimal'      % u, u_t and the delayed velocity only
        u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
        PDE_t = [diff(u,'t')==v1;
                 diff(v1,'t')==diff(u,s,2);
                 diff(z,'t')==-(1/tau)*diff(z,s);
                 subs(u,s,0)==0;  subs(v1,s,0)==0;
                 subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,s,1);
                 subs(z,s,0)==subs(v1,s,1)];
    case 'riemann'      % energy states w_s, w_t and the delayed velocity
        p = pde_var(s,[0,1],1);  q = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
        PDE_t = [diff(p,'t')==diff(q,s);
                 diff(q,'t')==diff(p,s);
                 diff(z,'t')==-(1/tau)*diff(z,s);
                 subs(q,s,0)==0;
                 subs(p,s,1)==-k*(1-mu)*subs(q,s,1)-k*mu*subs(z,s,1);
                 subs(z,s,0)==subs(q,s,1)];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

function pc_diag44f()
% PC_DIAG44F() follows DIAG44E (11/20/2026). The minimal form of PDE 44 (u,
% v = u_t, delayed velocity z; k = 1, tau = 1) certifies in 'Pd' at epneg
% 1e-2 (light) and 1e-3 (heavy); the X3 layout (z on s2) certifies at
% epneg 0, and at epneg 1e-3 passes the measure gate with eta 2-3e-7 (above
% 1e-7). Rows {variant, mu, form, tier, mods}:
%   'x3min' at mu 0.375 and 0.4: 'Pd' light epneg 1e-2, 'Pd' heavy epneg
%           1e-3 (the two exponential settings of 44min)
%   '44min' at mu 0.4: 'Pd' light epneg 3e-2 and 1e-1 (the decay rate the
%           dual certifies; spectral abscissa -0.289)
%   'dw1'   u_tt = u_ss, u(0) = 0, u_s(1) = -u_t(1) (states u of order 2,
%           v = u_t of order 1; exponentially stable): 'Q' light, 'Q' with
%           deg_Q{1} [8 6 6], 'P' light, 'Pd' light – does the direct form
%           certify a boundary-damped wave in the position state
%   '44min' 'Q' light with deg_Q{1} [8 6 6]
% results_DIAG44F.csv in this folder. A diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG44F.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,mu,form,tier,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'x3min',0.375,'Pd','light','st.epneg = 1e-2';
     'x3min',0.375,'Pd','heavy','st.epneg = 1e-3';
     'x3min',0.4,'Pd','light','st.epneg = 1e-2';
     'x3min',0.4,'Pd','heavy','st.epneg = 1e-3';
     '44min',0.4,'Pd','light','st.epneg = 3e-2';
     '44min',0.4,'Pd','light','st.epneg = 1e-1';
     'dw1',0,'Q','light','';
     'dw1',0,'Q','light','st.sop.deg_Q{1} = [8 6 6]';
     'dw1',0,'P','light','';
     'dw1',0,'Pd','light','';
     '44min',0.4,'Q','light','st.sop.deg_Q{1} = [8 6 6]'};
for i = 1:size(J,1)
    [v,mu,form,tier,mods] = J{i,:};
    try,  PIE = build44f(v,mu);
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,%g,%s,%s,"%s","build error: %s"\n',v,mu,form,tier,mods,strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG44F %s: build error %s\n',v,ME.message);  continue
    end
    if strcmp(form,'P')             % the primal 'Q' of the executor posed in form 'P'
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
    elseif any(strcmp(form,{'Pd','Qd'}))    % the dual test in that form
        setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
    elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
    end
    st = lpisettings_sop(tier,0,'','mosek');
    if ~isempty(mods),  eval([mods,';']);  end
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%g,%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v,mu,form,tier,mods,r.status, ...
        g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAG44F %s mu %g %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v,mu,form,tier,mods,r.status, ...
        g(r,'fr'),f(r,'exst'),g(r,'ndec'));
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG44F DONE %s\n',char(datetime('now')));
end

function PIE = build44f(v,mu)
k = 1;  tau = 1;
pvar s s2 t
clear stateNameGenerator
u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);
if strcmp(v,'dw1')                  % boundary-damped wave in the position state, no delay
    PDE_t = [diff(u,'t')==v1;
             diff(v1,'t')==diff(u,s,2);
             subs(u,s,0)==0;  subs(v1,s,0)==0;
             subs(diff(u,s),s,1)==-k*subs(v1,s,1)];
    PIE = convert(PDE_t,'pie');   return
end
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

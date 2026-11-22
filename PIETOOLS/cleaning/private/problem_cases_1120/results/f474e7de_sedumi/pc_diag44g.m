function pc_diag44g()
% PC_DIAG44G() sweeps the delay weight mu of PDE 44 (k = 1, tau = 1; stable
% whenever mu < 1/2, Xu 2006) across the stability boundary in the dual form
% 'Pd' at light, epneg 0: the minimal form ('44min': u, v = u_t, delayed v;
% certified at mu = 0.4 in DIAG44D-F) and the Riemann form ('riemann': w_s,
% w_t, delayed w_t), whose numerical witness (max Re of the Legendre
% pencil) is reliable for this system (-0.289 at mu = 0.4) and gives the
% reference; and the X3 layout ('x3min', the delay on s2) at mu 0.45, 0.55.
% A certificate past the boundary (mu > 1/2 with a positive reference) would
% make the minimal-form certificate unsound. results_DIAG44G.csv in this
% folder. A diagnostic.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
setappdata(0,'LB_FORM_MAP',struct('Q','Pd','only',{{'PIETOOLS_stability_sop'}}));   % the dual test
of = fullfile(here,'results_DIAG44G.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,mu,status,fr,numerr,exst,eta,ndec,m,sec,maxre\n');  fclose(fid);
J = {};
for mu = [0.45 0.48 0.49 0.55 0.6],  J(end+1,:) = {'44min',mu};  J(end+1,:) = {'riemann',mu};  end
J(end+1,:) = {'x3min',0.45};   J(end+1,:) = {'x3min',0.55};
for i = 1:size(J,1)
    [v,mu] = J{i,:};
    PIE = build44g(v,mu);
    st = lpisettings_sop('light',0,'','mosek');
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    mr = NaN;    % the witness of the executive, 1-D PIEs only
    try,  if strcmp(v,'riemann') || strcmp(v,'44min'),  W = pie_witness_sop(PIE,'rate');  if W.ok,  mr = W.maxre;  end,  end
    catch,  end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%g,%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f,%.6g\n',v,mu,r.status, ...
        g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'),mr);
    fclose(fid);
    fprintf('##### DIAG44G %s mu %g: %s fr %.3g ndec %g maxre %.4g\n',v,mu,r.status,g(r,'fr'),g(r,'ndec'),mr);
end
rmappdata(0,'LB_FORM_MAP');
fprintf('## DIAG44G DONE %s\n',char(datetime('now')));
end

function PIE = build44g(v,mu)
k = 1;  tau = 1;
pvar s s2 t
clear stateNameGenerator
switch v
    case 'riemann'      % energy states w_s, w_t and the delayed velocity (pc_diag44 'riemann')
        p = pde_var(s,[0,1],1);  q = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
        PDE_t = [diff(p,'t')==diff(q,s);
                 diff(q,'t')==diff(p,s);
                 diff(z,'t')==-(1/tau)*diff(z,s);
                 subs(q,s,0)==0;
                 subs(p,s,1)==-k*(1-mu)*subs(q,s,1)-k*mu*subs(z,s,1);
                 subs(z,s,0)==subs(q,s,1)];
    otherwise           % the minimal position form (pc_diag44e '44min', 'x3min')
        u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);
        if strcmp(v,'44min'),  z = pde_var(s,[0,1],1);   zs = s;
        else,                  z = pde_var(s2,[0,1],1);  zs = s2;
        end
        PDE_t = [diff(u,'t')==v1;
                 diff(v1,'t')==diff(u,s,2);
                 diff(z,'t')==-(1/tau)*diff(z,zs);
                 subs(u,s,0)==0;  subs(v1,s,0)==0;
                 subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,zs,1);
                 subs(z,zs,0)==subs(v1,s,1)];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

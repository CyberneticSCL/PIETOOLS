function pc_diag44h()
% PC_DIAG44H() continues DIAG44G: at tau = 1 the Riemann witness of PDE 44
% stays negative past mu = 1/2 (-0.117 at 0.55, -0.066 at 0.6), so the
% bound mu < 1/2 (for every delay) does not mark the boundary at tau = 1.
% For mu = 0.65, 0.7, 0.8, 0.9, 0.95: the Riemann witness at N = 24 and 48
% (Legendre-Galerkin pencil, reliable on this energy form) locates the
% boundary at tau = 1; '44min' and 'riemann' in 'Pd' light at epneg 0 test
% whether the dual certificate stops there. A '44min' certificate where the
% witness is positive would make it unsound. results_DIAG44H.csv. A
% diagnostic. MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
setappdata(0,'LB_FORM_MAP',struct('Q','Pd','only',{{'PIETOOLS_stability_sop'}}));   % the dual test
of = fullfile(here,'results_DIAG44H.csv');
fid = fopen(of,'w');  fprintf(fid,'mu,maxre_N24,maxre_N48,variant,status,fr,exst,eta,ndec\n');  fclose(fid);
for mu = [0.65 0.7 0.8 0.9 0.95]
    PR = build44h('riemann',mu);
    w = NaN(1,2);  NN = [24 48];
    for j = 1:2
        try,  W = pie_witness_sop(PR,'rate',struct('N_cheb',NN(j)));  if W.ok,  w(j) = W.maxre;  end
        catch,  end
    end
    for v = {'44min','riemann'}
        PIE = PR;   if strcmp(v{1},'44min'),  PIE = build44h('44min',mu);  end
        st = lpisettings_sop('light',0,'','mosek');
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%g,%.6g,%.6g,%s,%s,%.4g,%s,%.3g,%.0f\n',mu,w(1),w(2),v{1},r.status,g(r,'fr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'));
        fclose(fid);
        fprintf('##### DIAG44H mu %g (maxre %.4g / %.4g) %s: %s fr %.3g ndec %g\n',mu,w(1),w(2),v{1},r.status,g(r,'fr'),g(r,'ndec'));
    end
end
rmappdata(0,'LB_FORM_MAP');
fprintf('## DIAG44H DONE %s\n',char(datetime('now')));
end

function PIE = build44h(v,mu)
k = 1;  tau = 1;
pvar s t
clear stateNameGenerator
if strcmp(v,'riemann')
    p = pde_var(s,[0,1],1);  q = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
    PDE_t = [diff(p,'t')==diff(q,s);   diff(q,'t')==diff(p,s);   diff(z,'t')==-(1/tau)*diff(z,s);
             subs(q,s,0)==0;   subs(p,s,1)==-k*(1-mu)*subs(q,s,1)-k*mu*subs(z,s,1);   subs(z,s,0)==subs(q,s,1)];
else
    u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
    PDE_t = [diff(u,'t')==v1;   diff(v1,'t')==diff(u,s,2);   diff(z,'t')==-(1/tau)*diff(z,s);
             subs(u,s,0)==0;  subs(v1,s,0)==0;
             subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,s,1);   subs(z,s,0)==subs(v1,s,1)];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

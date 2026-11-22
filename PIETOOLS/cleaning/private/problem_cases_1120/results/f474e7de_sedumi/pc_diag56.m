function pc_diag56()
% PC_DIAG56() tests PDE 56 (2-D heat with an interior delay, delay as a 3-D
% transport state x2 of order (2,2,1); infeasible at every setting of
% 11/18-11/20) and its 1-D twin PDE 43 (x2 of order (2,1); certified in 'Q'
% light dD 0 with the cap at deg_Q{2} = [1 2 2], 132,609 variables) with the
% delay state of order 0 in the spatial variables: x2 transports in the
% delay variable only, so its only condition is x2 = x1 on the face where
% the delay variable is 0; the edge conditions x2 = 0 are dropped (they hold
% for x2 = x1(t - s) anyway). W56D: 56 at deg_Q{3} = [1 2 2] has 7,857,280
% variables (m 135,522). Variants:
%   '43lib'   the system of PDE 43 as the example declares it (control)
%   '43thin'  x2 of order (0,1)
%   '56thin'  x2 of order (0,0,1)
% each in 'Q' light dD 0 with the cap at deg_Q{nv} = [1 2 2] (nv = 2 for
% 43, 3 for 56) and at the default deg_Q; results_DIAG56.csv in this
% folder. A diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
of = fullfile(here,'results_DIAG56.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
for v = {'43lib','43thin','56thin'}
    nv = 2;   if startsWith(v{1},'56'),  nv = 3;  end
    J = {sprintf('st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''; st.sop.deg_Q{%d} = [1 2 2]',nv);
         'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto'''};
    try,  PIE = build56(v{1});
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v{1},strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG56 %s: build error %s\n',v{1},ME.message);  continue
    end
    for i = 1:numel(J)
        mods = J{i};
        st = lpisettings_sop('light',0,'','mosek');
        eval([mods,';']);
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,Q,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},strrep(mods,'''',''),r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAG56 %s Q [%s]: %s fr %.3g exst %s ndec %g\n',v{1},mods,r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
fprintf('## DIAG56 DONE %s\n',char(datetime('now')));
end

function PIE = build56(v)
pvar s1 s2 s3 t
clear stateNameGenerator
switch v
    case {'43lib','43thin'}             % PDE 43: c = 1, a0 = 1.9, a1 = 1, tau = 1 on [0,pi]
        tau = 1;
        x1 = pde_var(s1,[0,pi]);
        if strcmp(v,'43lib'),  x2 = pde_var([s1;s2],[0,pi;0,tau],[2;1]);
        else,                  x2 = pde_var([s1;s2],[0,pi;0,tau],[0;1]);
        end
        PDE_t = [diff(x1,'t')==diff(x1,s1,2)+1.9*x1-subs(x2,s2,tau);
                 diff(x2,'t')==-diff(x2,s2);
                 subs(x1,s1,0)==0;          subs(x1,s1,pi)==0;
                 subs(x2,s2,0)==x1];
        if strcmp(v,'43lib'),  PDE_t = [PDE_t; subs(x2,s1,0)==0; subs(x2,s1,pi)==0];  end
    case '56thin'                       % PDE 56: r = 2.9, tau = 0.9*1.034726 on [0,pi]^2
        r = 2.9;   tau = 0.9*1.034726;
        x1 = pde_var([s1;s2],[0,pi;0,pi]);
        x2 = pde_var([s1;s2;s3],[0,pi;0,pi;0,tau],[0;0;1]);
        PDE_t = [diff(x1,'t')==diff(x1,s1,2)+diff(x1,s2,2)+r*x1-subs(x2,s3,tau);
                 diff(x2,'t')==-diff(x2,s3);
                 subs(x1,s1,0)==0;   subs(x1,s1,pi)==0;   subs(x1,s2,0)==0;   subs(x1,s2,pi)==0;
                 subs(x2,s3,0)==x1];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

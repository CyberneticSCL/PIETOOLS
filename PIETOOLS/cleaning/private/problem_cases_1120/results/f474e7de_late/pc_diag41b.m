function pc_diag41b()
% PC_DIAG41B() settings on PDE 41 written with its delay as an explicit
% transport state (the branch expand_delays = true of the example file:
% x2_t = x2_s2 on [1, 1+tau], x1(1) = x2(1)). PC_DIAG41 (11/20/2026): this
% form gives a 'Q' program of 2,496 variables that MOSEK solves to optimal
% (eta 2e-10, measure 0.34), where the example as written (the delay
% expanded by the parser) ends 'unknown' at 19,776. Results
% results_DIAG41B.csv in this folder. MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG41B.csv');
fid = fopen(of,'w');  fprintf(fid,'form,tier,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'Q','light',''; 'Q','stripped',''; 'Q','heavy',''; 'Q','veryheavy','';
     'Q','light','st.epneg = 1e-3'; 'Q','light','st.epneg = 1e-2'; 'Q','light','st.eppos = 1e-2';
     'Q','light','st.lpi.ineq.dD = 0'; 'Q','light','st.lpi.ineq.dD = 2';
     'Q','light','st.sop.deg_Q{1} = [6 2 2]'; 'Q','light','st.sop.deg_Q{1} = [8 8 8]';
     'Q','light','st.sop.deg_Q{1} = [2 2 2]'; 'Q','light','st.epneg = 1e-3; st.lpi.ineq.dD = 0';
     'P','light',''; 'P','stripped',''; 'P','heavy',''; 'P','veryheavy','';
     'P','light','st.epneg = 1e-3'; 'P','light','st.eppos = 1e-2'; 'P','light','st.lpi.ineq.dD = 0'};
PIE = build41(1,1);
for i = 1:size(J,1)
    [form,tier,mods] = J{i,:};
    if strcmp(form,'P')
        setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
    elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
    end
    st = lpisettings_sop(tier,0,'','mosek');
    if ~isempty(mods),  eval([mods,';']);  end
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',form,tier,mods,r.status, ...
        g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAG41B %s %s [%s]: %s fr %.3g exst %s\n',form,tier,mods,r.status,g(r,'fr'),f(r,'exst'));
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG41B DONE %s\n',char(datetime('now')));
end

function PIE = build41(lam,tau)
% PDE 41 with its delay as the transport x2 on s2 in [1, 1+tau].
pvar s1 s2 t
clear stateNameGenerator
pde_var x1 x2 control u
x1.vars = s1;   x1.dom = [0,1];
x2.vars = s2;   x2.dom = [1,1+tau];
PDE_t = [diff(x1,t)==diff(x1,s1,2)+lam*x1;
         diff(x2,t)==diff(x2,s2);
         subs(x1,s1,0)==0;  subs(x1,s1,1)==subs(x2,s2,1);
         subs(x2,s2,1+tau)==u];
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

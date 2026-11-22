function pc_diag56t()
% PC_DIAG56T() measures how far below the delay threshold the default
% storage degree certifies PDE 56 (2-D heat with an interior delay, delay a
% 3-D transport state; the library delay tau = 0.9*tau*, tau* = 1.034726,
% infeasible at every setting within the local gate) and its 1-D twin PDE
% 43 (library tau = 1 = 0.966*tau*, certified only at deg_Q{2} = [1 2 2]).
% The delay tau = f*tau* is set through the parameters of the example files:
%   43 at f = 0.5, 0.75, 0.9: 'Q' light dD 0 with the cap (default deg_Q)
%   56 at f = 0.25, 0.5, 0.75: 'Q' light dD 0 (default deg_Q; at f = 0.9
%      840,156 variables, infeasible)
% results_DIAG56T.csv in this folder. A diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
of = fullfile(here,'results_DIAG56T.csv');
fid = fopen(of,'w');  fprintf(fid,'system,f,tau,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
ts = 1.034726;
J = {43,0.5,'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''';
     43,0.75,'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''';
     43,0.9,'st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''';
     56,0.25,'st.lpi.ineq.dD = 0';
     56,0.5,'st.lpi.ineq.dD = 0';
     56,0.75,'st.lpi.ineq.dD = 0'};
for i = 1:size(J,1)
    [sys,fr0,mods] = J{i,:};   tau = fr0*ts;
    try
        if sys==43,  PDE_t = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(0,{sprintf('tau=%.10g;',tau)});
        else,        PDE_t = PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman(0,{sprintf('tau=%.10g;',tau)});
        end
        PIE = convert(PDE_t,'pie');
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%d,%g,%.6g,"%s","build error: %s"\n',sys,fr0,tau,strrep(mods,'''',''),strrep(ME.message,',',';'));  fclose(fid);
        continue
    end
    st = lpisettings_sop('light',0,'','mosek');
    eval([mods,';']);
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    fid = fopen(of,'a');
    fprintf(fid,'%d,%g,%.6g,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',sys,fr0,tau,strrep(mods,'''',''),r.status, ...
        g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAG56T %d f %g [%s]: %s fr %.3g exst %s ndec %g\n',sys,fr0,mods,r.status, ...
        g(r,'fr'),f(r,'exst'),g(r,'ndec'));
end
fprintf('## DIAG56T DONE %s\n',char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

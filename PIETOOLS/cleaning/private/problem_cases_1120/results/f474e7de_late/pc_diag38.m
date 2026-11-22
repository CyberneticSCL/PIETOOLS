function pc_diag38()
% PC_DIAG38() tests whether the order of the velocity state alone decides
% the 2-D stability LPI of a wave in the position state. PDE 38 (u_tt =
% u_s1s1 + u_s2s2 on [0,1]^2, u = 0 on every edge; states u and v = u_t,
% both of order (2,2)) certifies in 'Q' light dD 0 at deg_Q{2} = [0 1 1];
% PDE 69 declares v of order (1,2) and PDE 70 of order (2,1), and neither
% certifies. Variants, without the couplings of 69 and 70:
%   '38'      the system of PDE 38: u, v of order (2,2); u = v = 0 on
%             every edge
%   '38v12'   v of order (1,2): v = 0 at s1 = 1, s2 = 0, s2 = 1 (as 69)
%   '38v21'   v of order (2,1): v = 0 at s1 = 0, s1 = 1, s2 = 1 (as 70)
% each in 'Q' light dD 0 with deg_Q{2} [0 1 1] and in 'Pd' light dD 0;
% results_DIAG38.csv in this folder. A diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG38.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]';
     'Pd','st.lpi.ineq.dD = 0'};
for v = {'38','38v12','38v21'}
    try,  PIE = build38(v{1});
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v{1},strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG38 %s: build error %s\n',v{1},ME.message);  continue
    end
    for i = 1:size(J,1)
        [form,mods] = J{i,:};
        if any(strcmp(form,{'Pd','Qd'}))    % the dual test in that form
            setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
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
        fprintf('##### DIAG38 %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v{1},form,mods,r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG38 DONE %s\n',char(datetime('now')));
end

function PIE = build38(v)
pvar s1 s2 t
clear stateNameGenerator
u = pde_var([s1;s2],[0,1;0,1],[2;2]);
switch v
    case '38',     w = pde_var([s1;s2],[0,1;0,1],[2;2]);
    case '38v12',  w = pde_var([s1;s2],[0,1;0,1],[1;2]);
    case '38v21',  w = pde_var([s1;s2],[0,1;0,1],[2;1]);
end
PDE_t = [diff(u,'t')==w;
         diff(w,'t')==diff(u,s1,2)+diff(u,s2,2);
         subs(u,s1,0)==0;   subs(u,s1,1)==0;   subs(u,s2,0)==0;   subs(u,s2,1)==0];
switch v
    case '38',     PDE_t = [PDE_t; subs(w,s1,0)==0; subs(w,s1,1)==0; subs(w,s2,0)==0; subs(w,s2,1)==0];
    case '38v12',  PDE_t = [PDE_t; subs(w,s1,1)==0; subs(w,s2,0)==0; subs(w,s2,1)==0];
    case '38v21',  PDE_t = [PDE_t; subs(w,s1,0)==0; subs(w,s1,1)==0; subs(w,s2,1)==0];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

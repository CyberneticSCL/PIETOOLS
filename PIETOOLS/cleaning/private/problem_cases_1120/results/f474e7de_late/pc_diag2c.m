function pc_diag2c()
% PC_DIAG2C() measures which coupling defeats the 2-D 'Q' test of PDE 69 and
% 70. The uncoupled wave of PDE 38 certifies in 'Q' light dD 0 at deg_Q{2}
% = [0 1 1] with the velocity of order (2,2), (1,2) or (2,1) (DIAG38); 69
% and 70 do not. Variants on [0,1]^2, states and orders as in 69 and 70:
%   '70neu'   ph (2,2), ps (2,1); ph = ps = 0 at s1 = 0, 1 and s2 = 1;
%             ph_s2(s1,0) = 0 (no edge states; conservative)
%   '70damp'  as '70neu' with ph_s2(s1,0) = ps(s1,0) (boundary damping at
%             s2 = 0 in place of the acoustic edge states de, et:
%             dE/dt = -int ps(s1,0)^2 ds1)
%   '69dec'   u (2,2), v (1,2), w (2,2) of 69 with the interface removed:
%             u_s1(0,s2) = 0, w(0,s2) = 0; u = v = w = 0 on the other edges
% each in 'Q' light dD 0 with deg_Q{2} [0 1 1] and 'Q' light dD 0;
% results_DIAG2C.csv in this folder. A diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
of = fullfile(here,'results_DIAG2C.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]';
     'Q','st.lpi.ineq.dD = 0'};
for v = {'70neu','70damp','69dec'}
    try,  PIE = build2c(v{1});
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v{1},strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG2C %s: build error %s\n',v{1},ME.message);  continue
    end
    for i = 1:size(J,1)
        [form,mods] = J{i,:};
        st = lpisettings_sop('light',0,'','mosek');
        if ~isempty(mods),  eval([mods,';']);  end
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},form,mods,r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAG2C %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v{1},form,mods,r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
fprintf('## DIAG2C DONE %s\n',char(datetime('now')));
end

function PIE = build2c(v)
pvar s1 s2 t
clear stateNameGenerator
switch v
    case {'70neu','70damp'}
        ph = pde_var([s1;s2],[0,1;0,1],[2;2]);
        ps = pde_var([s1;s2],[0,1;0,1],[2;1]);
        if strcmp(v,'70neu'),  bc = subs(diff(ph,s2),s2,0)==0;
        else,                  bc = subs(diff(ph,s2),s2,0)==subs(ps,s2,0);
        end
        PDE_t = [diff(ph,'t')==ps;
                 diff(ps,'t')==diff(ph,s1,2)+diff(ph,s2,2);
                 subs([ph;ps],s1,0)==0;     subs([ph;ps],s1,1)==0;
                 bc;
                 subs(ph,s2,1)==0;          subs(ps,s2,1)==0];
    case '69dec'
        u = pde_var([s1;s2],[0,1;0,1],[2;2]);
        w1 = pde_var([s1;s2],[0,1;0,1],[1;2]);
        w = pde_var([s1;s2],[0,1;0,1],[2;2]);
        PDE_t = [diff(u,'t')==w1;
                 diff(w1,'t')==diff(u,s1,2)+diff(u,s2,2);
                 diff(w,'t')==diff(w,s1,2)+diff(w,s2,2);
                 subs(u,s1,1)==0;           subs(diff(u,s1),s1,0)==0;
                 subs(w1,s1,1)==0;
                 subs(w,s1,1)==0;           subs(w,s1,0)==0;
                 subs([u;w1;w],s2,0)==0;    subs([u;w1;w],s2,1)==0];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

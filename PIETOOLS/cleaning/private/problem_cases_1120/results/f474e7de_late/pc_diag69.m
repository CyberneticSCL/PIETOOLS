function pc_diag69()
% PC_DIAG69() tests in one variable whether the states of PDE 69 and 70
% (2-D; infeasible with a certificate at every setting of 11/19-11/20)
% defeat the stability LPI. Both carry the position: u, v = u_t (69, wave
% coupled to heat at s1 = 0, Batty 2019) and ph, ps = ph_t (70, acoustic
% boundary at s2 = 0, Beale 1976). PDE 68, the 1-D twin of 69 in the
% energy states e = u_s, v, certifies. Variants on s in [0,1]:
%   '68'      the example file of PDE 68 (e, v, w)
%   '69pos'   u (order 2), v (order 1), w (order 2): u_t = v, v_t = u_ss,
%             w_t = w_ss; u(1) = v(1) = w(1) = 0, u_s(0) = -w_s(0),
%             w(0) = v(0) (the s1-direction of 69)
%   '70riem'  e = ph_s, v = ps (order 1) and ODE states de, et:
%             e_t = v_s, v_t = e_s, de_t = et, et_t = -et - de - v(0);
%             v(1) = 0, e(0) = -et
%   '70pos'   ph (order 2), ps (order 1), de, et: ph_t = ps, ps_t = ph_ss,
%             de_t = et, et_t = -et - de - ps(0); ph(1) = ps(1) = 0,
%             ph_s(0) = -et (the s2-direction of 70)
% The energy (1/2)(|e|^2 + |v|^2 + de^2 + et^2) has derivative -et^2 in
% '70riem'; 68 and '69pos' lose |w_s|^2. Each variant in 'Q', 'P', 'Pd' and
% 'Qd' at light and 'P', 'Pd' at stripped; results_DIAG69.csv in this
% folder. A diagnostic, not a setting.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
of = fullfile(here,'results_DIAG69.csv');
fid = fopen(of,'w');  fprintf(fid,'variant,form,tier,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
J = {'Q','light'; 'P','light'; 'Pd','light'; 'Qd','light'; 'P','stripped'; 'Pd','stripped'};
for v = {'68','69pos','70riem','70pos'}
    try,  PIE = build69(v{1});
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v{1},strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG69 %s: build error %s\n',v{1},ME.message);  continue
    end
    for i = 1:size(J,1)
        [form,tier] = J{i,:};
        if strcmp(form,'P')         % the primal 'Q' of the executor posed in form 'P'
            setappdata(0,'LB_FORM_MAP',struct('Q','P','Qd','Pd','only',{{'PIETOOLS_stability_sop'}}));
        elseif any(strcmp(form,{'Pd','Qd'}))    % the dual test in that form
            setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
        elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
        end
        st = lpisettings_sop(tier,0,'','mosek');
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,%s,%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},form,tier,r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAG69 %s %s %s: %s fr %.3g exst %s ndec %g\n',v{1},form,tier,r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG69 DONE %s\n',char(datetime('now')));
end

function PIE = build69(v)
if strcmp(v,'68')
    PDE_t = PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty(0,{});
    PIE = convert(PDE_t,'pie');   return
end
pvar s t
clear stateNameGenerator
switch v
    case '69pos'        % 69 along s1: position u, velocity v, heat w
        u = pde_var(s,[0,1],2);  vv = pde_var(s,[0,1],1);  w = pde_var(s,[0,1],2);
        PDE_t = [diff(u,'t')==vv;
                 diff(vv,'t')==diff(u,s,2);
                 diff(w,'t')==diff(w,s,2);
                 subs(u,s,1)==0;           subs(diff(u,s),s,0)==-subs(diff(w,s),s,0);
                 subs(vv,s,1)==0;
                 subs(w,s,1)==0;           subs(w,s,0)==subs(vv,s,0)];
    case '70riem'       % 70 along s2 in the energy states, edge states as ODEs
        e = pde_var(s,[0,1],1);  vv = pde_var(s,[0,1],1);  de = pde_var();  et = pde_var();
        PDE_t = [diff(e,'t')==diff(vv,s);
                 diff(vv,'t')==diff(e,s);
                 diff(de,'t')==et;
                 diff(et,'t')==-et-de-subs(vv,s,0);
                 subs(vv,s,1)==0;          subs(e,s,0)==-et];
    case '70pos'        % 70 along s2 in the position ph and velocity ps
        ph = pde_var(s,[0,1],2);  ps = pde_var(s,[0,1],1);  de = pde_var();  et = pde_var();
        PDE_t = [diff(ph,'t')==ps;
                 diff(ps,'t')==diff(ph,s,2);
                 diff(de,'t')==et;
                 diff(et,'t')==-et-de-subs(ps,s,0);
                 subs(ph,s,1)==0;          subs(diff(ph,s),s,0)==-et;
                 subs(ps,s,1)==0];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

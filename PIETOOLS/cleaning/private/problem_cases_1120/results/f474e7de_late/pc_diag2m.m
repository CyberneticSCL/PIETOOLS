% function pc_diag2m()                                                      % MMP, 11/20/2026 (12:37, was)
function pc_diag2m(J, tag, V)                                               % MMP, 11/20/2026 (12:37)
% PC_DIAG2M() tests PDE 69 and 70 in a mixed position-strain form. DIAG2C
% (11/20/2026): 69 without its heat interface and 70 without its acoustic
% boundary certify in 'Q'; 70 with a plain damping condition ph_s2 = ps at
% s2 = 0 does not; in one variable 'Q' and 'P' fail on a wave in the
% position state whose Neumann datum is another trace (dw1, DIAG44F), and
% the energy (Riemann) forms certify. Here the position keeps only the
% derivatives of the uncoupled direction and the strain of the coupled
% direction is a state, so the coupling becomes a Dirichlet condition on
% the strain:
%   '70mix'  ph (2,0), e2 = ph_s2 (0,1), ps (2,1), de, et (order 2 in s1):
%            ph_t = ps, ps_t = ph_s1s1 + e2_s2, e2_t = ps_s2, de_t = et,
%            et_t = -et - de - ps(s1,0); ph = ps = 0 at s1 = 0, 1;
%            ps(s1,1) = 0; e2(s1,0) = -et; de = et = 0 at s1 = 0, 1.
%            E = (|ps|^2 + |ph_s1|^2 + |e2|^2 + |de|^2 + |et|^2)/2,
%            dE/dt = -|et|^2.
%   '69mix'  u (0,2), e1 = u_s1 (1,0), v (1,2), w (2,2):
%            u_t = v, v_t = e1_s1 + u_s2s2, e1_t = v_s1, w_t = w_s1s1 +
%            w_s2s2; u = 0 at s2 = 0, 1; e1(0,s2) = -w_s1(0,s2);
%            v(1,s2) = 0, v = 0 at s2 = 0, 1; w(1,s2) = 0, w(0,s2) =
%            v(0,s2), w = 0 at s2 = 0, 1.
%            E = (|v|^2 + |e1|^2 + |u_s2|^2 + |w|^2)/2, dE/dt = -|grad w|^2.
% Data with e2 ~= ph_s2 (e1 ~= u_s1) keep the difference constant, so the
% extended systems are not asymptotically stable; the test is epneg = 0
% (a non-increasing storage), as for the examples. Each in 'Q' light dD 0
% with deg_Q{2} [0 1 1], 'Q' light dD 0 and 'Pd' light dD 0;
% results_DIAG2M.csv in this folder. A diagnostic, not a setting.
% Optional: J, rows {form, mods, tier} in place of those above (tier light
% when omitted); tag, the results file results_<tag>.csv (default DIAG2M);
% V, the variants (default both).
% MMP, 11/20/2026.
% MMP, 11/20/2026 (12:37): inputs J, tag, V; the tier per row (DIAG70N:
%                setting cuts on the '70mix' certificate).
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if nargin < 2 || isempty(tag),  tag = 'DIAG2M';  end                       % MMP, 11/20/2026 (12:37)
% of = fullfile(here,'results_DIAG2M.csv');                                 % MMP, 11/20/2026 (12:37, was)
of = fullfile(here,['results_',tag,'.csv']);                                % MMP, 11/20/2026 (12:37)
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
if nargin < 1 || isempty(J)                                                 % MMP, 11/20/2026 (12:37)
    J = {'Q','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]';
         'Q','st.lpi.ineq.dD = 0';
         'Pd','st.lpi.ineq.dD = 0'};
end                                                                         % MMP, 11/20/2026 (12:37)
if size(J,2) < 3,  J(:,3) = {'light'};  end     % tier per row             % MMP, 11/20/2026 (12:37)
if nargin < 3 || isempty(V),  V = {'70mix','69mix'};  end                   % MMP, 11/20/2026 (12:37)
% for v = {'70mix','69mix'}                                                 % MMP, 11/20/2026 (12:37, was)
for v = V                                                                   % MMP, 11/20/2026 (12:37)
    try,  PIE = build2m(v{1});
    catch ME
        fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v{1},strrep(ME.message,',',';'));  fclose(fid);
        fprintf('##### DIAG2M %s: build error %s\n',v{1},ME.message);  continue
    end
    for i = 1:size(J,1)
%         [form,mods] = J{i,:};                                             % MMP, 11/20/2026 (12:37, was)
        [form,mods,tier] = J{i,:};                                          % MMP, 11/20/2026 (12:37)
        if any(strcmp(form,{'Pd','Qd'}))    % the dual test in that form
            setappdata(0,'LB_FORM_MAP',struct('Q',form,'only',{{'PIETOOLS_stability_sop'}}));
        elseif isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');
        end
%         st = lpisettings_sop('light',0,'','mosek');                       % MMP, 11/20/2026 (12:37, was)
        st = lpisettings_sop(tier,0,'','mosek');                            % MMP, 11/20/2026 (12:37)
        if ~isempty(mods),  eval([mods,';']);  end
        try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
        catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
        end
        if ~strcmp(tier,'light'),  mods = [tier,': ',mods];  end   % tier in the record % MMP, 11/20/2026 (12:37)
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v{1},form,mods,r.status, ...
            g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
        fclose(fid);
        fprintf('##### DIAG2M %s %s [%s]: %s fr %.3g exst %s ndec %g\n',v{1},form,mods,r.status, ...
            g(r,'fr'),f(r,'exst'),g(r,'ndec'));
    end
end
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
fprintf('## DIAG2M DONE %s\n',char(datetime('now')));
end

function PIE = build2m(v)
pvar s1 s2 t
clear stateNameGenerator
switch v
    case '70mix'
        ph = pde_var([s1;s2],[0,1;0,1],[2;0]);
        e2 = pde_var([s1;s2],[0,1;0,1],[0;1]);
        ps = pde_var([s1;s2],[0,1;0,1],[2;1]);
        de = pde_var(s1,[0,1],2);
        et = pde_var(s1,[0,1],2);
        PDE_t = [diff(ph,'t')==ps;
                 diff(ps,'t')==diff(ph,s1,2)+diff(e2,s2);
                 diff(e2,'t')==diff(ps,s2);
                 diff(de,'t')==et;
                 diff(et,'t')==-et-de-subs(ps,s2,0);
                 subs([ph;ps],s1,0)==0;     subs([ph;ps],s1,1)==0;
                 subs(ps,s2,1)==0;
                 subs(e2,s2,0)==-et;
                 subs([de;et],s1,0)==0;     subs([de;et],s1,1)==0];
    case '69mix'
        u = pde_var([s1;s2],[0,1;0,1],[0;2]);
        e1 = pde_var([s1;s2],[0,1;0,1],[1;0]);
        vv = pde_var([s1;s2],[0,1;0,1],[1;2]);
        w = pde_var([s1;s2],[0,1;0,1],[2;2]);
        PDE_t = [diff(u,'t')==vv;
                 diff(vv,'t')==diff(e1,s1)+diff(u,s2,2);
                 diff(e1,'t')==diff(vv,s1);
                 diff(w,'t')==diff(w,s1,2)+diff(w,s2,2);
                 subs(u,s2,0)==0;           subs(u,s2,1)==0;
                 subs(e1,s1,0)==-subs(diff(w,s1),s1,0);
                 subs(vv,s1,1)==0;          subs(vv,s2,0)==0;    subs(vv,s2,1)==0;
                 subs(w,s1,1)==0;           subs(w,s1,0)==subs(vv,s1,0);
                 subs(w,s2,0)==0;           subs(w,s2,1)==0];
end
PIE = convert(PDE_t,'pie');
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

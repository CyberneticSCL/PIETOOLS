% function pc_diag69m()                                                     % MMP, 11/20/2026 (12:35, was)
function pc_diag69m(J, tag, T)                                              % MMP, 11/20/2026 (12:35)
% PC_DIAG69M() tests PDE 69 in the mixed position-strain form of DIAG2M
% with the strain of order (1,2). DIAG2M '69mix' (e1 = u_s1 of order (1,0))
% did not convert: "The boundary conditions along s1 are order
% incompatible" – the interface condition e1(0,s2) = -w_s1(0,s2) relates a
% state of order 0 in s2 to w of order 2 in s2. Since u = 0 at s2 = 0, 1,
% e1 = u_s1 = 0 there, so e1 of order (1,2) has these two conditions:
%   '69mix2'  u (0,2), e1 = u_s1 (1,2), v (1,2), w (2,2):
%             u_t = v, v_t = e1_s1 + u_s2s2, e1_t = v_s1, w_t = w_s1s1 +
%             w_s2s2; u = 0 at s2 = 0, 1; e1(0,s2) = -w_s1(0,s2), e1 = 0
%             at s2 = 0, 1; v(1,s2) = 0, v = 0 at s2 = 0, 1; w(1,s2) = 0,
%             w(0,s2) = v(0,s2), w = 0 at s2 = 0, 1.
%             E = (|v|^2 + |e1|^2 + |u_s2|^2 + |w|^2)/2, dE/dt = -|grad w|^2.
% Data with e1 ~= u_s1 keep the difference constant, so the test is epneg
% = 0, as for the example. In 'Q' light dD 0, 'Q' light dD 0 with deg_Q{2}
% [0 1 1]; results_DIAG69M.csv in this folder. A diagnostic, not a setting.
% Optional: J, mods per row in place of those above; tag, the results file
% results_<tag>.csv (default DIAG69M); T, the tier per row (default light;
% a tier other than light is written before the mods).
% MMP, 11/20/2026.
% MMP, 11/20/2026 (12:35): inputs J, tag and T (DIAG69N: settings around the
%                certificate at light dD 0, measure 0.722, near the gate).
here = fileparts(mfilename('fullpath'));    cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
if isappdata(0,'LB_FORM_MAP'),  rmappdata(0,'LB_FORM_MAP');  end
if nargin < 2 || isempty(tag),  tag = 'DIAG69M';  end                      % MMP, 11/20/2026 (12:35)
% of = fullfile(here,'results_DIAG69M.csv');                                % MMP, 11/20/2026 (12:35, was)
of = fullfile(here,['results_',tag,'.csv']);                                % MMP, 11/20/2026 (12:35)
fid = fopen(of,'w');  fprintf(fid,'variant,form,mods,status,fr,numerr,exst,eta,ndec,m,sec\n');  fclose(fid);
% J = {'st.lpi.ineq.dD = 0';                                                % MMP, 11/20/2026 (12:35, was)
%      'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]'};                    % MMP, 11/20/2026 (12:35, was)
if nargin < 1 || isempty(J)                                                 % MMP, 11/20/2026 (12:35)
    J = {'st.lpi.ineq.dD = 0';                                              % MMP, 11/20/2026 (12:35)
         'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]'};                  % MMP, 11/20/2026 (12:35)
end                                                                         % MMP, 11/20/2026 (12:35)
if nargin < 3 || isempty(T),  T = repmat({'light'},size(J));  end           % MMP, 11/20/2026 (12:35)
v = '69mix2';
try
    pvar s1 s2 t
    clear stateNameGenerator
    u = pde_var([s1;s2],[0,1;0,1],[0;2]);
    e1 = pde_var([s1;s2],[0,1;0,1],[1;2]);
    vv = pde_var([s1;s2],[0,1;0,1],[1;2]);
    w = pde_var([s1;s2],[0,1;0,1],[2;2]);
    PDE_t = [diff(u,'t')==vv;
             diff(vv,'t')==diff(e1,s1)+diff(u,s2,2);
             diff(e1,'t')==diff(vv,s1);
             diff(w,'t')==diff(w,s1,2)+diff(w,s2,2);
             subs(u,s2,0)==0;           subs(u,s2,1)==0;
             subs(e1,s1,0)==-subs(diff(w,s1),s1,0);
             subs(e1,s2,0)==0;          subs(e1,s2,1)==0;
             subs(vv,s1,1)==0;          subs(vv,s2,0)==0;    subs(vv,s2,1)==0;
             subs(w,s1,1)==0;           subs(w,s1,0)==subs(vv,s1,0);
             subs(w,s2,0)==0;           subs(w,s2,1)==0];
    PIE = convert(PDE_t,'pie');
catch ME
    fid = fopen(of,'a');  fprintf(fid,'%s,,,"build error: %s"\n',v,strrep(ME.message,',',';'));  fclose(fid);
    fprintf('##### DIAG69M %s: build error %s\n',v,ME.message);
    fprintf('## DIAG69M DONE %s\n',char(datetime('now')));   return
end
for i = 1:numel(J)
    mods = J{i};
%     st = lpisettings_sop('light',0,'','mosek');                           % MMP, 11/20/2026 (12:35, was)
    st = lpisettings_sop(T{i},0,'','mosek');                                % MMP, 11/20/2026 (12:35)
    eval([mods,';']);
    try,  r = lb_exec('sop:PIETOOLS_PIE2PDEstability', PIE, st);
    catch ME,  r = struct('status',['error: ',strrep(ME.message,',',';')]);
    end
    if ~strcmp(T{i},'light'),  mods = [T{i},': ',mods];  end    % tier in the record % MMP, 11/20/2026 (12:35)
    fid = fopen(of,'a');
    fprintf(fid,'%s,Q,"%s",%s,%.4g,%g,%s,%.3g,%.0f,%.0f,%.1f\n',v,mods,r.status, ...
        g(r,'fr'),g(r,'numerr'),f(r,'exst'),g(r,'eta'),g(r,'ndec'),g(r,'m'),g(r,'sec'));
    fclose(fid);
    fprintf('##### DIAG69M %s Q [%s]: %s fr %.3g exst %s ndec %g\n',v,mods,r.status, ...
        g(r,'fr'),f(r,'exst'),g(r,'ndec'));
end
fprintf('## DIAG69M DONE %s\n',char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end

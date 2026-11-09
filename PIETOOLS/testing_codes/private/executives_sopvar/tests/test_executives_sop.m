function test_executives_sop(quick)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_EXECUTIVES_SOP(QUICK) tests the container executives at
% lpisettings_sop('light') with MOSEK when installed, and errors at the first
% failed check:
%   (a) on the 1-D PIE  x_t = x_ss + pi^2/2 x + s(1-s) w + u, x(0) = x(1) = 0,
%       z = [int x; u], y = int x,  every form of PIETOOLS_Hinf_gain_sop and
%       PIETOOLS_H2_norm_sop and the four synthesis executives certify, each
%       bound between the value of the Legendre discretization and that
%       value plus 1e-3 (H-infinity) or 5e-3 (H2) relative, each synthesis
%       value at most the open-loop value;
%   (b) with settings.sop.slack = 'stock', PIETOOLS_Hinf_gain_sop and
%       PIETOOLS_H2_norm_sop (form 'c') return the recorded value of the
%       stock PIETOOLS_Hinf_gain and PIETOOLS_H2_norm_c to 1e-4 relative;
%   (c) on x_t = x_ss + pi^2/2 x, x(0) = x(1) = 0, the four forms of
%       PIETOOLS_stability_sop and PIETOOLS_well_posedness_sop certify, with
%       the largest real part of the discretized spectrum at -pi^2/2 to 1e-3;
%   (d) with QUICK false, PIETOOLS_stability_sop (form 'P') certifies the
%       2-D PIE  x_t = x_s1s1 + x_s2s2 + pi^2 x  with Dirichlet conditions.
%
% INPUT
% - quick:  (optional) logical, true (default) to skip check (d).
%
% OUTPUT
% none; one printed line per check.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_executives_sop
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 10/08/2026
% MMP, 10/28/2026: The executives called are the consolidated ones
%                (PIETOOLS_stability_sop, PIETOOLS_Hinf_gain_sop and
%                PIETOOLS_H2_norm_sop with a form, in place of the removed
%                per-form files); header to the header protocol, the licence
%                block added.
% MMP, 11/04/2026: The gains and solved operators are checked as 'copvar'
%                objects: the executives return containers.
% MMP, 11/07/2026: Part (b) compares with recorded values of the stock
%                executives: the PIEs hold 'copvar' operators, on which the
%                stock PIETOOLS_Hinf_gain and PIETOOLS_H2_norm_c fail at
%                their first poslpivar call (Top.dim), so the test ended
%                there, before parts (c) and (d).
% MMP, 11/08/2026: Settings lpisettings_sop('light') (was the stock
%                lpisettings('light'), retired with its folder).
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): the help names
%                lpisettings_sop; code unchanged.

% Settings, solver and messages; counter of passed checks
if nargin<1,    quick = true;   end
warning('off','sopvar:noncanonicalMultiplier');  warning('off','sdopvar:noncanonicalMultiplier');
% st = lpisettings('light');                                                % MMP, 11/08/2026 (was)
st = lpisettings_sop('light');                                              % MMP, 11/08/2026
if ~isempty(which('mosekopt')),     st.sos_opts.solver = 'mosek';   end
st.sop.verbose = false;
n = 0;
% ---- the plants
pvar s t
x = pde_var(1,s,[0,1]);     w = pde_var('in',1);     u = pde_var('control',1);
y = pde_var('sense',1);     z = pde_var('out',2);
PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x+(s*(1-s))*w+u;
       z==[int(x,s,[0,1]); u];     y==int(x,s,[0,1]);
       subs(x,s,0)==0;  subs(x,s,1)==0];
evalc('PIE = initialize(convert(PDE));');
x2 = pde_var(1,s,[0,1]);
PDE2 = [diff(x2,t,1)==diff(x2,s,2)+0.5*pi^2*x2;   subs(x2,s,0)==0;  subs(x2,s,1)==0];
evalc('PIEs = initialize(convert(PDE2));');
% numerical values
Wg = pie_witness_sop(PIE,'gain',struct('N_cheb',24));    g_num = Wg.gain;
Wh = pie_witness_sop(PIE,'h2',struct('N_cheb',24));      h_num = Wh.h2;
fprintf('numerical open-loop values: Hinf gain %.7g, H2 norm %.7g\n',g_num,h_num);
% ---- (a) gains, H2 norms, synthesis with the 'like' slack
% for nm = {'Hinf_gain','Hinf_gain_coercive','Hinf_gain_dual','Hinf_gain_dual_coercive'}
%     [~,~,g,info] = feval(['PIETOOLS_' nm{1} '_sop'],PIE,st);              % MMP, 10/28/2026 (was)
for fm = {'Q','P','Qd','Pd'}                                                % MMP, 10/28/2026
    [~,~,g,info] = PIETOOLS_Hinf_gain_sop(PIE,st,fm{1});                    % MMP, 10/28/2026
    n = check(n,info.certified && g>=g_num*(1-1e-6) && g<=g_num*(1+1e-3), ...
              sprintf('Hinf_gain_sop %s: %s, gamma %.7g (numerical %.7g)',fm{1},info.status,g,g_num)); % MMP, 10/28/2026
end
% for nm = {'H2_norm_c','H2_norm_o','H2_norm_c_coercive','H2_norm_o_coercive'}
%     out = cell(1,nargout(['PIETOOLS_' nm{1} '_sop']));
%     [out{:}] = feval(['PIETOOLS_' nm{1} '_sop'],PIE,st);
%     g = out{3};     info = out{end};                                      % MMP, 10/28/2026 (was)
for fm = {'c','o','cco','oco'}                                              % MMP, 10/28/2026
    [~,~,g,~,~,info] = PIETOOLS_H2_norm_sop(PIE,st,fm{1});                  % MMP, 10/28/2026
    n = check(n,info.certified && g>=h_num*(1-1e-6) && g<=h_num*(1+5e-3), ...
              sprintf('H2_norm_sop %s: %s, value %.7g (numerical %.7g)',fm{1},info.status,g,h_num)); % MMP, 10/28/2026
end
for nm = {'Hinf_control','Hinf_estimator','H2_control','H2_estimator'}
    out = cell(1,nargout(['PIETOOLS_' nm{1} '_sop']));
    [out{:}] = feval(['PIETOOLS_' nm{1} '_sop'],PIE,st);
    g = out{3};     info = out{end};
    ref = g_num;    if contains(nm{1},'H2'),    ref = h_num;    end
%   n = check(n,info.certified && isa(out{2},'opvar') && g<=ref*(1+1e-3), ... % MMP, 11/04/2026 (was)
    n = check(n,info.certified && isa(out{2},'copvar') && g<=ref*(1+1e-3), ...
              sprintf('%s_sop: %s, gamma %.7g (open loop %.7g), gain class %s',nm{1},info.status,g,ref,class(out{2}))); % MMP, 11/04/2026
end
% ---- (b) the stock slack reproduces the stock
st2 = st;   st2.sop.slack = 'stock';
% evalc('[~,~,gs] = PIETOOLS_Hinf_gain(PIE,st);');                          % MMP, 11/07/2026 (was)
% Stock values recorded 11/07/2026 at e6319ace, MOSEK 11.0.30: each         % MMP, 11/07/2026
% executive without its pie_struct check, on pie2opvar(PIE), settings st.   % MMP, 11/07/2026
gs = 0.03333235189;                         % stock PIETOOLS_Hinf_gain      % MMP, 11/07/2026
[~,~,g,info] = PIETOOLS_Hinf_gain_sop(PIE,st2);
n = check(n,info.certified && abs(g-gs)<=1e-4*gs,sprintf('Hinf_gain_sop stock slack %.8g vs stock %.8g',g,gs));
% evalc('[~,~,gs] = PIETOOLS_H2_norm_c(PIE,st);');                          % MMP, 11/07/2026 (was)
gs = 0.05271563553;                         % stock PIETOOLS_H2_norm_c      % MMP, 11/07/2026
% [~,~,g,~,~,info] = PIETOOLS_H2_norm_c_sop(PIE,st2);                      % MMP, 10/28/2026 (was)
[~,~,g,~,~,info] = PIETOOLS_H2_norm_sop(PIE,st2,'c');                       % MMP, 10/28/2026
% n = check(n,info.certified && abs(g-gs)<=1e-4*gs,sprintf('H2_norm_c_sop stock slack %.8g vs stock %.8g',g,gs)); % MMP, 10/28/2026 (was)
n = check(n,info.certified && abs(g-gs)<=1e-4*gs,sprintf('H2_norm_sop c stock slack %.8g vs stock %.8g',g,gs)); % MMP, 10/28/2026
% ---- (c) stability and well-posedness
sts = st;   sts.eppos = 1e-4;   sts.eppos2 = 1e-6;    sts.epneg = 0;
% for nm = {'PDEstability','PDEstability_dual','PIE2PDEstability','PIE2PDEstability_dual'}
%     [~,P,info] = feval(['PIETOOLS_' nm{1} '_sop'],PIEs,sts);              % MMP, 10/28/2026 (was)
for fm = {'P','Pd','Q','Qd'}                                                % MMP, 10/28/2026
    [~,P,info] = PIETOOLS_stability_sop(PIEs,sts,fm{1});                    % MMP, 10/28/2026
%   n = check(n,info.certified && isa(P,'opvar') && abs(info.maxre+pi^2/2)<=1e-3, ... % MMP, 11/04/2026 (was)
    n = check(n,info.certified && isa(P,'copvar') && abs(info.maxre+pi^2/2)<=1e-3, ...
              sprintf('stability_sop %s: %s, max Re eig(A,T) %.5g',fm{1},info.status,info.maxre)); % MMP, 10/28/2026; MMP, 11/04/2026
end
[~,P,R,om,info] = PIETOOLS_well_posedness_sop(PIEs,sts);
% n = check(n,info.certified && isa(P,'opvar') && isa(R,'opvar') && om==0, ... % MMP, 11/04/2026 (was)
n = check(n,info.certified && isa(P,'copvar') && isa(R,'copvar') && om==0, ...
          sprintf('well_posedness_sop: %s, omega %g',info.status,om));      % MMP, 11/04/2026
% ---- (d) 2-D stability
if ~quick
    pvar s1 s2
    x3 = pde_var(1,[s1;s2],[0,1;0,1]);
    PDE3 = [diff(x3,t,1)==diff(x3,s1,2)+diff(x3,s2,2)+0.5*2*pi^2*x3;
            subs(x3,s1,0)==0; subs(x3,s1,1)==0;  subs(x3,s2,0)==0; subs(x3,s2,1)==0];
    evalc('PIE3 = initialize(convert(PDE3));');
    st3 = st;   st3.sop.max_raise = 0;
%   [~,P,info] = PIETOOLS_stability_2D_sop(PIE3,st3);                       % MMP, 10/28/2026 (was)
    [~,P,info] = PIETOOLS_stability_sop(PIE3,st3,'P');                      % MMP, 10/28/2026
%   n = check(n,info.certified && isa(P,'opvar2d'),sprintf('stability_sop P (2-D): %s, nx %d',info.status,info.shape.nx)); % MMP, 11/04/2026 (was); MMP, 10/28/2026
    n = check(n,info.certified && isa(P,'copvar'),sprintf('stability_sop P (2-D): %s, nx %d',info.status,info.shape.nx)); % MMP, 11/04/2026
end
fprintf('test_executives_sop: %d checks passed\n',n);
end


function n = check(n,ok,msg)
% N = CHECK(N,OK,MSG) prints MSG and returns N+1 when OK is true, and
% errors with MSG otherwise.
% INPUT  n: double, the checks passed so far; ok: logical; msg: char.
% OUTPUT n: double, N+1.
if ok,  fprintf('  ok   %s\n',msg);   n = n+1;
else,   error('test_executives_sop:fail','FAILED %s',msg);
end
end

function test_executives_sop(quick)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_EXECUTIVES_SOP([QUICK]) asserts on the container executives
% (executives_sopvar), light settings, MOSEK when installed:
%   (a) on a 1-D plant with a finite-dimensional disturbance and input
%       (reaction-diffusion with the profile s(1-s) on w), every gain, H2
%       and synthesis executive certifies with the default 'like' slack, and
%       each certified bound sits between the numerical value of the
%       Legendre discretization (a lower bound) and that value plus a
%       tolerance (1e-3 relative for the H-infinity gains, 5e-3 for the H2
%       norms; the synthesis values at most the open-loop values);
%   (b) with settings.sop.slack = 'stock' the H-infinity gain and the H2 norm
%       (c form) reproduce the stock executive's value to 1e-4 relative;
%   (c) the four stability executives and well-posedness certify the
%       reaction-diffusion plant at half the critical coefficient, with the
%       numerical spectrum's largest real part at -pi^2/2 to 1e-3;
%   (d) QUICK false (default true) adds the 2-D stability executive on the
%       2-D reaction-diffusion plant at half the critical coefficient
%       (about 15 s with MOSEK).
% Each check prints its line; the function errors at the first failure.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1,    quick = true;   end
warning('off','sopvar:noncanonicalMultiplier');  warning('off','sdopvar:noncanonicalMultiplier');
st = lpisettings('light');
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
for nm = {'Hinf_gain','Hinf_gain_coercive','Hinf_gain_dual','Hinf_gain_dual_coercive'}
    [~,~,g,info] = feval(['PIETOOLS_' nm{1} '_sop'],PIE,st);
    n = check(n,info.certified && g>=g_num*(1-1e-6) && g<=g_num*(1+1e-3), ...
              sprintf('%s_sop: %s, gamma %.7g (numerical %.7g)',nm{1},info.status,g,g_num));
end
for nm = {'H2_norm_c','H2_norm_o','H2_norm_c_coercive','H2_norm_o_coercive'}
    out = cell(1,nargout(['PIETOOLS_' nm{1} '_sop']));
    [out{:}] = feval(['PIETOOLS_' nm{1} '_sop'],PIE,st);
    g = out{3};     info = out{end};
    n = check(n,info.certified && g>=h_num*(1-1e-6) && g<=h_num*(1+5e-3), ...
              sprintf('%s_sop: %s, value %.7g (numerical %.7g)',nm{1},info.status,g,h_num));
end
for nm = {'Hinf_control','Hinf_estimator','H2_control','H2_estimator'}
    out = cell(1,nargout(['PIETOOLS_' nm{1} '_sop']));
    [out{:}] = feval(['PIETOOLS_' nm{1} '_sop'],PIE,st);
    g = out{3};     info = out{end};
    ref = g_num;    if contains(nm{1},'H2'),    ref = h_num;    end
    n = check(n,info.certified && isa(out{2},'opvar') && g<=ref*(1+1e-3), ...
              sprintf('%s_sop: %s, gamma %.7g (open loop %.7g), gain class %s',nm{1},info.status,g,ref,class(out{2})));
end
% ---- (b) the stock slack reproduces the stock
st2 = st;   st2.sop.slack = 'stock';
evalc('[~,~,gs] = PIETOOLS_Hinf_gain(PIE,st);');
[~,~,g,info] = PIETOOLS_Hinf_gain_sop(PIE,st2);
n = check(n,info.certified && abs(g-gs)<=1e-4*gs,sprintf('Hinf_gain_sop stock slack %.8g vs stock %.8g',g,gs));
evalc('[~,~,gs] = PIETOOLS_H2_norm_c(PIE,st);');
[~,~,g,~,~,info] = PIETOOLS_H2_norm_c_sop(PIE,st2);
n = check(n,info.certified && abs(g-gs)<=1e-4*gs,sprintf('H2_norm_c_sop stock slack %.8g vs stock %.8g',g,gs));
% ---- (c) stability and well-posedness
sts = st;   sts.eppos = 1e-4;   sts.eppos2 = 1e-6;    sts.epneg = 0;
for nm = {'PDEstability','PDEstability_dual','PIE2PDEstability','PIE2PDEstability_dual'}
    [~,P,info] = feval(['PIETOOLS_' nm{1} '_sop'],PIEs,sts);
    n = check(n,info.certified && isa(P,'opvar') && abs(info.maxre+pi^2/2)<=1e-3, ...
              sprintf('%s_sop: %s, max Re eig(A,T) %.5g',nm{1},info.status,info.maxre));
end
[~,P,R,om,info] = PIETOOLS_well_posedness_sop(PIEs,sts);
n = check(n,info.certified && isa(P,'opvar') && isa(R,'opvar') && om==0, ...
          sprintf('well_posedness_sop: %s, omega %g',info.status,om));
% ---- (d) 2-D stability
if ~quick
    pvar s1 s2
    x3 = pde_var(1,[s1;s2],[0,1;0,1]);
    PDE3 = [diff(x3,t,1)==diff(x3,s1,2)+diff(x3,s2,2)+0.5*2*pi^2*x3;
            subs(x3,s1,0)==0; subs(x3,s1,1)==0;  subs(x3,s2,0)==0; subs(x3,s2,1)==0];
    evalc('PIE3 = initialize(convert(PDE3));');
    st3 = st;   st3.sop.max_raise = 0;
    [~,P,info] = PIETOOLS_stability_2D_sop(PIE3,st3);
    n = check(n,info.certified && isa(P,'opvar2d'),sprintf('stability_2D_sop: %s, nx %d',info.status,info.shape.nx));
end
fprintf('test_executives_sop: %d checks passed\n',n);
end


function n = check(n,ok,msg)
if ok,  fprintf('  ok   %s\n',msg);   n = n+1;
else,   error('test_executives_sop:fail','FAILED %s',msg);
end
end

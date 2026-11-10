% Test_PIESOS_sublevel.m -- tests the PIESOS_sublevel executive, which finds
% the largest sublevel set {V <= c} of a functional V in the semialgebraic
% set {g_i >= 0, h_j = 0}, against exact values (computed by hand, not by
% PIETOOLS). With a constant multiplier on V-c, the value is exact for
% quadratic V, g_i (S-lemma). Checked:
%   - two ODE states x1, x2 with V = x1^2 + x2^2 and
%       g = 1 - x1^2 - 2*x2^2:                 c = 1/2
%       the same with h = x2^2 (x2 = 0):        c = 1 (the equality is used)
%       g = {3 - x1^2, 1.5 - x2^2}:             c = min(3,1.5) = 1.5
%   - one ODE state, V = x^2 + x^4, g = 1 - x^2: c = 1, conservative, as the
%     exact value 2 would need a multiplier that is a polynomial in x
%   - the heat equation x_t = x_ss + x^2, x_s(0)=x(1)=0, with V = ||x_s||^2
%     and the ball g = 1 - ||x||^2 - ||x_s||^2: c = lambda_1/(1+lambda_1) =
%     0.7116, with the Poincare constant lambda_1 = pi^2/4 of x_s(0)=x(1)=0
%     (min ||x_s||^2/(||x||^2+||x_s||^2)); and c = 1 for V = ||x_s||^2 in
%     g = 1 - ||x_s||^2
%   - invalid inputs are refused: g(0)<0, V(0)~=0, empty g
%   - the Lyapunov functional V of PIESOS_stability_local for the heat
%     equation with integral feedback, x_t = x_ss + x^2 - k*int x with
%     k=2*pi^2, in the ball ||x_s||^2 <= R^2, R=14: since V >= ||x_s||^2 and
%     V <= C*||x_s||^2 on the ball, R^2*(1-1e-3) <= c <= C*R^2. This test
%     takes about a minute, and can be skipped by setting skip_executive=true
%     before running the script.
%
% DJ (CC), 11/10/2026: Initial coding;

clear stateNameGenerator
pvar s t
tol = 1e-5;
st = PIESOS_settings();
st.verbose = false;
if ~exist('skip_executive','var')
    skip_executive = false;
end

% % Cases 1-4: ODE states
x1 = pde_var('state',1,[],[]);  x2 = pde_var('state',1,[],[]);
evalc('PIE = convert(initialize([diff(x1,t)==-x1-x1^3; diff(x2,t)==-x2-x2^3]),{''silent''});');
u = PIE.T*PIE.x;
N1 = weighted_sobolev_norm(u,{1;0});  N2 = weighted_sobolev_norm(u,{0;1});
V = N1 + N2;
cases = {
    'g = 1-x1^2-2*x2^2',       1 - N1 - 2*N2,       [],     0.5;
    'same, h = x2^2',          1 - N1 - 2*N2,       N2,     1;
    'g = {3-x1^2, 1.5-x2^2}',  {3 - N1, 1.5 - N2},  [],     1.5};
for k=1:size(cases,1)
    [~,info] = PIESOS_sublevel(PIE,V,cases{k,2},cases{k,3},st);
    if ~info.feasible || abs(info.c-cases{k,4})>tol*cases{k,4}
        error("Test_PIESOS_sublevel: ODE case '%s': c=%g, expected %g.",cases{k,1},info.c,cases{k,4})
    end
end
clear stateNameGenerator
x = pde_var('state',1,[],[]);
evalc('PIE1 = convert(initialize(diff(x,t)==-x-x^3),{''silent''});');
u1 = PIE1.T*PIE1.x;
N = weighted_sobolev_norm(u1,{1});
[~,info] = PIESOS_sublevel(PIE1,N+N*N,1-N,[],st);
if ~info.feasible || abs(info.c-1)>tol
    error("Test_PIESOS_sublevel: quartic V: c=%g, expected 1.",info.c)
end

% % Cases 5-6: heat equation, Poincare constant
clear stateNameGenerator
x = pde_var(s,[0,1]);
evalc('PIE = convert(initialize([diff(x,t)==diff(x,s,2)+x^2; subs(diff(x,s),s,0)==0; subs(x,s,1)==0]),{''silent''});');
u = PIE.T*PIE.x;
Nxs = weighted_sobolev_norm(u,{[0,1]});
lam1 = pi^2/4;
[~,info] = PIESOS_sublevel(PIE,Nxs,1-weighted_sobolev_norm(u,{[1,1]}),[],st);
if ~info.feasible || abs(info.c-lam1/(1+lam1))>tol
    error("Test_PIESOS_sublevel: Poincare constant: c=%g, expected %g.",info.c,lam1/(1+lam1))
end
[~,info] = PIESOS_sublevel(PIE,Nxs,1-Nxs,[],st);
if ~info.feasible || abs(info.c-1)>tol
    error("Test_PIESOS_sublevel: V=N: c=%g, expected 1.",info.c)
end

% % Invalid inputs
bad = {@() PIESOS_sublevel(PIE,Nxs,-1-Nxs,[],st),  "g(0)<0";
       @() PIESOS_sublevel(PIE,Nxs+1,1-Nxs,[],st), "V(0)~=0";
       @() PIESOS_sublevel(PIE,Nxs,{},[],st),      "empty g"};
for k=1:size(bad,1)
    refused = false;
    try
        bad{k,1}();
    catch
        refused = true;
    end
    if ~refused
        error("Test_PIESOS_sublevel: invalid input (%s) was not refused.",bad{k,2})
    end
end

% % Case 7: V of PIESOS_stability_local, heat equation with integral feedback
if ~skip_executive
    clear stateNameGenerator
    k = 2*pi^2;
    x = pde_var(s,[0,1]);
    evalc('PIE = convert(initialize([diff(x,t)==diff(x,s,2)+x^2-k*int(x,s,0,1); subs(diff(x,s),s,0)==0; subs(x,s,1)==0]),{''silent''});');
    R = 14;
    st = PIESOS_settings();
    st.verbose = false;
    st.LF.mon_deg = 2;  st.mult.dist_deg = [1,1];  st.mult.mon_deg = [2,2];
    st.sos.psatz = [0,0,1];  st.simplify = false;
    [~,info_loc] = PIESOS_stability_local(PIE,R,[0,1],0,'pie',st);
    if ~info_loc.feasible
        error("Test_PIESOS_sublevel: the local executive did not verify R=%g.",R)
    end
    u = PIE.T*PIE.x;
    [~,info] = PIESOS_sublevel(PIE,info_loc.V,R^2-weighted_sobolev_norm(u,{[0,1]}),[],st);
    if ~info.feasible || info.c<R^2*(1-1e-3) || info.c>info_loc.C*R^2
        error("Test_PIESOS_sublevel: heat equation with feedback: c=%g, expected between %g and %g.",info.c,R^2*(1-1e-3),info_loc.C*R^2)
    end
end

disp("Test_PIESOS_sublevel passed.")

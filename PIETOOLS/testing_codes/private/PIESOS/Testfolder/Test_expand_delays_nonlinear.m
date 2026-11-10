% Test_expand_delays_nonlinear.m -- tests the conversion of PDEs with delayed
% factors in nonlinear terms, e.g. x(t)*x(t-tau) or x(t-tau)^3, by
% pde_struct/times_pdes, pde_struct/power and pde_struct/expand_delays.
%
% The PIE of such a PDE has the polynomial f(x) = (right-hand side of the
% PIE). At a state x of the PIE (random polynomials in the variables of each
% component), the value of f must equal the right-hand side of the PDE,
% evaluated with the value of each delayed factor x(t-tau), as given by the
% PIE of the linear PDE (d/dt)x = x(t-tau), at the same values of the
% history states, and of each non-delayed factor, given by the state itself
% or by T*x. So the products and powers of the PIE are checked against the
% linear conversion of the delay, and not against the converter of the
% products, which are tested for delayed factors only here. Also checked:
%   - delayed factors with the same delay share one history state, also in
%     linear terms (delays differing from each other do not);
%   - a PIE that does not change: values of f computed before the nonlinear
%     delay support was added (x(t-1) in a single term);
%   - delayed inputs in a product are still refused;
%   - the stability test of the PIESOS executive for Wright's equation
%       x'(t) = -a*x(t-1)*(1+x(t)),
%     which is locally stable iff a<pi/2: feasible for a=1.2 and infeasible
%     for a=1.7 (a certificate would be a bug). The stability tests can be
%     skipped by setting skip_executive=true before running the script.
%
% DJ (CC), 11/09/2026: Initial coding;

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'polyopvar')));    % eval_polyopvar
clear stateNameGenerator
pvar s t
tol = 1e-8;
ncases = 0;

% % Case 1: ODE with delayed factors, f compared with the linear functionals
% % ell_k = x(t-tau_k) of the linear PIEs, at states (c,h_1,...,h_m)
clear stateNameGenerator
x = pde_var('state',1,[],[]);
PIE_l1 = convert_pde(diff(x,t)==subs(x,t,t-1));
clear stateNameGenerator
x = pde_var('state',1,[],[]);
PIE_l2 = convert_pde(diff(x,t)==subs(x,t,t-2));
a = 0.8;
% name, PDE right-hand side (xd1=x(t-1), xd2=x(t-2)), number of states, the
% delays of the history states (in the order of their first occurrence in
% the PDE), and the expected right-hand side in terms of c=x(t), l1=x(t-1)
% and l2=x(t-2)
cases1 = {
    'Wright: -a*xd1-a*xd1*x',   @(x,xd1,xd2) -a*xd1-a*xd1*x,        2,  1,      @(c,l1,l2) -a*l1-a*l1*c;
    'factors reversed',         @(x,xd1,xd2) -a*xd1-a*x*xd1,        2,  1,      @(c,l1,l2) -a*l1-a*c*l1;
    'cubic: -x-xd1^3',          @(x,xd1,xd2) -x-xd1^3,              2,  1,      @(c,l1,l2) -c-l1^3;
    'square of a sum',          @(x,xd1,xd2) (x+xd1)^2,             2,  1,      @(c,l1,l2) (c+l1)^2;
    'repeated factors',         @(x,xd1,xd2) -x*xd1+0.5*xd1*xd1,    2,  1,      @(c,l1,l2) -c*l1+0.5*l1^2;
    'two delays: -x-xd1*xd2',   @(x,xd1,xd2) -x-xd1*xd2,            3,  [1,2],  @(c,l1,l2) -c-l1*l2;
    'two delays: x*xd2*xd1-xd2',@(x,xd1,xd2) x*xd2*xd1-xd2,         3,  [2,1],  @(c,l1,l2) c*l2*l1-l2};
for k=1:size(cases1,1)
    clear stateNameGenerator
    x = pde_var('state',1,[],[]);
    PIE = convert_pde(diff(x,t)==cases1{k,2}(x,subs(x,t,t-1),subs(x,t,t-2)));
    if numel(PIE.space.x.dims)~=cases1{k,3}
        error("Test_expand_delays_nonlinear: case '%s' has %d states, not %d.",cases1{k,1},numel(PIE.space.x.dims),cases1{k,3})
    end
    for seed=1:3
        xs = rand_state(PIE,seed);       % states (c,h_1,...), with the delays in cases1{k,4}
        c = double_val(xs{1});
        l1 = 0;     l2 = 0;
        if any(cases1{k,4}==1)
            l1 = lin_val(PIE_l1,{xs{1},xs{1+find(cases1{k,4}==1)}});
        end
        if any(cases1{k,4}==2)
            l2 = lin_val(PIE_l2,{xs{1},xs{1+find(cases1{k,4}==2)}});
        end
        Y = eval_polyopvar(PIE.f,xs);
        err = abs(double_val(Y(1))-cases1{k,5}(c,l1,l2));
        if err>tol*max(1,abs(c)+abs(l1)+abs(l2))^3
            error("Test_expand_delays_nonlinear: case '%s', seed %d: f differs from the expected value by %.2e.",cases1{k,1},seed,err)
        end
    end
    ncases = ncases + 1;
end

% % Case 2: PDE with a delayed factor, u_t = u_ss + 2u - u*u(t-1). The PIEs
% % of u_t = u_ss + 2u + c*u(t-1), c=1,2, have the same states and T as the
% % PIE of the nonlinear PDE, with f_c = f_ss + c*ell (ell the delayed
% % value, f_ss the rest), so that f = f_ss - u*ell, for ell = f_2 - f_1,
% % f_ss = 2*f_1 - f_2 and u=T*x. Also, ell must be the history state at the
% % lower end of its variable, the delay s_2 on [-1,1] (rescaled)
clear stateNameGenerator
dom = [0,1];
u = pde_var(s,dom);
PDE = [diff(u,t)==diff(u,s,2)+2*u-u*subs(u,t,t-1);    subs(u,s,dom(1))==0;    subs(u,s,dom(2))==0];
evalc('display(PDE)');          % the delayed nonlinear terms can be displayed
PIE = convert_pde(PDE);
PIE_c = cell(1,2);
for c=1:2
    clear stateNameGenerator
    u = pde_var(s,dom);
    PIE_c{c} = convert_pde([diff(u,t)==diff(u,s,2)+2*u+c*subs(u,t,t-1);    subs(u,s,dom(1))==0;    subs(u,s,dom(2))==0]);
end
if numel(PIE.space.x.dims)~=2 || numel(PIE_c{1}.space.x.dims)~=2 || numel(PIE_c{2}.space.x.dims)~=2
    error("Test_expand_delays_nonlinear: unexpected number of states in the PDE case.")
end
vnames = cellstr(PIE.space.vars(:));
s_delay = polynomial({vnames{2}});
for seed=1:3
    xs = rand_state(PIE,seed);
    Y = eval_polyopvar(PIE.f,xs);
    Y1 = eval_polyopvar(PIE_c{1}.f,xs);
    Y2 = eval_polyopvar(PIE_c{2}.f,xs);
    Uv = eval_polyopvar(PIE.T*PIE.x,xs);
    ell = Y2(1) - Y1(1);
    f_ss = 2*Y1(1) - Y2(1);
    ell_hist = subs(Uv(2),s_delay,-1);         % history of u at the lower end
    for ell_k = {ell,ell_hist}
        err = poly_max(Y(1) - (f_ss - Uv(1)*ell_k{1}));
        if err>tol*max(1,poly_max(Y(1)))
            error("Test_expand_delays_nonlinear: PDE case, seed %d: f differs from the expected value by %.2e.",seed,err)
        end
    end
end
ncases = ncases + 1;
% % Case 3: delayed factors with the same delay share a history state, in
% % linear terms too; a different delay has its own
clear stateNameGenerator
x = pde_var('state',1,[],[]);
xd = subs(x,t,t-1);
PIE_rep = convert_pde(diff(x,t)==-x+0.2*xd+0.1*xd);
clear stateNameGenerator
x = pde_var('state',1,[],[]);
PIE_one = convert_pde(diff(x,t)==-x+0.3*subs(x,t,t-1));
clear stateNameGenerator
x = pde_var('state',1,[],[]);
PIE_two = convert_pde(diff(x,t)==-x+0.2*subs(x,t,t-1)+0.3*subs(x,t,t-2)+0.1*subs(x,t,t-1));
if numel(PIE_rep.space.x.dims)~=2 || numel(PIE_one.space.x.dims)~=2 || numel(PIE_two.space.x.dims)~=3
    error("Test_expand_delays_nonlinear: repeated delays do not share a history state (%d, %d, %d states, not 2, 2, 3).",...
            numel(PIE_rep.space.x.dims),numel(PIE_one.space.x.dims),numel(PIE_two.space.x.dims))
end
for seed=1:3
    xs = rand_state(PIE_rep,seed);
    Y_rep = eval_polyopvar(PIE_rep.f,xs);       Y_one = eval_polyopvar(PIE_one.f,xs);
    if abs(double_val(Y_rep(1))-double_val(Y_one(1)))>tol
        error("Test_expand_delays_nonlinear: repeated linear delay, seed %d: f differs from that of the combined term.",seed)
    end
    xs3 = rand_state(PIE_two,seed);
    Y_two = eval_polyopvar(PIE_two.f,xs3);
    l1 = lin_val(PIE_l1,{xs3{1},xs3{2}});       l2 = lin_val(PIE_l2,{xs3{1},xs3{3}});
    if abs(double_val(Y_two(1))-(-double_val(xs3{1})+0.3*l1+0.3*l2))>tol*max(1,abs(l1)+abs(l2))
        error("Test_expand_delays_nonlinear: linear terms with two delays, seed %d: f differs from the expected value.",seed)
    end
end
ncases = ncases + 1;

% % Case 4: a PIE that has not changed: x'=-x+0.5*x(t-1)-x^3 (values of f
% % at the states below, computed before the support for nonlinear delays)
clear stateNameGenerator
x = pde_var('state',1,[],[]);
PIE = convert_pde(diff(x,t)==-x+0.5*subs(x,t,t-1)-x^3);
f_ref = [-0.2134113034,2.670880519];
for seed=1:2
    rng(seed);
    xs = {randn(1), polynomial(randn(1,4)*[1;pvar('s1');pvar('s1')^2;pvar('s1')^3])};
    Y = eval_polyopvar(PIE.f,xs);
    if abs(double_val(Y(1))-f_ref(seed))>1e-8
        error("Test_expand_delays_nonlinear: the PIE of x'=-x+0.5*x(t-1)-x^3 has changed (seed %d).",seed)
    end
end
ncases = ncases + 1;

% % Case 5: delayed inputs cannot appear in a product (inputs only appear
% % linearly in a PIE)
clear stateNameGenerator
x = pde_var('state',1,[],[]);
w = pde_var('in',1);
refused = false;
try
    PDE = diff(x,t)==-x*subs(w,t,t-1);
    convert(PDE,'silent');
catch ME
    refused = contains(ME.message,"Multiplication with input or output signals");
end
if ~refused
    error("Test_expand_delays_nonlinear: no error for a delayed input in a product.")
end
ncases = ncases + 1;

% % Case 6: stability test of the executive, Wright's equation
if ~exist('skip_executive','var') || ~skip_executive
    stts = PIESOS_settings();       stts.verbose = false;
    for a = [1.2,1.7]
        clear stateNameGenerator
        x = pde_var('state',1,[],[]);
        xd = subs(x,t,t-1);
        PIE = convert_pde(diff(x,t)==-a*xd-a*xd*x);
        [~,info] = PIESOS_stability_local(PIE,0.1,1,0,'pie',stts);
        if info.feasible ~= (a<pi/2)
            error("Test_expand_delays_nonlinear: Wright's equation with a=%g: local stability test feasible=%d, expected %d.",a,info.feasible,a<pi/2)
        end
        ncases = ncases + 1;
    end
end

fprintf('Test_expand_delays_nonlinear passed (%d cases).\n',ncases)



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function PIE = convert_pde(PDE)
% PIE = CONVERT_PDE(PDE) converts the PDE to a PIE in the nonlinear format,
% without display.
%
% INPUT
% - PDE:    'pde_struct' object;
%
% OUTPUT
% - PIE:    'npie_struct' object;

evalc('PIE = convert(initialize(PDE),{''silent'',''nonlinear''});');
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function xs = rand_state(PIE,seed)
% XS = RAND_STATE(PIE,SEED) returns random polynomial values of the
% components of the (scalar) state of the PIE, each a polynomial in the
% variables of the component, with nonzero coefficients for each variable.
%
% INPUT
% - PIE:    'npie_struct' object with scalar state components;
% - seed:   integer, the seed of the random number generator;
%
% OUTPUT
% - xs:     1 x n cell with a 'polynomial' for each state component;

rng(seed);
sp = PIE.space.x;
vnames = cellstr(PIE.space.vars(:));
xs = cell(1,numel(sp.dims));
for i=1:numel(xs)
    vi = vnames(sp.depmat(i,:));
    mons = polynomial(1);
    for k=1:numel(vi)
        pk = polynomial({vi{k}});
        mons = [mons; pk; pk^2];                                                %#ok<AGROW>
    end
    if numel(vi)==2
        mons = [mons; polynomial({vi{1}})*polynomial({vi{2}})];
    end
    xs{i} = randn(1,numel(mons))*mons;
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function v = lin_val(PIE_lin,xs)
% V = LIN_VAL(PIE_LIN,XS) returns the value of the delayed state x(t-tau)
% as a function of the state, given by the right-hand side f of the PIE of
% the PDE x'(t) = x(t-tau), at the state values xs.
%
% INPUT
% - PIE_lin:    'npie_struct' of an ODE with a single delayed term;
% - xs:         1 x 2 cell with the values of the ODE state and of its
%               history;
%
% OUTPUT
% - v:          scalar, the value;

Y = eval_polyopvar(PIE_lin.f,xs);
v = double_val(Y(1));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function v = double_val(P)
% V = DOUBLE_VAL(P) returns the value of a constant polynomial, as double.
%
% INPUT
% - P:      constant 'polynomial' or double;
%
% OUTPUT
% - v:      scalar, the value;

v = double(P);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function m = poly_max(P)
% M = POLY_MAX(P) returns the largest absolute value of the coefficients of
% a polynomial.
%
% INPUT
% - P:      'polynomial' object;
%
% OUTPUT
% - m:      scalar, the maximal absolute coefficient;

m = max([abs(full(P.coeff(:)));0]);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

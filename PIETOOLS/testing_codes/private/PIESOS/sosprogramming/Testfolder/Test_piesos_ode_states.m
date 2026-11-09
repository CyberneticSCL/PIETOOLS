function Test_piesos_ode_states()
% TEST_PIESOS_ODE_STATES tests, by solving small programs, the support of
% states without spatial variables (ODE states) in PIESOS programs, alone and
% together with PDE states, including products of ODE states in the monomial
% bases and in the polynomials. The programs have a known answer, so that both
% sides of the threshold are tested:
%   (1) product of ODE states, as functionals: the product of the functionals
%       (y'*y)^2 evaluated at a given y, for a vector ODE state y, vs. the value
%       computed directly;
%   (2) positivity of the polynomial (y'*y)^2-2*y'*y+c of a vector ODE state y
%       in R^2, with c>=1 (verified) and c<1 (not verified): the basis [1;y o y];
%   (3) the 2D ODE x' = [-a w;-w -a]*x-|x|^2*x in two scalar states, with
%       V=|x|^2, dV=-2*a*|x|^2-2*|x|^4: verified iff a>=0 (products x1*x2);
%   (4) the delay ODE x' = -alp*x+bet*x(t-tau)-x^3, with a transport state
%       y(t,s) = x(t-s), V=x^2+|bet|*int y^2: verified iff alp>=|bet|, for
%       several delays (ODE state times PDE state);
%   (5) the cascade of Burgers' equation u_t=u_ss-u*u_s and the ODE
%       x' = -a*x+int u-x^3, with V=x^2+c*<u,u>: dV = -2*a*x^2+2*x*int u
%       -2*c*||u_s||^2-2*x^4, which is negative iff a*c>=1/48 (as
%       sup (int u)^2/||u_s||^2 = 1/12); and with the free functional of
%       piesos_LF, for which a<0 is not verified.
% A program is verified if it is feasible: pinf, dinf and numerr are zero
% and the feasibility ratio is not far from 1.
%
% DJ (CC), 11/08/2026: Initial coding;

tests = {@case_functional_product, @case_vector_ode_positivity, @case_two_ode_states, ...
         @case_delay_ode, @case_ode_pde_cascade, @case_ode_pde_lf};
nfail = 0;
for i=1:numel(tests)
    tstart = tic;
    try
        tests{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(tests{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_piesos_ode_states: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_ode_states: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function ok = is_verified(prog)
% Solve the program, and check that it is feasible
so.simplify = true;
so.params.fid = 0;
ps = piesos_solve(prog,so);
i = ps.solinfo.info;
ok = ~(i.pinf || i.dinf || i.numerr || abs(i.feasratio-1)>0.5);
end

function ok = verify_neg(PIE,Pw,degs,psatz)
% Verify <w,Pw*f> <= 0 for the PDE states w=T*x: the derivative of V=<w,Pw*w>
prog = piesos_program(PIE.x);
w = PIE.T*PIE.x;
dV = 2*innerprod(w,Pw,PIE.f);
o.degs = degs;      o.psatz = psatz;
prog = piesos_ineq(prog,-dV,o);
ok = is_verified(prog);
end

function check(ok,expected,label)
assert(ok==expected,'%s: verified=%d, expected %d.',label,ok,expected);
end


%% ---- cases ------------------------------------------------------------
function case_functional_product()
% (y'*y)^2 for y in R^2 and (x'*x)^2 for scalar x, from the product of functionals
% without spatial variables, vs. the value at given states
yo = polyopvar({'yo'},2,{'s1'},[0,1],false);
q = yo'*yo;
q2 = opower(q,2);
y0 = [1.1;-0.7];
val = apply_poly(q2,{polynomial(y0)});
assert(abs(double(val)-(y0'*y0)^2)<1e-10,'(y''*y)^2 differs from the direct value.');
assert(isequal(full(q2.Z),[1 1 1 1]),'(y''*y)^2 should be a monomial of degree 4.');
xo = polyopvar({'xo'},1,{'s1'},[0,1],false);
val = apply_poly(opower(xo'*xo,2),{polynomial(1.3)});
assert(abs(double(val)-1.3^4)<1e-10,'(x''*x)^2 differs from the direct value.');
end

function case_vector_ode_positivity()
% (y'*y)^2-2*y'*y+c = (y'*y-1)^2+c-1, with y in R^2: nonnegative iff c>=1
yo = polyopvar({'yo'},2,{'s1'},[0,1],false);
q = yo'*yo;
for c = [1.01,1,0.99,0.9]
    prog = piesos_program(yo);
    o.degs = 0;
    prog = piesos_ineq(prog,opower(q,2)-2*q+c,o);
    check(is_verified(prog),c>=1,sprintf('(y''y)^2-2y''y+%g',c));
end
end

function case_two_ode_states()
% Two scalar ODE states: V=|x|^2 verifies stability iff a>=0
w = 3;
for a = [0.5,0.001,0,-0.001,-0.05]
    clear stateNameGenerator
    pvar t
    x1 = pde_var('state',1,[],[]);      x2 = pde_var('state',1,[],[]);
    PDE = [diff(x1,t)==-a*x1+w*x2-x1*(x1^2+x2^2);
           diff(x2,t)==-w*x1-a*x2-x2*(x1^2+x2^2)];
    PIE = convert(PDE,'silent');
    check(verify_neg(PIE,eye(2),0,0),a>=0,sprintf('two ODE states, a=%g',a));
end
end

function case_delay_ode()
% x' = -alp*x+bet*x(t-tau)-x^3, transport state z(t,s)=x(t-s), s in [0,tau]:
% verified iff alp>=|bet|, for any tau
cfgs = [1 0.95 1; 1 1 1; 1 1.05 1; 2 1.9 3; 2 1.9 0.2; 1 -1 1; 1 -1.05 1; 1 -0.5 1];
for k=1:size(cfgs,1)
    alp = cfgs(k,1);    bet = cfgs(k,2);    tau = cfgs(k,3);
    clear stateNameGenerator
    pvar t s1
    x = pde_var('state',1,[],[]);
    z = pde_var(s1,[0,tau]);
    PDE = [diff(x,t)==-alp*x+bet*subs(z,s1,tau)-x^3;
           diff(z,t)==-diff(z,s1);
           subs(z,s1,0)==x];
    PIE = convert(PDE,'silent');
    check(verify_neg(PIE,diag([1,abs(bet)]),2,1),alp>=abs(bet),sprintf('delay ODE, (alp,bet,tau)=(%g,%g,%g)',alp,bet,tau));
end
end

function case_ode_pde_cascade()
% x' = -a*x+int u-x^3, u_t=u_ss-u*u_s, V=x^2+c<u,u>: verified iff a*c>=1/48
c = 1;
for a = [0.1,0.0215,0.0203,0.01]
    clear stateNameGenerator
    pvar t s
    x = pde_var('state',1,[],[]);
    u = pde_var(s,[0,1]);
    PDE = [diff(x,t)==-a*x+int(u,s,[0,1])-x^3;
           diff(u,t)==diff(u,s,2)-u*diff(u,s);
           subs(u,s,0)==0;  subs(u,s,1)==0];
    PIE = convert(PDE,'silent');
    check(verify_neg(PIE,diag([1,c]),2,1),a*c>=1/48,sprintf('ODE-PDE cascade, a=%g',a));
end
end

function case_ode_pde_lf()
% The free quadratic functional of piesos_LF in the ODE and PDE states: verified for
% a>0 where the cascade is stable (more than with c=1 fixed), and not for a<0
for a = [0.1,0.005,-0.01]
    clear stateNameGenerator
    pvar t s
    x = pde_var('state',1,[],[]);
    u = pde_var(s,[0,1]);
    PDE = [diff(x,t)==-a*x+int(u,s,[0,1])-x^3;
           diff(u,t)==diff(u,s,2)-u*diff(u,s);
           subs(u,s,0)==0;  subs(u,s,1)==0];
    PIE = convert(PDE,'silent');
    prog = piesos_program(PIE.x);
    Popts.degs = 2;     Popts.exclude = [0;1;1];    Popts.pos = true;
    [prog,~,dV] = piesos_LF(prog,PIE,'pde',1,Popts);
    w = PIE.T*PIE.x;
    dV = dV + 2e-3*innerprod(w,PIE.f);
    o.degs = 2;     o.psatz = 1;
    prog = piesos_ineq(prog,-dV,o);
    check(is_verified(prog),a>0,sprintf('piesos_LF, ODE-PDE cascade, a=%g',a));
end
end

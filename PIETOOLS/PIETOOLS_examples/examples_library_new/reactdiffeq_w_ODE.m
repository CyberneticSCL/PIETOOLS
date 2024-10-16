% % ODE        xo_{t} = A * xo + Bxr * x_{s}(s=a)
% % PDE        x_{t} = lam * x + x_{ss} + Bpv * xo
% % with BCs   x(s=a) = 0
% %            x(s=b) = 0
function obj = reactdiffeq_w_ODE()
pvar t s1;
odevar x; pdevar X;
X.len = 4;

A = [-1.2142,  1.9649,  0.2232,  0.5616;
    -1.8042, -0.7260, -0.3479,  5.4355;
    -0.2898,  0.7381, -1.7606,  0.8294;
    -0.9417, -5.3399, -1.0704, -0.7590];
Bxr = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];
Bpv = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
lam = pi^2-1;
eqns = [diff(x,t)==A*x+Bxr*diff(X,s1,0); 
             diff(X,t)==lam*X+diff(X,s1,2)+Bpv*x
             subs(X,s1,0)==0;
             subs(X,s1,1)==0];
obj = sys('pde',eqns);
end
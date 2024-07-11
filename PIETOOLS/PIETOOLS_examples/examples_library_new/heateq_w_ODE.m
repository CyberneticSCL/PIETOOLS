% % ODE        xo_{t} = k * xo
% % PDE        x_{t} = x_{ss} 
% % with BCs   x(s=0) = 0
% %            x(s=1) = xo
function obj = heateq_w_ODE()
pvar t s1;
odevar x; pdevar X;
k=-1;
eqns = [diff(x,t)==k*x; 
             diff(X,t)==diff(X,s1,2)
             subs(X,s1,0)==0;
             subs(X,s1,1)==x];
obj = sys('pde',eqns);
end
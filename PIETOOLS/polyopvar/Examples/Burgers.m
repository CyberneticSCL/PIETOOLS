

pvar s t
r = 9.0;
u = pde_var(s,[0,1],2);
pde1 = u;
pde2 = diff(u,s,1);
pde2.free = {diff(u,s,1)}; % class(pde2.free)



%diff(u,s,2) + r*u - u*diff(u,s,1);



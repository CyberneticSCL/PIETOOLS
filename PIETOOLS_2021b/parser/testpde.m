x = state('ode',2); X = state('pde',2); w = state('in',2); u = state('in',2);
z = state('out',2); y = state('out',2);
pvar s theta t;
ssobj = sys(); ssobj.type = 'pde';

eqns = [diff(x,t)+5*x-3*u; % dot x = -5 x + 3 u
        diff(X,t)-3*X+diff(X,s,2)+w;
        z-x-int(X,s,[0,1])-u;
        y-diff(subs(X,s,0),s)-w;
        subs(X,s,0);
        subs(diff(X,s),s,1)-x];
ssobj = addequation(ssobj,eqns);
ssobj = setControl(ssobj,u);
ssobj = setObserve(ssobj,y);

ssobj = getParams(ssobj);
%%
initialize_PIETOOLS_PDE_terms(ssobj.params);
display_PDE(ssobj.params);

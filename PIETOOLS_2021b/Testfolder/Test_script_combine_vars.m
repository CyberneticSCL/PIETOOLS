

nvars = 8;

pvar aa bb cc dd ee ff gg hh aa2 bb2 cc2 dd2 ee2 ff2 gg2 hh2
old_vars = [aa,aa2; bb,bb2; cc,cc2; dd,dd2; ee,ee2; ff,ff2; gg,gg2; hh,hh2];
nvars = size(old_vars);

old_dom = randi(5,[8,1])-3;
old_dom = [old_dom,old_dom + randi(4,[8,1])];

dep_tab = logical(randi([0,1],5,8))
dep_tab_old = dep_tab;



pvar x y z q r tt nu mu
pvar tau1 tau2 tau3
a = 0;      b = 1;
c = 3;      d = 5;
e = -1;     f = 2;
g = 0;      h = 1;
i = -1;     j = 1;
PDE = pde_struct();
PDE.tau = [tau1,1; tau2, 3];

PDE.x{1}.vars = x;
PDE.x{1}.dom = [a,b];

PDE.x{2}.vars = [x,nu;z,mu];
PDE.x{2}.dom = [a,b;e,f];

PDE.w{1}.vars = y;
PDE.w{1}.dom = [c,d];

PDE.u{1}.vars = q;
PDE.u{1}.dom = [g,h];

PDE.y{1}.vars = r;
PDE.y{1}.dom = [a,b];

PDE.x{1}.term{1}.x = 1;
PDE.x{1}.term{1}.D = [2;0];
PDE.x{1}.term{1}.delay = tau1;

PDE.x{1}.term{2}.w = 1;
PDE.x{1}.term{2}.I{1} = [c,d];
PDE.x{1}.term{2}.C = x*y;

PDE.x{2}.term{1}.x = 2;
PDE.x{2}.term{1}.D = [0,1];
PDE.x{2}.term{1}.I{1} = [a,x];
PDE.x{2}.term{1}.C = (x-nu)*tau2;
PDE.x{2}.term{1}.delay = tau2;

PDE.x{2}.term{2}.u = 1;
PDE.x{2}.term{2}.I{1} = [g,h];
PDE.x{2}.term{2}.delay = tau1;

PDE.y{1}.term{1}.x = 2;
PDE.y{1}.term{1}.I{1} = [a,b];
PDE.y{1}.term{1}.I{2} = [e,f];
PDE.y{1}.term{1}.C = r;

PDE.BC{1}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = a;
PDE.BC{2}.term{1}.x = 1;
PDE.BC{2}.term{1}.I{1} = [a,b];

PDE.BC{3}.term{1}.x = 2;
PDE.BC{3}.term{1}.loc = [x,f];
PDE.BC{3}.term{1}.delay = 2;

PDE = initialize(PDE);


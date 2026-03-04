

pvar s1 s1_dum
dom = [0,1];

Aop = rand_ndopvar([1,1],2,[0,1],s1,s1_dum);
Bop = rand_ndopvar([1,1],2,[0,1],s1,s1_dum);
Cop = rand_ndopvar([1,1],2,[0,1],s1,s1_dum);
Dop = rand_ndopvar([1,1],2,[0,1],s1,s1_dum);

C1 = tensopvar();   C2 = tensopvar();
C1.ops{1} = Aop;
C1.ops{2} = {Aop,Bop};
C2.ops{1} = Cop;
C2.ops{2} = {Cop,Dop};

p1 = polyopvar();
p1.varname = {'x';'y'};
p1.pvarname = {'s1'};
p1.varmat = [1;1];
p1.dom = dom;
p1.degmat = [1,0;1,1];
p1.C = C1;

p2 = polyopvar();
p2.varname = {'y';'x'};
p2.pvarname = {'s1'};
p2.varmat = [1;1];
p2.dom = dom;
p2.degmat = [1,0;1,1];
p2.C = C2;

p3 = p1+p2;
C3 = p3.C;


p4 = p3+p2;
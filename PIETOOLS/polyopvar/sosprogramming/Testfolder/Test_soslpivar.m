clear
pvar s1
dom = [0,1];

% Declare the monomial basis in distributed state x1
Z = polyopvar();
Z.varname = {'x1'};
Z.pvarname = s1.varname;
Z.dom = dom;
Z.degmat = [1;3];

% Declare an LPI optimization program
prog = lpiprogram(s1,dom);

% Declare the positive operator variable
pdeg = 1;
[prog,Pmat,Zop] = soslpivar(prog,Z,pdeg);

%Zop1 = ndopvar2dopvar(Zop{1});
%Zop21 = ndopvar2dopvar(Zop{2}{1});
%Zop22 = ndopvar2dopvar(Zop{2}{2});

Zop31 = ndopvar2dopvar(Zop.ops{2}{1});
Zop32 = ndopvar2dopvar(Zop.ops{2}{2});
Zop33 = ndopvar2dopvar(Zop.ops{2}{3});
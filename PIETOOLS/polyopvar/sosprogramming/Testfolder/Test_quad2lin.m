clear
pvar s1 s1_dum s1_dum_1 s1_dum_2 s1_dum_3 s1_dum_4
dom = [0,1];
nvars = 1;

% Declare the monomial basis in distributed state x1
Zx = polyopvar();
Zx.varname = {'x1'};
Zx.pvarname = s1.varname;
Zx.dom = dom;
Zx.degmat = [1;2];
Zx.varmat = 1;

% Declare an LPI optimization program
prog = lpiprogram(s1,dom);

% Declare the positive operator variable
pdeg = 1;
[prog,Pmat,Zop] = soslpivar(prog,Zx,pdeg);

% Construct the associated polynomial functional
Vx = quad2lin(Pmat,Zop,Zx);

% Convert Zop to opvar objects
Zop1 = ndopvar2dopvar(Zop.ops{1});
Zop2 = cell(1,2);
Zop2{1} = ndopvar2dopvar(Zop.ops{2,2}{1});
Zop2{2} = ndopvar2dopvar(Zop.ops{2,2}{2});

% Generate random functions xi(s)
x_tst = rand_poly([1,1],s1,3);
z1_tst = apply_opvar(Zop1,x_tst);
z2_tst = apply_opvar(Zop2{1},x_tst).*apply_opvar(Zop2{2},x_tst);
P11 = Pmat(1:numel(z1_tst),1:numel(z1_tst));
P12 = Pmat(1:numel(z1_tst),numel(z1_tst)+1:end);
P21 = Pmat(numel(z1_tst)+1:end,1:numel(z1_tst));
P22 = Pmat(numel(z1_tst)+1:end,numel(z1_tst)+1:end);

fval11 = int(z1_tst'*P11*z1_tst,s1,dom(1),dom(2));
fval12 = int(z1_tst'*P12*z2_tst,s1,dom(1),dom(2));
fval21 = int(z2_tst'*P21*z1_tst,s1,dom(1),dom(2));
fval22 = int(z2_tst'*P22*z2_tst,s1,dom(1),dom(2));

K1 = Vx.C.ops{1}.params;    
idx1 = Vx.C.ops{1}.omat;
vars1 = Vx.C.ops{1}.vars;
fval_alt11 = apply_functional(Vx.C.ops{1},x_tst,Vx.degmat(1,:));
K2 = Vx.C.ops{2}.params;    
idx2 = Vx.C.ops{2}.omat;
vars2 = Vx.C.ops{2}.vars;
fval_alt12 = apply_functional(Vx.C.ops{2},x_tst,Vx.degmat(2,:));
K3 = Vx.C.ops{3}.params;    
idx3 = Vx.C.ops{3}.omat;
vars3 = Vx.C.ops{3}.vars;
fval_alt22 = apply_functional(Vx.C.ops{3},x_tst,Vx.degmat(3,:));
  
err_11 = fval11 - fval_alt11;
err_12 = fval12 + fval21 - fval_alt12;
err_22 = fval22 - fval_alt22;

max(max(abs(err_11.C)))
max(max(abs(err_12.C)))
max(max(abs(err_22.C)))
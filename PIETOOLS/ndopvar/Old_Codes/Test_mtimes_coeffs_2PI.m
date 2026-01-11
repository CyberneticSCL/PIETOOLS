
pvar s1 s1_dum
dom = [-1,1];
d = 5;
m = 2;
q = 4;
n = 3;

% Generate random opvar objects
Qop = rand_opvar([0,0;m,q],d,s1,s1_dum,dom);
Qop.R.R0 = zeros(m,q);
Rop = rand_opvar([0,0;q,n],d,s1,s1_dum,dom);
Rop.R.R0 = zeros(q,n);

% Compute coefficients representing kernels in quadratic form
C = opvar2coeffs(Qop,d);
D = opvar2coeffs(Rop,d);

Qop_alt = coeffs2opvar(C,[m,q],dom,s1,s1_dum);
Rop_alt = coeffs2opvar(D,[q,n],dom,s1,s1_dum);

% Compute coefficients defining the composition
C_str = struct();
C_str.C = C;    C_str.dim = [m,q];  C_str.dom = dom;
D_str = struct();
D_str.C = D;    D_str.dim = [q,n];  D_str.dom = dom;
tic
B = mtimes_coeffs(C_str,D_str);
t_coeffs = toc;

% Also compute composition using standard routing
tic
Pop_alt = Qop*Rop;
t_opvar = toc;

% Check that composition is correct
Pop = coeffs2opvar(B,[m,n],dom,s1,s1_dum);
Pop_err = clean_opvar(Pop-Pop_alt,1e-10);

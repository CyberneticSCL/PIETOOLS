
pvar s1 s1_dum
dom = [-1,1];
d = 5;

% Generate random opvar objects
Qop = rand_opvar([0,0;1,1],d,s1,s1_dum,dom);
Qop.R.R0 = 0;
Rop = rand_opvar([0,0;1,1],d,s1,s1_dum,dom);
Rop.R.R0 = 0;

% Compute coefficients representing kernels in quadratic form
C = opvar2coeffs(Qop,d);
D = opvar2coeffs(Rop,d);

% Compute coefficients defining the composition
tic
B = mtimes_coeffs(C,D,dom);
t_coeffs = toc;

% Also compute composition using standard routing
tic
Pop_alt = Qop*Rop;
t_opvar = toc;

% Check that composition is correct
Pop = coeffs2opvar(B,dom,s1,s1_dum);
Pop_err = clean_opvar(Pop-Pop_alt,1e-10);

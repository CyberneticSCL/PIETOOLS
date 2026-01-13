
pvar s1 s1_dum
dom = [-1,1];
d = 4;
m = 2;
q = 4;
n = 3;

% Generate random opvar objects
Qop = rand_opvar([0,0;m,q],d,s1,s1_dum,dom);
%Qop.R.R0 = zeros(m,q);
Rop = rand_opvar([0,0;q,n],d,s1,s1_dum,dom);
%Rop.R.R0 = zeros(q,n);

% Compute coefficients representing kernels in quadratic form
C_str = opvar2coeffs(Qop,d);
D_str = opvar2coeffs(Rop,d);

Qop_alt = coeffs2opvar(C_str);
Rop_alt = coeffs2opvar(D_str);

% Compute coefficients defining the composition
tic
B_str = mtimes_coeffs(C_str,D_str);
t_coeffs = toc;

% Also compute composition using standard routing
tic
Pop_alt = Qop*Rop;
t_opvar = toc;

% Check that composition is correct
Pop = coeffs2opvar(B_str);
Pop_err = clean_opvar(Pop-Pop_alt,1e-10);

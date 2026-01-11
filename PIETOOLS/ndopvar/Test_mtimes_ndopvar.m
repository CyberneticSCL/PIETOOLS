
pvar s1 s1_dum
dom = [-1,1];
d = 5;
m = 1;
q = 1;
n = 1;

% Generate random opvar objects
Qop = rand_opvar([0,0;m,q],d,s1,s1_dum,dom);
%Qop.R.R0 = zeros(m,q);
Rop = rand_opvar([0,0;q,n],d,s1,s1_dum,dom);
%Rop.R.R0 = zeros(q,n);

% Compute coefficients representing kernels in quadratic form
Qop_nd = opvar2ndopvar(Qop,d);
Rop_nd = opvar2ndopvar(Rop,d);

Qop_alt = ndopvar2opvar(Qop_nd);
Rop_alt = ndopvar2opvar(Rop_nd);

% Compute coefficients defining the composition
tic
Pop_nd = Qop_nd*Rop_nd;
t_coeffs = toc;

% Also compute composition using standard routing
tic
Pop_alt = Qop*Rop;
t_opvar = toc;

% Check that composition is correct
Pop = ndopvar2opvar(Pop_nd);
Pop_err = clean_opvar(Pop-Pop_alt,1e-10);

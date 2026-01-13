
pvar s1 s1_dum s2 s2_dum
dom = [0,1;-1,1];
d = 5;
m = 2;
q = 4;
n = 3;

% Generate random opvar objects
Qop = opvar2d();                Rop = opvar2d();
Qop.I = dom;                    Rop.I = dom;
Qop.dim = [0,0;0,0;0,0;m,q];    Rop.dim = [0,0;0,0;0,0;q,n];
Qop.R22{1,1} = rand_poly([m,q],[s1;s2],d);
Rop.R22{1,1} = rand_poly([q,n],[s1;s2],d);
for ii=2:3
    Qop.R22{ii,1} = rand_poly([m,q],[s1;s2;s1_dum],d);
    Rop.R22{ii,1} = rand_poly([q,n],[s1;s2;s1_dum],d);
    Qop.R22{1,ii} = rand_poly([m,q],[s1;s2;s2_dum],d);
    Rop.R22{1,ii} = rand_poly([q,n],[s1;s2;s2_dum],d);
    for jj=2:3
        Qop.R22{ii,jj} = rand_poly([m,q],[s1;s2;s1_dum;s2_dum],d);
        Rop.R22{ii,jj} = rand_poly([q,n],[s1;s2;s1_dum;s2_dum],d);
    end
end

% Qop = rand_opvar([0,0;m,q],d,s1,s1_dum,dom);
% %Qop.R.R0 = zeros(m,q);
% Rop = rand_opvar([0,0;q,n],d,s1,s1_dum,dom);
% %Rop.R.R0 = zeros(q,n);

% Compute coefficients representing kernels in quadratic form
C_str = opvar2d2coeffs(Qop,d);
D_str = opvar2d2coeffs(Rop,d);

Qop_alt = coeffs2opvar2d(C_str);
Rop_alt = coeffs2opvar2d(D_str);

% Compute coefficients defining the composition
tic
B_str = mtimes_coeffs2d(C_str,D_str);
t_coeffs = toc;

% Also compute composition using standard routing
tic
Pop_alt = Qop*Rop;
t_opvar = toc;

% Check that composition is correct
Pop = coeffs2opvar2d(B_str);
Pop_err = clean_opvar(Pop-Pop_alt,1e-10);

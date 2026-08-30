
% This script is to test that composition of 'nopvar' or 'ndopvar' objects
% in 1D produces the same result as using the 'opvar' and 'dopvar'
% composition functions.
clear
pvar s1 s1_dum
dom = [-1,1];
d = 4;
m = 2;
q = 2;
n = 2;

use_dpvar = true;

% Generate random opvar/dopvar objects
if use_dpvar
    %opts.exclude = [1,0,0,0];
    prog = lpiprogram([s1,s1_dum],dom);
    [prog,Qop] = lpivar(prog,[0,0;m,q],d);
else
    Qop = rand_opvar([0,0;m,q],d,s1,s1_dum,dom);
end
%Qop.R.R0 = zeros(m,q);
Rop = rand_opvar([0,0;q,n],d,s1,s1_dum,dom);
%Rop.R.R0 = zeros(q,n);

% Compute coefficients representing kernels in quadratic form
Qop_nd = dopvar2ndopvar(Qop,d);
Rop_nd = dopvar2ndopvar(Rop,d);

Qop_alt = ndopvar2dopvar(Qop_nd);
Rop_alt = ndopvar2dopvar(Rop_nd);

% Compute coefficients defining the composition
tic
Pop_nd = Qop_nd*Rop_nd;
t_coeffs = toc;

% Also compute composition using standard routing
tic
Pop_alt = Qop*Rop;
t_opvar = toc;

% Check that composition is correct
if use_dpvar
    Pop_alt_nd = dopvar2ndopvar(Pop_alt,Pop_nd.deg,Pop_nd.dvarname);
else
    Pop_alt_nd = dopvar2ndopvar(Pop_alt,Pop_nd.deg);
end
Pop_nd_err = Pop_alt_nd-Pop_nd;
C_err = zeros(size(Pop_nd_err.C));
for ii=1:numel(Pop_nd_err.C)
    C_err(ii) = max(max(abs(Pop_nd_err.C{ii})));
end
C_err
%Pop = ndopvar2dopvar(Pop_nd);
%Pop_err = clean_opvar(Pop-Pop_alt,1e-10);

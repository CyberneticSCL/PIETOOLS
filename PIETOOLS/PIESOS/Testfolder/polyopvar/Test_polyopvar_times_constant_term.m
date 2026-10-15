% Test_polyopvar_times_constant_term.m -- product of distributed
% polynomials, one of which has a constant (degree-0) term, e.g.
% lam*(R^2-<u,u>). The constant is stored as a tensopvar without
% variables, so the product needs otimes(intop,tensopvar) and
% otimes(tensopvar,intop).
% The value of the product at polynomial states is compared with that
% computed by direct integration (ground truth independent of the product
% code).
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-9;
pvar s
dom = [0,1];
ncases = 0;

% States a, b (functions of s), functionals built from innerprod
x = polyopvar({'a';'b'},[1;1],{'s'},dom,true(2,1));
a_x = polyopvar({'a'},1,{'s'},dom,true);
b_x = polyopvar({'b'},1,{'s'},dom,true);
rng(7);
for icase = 1:3
    a = rand_poly([1,1],s,2,5);         b = rand_poly([1,1],s,2,5);
    vals = struct('a',a,'b',b);
    R2 = 3+icase;
    c1 = 0.5*icase;
    switch icase
        case 1
            lam = innerprod(a_x,a_x);                    lam_val = int(a^2,s,0,1);
            ipg = innerprod(a_x,a_x);                    ipg_val = int(a^2,s,0,1);
        case 2
            lam = innerprod(a_x,b_x)+2*innerprod(b_x,b_x);  lam_val = int(a*b,s,0,1)+2*int(b^2,s,0,1);
            ipg = innerprod(a_x,a_x);                    ipg_val = int(a^2,s,0,1);
        case 3
            lam = c1*innerprod(b_x,b_x);                 lam_val = c1*int(b^2,s,0,1);
            ipg = innerprod(a_x,b_x);                    ipg_val = int(a*b,s,0,1);
    end
    g = R2-ipg;                                          g_val = R2-ipg_val;
    expected = double(lam_val*g_val);
    if abs(expected)<1e-6
        error("Test_polyopvar_times_constant_term (case %d): value is zero, the test is vacuous.",icase)
    end
    % lam*g, g*lam, lam.*g, g.*lam
    prods = {lam*g, g*lam, lam.*g, g.*lam};
    for k=1:numel(prods)
        val = double(eval_V(prods{k},vals));
        if abs(val-expected)>tol*max(1,abs(expected))
            error("Test_polyopvar_times_constant_term failed (case %d, product %d): value %g, expected %g.",icase,k,val,expected)
        end
        ncases = ncases+1;
    end
end

% Product of two functionals which both have constant terms: (R2-<a,a>)*(R3-<b,b>)
a = rand_poly([1,1],s,2,5);             b = rand_poly([1,1],s,2,5);
vals = struct('a',a,'b',b);
g1 = 2-innerprod(a_x,a_x);              g2 = 5-innerprod(b_x,b_x);
expected = double((2-int(a^2,s,0,1))*(5-int(b^2,s,0,1)));
val = double(eval_V(g1*g2,vals));
if abs(val-expected)>tol*max(1,abs(expected))
    error("Test_polyopvar_times_constant_term failed (two constant terms): value %g, expected %g.",val,expected)
end
ncases = ncases+1;

fprintf('Test_polyopvar_times_constant_term passed (%d cases).\n',ncases)


function Y = eval_V(V,vals)
% Evaluate the functional V (polyopvar of intop kernels, or constant tensopvar) at the
% states vals.(name), polynomials in s. The dummy variables of s are numbered over the
% factors of the monomial, and named s_j if the monomial has more than 1 factor.
Y = 0;
Zf = full(V.Z);
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op), continue, end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    if d==0
        % Constant term: tensopvar with a sopvar without variables
        Y = Y + op.ops{1}.params{1};
        continue
    end
    v = 1;
    for j=1:d
        xj = vals.(V.varname{idcs(j)});
        if d>1
            xj.varname(strcmp(xj.varname,'s')) = {['s_',num2str(j)]};
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

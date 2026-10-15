% Test_piesos_eq_canonical.m -- tests that piesos_eq does not
% over-constrain the coefficients of equal functionals. The equation
%   Kop*(x o x) = T*(x o x),
% for T a given functional on the full square [0,3,0,3] with a kernel that
% is not symmetric, and Kop the unknown functional of piesos_intopvar with
% all integrals (multiplier, lower [0,3,0,1], and upper [0,2,0,3]), must
% be feasible, and the solved functional must equal T.
%
% The unknown functional has no integral over the full square, only the
% two chains: the equation is feasible only if combine_terms writes the
% square as the chain [0,3,0,1] with the kernel k(s1,s2)+k(s2,s1). This is
% also tested for T given by the lower chain and the kernel
% k(s1,s2)+k(s2,s1), and by the upper chain, which give the same T.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Only the solved functional is compared with T, on random products
%   p(s_1)*p(s_2): the kernels of Kop are not unique.

pvar s s_dum1 s_dum2
rng(3);
kap = 1+2*s_dum1-s_dum2+0.5*s_dum1*s_dum2;          % not symmetric
ksw = subs(kap,[s_dum1;s_dum2],[s_dum2;s_dum1]);
kmul = 2-s_dum1+s_dum1^2;

% T given in three ways: the square, the lower chain, and the upper chain
Tarr = cell(1,3);
Tarr{1} = intop({kmul;kap},{[0,3,1,1;0,3,0,3]},{'s'},[0,1]);
Tarr{2} = intop({kmul;kap+ksw},{[0,3,1,1;0,3,0,1]},{'s'},[0,1]);
Tarr{3} = intop({kmul;kap+ksw},{[0,3,1,1;0,2,0,3]},{'s'},[0,1]);

for j=1:3
    x1 = polyopvar('x1',polynomial({'s'}),[0,1]);
    prog = piesos_program(x1);
    [prog,Kop] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',2);
    prog = piesos_eq(prog,make_V(x1,Kop-Tarr{j},2));
    [prog,~] = sosposmatr(prog,1);      % an equality-only program cannot be solved by SeDuMi
    out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
    if ~(info.pinf==0 && info.dinf==0)
        error("Test_piesos_eq_canonical failed: form %d of T, the equation Kop = T is infeasible.",j)
    end
    Ks = piesos_getsol(prog,Kop);
    for t=1:3
        c = randn(3,1);
        pvar s_1 s_2
        xval = (c(1)+c(2)*s_1+c(3)*s_1^2)*(c(1)+c(2)*s_2+c(3)*s_2^2);
        val = double(apply_op(Ks,xval));
        ref = double(apply_op(Tarr{j},xval));
        if abs(val-ref)>1e-6*max(1,abs(ref))
            error("Test_piesos_eq_canonical failed: form %d of T, the solved functional differs from T: %g vs %g.",j,val,ref)
        end
    end
end

disp('Test_piesos_eq_canonical passed.')


%% ---- local functions ----------------------------------------------------
function V = make_V(x1,K,d)
% d-th power of the state, with coefficient functional K, as in piesos_sosvar
V = x1;
for i=2:d
    V = V*x1;
end
V.C = tensopmat({K});
V.dims.out = V.C.dims.out;      V.vars.out = V.C.vars.out;
V.dom.out = V.C.dom.out;        V.depmat.out = V.C.depmat.out;
V.chkval = true;
end

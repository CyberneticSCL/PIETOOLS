function Test_piesos_tensor_gram()
% TEST_PIESOS_TENSOR_GRAM tests, using assertions, the SOS functional with a
% tensor-product basis operator (opts.otype=[2,2], constant Pmat),
%   [prog,V,Pmat,Zop] = piesos_sosvar(prog,Z,opts),
%   V = <Zop*Z, Pmat*Zop*Z>,   Zop*Z a function of d copies of the primary
%                              variables for a monomial of d factors,
% which represents products of independent integrals, e.g. (int x^2)^2 =
% <x o x,x o x>, and the recipe to certify them in piesos_ineq,
%   [prog,gt] = piesos_sosvar(prog,dmonomials(x,2),tensor_opts);
%   prog = piesos_ineq(prog,f-gt);
%
% Ground truths are independent of the function under test:
%  - values (case_values_*): V at fixed decision variables, chosen such that
%    Pmat is the identity I or the matrix of ones (a function of the entries
%    of Pmat, which does not depend on the order of the rows of Zop), equals
%    the integrals of the squares of the rows of Zop*Z, computed by
%    polynomial integration: for Z=x o x, the rows are products
%    (A_a x)(s1)*(A_b x)(s2) of the multiplier, lower and upper integral of
%    x, so that V(I) = (sum_a int (A_a x)^2)^2 and V(ones) = (int (sum_a A_a
%    x)^2)^2; likewise for N=2 spatial variables, two states, a vector state,
%    and a mixed basis [x; x o x], for which Zop*x is constant in the second
%    copy; with Pmat a polynomial in the copies s_1, s_2 (degs.out, psatz),
%    V at random decision variable values equals <Zop*z,Pmat*Zop*z>, by
%    apply_op on the monomials z and polynomial integration over the copies
%    (case_values_pmat);
%  - feasibility (case_certification): f=(int x^2)^2 + 0.1 int x^4 is not in
%    the cone of piesos_ineq (pinf=dinf=1), but f-gt is; a negative f, or
%    f=(int x^2)^2 - int x^4 <= 0 (Cauchy-Schwarz), is infeasible; a
%    multiplier p1 with a Gram matrix [1 1;1 1] (moment, (int x)^2) gives a
%    product p1*<x,x> in the cone of piesos_ineq, one with [1 -1;-1 1]
%    ((int (L-U)x)^2, piecewise kernel) only with gt; weighted products such
%    as (int s^2 x^2)(int x^2) only with monomials s_1 in Zs, and with
%    weights that are not squares, such as s and s(1-s), only with psatz
%    multipliers, Pmat being positive semidefinite on the box
%    (case_certification_pmat);
%  - the guards of piesos_quadopvar for otype(1)=2.
% The programs of case_certification are solved by SeDuMi: for exact squares
% (boundary of the PSD cone) it may report numerr=1, with pinf=dinf=0, so
% only pinf, dinf are used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_tensor_gram
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ (CC), 11/10/2026: Initial coding;
% DJ (CC), 11/10/2026: Add Pmat in the copies s_1, s_2 with psatz:
%                        case_values_pmat, case_certification_pmat;

% NOTES
% - A scalar degs.in applies to the first factor of a monomial only: the
%   degrees are given per factor, e.g. degs.in=[d;d] for N=1. A scalar
%   degs.out is the degree of Zs in each copy, and [w1;w2] a degree per copy;
% - The mixed basis [x; x o x] (case_values_mixed) is accepted: block 1
%   maps to one copy, block 2 to two, and the product is taken with block 1
%   constant in the second copy. It is not used by the recipe;
% - int x^4 - (int x^2)^2 >= 0, true by Cauchy-Schwarz, has the form of a
%   variance and is not in the sum of the cones: case_certification pins
%   this, and must be updated if it becomes certified;
% - the test needs the helpers eval_dpvar_at, eval_intop_dpvar and
%   eval_polyopvar (Testfolder/intop, Testfolder/polyopvar).

tests = {@case_values_scalar, @case_values_pointwise, @case_values_N2, ...
         @case_values_two_states, @case_values_vector, @case_values_mixed, ...
         @case_values_pmat, @case_certification, @case_certification_pmat, ...
         @case_guards};
nfail = 0;
for i=1:numel(tests)
    tstart = tic;
    try
        tests{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(tests{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_piesos_tensor_gram: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_tensor_gram: all %d cases passed.\n',numel(tests));
end


%% ---- cases ------------------------------------------------------------

function case_values_scalar()
% Z = x o x, N=1, all types (multiplier, lower, upper), degs.in=0: V(I), V(ones)
pvar s
x = polyopvar('x',s,[0,1]);
xv = 1+2*s-s^2;
Lx = int(xv,s,0,s);     Ux = int(xv,s,s,1);
ints = [int_sq(xv,s),int_sq(Lx,s),int_sq(Ux,s)];
Zx = dmonomials(x,2);
[V,Pmat] = declare(x,Zx,struct('degs',struct('in',[0;0],'out',0),'otype',[2,2]));
check_V(V,Pmat,xv,@(b,c,i,j) double(i==j),sum(ints)^2,'N=1, d=0, P=I');
check_V(V,Pmat,xv,@(b,c,i,j) 1,int_sq(xv+Lx+Ux,s)^2,'N=1, d=0, P=ones');
% degrees per factor, integral rows only: rows {lower, upper} x {eta^k}
pvar eta
xe = subs(xv,s,eta);
for dd = {[1;1],[2;1],[0;2]}
    d = dd{1};
    [V,Pmat] = declare(x,Zx,struct('degs',struct('in',d,'out',0),'otype',[2,2],'exclude',[1;0;0]));
    Ifac = zeros(2,1);
    for f=1:2
        for k=0:d(f)
            r1 = int(eta^k*xe,eta,0,s);     r2 = int(eta^k*xe,eta,s,1);
            Ifac(f) = Ifac(f)+int_sq(r1,s)+int_sq(r2,s);
        end
    end
    check_V(V,Pmat,xv,@(b,c,i,j) double(i==j),prod(Ifac),sprintf('N=1, degs.in=%s, integral rows, P=I',mat2str(d)));
end
% degs.in scalar: the first factor only
[V,Pmat] = declare(x,Zx,struct('degs',struct('in',1,'out',0),'otype',[2,2],'exclude',[1;0;0]));
Ifac = zeros(2,1);
for f=1:2
    for k=0:(f==1)
        r1 = int(eta^k*xe,eta,0,s);     r2 = int(eta^k*xe,eta,s,1);
        Ifac(f) = Ifac(f)+int_sq(r1,s)+int_sq(r2,s);
    end
end
check_V(V,Pmat,xv,@(b,c,i,j) double(i==j),prod(Ifac),'N=1, scalar degs.in=1 (first factor), P=I');
end

function case_values_pointwise()
% The default (pointwise) basis is unchanged: rows (A_a x)(s)*(A_b x)(s), P=I:
% V = sum_ab int ((A_a x)(A_b x))^2 ds, a single integral
pvar s
x = polyopvar('x',s,[0,1]);
xv = 1+2*s-s^2;
A = {xv,int(xv,s,0,s),int(xv,s,s,1)};
orc = 0;
for a=1:3
    for b=1:3
        orc = orc+int_sq(A{a}*A{b},s);
    end
end
[V,Pmat] = declare(x,dmonomials(x,2),struct('degs',struct('in',[0;0],'out',0)));
check_V(V,Pmat,xv,@(b,c,i,j) double(i==j),orc,'pointwise, d=0, P=I');
% and it differs from the tensor-product value (the squares of the integrals)
ints = cellfun(@(a) int_sq(a,s),A);
assert(abs(orc-sum(ints)^2)>1e-3*orc,'pointwise and tensor values must differ');
end

function case_values_N2()
% N=2 spatial variables, integral rows only: 4 rows per factor, P=I
pvar s1 s2
x = polyopvar('x',[s1;s2],[0,1;0,1]);
xv = 1+s1+2*s2-s1*s2+s2^2;
% degs.in is specified per factor and per variable: 2 x 2
[V,Pmat] = declare(x,dmonomials(x,2),struct('degs',struct('in',zeros(2,2),'out',0),'otype',[2,2],'exclude',[1;0;0]));
ops = {{[0,s1],[0,s2]},{[0,s1],[s2,1]},{[s1,1],[0,s2]},{[s1,1],[s2,1]}};
I = 0;
for k=1:4
    Ax = int(xv,s1,ops{k}{1}(1),ops{k}{1}(2));
    Ax = int(Ax,s2,ops{k}{2}(1),ops{k}{2}(2));
    I = I+double(int(int(Ax^2,s1,0,1),s2,0,1));
end
check_V(V,Pmat,xv,@(b,c,i,j) double(i==j),I^2,'N=2, integral rows, P=I');
end

function case_values_two_states()
% Two scalar states: Z = [x1 o x1; x1 o x2; x2 o x2], block-diagonal P=I:
% V = I1^2 + I1*I2 + I2^2, I_a = sum_types int (A x_a)^2
pvar s
xs = polyopvar({'x1','x2'},[1;1],{'s'},[0,1],true(2,1));
x1 = 1+2*s-s^2;     x2 = s+3*s^3;
[V,Pmat] = declare(xs,dmonomials(xs,[2 0;1 1;0 2]),struct('degs',struct('in',[0;0],'out',0),'otype',[2,2]));
Ia = @(a) int_sq(a,s)+int_sq(int(a,s,0,s),s)+int_sq(int(a,s,s,1),s);
I1 = Ia(x1);    I2 = Ia(x2);
check_V(V,Pmat,{x1;x2},@(b,c,i,j) double(b==c && i==j),I1^2+I1*I2+I2^2,'two states, block-diagonal P=I');
end

function case_values_vector()
% Vector state y in R^2: Z = y o y (4 entries), P=I: V = (sum_i sum_a int
% (A_a y_i)^2)^2, the rows of Zop acting on each component
pvar s
y = polyopvar({'y'},2,{'s'},[0,1],true(1,1));
yv = [1+2*s-s^2; s+3*s^3];
Ia = @(a) int_sq(a,s)+int_sq(int(a,s,0,s),s)+int_sq(int(a,s,s,1),s);
[V,Pmat] = declare(y,dmonomials(y,2),struct('degs',struct('in',[0;0],'out',0),'otype',[2,2]));
check_V(V,Pmat,yv,@(b,c,i,j) double(i==j),(Ia(yv(1))+Ia(yv(2)))^2,'vector y, P=I');
end

function case_values_pmat()
% Pmat in the copies s_1, s_2 of the primary variable: Pmat=Zs'*Q*Zs with Zs a
% basis in the copies (degs.out), with psatz multipliers g_k(s_1,s_2). V at
% random decision variable values against the direct evaluation of
% <Zop*z,Pmat*Zop*z>, by apply_op on the monomials z and polynomial integration
% over the copies, independent of the composition in piesos_quadvar and of int
pvar s s1 s2
x = polyopvar('x',s,[0,1]);
X = struct('x',1+2*s-s^2);
Z2 = dmonomials(x,2);
% degs.in, degs.out, exclude, psatz: out is the degree in each copy (scalar, or a
% row per copy), [1;0] has degree 1 in copy 1 only
specs = {[0;0],1,[],0; [0;0],[2;1],[],0; [1;0],[1;1],[],0; [0;0],1,[],1; ...
         [0;0],2,[],1; [1;1],1,[1;0;0],2};
for c=1:size(specs,1)
    opts = struct('degs',struct('in',specs{c,1},'out',specs{c,2}),'otype',[2,2]);
    if ~isempty(specs{c,3}),    opts.exclude = specs{c,3};  end
    if specs{c,4}>0,            opts.psatz = specs{c,4};    end
    label = sprintf('N=1, degs.in=%s, degs.out=%s, psatz=%d',mat2str(specs{c,1}),mat2str(specs{c,2}),specs{c,4});
    [~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x),Z2,opts);
    assert(all(ismember(Pmat{1,1}.varname,{'s_1','s_2'})),'%s: Pmat must depend on the copies s_1, s_2 only',label);
    check_V_random(V,Pmat,Zop,Z2,X,label,c);
end
% mixed basis [x; x o x]: block 1 depends on the copy s_1 only
Zm = dmonomials(x,[1;2]);
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x),Zm,struct('degs',struct('in',[0;0],'out',1),'otype',[2,2]));
assert(isequal(Pmat{1,1}.varname,{'s_1'}),'block 1 of the mixed basis depends on s_1 only');
check_V_random(V,Pmat,Zop,Zm,X,'mixed [x; x o x], degs.out=1',10);
% two states
xs = polyopvar({'x1','x2'},[1;1],{'s'},[0,1],true(2,1));
Xs = struct('x1',1+2*s-s^2,'x2',s+3*s^3);
Zs2 = dmonomials(xs,[2 0;1 1;0 2]);
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(xs),Zs2,struct('degs',struct('in',[0;0],'out',1),'otype',[2,2]));
check_V_random(V,Pmat,Zop,Zs2,Xs,'two states, degs.out=1',11);
% N=2: copies s1_1, s1_2, s2_1, s2_2; degree 1 in the copies of s1, 0 in those of s2
x2d = polyopvar('x',[s1;s2],[0,1;0,1]);
X2 = struct('x',1+s1+2*s2-s1*s2+s2^2);
Z22 = dmonomials(x2d,2);
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x2d),Z22,struct('degs',struct('in',zeros(2,2),'out',[1 0]),'otype',[2,2],'exclude',[1;0;0]));
assert(all(ismember(Pmat{1,1}.varname,{'s1_1','s1_2','s2_1','s2_2'})),'N=2: Pmat must depend on the copies si_j');
check_V_random(V,Pmat,Zop,Z22,X2,'N=2, integral rows, degs.out=[1 0]',12);
end

function case_values_mixed()
% Mixed basis [x; x o x], P=ones: with S=(I+L+U)x, block 1 is S(s1), constant in
% the second copy, block 2 is S(s1)*S(s2): V = int int (S(s1)+S(s1)S(s2))^2
pvar s
x = polyopvar('x',s,[0,1]);
xv = 1+2*s-s^2;
S = xv+int(xv,s,0,s)+int(xv,s,s,1);
[V,Pmat] = declare(x,dmonomials(x,[1;2]),struct('degs',struct('in',[0;0],'out',0),'otype',[2,2]));
check_V(V,Pmat,xv,@(b,c,i,j) 1,int_sq(S,s)*int_sq(1+S,s),'mixed [x; x o x], P=ones');
end

function case_certification()
pvar s
x = polyopvar('x',s,[0,1]);
q = innerprod(x,x);                         % <x,x>
x2 = x.*x;      x4 = innerprod(x2,x2);      % int x^4
Z2 = dmonomials(x,2);
topts = struct('degs',struct('in',[0;0],'out',0),'otype',[2,2]);
% f = (int x^2)^2 + 0.1 int x^4: not in the cone of piesos_ineq, but f-gt is
f = q*q+0.1*x4;
assert(~is_feasible(piesos_ineq(piesos_program(x),f)),'(int x^2)^2 + 0.1 int x^4 must be infeasible in piesos_ineq');
[prog,gt] = piesos_sosvar(piesos_program(x),Z2,topts);
assert(is_feasible(piesos_ineq(prog,f-gt)),'(int x^2)^2 + 0.1 int x^4 - gt must be feasible');
% negative functionals
[prog,gt] = piesos_sosvar(piesos_program(x),Z2,topts);
assert(~is_feasible(piesos_ineq(prog,-1*(q*q)+0.1*x4-gt)),'-(int x^2)^2 + 0.1 int x^4 must be infeasible');
[prog,gt] = piesos_sosvar(piesos_program(x),Z2,topts);
assert(~is_feasible(piesos_ineq(prog,q*q-x4-gt)),'(int x^2)^2 - int x^4 <= 0 must be infeasible');
% int x^4 - (int x^2)^2 >= 0 (Cauchy-Schwarz) is not certified (conservative)
[prog,gt] = piesos_sosvar(piesos_program(x),Z2,topts);
assert(~is_feasible(piesos_ineq(prog,x4-q*q-gt)),'int x^4 - (int x^2)^2 is not certified (see NOTES)');
% vector state y in R^2: (int |y|^2)^2 + 0.1 int |y|^4, with y o y as basis
y = polyopvar({'y'},2,{'s'},[0,1],true(1,1));
yy = otimes(y,y,2);     y2 = y.*y;
fy = innerprod(yy,yy)+0.1*innerprod(y2,ones(2),y2);
assert(~is_feasible(piesos_ineq(piesos_program(y),fy)),'vector: (int |y|^2)^2 + 0.1 int |y|^4 must be infeasible in piesos_ineq');
[prog,gt] = piesos_sosvar(piesos_program(y),dmonomials(y,2),topts);
assert(is_feasible(piesos_ineq(prog,fy-gt)),'vector: (int |y|^2)^2 + 0.1 int |y|^4 - gt must be feasible');
% multiplier p1 = <Zop*x,P*Zop*x> with lower and upper integral rows, P fixed
for P = {[1 1;1 1],[1 -1;-1 1]}
    Pv = P{1};
    is_moment = isequal(Pv,[1 1;1 1]);
    for use_gt = [false,true]
        prog = piesos_program(x);
        [prog,p1,Pm] = piesos_sosvar(prog,dmonomials(x,1),struct('degs',0,'exclude',[1;0;0]));
        for i=1:2
            for j=1:2
                prog = soseq(prog,Pm{1,1}(i,j)-Pv(i,j));
            end
        end
        F = p1*q+0.1*x4;
        if use_gt
            [prog,gt] = piesos_sosvar(prog,Z2,topts);
            F = F-gt;
        end
        ok = is_feasible(piesos_ineq(prog,F));
        % p1*<x,x> in the pointwise cone only for the moment (int x)^2, in the sum of the cones always
        assert(ok==(is_moment || use_gt),sprintf('p1*<x,x> + 0.1 int x^4 with P=%s, gt=%d: feasible=%d',mat2str(Pv),use_gt,ok));
    end
end
end

function case_certification_pmat()
% Products with weights in Pmat(s_1,s_2). f1 = (int s^2 x^2)(int x^2) + 0.1 int
% x^4 needs the monomial s_1 in Zs: infeasible for degs.out=0, feasible for
% degs.out=[1;0] (copy 1 only). f2 = (int s^2 x^2)(int (1-s)^2 x^2) needs s_1 and
% s_2: infeasible for [1;0], feasible for [1;1]. The weights s and s(1-s), that
% are nonnegative on [0,1] but not squares, are certified only with psatz
% multipliers; s-2, negative on [0,1], never. With psatz, Pmat is positive
% semidefinite at points of [0,1]^2
pvar s
x = polyopvar('x',s,[0,1]);
Z2 = dmonomials(x,2);
x2 = x.*x;      x4 = innerprod(x2,x2);      q = innerprod(x,x);
qs2 = innerprod(x,s^2,x);       qc2 = innerprod(x,(1-s)^2,x);
f1 = qs2*q+0.1*x4;      f2 = qs2*qc2+0.1*x4;
fs = innerprod(x,s,x)*q+0.1*x4;
fb = innerprod(x,s*(1-s),x)*q+0.1*x4;
fn = innerprod(x,s-2,x)*q+0.1*x4;
% f, degs.out, psatz, expected feasibility
specs = {f1,0,0,false; f1,[1;0],0,true; f2,[1;0],0,false; f2,[1;1],0,true; ...
         fs,1,0,false; fs,0,1,true; fs,1,1,true; fb,2,0,false; fb,1,2,true; ...
         fn,1,1,false};
for c=1:size(specs,1)
    opts = struct('degs',struct('in',[0;0],'out',specs{c,2}),'otype',[2,2]);
    if specs{c,3}>0,    opts.psatz = specs{c,3};    end
    [prog,gt,Pm] = piesos_sosvar(piesos_program(x),Z2,opts);
    [ok,prog_sol] = is_feasible(piesos_ineq(prog,specs{c,1}-gt));
    msg = sprintf('weights case %d (degs.out=%s, psatz=%d): feasible=%d, expected %d',c,mat2str(specs{c,2}),specs{c,3},ok,specs{c,4});
    assert(ok==specs{c,4},msg);
    if ok && specs{c,3}>0
        % Pmat is PSD on the box where the multipliers are nonnegative
        Pn = piesos_getsol(prog_sol,Pm{1,1});
        rng(c);
        for t=1:20
            A = eval_poly_at(Pn,containers.Map({'s_1','s_2'},{rand,rand}));
            assert(min(eig((A+A')/2))>-1e-5*max(1,norm(A)),sprintf('%s: Pmat is not PSD in [0,1]^2',msg));
        end
    end
end
% the products are not in the pointwise cone
assert(~is_feasible(piesos_ineq(piesos_program(x),f1)),'(int s^2 x^2)(int x^2) + 0.1 int x^4 must be infeasible in piesos_ineq');
end

function case_guards()
pvar s
x = polyopvar('x',s,[0,1]);
Z2 = dmonomials(x,2);
prog = piesos_program(x);
d0 = struct('in',[0;0],'out',0);
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',d0,'otype',[2,2],'psatz',true(2,1))),'psatz must be specified as numeric');
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',struct('in',[0;0],'out',[1 2]),'otype',[2,2])),'Degrees of the output variables');
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',struct('in',[0;0],'out',-1),'otype',[2,2])),'nonnegative integers');
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',d0,'L',struct('otype',[2,2]),'R',struct())),'both ZopL and ZopR');
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',d0,'otype',[1,1])),'Pointwise product of the input');
% otype(1)=1, and otype(1)=2 with psatz (numeric), with degs.out a scalar or a row
% per copy, are accepted
piesos_sosvar(prog,Z2,struct('degs',d0,'otype',[1,2]));
piesos_sosvar(prog,Z2,struct('degs',d0,'otype',[2,2],'psatz',0));
piesos_sosvar(prog,Z2,struct('degs',d0,'otype',[2,2],'psatz',2));
piesos_sosvar(prog,Z2,struct('degs',struct('in',[0;0],'out',[1;2]),'otype',[2,2]));
piesos_sosvar(prog,Z2,struct('degs',struct('in',[0;0],'out',[1;2;3]),'otype',[2,2]));
piesos_sosvar(prog,Z2,struct('degs',1,'otype',[2,2]));
end


%% ---- helpers ----------------------------------------------------------

function I = int_sq(f,s)
% int_0^1 f(s)^2 ds
I = double(int(f^2,s,0,1));
end

function [V,Pmat] = declare(x,Z,opts)
% V = <Zop*Z,Pmat*Zop*Z> as declared by piesos_sosvar
[~,V,Pmat] = piesos_sosvar(piesos_program(x),Z,opts);
end

function check_V(V,Pmat,xvals,target,oracle,label)
% Value of V at the state xvals (polynomial, or a cell of them), for decision
% variables such that Pmat(b,c)(i,j) = target(b,c,i,j), against the oracle
if ~iscell(xvals),  xvals = {xvals};    end
d0map = gram_values(Pmat,target);
Vn = V;
ops = V.C.ops;
for k=1:numel(ops)
    if isa(ops{1,k},'intop')
        ops{1,k} = eval_intop_dpvar(ops{1,k},d0map);
    end
end
Vn.C.ops = ops;     Vn.C.chkval = true;
val = double(eval_polyopvar(Vn,xvals));
relerr = abs(val-oracle)/abs(oracle);
assert(relerr<1e-8,sprintf('%s: V=%.10g, oracle=%.10g (relative error %.2e)',label,val,oracle,relerr));
end

function d0map = gram_values(Pmat,target)
% Values d0map of the decision variables in the constant blocks Pmat{b,c} (dpvar
% layout: matrix row outer, decision variable inner, constant term first) such
% that Pmat{b,c}(i,j) = target(b,c,i,j), by least squares
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)];                                    %#ok<AGROW>
end
names = unique(names,'stable');
A = zeros(0,numel(names));      rhs = zeros(0,1);
[nb,nc] = size(Pmat);
for b=1:nb
    for c=1:nc
        P = Pmat{b,c};
        [m,n] = size(P);
        nd = numel(P.dvarname);
        [~,loc] = ismember(P.dvarname,names);
        for i=1:m
            for j=1:n
                row = zeros(1,numel(names));
                row(loc) = full(P.C((i-1)*(nd+1)+2:i*(nd+1),j)).';
                const = full(P.C((i-1)*(nd+1)+1,j));
                A(end+1,:) = row;                                           %#ok<AGROW>
                rhs(end+1,1) = target(b,c,i,j)-const;                       %#ok<AGROW>
            end
        end
    end
end
v = A\rhs;
assert(norm(A*v-rhs)<1e-9,'the target Gram matrix is not in the span of the decision variables');
d0map = containers.Map(names,num2cell(v));
end

function check_V_random(V,Pmat,Zop,Z,X,label,seed)
% V at the states X (struct: state name -> polynomial) and random values of the
% decision variables of Pmat, against the direct evaluation of <Zop*z,Pmat*Zop*z>:
% z the monomials of Z, Zop applied by apply_op, Pmat evaluated at the decision
% variables, and the integral over all variables of the product (the copies
% si_j, on [0,1]) by polynomial integration
rng(seed);
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)];                                    %#ok<AGROW>
end
names = unique(names,'stable');
d0map = containers.Map(names,num2cell(0.5+rand(numel(names),1)));
% V through its kernels
Vn = V;
ops = V.C.ops;
for k=1:numel(ops)
    if isa(ops{1,k},'intop')
        ops{1,k} = eval_intop_dpvar(ops{1,k},d0map);
    end
end
Vn.C.ops = ops;     Vn.C.chkval = true;
xv = cell(numel(Vn.varname),1);
for k=1:numel(xv),  xv{k} = X.(Vn.varname{k});     end
val = double(eval_polyopvar(Vn,xv));
% direct evaluation
Y = apply_side(Zop,Z,X);
[q1,q2] = size(Pmat);
rows = cell(q1,1);
for i=1:q1
    blks = cell(1,q2);
    for j=1:q2
        blks{j} = eval_dpvar_at(Pmat{i,j},d0map);
    end
    rows{i} = horzcat(blks{:});
end
F = Y.'*vertcat(rows{:})*Y;
vn = F.varname;
for k=1:numel(vn)
    sk = polynomial(vn(k));
    F = int(F,sk(1),0,1);
end
orc = double(F);
assert(abs(orc)>1e-6,'%s: the value %g is too small to test.',label,orc);
assert(abs(val-orc)<=1e-8*abs(orc),'%s: V=%.10g, direct evaluation=%.10g (relative error %.2e)',label,val,orc,abs(val-orc)/abs(orc));
end

function Y = apply_side(Zop,Z,X)
% Zop applied to the entries of a monomial basis Z (identity coefficients): block
% i acts on the monomial of row i of Z.C, the product of the states in X listed in
% the row of Z.Z that has a block in row i. As in Test_piesos_quadvar
Zf = full(Z.Z);
nr = size(Z.C.ops,1);
Xc = cell(1,nr);
for i=1:nr
    k = find(~cellfun(@isempty,Z.C.ops(i,:)));
    assert(isscalar(k),'Row %d of the basis does not have a single monomial.',i);
    idc = Zf(k,Zf(k,:)~=0);
    Xc{i} = cellfun(@(nm) X.(nm),reshape(Z.varname(idc),1,[]),'UniformOutput',false);
end
Y = apply_op(fix_empty_params(Zop),Xc);
end

function T = fix_empty_params(T)
% Zeros for the empty parameters params{k}=[] of the excluded terms in the sopvar
% factors of the tensor-PI operators of T (tensopmat), as apply_op cannot handle them
if isa(T,'tensopmat')
    for k=1:numel(T.ops)
        if ~isempty(T.ops{k}),  T.ops{k} = fix_empty_params(T.ops{k});  end
    end
elseif isa(T,'tensopvar')
    for k=1:numel(T.ops)
        T.ops{k} = fix_empty_params(T.ops{k});
    end
elseif isa(T,'sopvar')
    params = T.params;
    nrow = T.dims(1)*prod(cellfun(@numel,T.ZL));    ncol = T.dims(2)*prod(cellfun(@numel,T.ZR));
    for k=1:numel(params)
        if isempty(params{k}),  params{k} = sparse(nrow,ncol);  end
    end
    T = sopvar(params,T.vars,T.ZL,T.ZR,T.dom,T.dims);
end
end

function [ok,prog_sol] = is_feasible(prog)
% Solve with SeDuMi; feasible if neither pinf nor dinf is reported. A program
% without PSD variable fails in SeDuMi: sosposmatr adds a trivial one
[prog,~] = sosposmatr(prog,1);
evalc('prog_sol = piesos_solve(prog,struct(''solver'',''sedumi'',''simplify'',true));');
info = prog_sol.solinfo.info;
ok = ~info.pinf && ~info.dinf;
end

function assert_err(fn,pattern)
try
    fn();
catch err
    assert(~isempty(strfind(err.message,pattern)),sprintf('Unexpected error message: "%s" (expected to contain "%s")',err.message,pattern)); %#ok<STREMP>
    return
end
error('No error thrown (expected "%s").',pattern);
end

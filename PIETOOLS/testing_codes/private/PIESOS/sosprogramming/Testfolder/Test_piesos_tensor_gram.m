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
%  - the guards of piesos_quadopvar for otype(1)=2;
%  - otype per block in monopmat (case_monopmat_otype): the block rows of a
%    q x 1 otype equal those of separate calls with each type, and otype=3
%    gives the rows of type 1 and 2 of a block acting on its input, one
%    row for a single factor;
%  - mixed otype (case_values_mixed_otype): a block row of each type for
%    otype=3, with block counts (no second row for a single factor:
%    [x; x o x] gives 3 x 3 blocks, x alone 1 x 1), and ZopL, ZopR of
%    different type (general Pmat): V at random decision variable values
%    against the direct evaluation over all copies, also on [0,2] and with
%    psatz, where the otype-2 basis [x; x o x] was off by the domain length
%    in the (x,x) term and failed with psatz before 11/19/2026
%    (int_other_copies in piesos_quadvar); and certification
%    (case_certification_mixed): (int x^2)^2 + 0.1 int x^4 is certified by
%    piesos_ineq with otype=3 alone, negative functionals and
%    int x^4 - (int x^2)^2 are not.
%  - regions of integration (opts.intmat, case_*_regions): V at random
%    decision variables against the sum over the regions of the integrand
%    integrated by polynomial int/subs over each region (explicit rows,
%    with a delta, N = 1 and 2); the presets with the same P in each
%    region against the whole domain plus the orderings with ties, built
%    from the ranks of the copies (not from piesos_quadvar); errors; the
%    option otype used for ZopL and ZopR; and the embedding inequality
%    (int z^2)(int x^2) - int z^4 >= 0, z(s) = int_0^s x, certified with
%    otype 3 and the triangles ('all'), not over the square ('full'); psatz
%    on the faces of a region (case_psatz_regions): int int_{s_2<s_1}
%    (s_1-s_2) x(s_1)^2 x(s_2)^2 certified with psatz, not without.
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
% DJ (CC), 11/19/2026: Mixed otype: case_monopmat_otype,
%                        case_values_mixed_otype,
%                        case_certification_mixed; check_V_random with
%                        optional ZopR and domain, integrating over all
%                        copies; case_guards: ZopL, ZopR of different type
%                        are accepted, an SOS Pmat then fails on degrees;
% DJ (CC), 11/19/2026: Regions (opts.intmat) and option otype:
%                        case_values_regions, case_values_presets_regions,
%                        case_errors_regions, case_otype_regions,
%                        case_certification_regions;
% DJ (CC), 11/20/2026: Default otype 3 of piesos_sosvar: otype 1 explicit
%                        where the pointwise cone is meant, and the products
%                        certified by default (case_certification);
% DJ (CC), 11/20/2026: case_values_factor_orders: both orders of a monomial
%                        (dmonomials 'full'), pointwise rows of one order
%                        (opts.unique of piesos_quadvar);
% DJ (CC), 11/20/2026: Default otype of piesos_ineq 'auto' (was 3): the
%                        checks of otype 3 in case_certification pass it
%                        (for p1*<x,x> + 0.1 int x^4 'auto' gives otype 1);

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
         @case_guards, @case_monopmat_otype, @case_values_mixed_otype, ...      % DJ (CC), 11/19/2026
         @case_certification_mixed, @case_values_regions, ...
         @case_values_presets_regions, @case_errors_regions, ...
         @case_otype_regions, @case_certification_regions, @case_psatz_regions, ... % DJ (CC), 11/19/2026
         @case_values_factor_orders};                                       % DJ (CC), 11/20/2026
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
% The pointwise basis (otype 1, the default before 11/20/2026) is unchanged:
% rows (A_a x)(s)*(A_b x)(s), P=I: V = sum_ab int ((A_a x)(A_b x))^2 ds, a
% single integral
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
[V,Pmat] = declare(x,dmonomials(x,2),struct('degs',struct('in',[0;0],'out',0),'otype',1));    % DJ (CC), 11/20/2026
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
% f = (int x^2)^2 + 0.1 int x^4: not in the cone of piesos_ineq with otype 1,
% but f-gt is, and f is in the cone of otype 3 (also chosen by the default
% otype 'auto' of piesos_ineq for this f)
f = q*q+0.1*x4;
o1 = struct('otype',1);                                                     % DJ (CC), 11/20/2026
o3 = struct('otype',3);                                                     % DJ (CC), 11/20/2026
assert(~is_feasible(piesos_ineq(piesos_program(x),f,o1)),'(int x^2)^2 + 0.1 int x^4 must be infeasible in piesos_ineq, otype 1');   % DJ (CC), 11/20/2026
assert(is_feasible(piesos_ineq(piesos_program(x),f,o3)),'(int x^2)^2 + 0.1 int x^4 must be feasible in piesos_ineq, otype 3');     % DJ (CC), 11/20/2026
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
assert(~is_feasible(piesos_ineq(piesos_program(y),fy,o1)),'vector: (int |y|^2)^2 + 0.1 int |y|^4 must be infeasible in piesos_ineq, otype 1');   % DJ (CC), 11/20/2026
assert(is_feasible(piesos_ineq(piesos_program(y),fy,o3)),'vector: (int |y|^2)^2 + 0.1 int |y|^4 must be feasible in piesos_ineq, otype 3');     % DJ (CC), 11/20/2026
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
        ok = is_feasible(piesos_ineq(prog,F,o1));                          % DJ (CC), 11/20/2026
        % p1*<x,x> in the pointwise cone only for the moment (int x)^2, in the sum of the cones always
        assert(ok==(is_moment || use_gt),sprintf('p1*<x,x> + 0.1 int x^4 with P=%s, gt=%d: feasible=%d',mat2str(Pv),use_gt,ok));
        % and always in the cone of otype 3 (tensor rows (A x)(s_1)*x(s_2))  % DJ (CC), 11/20/2026
        ok = is_feasible(piesos_ineq(prog,F,o3));                           % DJ (CC), 11/20/2026
        assert(ok,sprintf('p1*<x,x> + 0.1 int x^4 with P=%s, gt=%d, otype 3: infeasible',mat2str(Pv),use_gt));    % DJ (CC), 11/20/2026
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
    % otype 1 in piesos_ineq, so that only gt can certify the products
    [ok,prog_sol] = is_feasible(piesos_ineq(prog,specs{c,1}-gt,struct('otype',1)));    % DJ (CC), 11/20/2026
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
assert(~is_feasible(piesos_ineq(piesos_program(x),f1,struct('otype',1))),'(int s^2 x^2)(int x^2) + 0.1 int x^4 must be infeasible in piesos_ineq');    % DJ (CC), 11/20/2026
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
%assert_err(@() piesos_sosvar(prog,Z2,struct('degs',d0,'L',struct('otype',[2,2]),'R',struct())),'both ZopL and ZopR');   % DJ (CC), 11/19/2026 (was)
% ZopL, ZopR of different type are accepted (mixed otype), but an SOS Pmat
% needs the same basis on both sides: rows in two copies on the left, one
% on the right
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',d0,'L',struct('otype',[2,2]),'R',struct('otype',1))),'Degrees of the primary variables');    % DJ (CC), 11/20/2026
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',d0,'otype',[2,2;2,2])),'q x 1 or q x 2 array');    % DJ (CC), 11/19/2026
assert_err(@() piesos_sosvar(prog,Z2,struct('degs',d0,'otype',4)),'or 3 for both');    % DJ (CC), 11/19/2026
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

function case_monopmat_otype()                                              % DJ (CC), 11/19/2026
% otype per block in monopmat, for the blocks x (one factor) and x o x (two):
% the block rows of a q x 1 otype are those of separate calls with a single
% type, and otype=3 gives the rows of type 1 and 2 of a block, both acting on
% its input (block column), one row for the single factor
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,1],'out',[0,1]);
depmat = struct('in',[1;2],'out',[true;true]);
degs = struct('in',[1;1],'out',0);
mono = @(otype) monopmat([1;1],vars,dom,depmat,degs,struct('otype',otype));
Zop1 = mono(1);     Zop2 = mono(2);
% a type per block
Zop = mono([1;2]);
assert(isequal(size(Zop.ops),[2,2]),'otype [1;2]: %s blocks, expected [2 2]',mat2str(size(Zop.ops)));
assert(isequal(Zop.ops{1,1},Zop1.ops{1,1}) && isequal(Zop.ops{2,2},Zop2.ops{2,2}),'otype [1;2]: blocks differ from those of a single type');
% both types: rows (x, type 1), (x o x, type 1), (x o x, type 2)
for otype = {3,[1;3],[3;3]}
    Zop = mono(otype{1});
    label = sprintf('otype %s',mat2str(otype{1}));
    assert(isequal(~cellfun(@isempty,Zop.ops),logical([1 0;0 1;0 1])),'%s: block pattern %s',label,mat2str(~cellfun(@isempty,Zop.ops)));
    assert(isequal(full(Zop.depmat.out(:)),[1;1;2]),'%s: copies per block row %s, expected [1;1;2]',label,mat2str(full(Zop.depmat.out(:))));
    assert(isequal(Zop.ops{2,2},Zop1.ops{2,2}) && isequal(Zop.ops{3,2},Zop2.ops{2,2}),'%s: blocks differ from those of a single type',label);
end
% otype 3 for the single factor only: no second row
Zop = mono([3;2]);
assert(isequal(size(Zop.ops),[2,2]),'otype [3;2]: %s blocks, expected [2 2]',mat2str(size(Zop.ops)));
assert(isequal(Zop.ops{2,2},Zop2.ops{2,2}),'otype [3;2]: block 2 differs from that of type 2');
end

function case_values_mixed_otype()                                          % DJ (CC), 11/19/2026
% Mixed otype: a pointwise and a tensor-product block row for otype=3, the
% pointwise rows in copy s_1 only. Block counts: a second (otype-2) row only
% for the monomials with several copies. V at random decision variable
% values against the direct evaluation of <ZopL*z,Pmat*ZopR*z>, integrated
% over all copies
pvar s
mixed = 3;
x = polyopvar('x',s,[0,1]);
X = struct('x',1+2*s-s^2);
% Z, degs.in, degs.out, exclude, psatz, expected blocks of Pmat
specs = {dmonomials(x,2),[0;0],1,[],0,[2,2]; dmonomials(x,2),[0;0],1,[],1,[2,2]; ...
         dmonomials(x,[1;2]),[0;0],1,[],0,[3,3]; dmonomials(x,1),0,1,[],0,[1,1]; ...
         dmonomials(x,[1;2]),[1;1],1,[1;0;0],0,[3,3]};
for c=1:size(specs,1)
    opts = struct('degs',struct('in',specs{c,2},'out',specs{c,3}),'otype',mixed);
    if ~isempty(specs{c,4}),    opts.exclude = specs{c,4};  end
    if specs{c,5}>0,            opts.psatz = specs{c,5};    end
    label = sprintf('mixed otype, case %d',c);
    [~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x),specs{c,1},opts);
    assert(isequal(size(Pmat),specs{c,6}),'%s: Pmat has %s blocks, expected %s',label,mat2str(size(Pmat)),mat2str(specs{c,6}));
    check_V_random(V,Pmat,Zop,specs{c,1},X,label,20+c);
end
% two states, degrees 1 and 2: 2 + 3 pointwise rows, 3 tensor-product rows
xs = polyopvar({'x1','x2'},[1;1],{'s'},[0,1],true(2,1));
Xs = struct('x1',1+2*s-s^2,'x2',s+3*s^3);
Zs12 = dmonomials(xs,[1 0;0 1;2 0;1 1;0 2]);
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(xs),Zs12,struct('degs',struct('in',[0;0],'out',0),'otype',mixed));
assert(isequal(size(Pmat),[8,8]),'two states: Pmat has %s blocks, expected [8 8]',mat2str(size(Pmat)));
check_V_random(V,Pmat,Zop,Zs12,Xs,'mixed otype, two states',30);
% the domain [0,2] (V over both copies), mixed and otype [2,2], with psatz:
% the (x,x) term of [x; x o x] depends on copy 1 only
x2 = polyopvar('x',s,[0,2]);
Zm2 = dmonomials(x2,[1;2]);
for ot = {mixed,[2,2]}
    for ps = [0,1]
        opts = struct('degs',struct('in',[0;0],'out',1),'otype',ot{1});
        if ps>0,    opts.psatz = ps;    end
        label = sprintf('[x; x o x] on [0,2], otype %s, psatz %d',mat2str(ot{1}),ps);
        [~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x2),Zm2,opts);
        check_V_random(V,Pmat,Zop,Zm2,X,label,40+ps,Zop,[0,2]);
    end
end
% ZopL and ZopR of different type (general Pmat, piesos_quadvar)
for LR = {{[1,2],[2,2],dmonomials(x,2)},{mixed,[2,2],dmonomials(x,[1;2])}}
    opts = struct('degs',struct('in',[0;0],'out',1),'L',struct('otype',LR{1}{1}),'R',struct('otype',LR{1}{2}));
    [~,V,Pmat,ZopL,ZopR] = piesos_quadvar(piesos_program(x),LR{1}{3},[],opts);
    check_V_random(V,Pmat,ZopL,LR{1}{3},X,sprintf('ZopL otype %s, ZopR otype %s',mat2str(LR{1}{1}),mat2str(LR{1}{2})),50,ZopR);
end
end

function case_certification_mixed()                                         % DJ (CC), 11/19/2026
% One Gram matrix of mixed otype in piesos_ineq contains both cones:
% (int x^2)^2 + 0.1 int x^4 is certified without a separate gt (not with
% otype 1 alone); negative functionals, and int x^4 - (int x^2)^2 >= 0
% (Cauchy-Schwarz), are not
pvar s
x = polyopvar('x',s,[0,1]);
q = innerprod(x,x);     x2 = x.*x;      x4 = innerprod(x2,x2);
d0 = struct('in',[0;0],'out',0);
% f, otype, expected feasibility
specs = {q*q+0.1*x4,1,false; q*q+0.1*x4,3,true; ...
         -1*(q*q)+0.1*x4,3,false; x4-q*q,3,false; q*q-x4,3,false};
for c=1:size(specs,1)
    ok = is_feasible(piesos_ineq(piesos_program(x),specs{c,1},struct('otype',specs{c,2},'degs',d0)));
    assert(ok==specs{c,3},'mixed certification case %d (otype %s): feasible=%d, expected %d',c,mat2str(specs{c,2}),ok,specs{c,3});
end
end

% DJ (CC), 11/19/2026: start. Cases for regions of integration (opts.intmat).
function case_values_regions()
% Explicit regions, a Pmat each, at random decision variable values: the
% lower triangle, the diagonal (delta) and the upper triangle for N = 1;
% for N = 2, two rows per variable (4 regions, k the linear index of the
% rows (k1,k2)), with a delta in r. otype 3, so that the pointwise rows (in
% copy 1) meet each region differently
pvar s r
x = polyopvar('x',s,[0,1]);
X = struct('x',1+2*s-s^2);
im = [0,3,0,1; 0,3,1,1; 0,3,1,3];
opts = struct('degs',struct('in',[0;0],'out',1),'otype',3,'intmat',{{im}});
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x),dmonomials(x,2),opts);
assert(isequal(size(regions_of(Pmat)),[3,1]),'N=1: %d regions, expected 3',numel(regions_of(Pmat)));
check_V_regions(V,Pmat,Zop,dmonomials(x,2),X,{im},'N=1, lower, diagonal, upper',61);
x2 = polyopvar('x',[s;r],[0,1;0,2]);
X2 = struct('x',1+s+2*s*r-r^2);
ims = {[0,3,0,1; 0,3,1,3], [0,3,1,1; 0,3,0,1]};
opts = struct('degs',struct('in',zeros(2,2),'out',[1 0]),'otype',3,'exclude',[0;1;1],'intmat',{ims});
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x2),dmonomials(x2,2),opts);
assert(isequal(size(regions_of(Pmat)),[2,2]),'N=2: %s regions, expected [2 2]',mat2str(size(regions_of(Pmat))));
check_V_regions(V,Pmat,Zop,dmonomials(x2,2),X2,ims,'N=2, 2 x 2 regions',62);
end

function case_values_presets_regions()
% The presets, with the same constant P in each region: V equals the
% integral over the whole domain for 'full' and 'steps' (the strict
% orderings partition it), plus the orderings with ties but not all copies
% tied for 'all' (int_cells, from the ranks of the copies). Region counts
% d!, F(d)-1 (F(3) = 13). Degree 2 with integral rows, degree 3 with
% multiplier rows, and N = 2 with a preset per variable
pvar s r
x = polyopvar('x',s,[0,1]);
X = struct('x',1+2*s-s^2);
specs = {2,[0;0],[0;0;0],{'full','steps','all'},[1,2,2]; ...
         3,[0;0;0],[0;1;1],{'full','steps','all'},[1,6,12]};
for c=1:size(specs,1)
    Z = dmonomials(x,specs{c,1});
    for p=1:3
        opts = struct('degs',struct('in',specs{c,2},'out',0),'otype',3,'exclude',specs{c,3},'intmat',specs{c,4}{p});
        [~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x),Z,opts);
        label = sprintf('degree %d, %s',specs{c,1},specs{c,4}{p});
        assert(numel(regions_of(Pmat))==specs{c,5}(p),'%s: %d regions, expected %d',label,numel(regions_of(Pmat)),specs{c,5}(p));
        assert(iscell(Pmat{1})==(p>1),'%s: a single Pmat only for the whole domain',label);
        check_V_presets(V,Pmat,Zop,Z,X,specs{c,4}(p),label,70+10*c+p);
    end
end
x2 = polyopvar('x',[s;r],[0,1;0,2]);
X2 = struct('x',1+s+2*s*r-r^2);
pre = {'all','steps'};
opts = struct('degs',struct('in',zeros(2,2),'out',0),'otype',3,'exclude',[0;1;1],'intmat',{pre});
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x2),dmonomials(x2,2),opts);
assert(isequal(size(regions_of(Pmat)),[2,2]),'N=2: %s regions, expected [2 2]',mat2str(size(regions_of(Pmat))));
check_V_presets(V,Pmat,Zop,dmonomials(x2,2),X2,pre,'N=2, all and steps',90);
end

function case_errors_regions()
pvar s
x = polyopvar('x',s,[0,1]);
Z = dmonomials(x,2);        prog = piesos_program(x);
o = @(im) struct('degs',struct('in',[0;0],'out',0),'otype',3,'intmat',{im});
assert_err(@() piesos_sosvar(prog,Z,o('tri')),'Presets');
assert_err(@() piesos_sosvar(prog,Z,o({'all','all'})),'1 x N cell');
assert_err(@() piesos_sosvar(prog,Z,o([2,3,0,3])),'copies 1,...,j-1');    % s_1 bounded by s_2
assert_err(@() piesos_sosvar(prog,Z,o([0,3])),'q x 4');
assert_err(@() piesos_sosvar(prog,Z,o(zeros(0,4))),'q x 4');
assert_err(@() piesos_sosvar(prog,Z,o([0,3,0,4])),'q x 4');
end

function case_otype_regions()
% opts.otype is used for ZopL and ZopR, unless opts.L or opts.R specify it
pvar s
x = polyopvar('x',s,[0,1]);
Z = dmonomials(x,2);
d0 = struct('in',[0;0],'out',0);
[~,~,~,ZopL,ZopR] = piesos_quadvar(piesos_program(x),Z,[],struct('degs',d0,'otype',3));
assert(isequal(full(ZopL.depmat.out(:)),[1;2]) && isequal(full(ZopR.depmat.out(:)),[1;2]),'otype 3 for both');
[~,~,~,ZopL,ZopR] = piesos_quadvar(piesos_program(x),Z,[],struct('degs',d0,'otype',3,'L',struct('otype',2),'R',struct()));
assert(isequal(full(ZopL.depmat.out(:)),2) && isequal(full(ZopR.depmat.out(:)),[1;2]),'otype 2 for L, 3 for R');
end

function case_certification_regions()
% (int z^2)(int x^2) - int z^4 >= 0 for x = z_s, z(0) = 0 (the embedding step
% of README_TODOs): a lower and an upper triangle Gram of otype 3 certify it,
% the square does not; also through piesos_ineq
pvar s t
zv = pde_var(s,[0,1]);
PDE = [diff(zv,t)==diff(zv,s)+zv*zv; subs(zv,s,0)==0];     % nonlinear: PIESOS form
PIE = convert(PDE,'silent');
x = PIE.x;      z = PIE.T*x;
f = innerprod(z,z)*innerprod(x,x) - innerprod(z.*z,z.*z);
d0 = struct('in',[0;0],'out',0);
for im = {'full','all'}
    prog = piesos_program(x);
    [prog,V] = piesos_sosvar(prog,dmonomials(x,2),struct('degs',d0,'otype',3,'intmat',im{1}));
    ok = is_feasible(piesos_eq(prog,f-V));
    assert(ok==strcmp(im{1},'all'),'embedding inequality, intmat %s: feasible=%d',im{1},ok);
end
ok = is_feasible(piesos_ineq(piesos_program(x),f,struct('otype',3,'intmat','all','degs',d0)));
assert(ok,'embedding inequality through piesos_ineq with intmat all: infeasible');
end

function case_psatz_regions()
% psatz on the faces of a region: on the triangle s_2<=s_1, P = s_1-s_2 is
% nonnegative (not on the square), so int int_{s_2<s_1} (s_1-s_2)
% x(s_1)^2 x(s_2)^2 is certified with psatz 1 (s_1-s_2 is a face), not
% without; a constant P stays certified with psatz (the products of the
% faces sum to (b-a)^2). Multiplier rows, otype 2. Also V at random
% decision variables against the direct evaluation, with psatz 1
pvar s s_1 s_2
x = polyopvar('x',s,[0,1]);
X = struct('x',1+2*s-s^2);
Z = dmonomials(x,2);
u = otimes(x,x,2);      Cop4 = otimes(u.C.ops{1,1},u.C.ops{1,1},[1,2]);
f1 = innerprod(u,u);    f1.C.ops{1} = int(Cop4,s_1-s_2,[s_1;s_2],[0,1;0,s_1]);
f0 = innerprod(u,u);    f0.C.ops{1} = int(Cop4,[],[s_1;s_2],[0,1;0,s_1]);
o = struct('degs',struct('in',[0;0],'out',1),'otype',[2,2],'exclude',[0;1;1],'intmat',[0,3,0,1]);
% f, psatz, expected feasibility
specs = {f1,0,false; f1,1,true; f0,1,true};
for c=1:size(specs,1)
    oc = o;
    if specs{c,2}>0,    oc.psatz = specs{c,2};      end
    [prog,V] = piesos_sosvar(piesos_program(x),Z,oc);
    ok = is_feasible(piesos_eq(prog,specs{c,1}-V));
    assert(ok==specs{c,3},'psatz on the triangle, case %d: feasible=%d, expected %d',c,ok,specs{c,3});
end
o.psatz = 1;    o.intmat = 'all';
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(x),Z,o);
check_V_regions(V,Pmat,Zop,Z,X,{[0,3,0,1;0,3,1,3]},'triangles with psatz 1',95);
end
% DJ (CC), 11/19/2026: end.


%% ---- helpers ----------------------------------------------------------

function case_values_factor_orders()                                        % DJ (CC), 11/20/2026
% A basis with both orders of a monomial of two states, x1 o x2 and
% x2 o x1 (dmonomials with 'full'), otype 3: by default (opts.unique) the
% pointwise rows of x1 o x2 only (blocks: pointwise and tensor of x1 o x2,
% tensor of x2 o x1), with unique=false those of both (4 blocks). V at
% random decision variables against <Zop*z,Pmat*Zop*z>, each block applied
% to the factors of its entry in their order (apply_side), over the square
% and over the regions 'all'. With otype 1, the entry x2 o x1 is left out;
% a monomial of a single state has one order only
pvar s
xs = polyopvar({'x1','x2'},[1;1],{'s'},[0,1],true(2,1));
X = struct('x1',1+2*s-s^2,'x2',s+3*s^3);
Z = dmonomials(xs,[1 1],'full');
assert(isequal(sortrows(full(Z.Z)),[1 2;2 1]),'dmonomials full: the basis must be x1 o x2 and x2 o x1.');
for c = {{[0;0],[0;0;0]},{[1;1],[1;0;1]}}
    opts = struct('degs',struct('in',c{1}{1},'out',1),'otype',3,'exclude',c{1}{2});
    label = sprintf('x1 o x2, x2 o x1, otype 3, degs.in %s',mat2str(c{1}{1}));
    [~,V,Pmat,Zop] = piesos_sosvar(piesos_program(xs),Z,opts);
    assert(isequal(size(Pmat),[3,3]),'%s: Pmat has %s blocks, expected [3 3].',label,mat2str(size(Pmat)));
    check_V_random(V,Pmat,Zop,Z,X,[label,', unique'],81);
    opts.unique = false;
    [~,V,Pmat,Zop] = piesos_sosvar(piesos_program(xs),Z,opts);
    assert(isequal(size(Pmat),[4,4]),'%s, not unique: Pmat has %s blocks, expected [4 4].',label,mat2str(size(Pmat)));
    check_V_random(V,Pmat,Zop,Z,X,[label,', not unique'],82);
end
opts = struct('degs',struct('in',[0;0],'out',0),'otype',3,'intmat','all');
[~,V,Pmat,Zop] = piesos_sosvar(piesos_program(xs),Z,opts);
check_V_presets(V,Pmat,Zop,Z,X,'all','x1 o x2, x2 o x1, otype 3, regions all',83);
% otype 1: x2 o x1 has no rows left, and is left out
[~,~,Pmat,Zop] = piesos_sosvar(piesos_program(xs),Z,struct('degs',0,'otype',1));
assert(isequal(size(Pmat),[1,1]) && size(Zop.ops,2)==1,'otype 1: Pmat must have 1 block, for x1 o x2 only.');
% a single state: one order
Z11 = dmonomials(xs,[2 0],'full');
assert(size(Z11.Z,1)==1,'x1 o x1 has a single order.');
end

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

function check_V_random(V,Pmat,Zop,Z,X,label,seed,ZopR,dom)                 % DJ (CC), 11/19/2026
% V at the states X (struct: state name -> polynomial) and random values of the
% decision variables of Pmat, against the direct evaluation of <Zop*z,Pmat*ZopR*z>:
% z the monomials of Z, Zop applied by apply_op, Pmat evaluated at the decision
% variables, and the integral over all copies si_j of the output variables (on
% DOM, default [0,1]; the output si of a pointwise block row being copy si_1),
% also those on which the product does not depend, by polynomial integration.
% ZopR defaults to Zop
if nargin<8 || isempty(ZopR),   ZopR = Zop;     end                         % DJ (CC), 11/19/2026
if nargin<9,    dom = [0,1];    end                                         % DJ (CC), 11/19/2026
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
YR = apply_side(ZopR,Z,X);                                                  % DJ (CC), 11/19/2026
[q1,q2] = size(Pmat);
rows = cell(q1,1);
for i=1:q1
    blks = cell(1,q2);
    for j=1:q2
        blks{j} = eval_dpvar_at(Pmat{i,j},d0map);
    end
    rows{i} = horzcat(blks{:});
end
%F = Y.'*vertcat(rows{:})*Y;                                                % DJ (CC), 11/19/2026 (was)
%vn = F.varname;                                                            % DJ (CC), 11/19/2026 (was)
%for k=1:numel(vn)                                                          % DJ (CC), 11/19/2026 (was)
%    sk = polynomial(vn(k));                                                % DJ (CC), 11/19/2026 (was)
%    F = int(F,sk(1),0,1);                                                  % DJ (CC), 11/19/2026 (was)
%end                                                                        % DJ (CC), 11/19/2026 (was)
F = Y.'*vertcat(rows{:})*YR;                                                % DJ (CC), 11/19/2026
% DJ (CC), 11/19/2026: start. All copies, also those absent from F (a factor
% of the domain length): the number of copies of each output variable is the
% maximum over the block rows of Zop, ZopR
vars = Zop.vars.out;
ncop = max([Zop.depmat.out;ZopR.depmat.out],[],1);
for i=1:numel(vars)
    if ismember(vars{i},F.varname)
        F = subs(F,polynomial(vars(i)),polynomial({[vars{i},'_1']}));
    end
    for j=1:ncop(i)
        F = int(F,polynomial({sprintf('%s_%d',vars{i},j)}),dom(1),dom(2));
    end
end
assert(isempty(F.varname),'%s: the direct evaluation depends on %s',label,strjoin(F.varname,','));
% DJ (CC), 11/19/2026: end.
orc = double(F);
assert(abs(orc)>1e-6,'%s: the value %g is too small to test.',label,orc);
assert(abs(val-orc)<=1e-8*abs(orc),'%s: V=%.10g, direct evaluation=%.10g (relative error %.2e)',label,val,orc,abs(val-orc)/abs(orc));
end

% DJ (CC), 11/19/2026: start. Helpers for regions of integration.
function check_V_regions(V,Pmat,Zop,Z,X,regions,label,seed)
% V at random decision variable values against the sum over the regions k
% of <Zop*z,Pmat{k}*Zop*z>, integrated over region k (rows (k1,...,kN) of
% regions{1},...,regions{N}, k the linear index) by int_rows
rng(seed);
Preg = regions_of(Pmat);
d0map = rand_dvars([Preg{:}]);
nreg = [cellfun(@(A) size(A,1),regions),1];
rows = cell(1,numel(nreg));
orc = 0;
Y = apply_side(Zop,Z,X);
for k=1:prod(nreg)
    [rows{:}] = ind2sub(nreg,k);
    P = cellfun(@(Pb) eval_dpvar_at(Pb,d0map),Preg{k},'UniformOutput',false);
    Prows = cell(size(P,1),1);
    for b=1:size(P,1),  Prows{b} = horzcat(P{b,:});    end
    F = Y.'*vertcat(Prows{:})*Y;
    [F,cps,dom] = copies_of(F,Zop);
    for i=1:numel(cps)
        F = int_rows(F,cps{i},dom(i,:),regions{i}(rows{i},:));
    end
    orc = orc + double(F);
end
compare(eval_V(V,d0map,X),orc,label);
end

function check_V_presets(V,Pmat,Zop,Z,X,presets,label,seed)
% V with the same random constant P in each region against <Zop*z,P*Zop*z>
% integrated over the cells of the preset of each variable (int_cells)
rng(seed);
Preg = regions_of(Pmat);
n = sum(cellfun(@(Pb) size(Pb,1),Preg{1}(:,1)));
T = randn(n);       T = T+T.';
d0map = containers.Map('KeyType','char','ValueType','any');
for k=1:numel(Preg)
    m = gram_values(Preg{k},blocks4(T,Preg{k}));
    for key=keys(m),    d0map(key{1}) = m(key{1});  end
end
Y = apply_side(Zop,Z,X);
[F,cps,dom] = copies_of(Y.'*T*Y,Zop);
if ~iscell(presets),    presets = repmat({presets},1,numel(cps));   end
for i=1:numel(cps)
    F = int_cells(F,cps{i},dom(i,:),presets{i});
end
compare(eval_V(V,d0map,X),double(F),label);
end

function F = int_rows(F,cps,dom,row)
% F integrated over the region ROW (intmat of an 'intop') of the copies CPS
% on DOM, the last copy first: copy j at copy p for [p,p], else over
% [L,U], with 0 for dom(1), d+1 for dom(2) and p for copy p
d = numel(cps);
lim = [{dom(1)},arrayfun(@(j) cps(j),1:d,'UniformOutput',false),{dom(2)}];
for j=d:-1:1
    if row(2*j-1)==row(2*j)
        F = subs(F,cps(j),lim{row(2*j)+1});
    else
        F = int(F,cps(j),lim{row(2*j-1)+1},lim{row(2*j)+1});
    end
end
end

function G = int_cells(F,cps,dom,preset)
% F integrated over the regions of a preset, built independently of
% piesos_quadvar: the whole domain (the strict orderings partition it),
% plus for 'all' each ordering with ties but not all copies tied, from the
% dense ranks r of the copies (copy j in group r(j), groups increasing):
% the copies of group m set to u_m, over dom(1) <= u_1 <= ... <= u_k <= dom(2)
d = numel(cps);
G = F;
for j=1:d,  G = int(G,cps(j),dom(1),dom(2));    end
if ~strcmp(preset,'all'),   return;     end
idx = cell(1,d);    [idx{:}] = ndgrid(1:d);
R = reshape(cat(d+1,idx{:}),[],d);
for t=1:size(R,1)
    k = max(R(t,:));
    if k<2 || k>d-1 || ~isequal(unique(R(t,:)),1:k),    continue;   end
    u = polynomial(arrayfun(@(m) sprintf('u_%d',m),1:k,'UniformOutput',false));
    H = subs(F,cps,reshape(u(R(t,:)),d,1));
    for m=1:k-1,    H = int(H,u(m),dom(1),u(m+1));  end
    G = G + int(H,u(k),dom(1),dom(2));
end
end

function [F,cps,dom] = copies_of(F,Zop)
% The copies of each output variable of Zop (the copies in all block rows),
% with the bare name of a pointwise output in F renamed to its copy 1
vars = Zop.vars.out;        dom = Zop.dom.out;
ncop = full(max(Zop.depmat.out,[],1));
cps = cell(1,numel(vars));
for i=1:numel(vars)
    cps{i} = get_tensorvars(vars{i},ncop(i),{'full','pvar'});
    if ismember(vars{i},F.varname)
        F = subs(F,polynomial(vars(i)),cps{i}(1));
    end
end
end

function P = regions_of(Pmat)
% The Pmat of each region, as a cell, also for a single Pmat (whole domain)
if iscell(Pmat{1}),     P = Pmat;
else,                   P = {Pmat};
end
end

function d0map = rand_dvars(Pmat)
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)];                                    %#ok<AGROW>
end
names = unique(names,'stable');
d0map = containers.Map(names,num2cell(0.5+rand(numel(names),1)));
end

function target = blocks4(T,Pmat)
% the blocks of T as target(b,c,i,j), for gram_values
szs = cellfun(@(Pb) size(Pb,1),Pmat(:,1));      cs = [0;cumsum(szs)];
nb = numel(szs);    target = zeros(nb,nb,max(szs),max(szs));
for b=1:nb
    for c=1:nb
        target(b,c,1:szs(b),1:szs(c)) = T(cs(b)+1:cs(b+1),cs(c)+1:cs(c+1));
    end
end
end

function val = eval_V(V,d0map,X)
% V at the decision variable values d0map and the states X
Vn = V;     ops = V.C.ops;
for k=1:numel(ops)
    if isa(ops{1,k},'intop'),   ops{1,k} = eval_intop_dpvar(ops{1,k},d0map);   end
end
Vn.C.ops = ops;     Vn.C.chkval = true;
xv = cell(numel(Vn.varname),1);
for k=1:numel(xv),  xv{k} = X.(Vn.varname{k});     end
val = double(eval_polyopvar(Vn,xv));
end

function compare(val,orc,label)
assert(abs(orc)>1e-6,'%s: the value %g is too small to test.',label,orc);
assert(abs(val-orc)<=1e-8*abs(orc),'%s: V=%.10g, direct evaluation=%.10g (relative error %.2e)',label,val,orc,abs(val-orc)/abs(orc));
end
% DJ (CC), 11/19/2026: end.

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

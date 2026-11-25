function Test_tensopvar_int_weights()
% TEST_TENSOPVAR_INT_WEIGHTS tests, using assertions, the weights W(s_1,...)
% of @tensopvar/int.m that depend on the copies of the output variables,
%   KCop = int(Cop,W) = int_Omega W*Cop dOmega,
% for an output of degree d>1 (a tensor product of d copies, or a pointwise
% product of such outputs), with W in the variables 'si_j', the jth copy of
% the output variable si, and for the other names rejected.
%
% Ground truth, independent of int.m: the functional
%   V = innerprod(u,P,u),   u = x1 o x2 (the tensor product x1(s_1)*x2(s_2)),
% which passes P to int(otimes(C1,C2,[1,2]),P), equals
%   V = int int P(s_1,s_2) x1(s_1)^2 x2(s_2)^2 ds_1 ds_2
% (copy j of the output is factor j of u); evaluated at polynomial states
% by eval_polyopvar (apply_op), and compared with the polynomial integral of
% the integrand. States x1, x2 differ, so that a swap of the copies is
% detected: the weights depend on the copies asymmetrically. Also for two
% spatial variables, copies si_1, si_2 of each of s1, s2 (4 integrals).
% No programs are solved.
%
% Variable limits in dom (case_var_limits*): KCop = int(Cop,W,vars,dom)
% applied to polynomial factors by apply_op of 'intop', against the
% iterated polynomial integral of W*apply_op(Cop,X) with the same limits.
% Neither oracle uses int.m, intop.canonicalize_intmat or the composition
% mtimes_tensopvar.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_tensopvar_int_weights
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
% DJ (CC), 11/10/2026: The only copy of a pointwise output is also named
%                        si_1: case_names;
% DJ (CC), 11/19/2026: The inputs vars, dom of int: case_vars_dom;
% DJ (CC), 11/19/2026: Variable limits in dom: case_var_limits,
%                        case_var_limits_PI, case_var_limits_other,
%                        case_var_limits_errors; case_vars_dom expects the
%                        error of a variable in its own limits;
% DJ (CC), 11/19/2026: Equal variable limits (Dirac delta):
%                        case_var_limits_delta, and in case_var_limits_PI;
%                        edge cases: case_var_limits_edge; a limit must be
%                        named as in vars;
% DJ (CC), 11/22/2026: The form int(Cop,W,intmat): case_intmat_input;

% NOTES
% - The test needs eval_polyopvar (Testfolder/polyopvar), and
%   eval_intop_dpvar (Testfolder/intop);
% - The copies are named 'si_j'; the only copy of a pointwise output (degree
%   1, otype(1)=1) is also named 'si', as before. A single-factor tensopvar
%   has otype [2,2] by default, so its copy is named 'si_1' only.

tests = {@case_values_1D, @case_values_2D, @case_names, @case_no_weight_unchanged, @case_vars_dom, ...
         @case_var_limits, @case_var_limits_PI, @case_var_limits_other, @case_var_limits_delta, ...
         @case_var_limits_edge, @case_var_limits_random, @case_var_limits_errors, ...  % DJ (CC), 11/19/2026
         @case_intmat_input};                                               % DJ (CC), 11/22/2026
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
    error('Test_tensopvar_int_weights: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_tensopvar_int_weights: all %d cases passed.\n',numel(tests));
end


%% ---- cases ------------------------------------------------------------

function case_values_1D()
% x1(s), x2(s) on [0,1]; weights in the copies s_1, s_2
pvar s s_1 s_2
x1 = polyopvar('x1',s,[0,1]);       x2 = polyopvar('x2',s,[0,1]);
u = otimes(x1,x2,2);
X = struct('x1',1+2*s-s^2,'x2',s+3*s^3);
x1c = subs(X.x1,s,s_1);             x2c = subs(X.x2,s,s_2);
Ps = {s_1, s_2, s_1*s_2^2+3*s_2, 2-s_1^2*s_2+s_2^3, polynomial(3)};
for k=1:numel(Ps)
    P = Ps{k};
    V = innerprod(u,P,u);
    val = eval_at(V,X);
    orc = double(int(int(P*x1c^2*x2c^2,s_1,0,1),s_2,0,1));
    check(val,orc,sprintf('N=1, weight %d',k));
end
% a weight in the copy of the other order is a different functional
V = innerprod(u,s_1,u);
val = eval_at(V,X);
orc_swapped = double(int(int(s_2*x1c^2*x2c^2,s_1,0,1),s_2,0,1));
assert(abs(val-orc_swapped)>1e-6*abs(orc_swapped),'weights s_1 and s_2 must give different functionals');
end

function case_values_2D()
% x1(s1,s2), x2(s1,s2): factor j has the copies (s1_j,s2_j) of (s1,s2)
pvar s1 s2 s1_1 s1_2 s2_1 s2_2
dom = [0,1;0,1];
x1 = polyopvar('x1',[s1;s2],dom);   x2 = polyopvar('x2',[s1;s2],dom);
u = otimes(x1,x2,2);
X = struct('x1',1+s1+2*s1*s2,'x2',s2+s1^2);
x1c = subs(X.x1,[s1;s2],[s1_1;s2_1]);       x2c = subs(X.x2,[s1;s2],[s1_2;s2_2]);
Ps = {s1_1*s2_2+s1_2^2*s2_1, s2_1, s1_2};
for k=1:numel(Ps)
    P = Ps{k};
    V = innerprod(u,P,u);
    val = eval_at(V,X);
    F = P*x1c^2*x2c^2;
    for v = {s1_1,s1_2,s2_1,s2_2}
        F = int(F,v{1},0,1);
    end
    check(val,double(F),sprintf('N=2, weight %d',k));
end
end

function case_names()
pvar s s_1 s_2 s_3 t_1
x1 = polyopvar('x1',s,[0,1]);       x2 = polyopvar('x2',s,[0,1]);
u = otimes(x1,x2,2);
w = @(p) polynomial(p);
% Kronecker product of two copies (otype(1)=2, degree 2): copies s_1, s_2
Cop2 = u.C.ops{1,1};
assert(Cop2.otype(1)==2 && Cop2.deps.out(1)==2,'unexpected output of x1 o x2');
int(Cop2,s_1);      int(Cop2,s_2);      int(Cop2,s_1+s_2^2);    int(Cop2,w(2));
assert_err(@() int(Cop2,s),'si_j');
assert_err(@() int(Cop2,s_3),'si_j');
assert_err(@() int(Cop2,s_1+t_1),'si_j');
% pointwise product of such operators (degree 2, otype(1)=1), as in innerprod
Cop3 = otimes(Cop2,Cop2,[1,2]);
assert(Cop3.otype(1)==1 && Cop3.deps.out(1)==2,'unexpected output of the product');
int(Cop3,s_1*s_2);
assert_err(@() int(Cop3,s),'si_j');
% one copy of a pointwise output (degree 1, otype(1)=1): 's' or 's_1', the same
% variable, but not both
Cop1 = x1.C.ops{1,1};
Cop1p = otimes(Cop1,Cop1,[1,2]);
assert(Cop1p.otype(1)==1 && Cop1p.deps.out(1)==1,'unexpected output of the pointwise product');
int(Cop1p,s);       int(Cop1p,s_1);     int(Cop1p,s+s^2);       int(Cop1p,s_1+s_1^2);
assert_err(@() int(Cop1p,s+s_1),'both');
assert_err(@() int(Cop1p,s_2),'si_j');
% a single factor of otype(1)=2 (the default of a tensopvar): 's_1' only
if Cop1.otype(1)==2
    int(Cop1,s_1);
    assert_err(@() int(Cop1,s),'si_j');
else
    int(Cop1,s);        int(Cop1,s_1);
end
end

function case_no_weight_unchanged()
% Without weights, or with a constant, int is unchanged: (int x1^2)(int x2^2)
pvar s s_1
x1 = polyopvar('x1',s,[0,1]);       x2 = polyopvar('x2',s,[0,1]);
u = otimes(x1,x2,2);
X = struct('x1',1+2*s-s^2,'x2',s+3*s^3);
orc = double(int(X.x1^2,s,0,1))*double(int(X.x2^2,s,0,1));
check(eval_at(innerprod(u,u),X),orc,'no weight');
check(eval_at(innerprod(u,2,u),X),2*orc,'constant weight');
% a pointwise weight on a degree-1 state is named 's', or 's_1'
V = innerprod(x1,s,x1);
check(eval_at(V,X),double(int(s*X.x1^2,s,0,1)),'pointwise weight s');
V = innerprod(x1,s_1^2+s_1,x1);
check(eval_at(V,X),double(int((s^2+s)*X.x1^2,s,0,1)),'pointwise weight s_1');
end

function case_vars_dom()
% int(Cop,W,vars,dom): W is integrated over the variables of vars on which
% the output of Cop does not depend; each copy of the output must be in
% vars once, over its full domain. KCop is evaluated as the coefficient of
% V = innerprod(.,.) (a single term), against the polynomial integral of
% the integrand. Domain [0,2], so that a missing length is detected, and
% x1 not symmetric on it, so that a dropped weight s_1 is detected
pvar s s_1 s_2 s_3 t
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
X = struct('x1',1+s+2*s^2,'x2',s+3*s^3);
x1c = subs(X.x1,s,s_1);
% one copy of a pointwise output (x1^2), with the variables s_2, s_3, t
% on which it does not depend
Cop1 = x1.C.ops{1,1};               Cop1p = otimes(Cop1,Cop1,[1,2]);
V = innerprod(x1,x1);
assert(numel(V.C.ops)==1,'innerprod(x1,x1) must have a single term');
W = 1+s_1*s_2^2+s_3*t;
vars = {'s_1','s_2','s_3','t'};     dom = [0,2;0,2;0,2;-1,2];
F = W*x1c^2;
for k=1:4
    F = int(F,polynomial(vars(k)),dom(k,1),dom(k,2));
end
V.C.ops{1} = int(Cop1p,W,vars,dom);
check(eval_at(V,X),double(F),'pointwise output, weight in other variables');
% the same with pvar input, and the name s of the only copy
V.C.ops{1} = int(Cop1p,subs(W,s_1,s),[s;s_2;s_3;t],dom);
check(eval_at(V,X),double(F),'pvar input, copy named s');
% no weight, or a constant: multiplied by the length 3 of the domain of t
orc = double(int(X.x1^2,s,0,2));
V.C.ops{1} = int(Cop1p,[],{'s','t'},[0,2;-1,2]);
check(eval_at(V,X),3*orc,'no weight, other variable t');
V.C.ops{1} = int(Cop1p,2,{'t','s_1'},[-1,2;0,2]);
check(eval_at(V,X),6*orc,'constant weight, other variable t');
% default domain of the copies of the output
V.C.ops{1} = int(Cop1p,s_1,'s_1');
check(eval_at(V,X),double(int(s*X.x1^2,s,0,2)),'default domain');
% tensor product of two copies (x1 o x2)^2, with the copy s_3: the same as
% the weight integrated over s_3 first, without vars (see case_values_1D).
% V.C of innerprod(u,u) orders the factors x1,x1,x2,x2, not as Cop3, so V
% is not compared with the integral of the integrand
u = otimes(x1,x2,2);
Cop3 = otimes(u.C.ops{1,1},u.C.ops{1,1},[1,2]);
V = innerprod(u,u);
assert(numel(V.C.ops)==1,'innerprod(u,u) must have a single term');
W = s_1+s_2*s_3^2;
V.C.ops{1} = int(Cop3,W,{'s_1','s_2','s_3'},[0,2;0,2;0,2]);
val = eval_at(V,X);
V.C.ops{1} = int(Cop3,int(W,s_3,0,2));
check(val,eval_at(V,X),'tensor product, other copy s_3');
% W in s_3 is rejected if s_3 is not in vars
assert_err(@() int(Cop3,W,{'s_1','s_2'},[0,2;0,2]),'si_j');
% errors
assert_err(@() int(Cop1p,[],{'s_2','t'},[0,2;-1,2]),'each copy');
assert_err(@() int(Cop1p,[],{'s','s_1'},[0,2;0,2]),'each copy');
assert_err(@() int(Cop3,[],{'s_1'},[0,2]),'each copy');
assert_err(@() int(Cop1p,[],{'s'},[0,1]),'full domain');
assert_err(@() int(Cop1p,[],{'s','t'},{0,2;0,'t'}),'own limits');          % DJ (CC), 11/19/2026
assert_err(@() int(Cop1p,[],{'s','t'}),'domain of variable');
assert_err(@() int(Cop1p,t,{'s'},[0,2]),'si_j');
end

% DJ (CC), 11/19/2026: start. Cases for variable limits in dom.
function case_var_limits()
% Limits of the copies of the output in other copies (the intmat of Kop),
% for factors without integrals (x1 o x2): the triangles s_2<=s_1 and
% s_2>=s_1, in three syntaxes of dom and also with the limit on s_1; three
% copies (a chain, and a copy between two copies); two spatial variables.
% States and weights are not symmetric in the copies, and the domain is
% [0,2], so that a limit on the wrong copy or side, or a dropped length,
% is detected
pvar s s_1 s_2 s_3
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
u = otimes(x1,x2,2);                Cop2 = u.C.ops{1,1};
X = {1+2*s-s^2, s+3*s^3};
Y = apply_op(Cop2,X);               % x1(s_1)*x2(s_2)
W = 1+s_1+2*s_2^2;
lo = double(int(int(W*Y,s_2,0,s_1),s_1,0,2));
up = double(int(int(W*Y,s_2,s_1,2),s_1,0,2));
sq = double(int(int(W*Y,s_2,0,2),s_1,0,2));
check(apply_val(int(Cop2,W,{'s_1','s_2'},{0,2;0,'s_1'}),X),lo,'lower triangle, names in a cell');
check(apply_val(int(Cop2,W,[s_1;s_2],[0,2;0,s_1]),X),lo,'lower triangle, polynomial dom');
check(apply_val(int(Cop2,W,[s_1;s_2],{0,2;s_1,2}),X),up,'upper triangle, pvar in a cell');
check(apply_val(int(Cop2,W,[s_1;s_2],[s_2,2;0,2]),X),lo,'lower triangle as s_1 in [s_2,2]');
check(apply_val(int(Cop2,W,[s_2;s_1],[0,s_1;0,2]),X),lo,'lower triangle, vars in other order');
check(lo+up,sq,'oracle: the triangles add up to the square');
check(apply_val(int(Cop2,W,[s_1;s_2],[0,2;0,s_1]),X)+apply_val(int(Cop2,W,[s_1;s_2],[0,2;s_1,2]),X),sq,'triangles add up to the square');
% no weight, and a constant weight
lo0 = double(int(int(Y,s_2,0,s_1),s_1,0,2));
check(apply_val(int(Cop2,[],[s_1;s_2],[0,2;0,s_1]),X),lo0,'lower triangle, no weight');
check(apply_val(int(Cop2,3,[s_1;s_2],[0,2;0,s_1]),X),3*lo0,'lower triangle, constant weight');
% three copies: the chain s_3<=s_2<=s_1, and s_1<=s_3<=s_2 (s_3 between two
% copies, ordered as s_2 has lower limit s_1)
x3 = polyopvar('x3',s,[0,2]);
u3 = otimes(u,x3,2);                Cop3 = u3.C.ops{1,1};
assert(Cop3.deps.out(1)==3,'unexpected output of x1 o x2 o x3');
X3 = [X,{2-s+s^2}];
Y3 = apply_op(Cop3,X3);
W3 = 1+s_1+s_2^2+2*s_3;
orc = int(int(int(W3*Y3,s_3,0,s_2),s_2,0,s_1),s_1,0,2);
check(apply_val(int(Cop3,W3,[s_1;s_2;s_3],[0,2;0,s_1;0,s_2]),X3),double(orc),'three copies, chain');
orc = int(int(int(W3*Y3,s_3,s_1,s_2),s_2,s_1,2),s_1,0,2);
check(apply_val(int(Cop3,W3,[s_1;s_2;s_3],[0,2;s_1,2;s_1,s_2]),X3),double(orc),'three copies, s_3 between s_1 and s_2');
% two spatial variables r1, r2, with copies r1_1, r1_2, r2_1, r2_2: a
% triangle in each, on different sides
pvar r1 r2 r1_1 r1_2 r2_1 r2_2
y1 = polyopvar('y1',[r1;r2],[0,1;0,2]);     y2 = polyopvar('y2',[r1;r2],[0,1;0,2]);
uy = otimes(y1,y2,2);               CopY = uy.C.ops{1,1};
XY = {1+r1+2*r1*r2, r2+r1^2};
YY = apply_op(CopY,XY);
WY = 1+r1_1*r2_2+r1_2^2;
F = int(int(WY*YY,r1_2,0,r1_1),r1_1,0,1);
F = int(int(F,r2_1,r2_2,2),r2_2,0,2);
KCop = int(CopY,WY,[r1_1;r1_2;r2_1;r2_2],[0,1;0,r1_1;r2_2,2;0,2]);
check(apply_val(KCop,XY),double(F),'two spatial variables, a triangle in each');
end

function case_var_limits_PI()
% The triangles for the products of the block rows of a basis of otype 3
% with integral terms (degs.in [1;1]), as piesos_quadvar forms them: a
% pointwise row (output in copy s_1 only) and a tensor-product row, in each
% order. The integrals of Kop then compose with those of the factors
% (switch_ints in mtimes_tensopvar). Also the delta s_2=s_1 (the diagonal).
% The output s of a pointwise factor is copy s_1
pvar s s_1 s_2
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,2],'out',[0,2]);
Zop = monopmat(1,vars,dom,struct('in',2,'out',true),struct('in',[1;1],'out',0),struct('otype',3));
assert(isequal(full(Zop.depmat.out(:)),[1;2]),'expected a pointwise and a tensor-product row');
ZP = fix_empty_params(Zop.ops{1,1});        ZT = fix_empty_params(Zop.ops{2,1});
X = {1+2*s-s^2, s+3*s^3, 2-s, 1+s^2};       % a different state per factor
pairs = {ZP,ZP,'pointwise x pointwise'; ZP,ZT,'pointwise x tensor'; ...
         ZT,ZP,'tensor x pointwise'; ZT,ZT,'tensor x tensor'};
for k=1:size(pairs,1)
    A = pairs{k,1};     B = pairs{k,2};
    m = size(A,1);      n = size(B,1);
    Cop = otimes(repmat(A,n,1),repelem(B,m,1),[1,2]);
    Y = apply_op(Cop,X);
    if ismember('s',Y.varname),     Y = subs(Y,s,s_1);  end
    p = size(Cop,1);
    W = (1+s_1^2+3*s_2)*ones(1,p) + (s_1*s_2)*(1:p)/p;
    lo = int(int(W*Y,s_2,0,s_1),s_1,0,2);
    up = int(int(W*Y,s_2,s_1,2),s_1,0,2);
    check(apply_val(int(Cop,W,[s_1;s_2],[0,2;0,s_1]),X),double(lo),[pairs{k,3},', lower triangle']);
    check(apply_val(int(Cop,W,[s_1;s_2],[0,2;s_1,2]),X),double(up),[pairs{k,3},', upper triangle']);
    dg = int(subs(W*Y,s_2,s_1),s_1,0,2);
    check(apply_val(int(Cop,W,[s_1;s_2],[0,2;s_1,s_1]),X),double(dg),[pairs{k,3},', diagonal']);
    % the lower triangle as s_1 in [s_2,2]: for a pointwise x pointwise
    % product, s_2 is not a copy of the output, so it cannot bound s_1
    if k==1
        assert_err(@() int(Cop,W,[s_1;s_2],[s_2,2;0,2]),'full domain');
    else
        check(apply_val(int(Cop,W,[s_1;s_2],[s_2,2;0,2]),X),double(lo),[pairs{k,3},', lower triangle as s_1 in [s_2,2]']);
    end
end
end

function case_var_limits_other()
% Variable limits of the other variables, over which W is integrated first:
% a copy as limit, no weight, a chain of other variables listed in another
% order, a dpvar weight, and a variable between two copies of a
% tensor-product output (a signed integral)
pvar s s_1 s_2 t r
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
X1 = 1+s+2*s^2;                     x1c = subs(X1,s,s_1);
Cop1 = x1.C.ops{1,1};               Cop1p = otimes(Cop1,Cop1,[1,2]);    % x1(s)^2
W = 1+s_1*s_2^2+2*s_2;
orc = double(int(int(W,s_2,0,s_1)*x1c^2,s_1,0,2));
check(apply_val(int(Cop1p,W,[s_1;s_2],[0,2;0,s_1]),{X1,X1}),orc,'pointwise output, s_2 in [0,s_1]');
check(apply_val(int(Cop1p,[],[s_1;s_2],[0,2;0,s_1]),{X1,X1}),double(int(s*X1^2,s,0,2)),'no weight, s_2 in [0,s_1]');
% t in [0,s] and r in [t,2]: r is integrated first, then t
Wc = 1+s*t+r^2*t;
orc = double(int(int(int(Wc,r,t,2),t,0,s)*X1^2,s,0,2));
check(apply_val(int(Cop1p,Wc,[r;s;t],[t,2;0,2;0,s]),{X1,X1}),orc,'chain of other variables');
% dpvar weight, evaluated at values of the decision variables
c1 = dpvar('c1');       c2 = dpvar('c2');
KCop = int(Cop1p,c1*s_2+c2*s_1*s_2,[s_1;s_2],[0,2;s_1,2]);
val = apply_val(eval_intop_dpvar(KCop,containers.Map({'c1','c2'},{0.7,-1.3})),{X1,X1});
orc = double(int(int(0.7*s_2-1.3*s_1*s_2,s_2,s_1,2)*x1c^2,s_1,0,2));
check(val,orc,'dpvar weight, s_2 in [s_1,2]');
% t between the copies s_2 and s_1 of x1 o x2, over the square
u = otimes(x1,x2,2);                Cop2 = u.C.ops{1,1};
X = {X1, s+3*s^3};                  Y = apply_op(Cop2,X);
Wt = 1+t*s_1;
orc = double(int(int(int(Wt,t,s_2,s_1)*Y,s_2,0,2),s_1,0,2));
check(apply_val(int(Cop2,Wt,[s_1;s_2;t],[0,2;0,2;s_2,s_1]),X),orc,'other variable between two copies');
end

function case_var_limits_delta()
% Equal variable limits [t,t] evaluate at t (Dirac delta): a copy at another
% copy (the diagonal), in a chain with an integral, for two spatial
% variables, and an other variable at a copy (W evaluated, also a dpvar W),
% against the polynomial W*apply_op(Cop,X) evaluated there and integrated
% over the remaining copies. Equal numbers are not an evaluation: a copy at
% a or b is rejected, and an other variable over [1,1] gives zero
pvar s s_1 s_2 s_3 t
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
x3 = polyopvar('x3',s,[0,2]);
u = otimes(x1,x2,2);                Cop2 = u.C.ops{1,1};
X = {1+2*s-s^2, s+3*s^3};
Y = apply_op(Cop2,X);
W = 1+s_1+2*s_2^2;
F = W*Y;
check(apply_val(int(Cop2,W,[s_1;s_2],[0,2;s_1,s_1]),X),double(int(subs(F,s_2,s_1),s_1,0,2)),'s_2 at s_1');
check(apply_val(int(Cop2,W,[s_1;s_2],{'s_2','s_2';0,2}),X),double(int(subs(F,s_1,s_2),s_2,0,2)),'s_1 at s_2, names in a cell');
assert_err(@() int(Cop2,W,[s_1;s_2],[0,2;0,0]),'full domain');
assert_err(@() int(Cop2,W,[s_1;s_2],[2,2;0,2]),'full domain');
% three copies: s_3 at s_2, with s_2 in [0,s_1]
u3 = otimes(u,x3,2);                Cop3 = u3.C.ops{1,1};
X3 = [X,{2-s+s^2}];                 Y3 = apply_op(Cop3,X3);
W3 = 1+s_1+s_2^2+2*s_3;
orc = int(int(subs(W3*Y3,s_3,s_2),s_2,0,s_1),s_1,0,2);
check(apply_val(int(Cop3,W3,[s_1;s_2;s_3],[0,2;0,s_1;s_2,s_2]),X3),double(orc),'three copies, s_3 at s_2 in a triangle');
% two spatial variables: r1_2 at r1_1, and r2_1 in [r2_2,2]
pvar r1 r2 r1_1 r1_2 r2_1 r2_2
y1 = polyopvar('y1',[r1;r2],[0,1;0,2]);     y2 = polyopvar('y2',[r1;r2],[0,1;0,2]);
uy = otimes(y1,y2,2);               CopY = uy.C.ops{1,1};
XY = {1+r1+2*r1*r2, r2+r1^2};
FY = (1+r1_1*r2_2+r1_2^2)*apply_op(CopY,XY);
FY = int(subs(FY,r1_2,r1_1),r1_1,0,1);
FY = int(int(FY,r2_1,r2_2,2),r2_2,0,2);
KCop = int(CopY,1+r1_1*r2_2+r1_2^2,[r1_1;r1_2;r2_1;r2_2],[0,1;r1_1,r1_1;r2_2,2;0,2]);
check(apply_val(KCop,XY),double(FY),'two spatial variables, a delta and a triangle');
% other variables: t at s_1 (W evaluated), also with W=[] and a dpvar W,
% and t over [1,1] (zero)
X1 = X{1};                          x1c = subs(X1,s,s_1);
Cop1 = x1.C.ops{1,1};               Cop1p = otimes(Cop1,Cop1,[1,2]);    % x1(s)^2
Wt = 1+t^2*s_1+t;
check(apply_val(int(Cop1p,Wt,[s_1;t],[0,2;s_1,s_1]),{X1,X1}),double(int(subs(Wt,t,s_1)*x1c^2,s_1,0,2)),'other variable at a copy');
check(apply_val(int(Cop1p,[],[s_1;t],[0,2;s_1,s_1]),{X1,X1}),double(int(x1c^2,s_1,0,2)),'no weight, other variable at a copy');
check(apply_val(int(Cop1p,Wt,[s_1;t],[0,2;1,1]),{X1,X1}),0,'other variable over [1,1]');
check(apply_val(int(Cop1p,3,[s_1;t],[0,2;1,1]),{X1,X1}),0,'constant weight, other variable over [1,1]');
c1 = dpvar('c1');
KCop = int(Cop1p,c1*t*s_1,[s_1;t],[0,2;s_1,s_1]);
val = apply_val(eval_intop_dpvar(KCop,containers.Map({'c1'},{0.7})),{X1,X1});
check(val,double(int(0.7*s_1^2*x1c^2,s_1,0,2)),'dpvar weight, other variable at a copy');
end

function case_var_limits_edge()
% Edge cases: dom = the full square equals the default; a string limit; a
% matrix weight (2 rows); no weight on a basis operator with integral terms
% and 9 rows (vector output); a dpvar weight in the copies over a triangle;
% three copies in a reversed chain, with a delta chain, and with a delta and
% a triangle; an other variable at another other variable; and two spatial
% variables with the variables of integration shuffled and an other
% variable between copies of both
pvar s s_1 s_2 s_3 t r
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
x3 = polyopvar('x3',s,[0,2]);
u = otimes(x1,x2,2);                Cop2 = u.C.ops{1,1};
X = {1+2*s-s^2, s+3*s^3};           Y = apply_op(Cop2,X);
W = 1+s_1+2*s_2^2;
lo = double(int(int(W*Y,s_2,0,s_1),s_1,0,2));
check(apply_val(int(Cop2,W,[s_1;s_2],[0,2;0,2]),X),apply_val(int(Cop2,W),X),'full square equals the default');
check(apply_val(int(Cop2,W,[s_1;s_2],{0,2;0,"s_1"}),X),lo,'string limit');
W2 = [1+s_1; s_2^2-s_1];
orc = [double(int(int(W2(1)*Y,s_2,0,s_1),s_1,0,2)); double(int(int(W2(2)*Y,s_2,0,s_1),s_1,0,2))];
check(apply_val(int(Cop2,W2,[s_1;s_2],[0,2;0,s_1]),X),orc,'matrix weight');
% no weight: a tensor-product basis operator with integral terms of degree 0
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,2],'out',[0,2]);
Zop = monopmat(1,vars,dom,struct('in',2,'out',true),struct('in',[0;0],'out',0),struct('otype',2));
ZT = fix_empty_params(Zop.ops{1,1});
YT = apply_op(ZT,X);
assert(size(YT,1)>1,'expected a basis operator with several rows');
check(apply_val(int(ZT,[],[s_1;s_2],[0,2;0,s_1]),X),double(int(int(YT,s_2,0,s_1),s_1,0,2)),'no weight, vector output, lower triangle');
check(apply_val(int(ZT,[],[s_1;s_2],[0,2;s_1,s_1]),X),double(int(subs(YT,s_2,s_1),s_1,0,2)),'no weight, vector output, diagonal');
% dpvar weight in the copies, over the upper triangle
c1 = dpvar('c1');       c2 = dpvar('c2');
KCop = int(Cop2,c1*s_1+c2*s_2^2,[s_1;s_2],[0,2;s_1,2]);
val = apply_val(eval_intop_dpvar(KCop,containers.Map({'c1','c2'},{0.7,-1.3})),X);
check(val,double(int(int((0.7*s_1-1.3*s_2^2)*Y,s_2,s_1,2),s_1,0,2)),'dpvar weight, upper triangle');
% three copies
u3 = otimes(u,x3,2);                Cop3 = u3.C.ops{1,1};
X3 = [X,{2-s+s^2}];                 Y3 = apply_op(Cop3,X3);
W3 = 1+s_1+s_2^2+2*s_3;             F3 = W3*Y3;
orc = int(int(int(F3,s_1,s_2,2),s_2,s_3,2),s_3,0,2);
check(apply_val(int(Cop3,W3,[s_1;s_2;s_3],[s_2,2;s_3,2;0,2]),X3),double(orc),'three copies, reversed chain');
orc = int(subs(subs(F3,s_3,s_2),s_2,s_1),s_1,0,2);
check(apply_val(int(Cop3,W3,[s_1;s_2;s_3],[0,2;s_1,s_1;s_2,s_2]),X3),double(orc),'three copies, delta chain');
orc = int(int(subs(F3,s_3,s_1),s_2,s_1,2),s_1,0,2);
check(apply_val(int(Cop3,W3,[s_1;s_2;s_3],[0,2;s_1,2;s_1,s_1]),X3),double(orc),'three copies, a delta and a triangle');
% an other variable at another other variable: t at r, r in [0,s]
Cop1 = x1.C.ops{1,1};               Cop1p = otimes(Cop1,Cop1,[1,2]);
Wr = 1+t*r+s*t^2;
orc = double(int(int(subs(Wr,t,r),r,0,s)*X{1}^2,s,0,2));
check(apply_val(int(Cop1p,Wr,[t;r;s],[r,r;0,s;0,2]),{X{1},X{1}}),orc,'other variable at another');
% two spatial variables, variables of integration shuffled, an other
% variable t between r1_1 and r2_2
pvar r1 r2 r1_1 r1_2 r2_1 r2_2
y1 = polyopvar('y1',[r1;r2],[0,1;0,2]);     y2 = polyopvar('y2',[r1;r2],[0,1;0,2]);
uy = otimes(y1,y2,2);               CopY = uy.C.ops{1,1};
XY = {1+r1+2*r1*r2, r2+r1^2};
WY = 1+t*r1_2+r2_1^2;
FY = int(WY,t,r1_1,r2_2)*apply_op(CopY,XY);
FY = int(int(FY,r1_2,r1_1,1),r1_1,0,1);
FY = int(int(FY,r2_1,0,2),r2_2,0,2);
KCop = int(CopY,WY,[r2_2;t;r1_2;r2_1;r1_1],[0,2;r1_1,r2_2;r1_1,1;0,2;0,1]);
check(apply_val(KCop,XY),double(FY),'two spatial variables, shuffled, other variable between both');
% a copy bounded by a copy of the other spatial variable
assert_err(@() int(CopY,WY,[r1_1;r1_2;r2_1;r2_2;t],[0,1;0,r2_1;0,2;0,2;0,1]),'full domain');
end

function case_var_limits_random()
% Random regions, against the iterated polynomial integral: the copies in a
% random order, the first over [0,2], each other copy below ([0,q]), above
% ([q,2]) or at ([q,q]) a random earlier copy q; an other variable t over
% [c1,c2], [c1,q], [q,c2] or at q, for a random copy q and numbers c1<c2;
% random polynomial weights in the copies and t. For factors without
% integrals (2 and 3 copies), and a basis operator with integral terms
% (degs.in [1;1], 25 rows, 2 copies)
rng(20261119);
pvar s s_1 s_2 s_3 t
cps = [s_1;s_2;s_3];
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
x3 = polyopvar('x3',s,[0,2]);
u = otimes(x1,x2,2);                u3 = otimes(u,x3,2);
Xs = {1+2*s-s^2, s+3*s^3, 2-s+s^2};
Zop = monopmat(1,struct('in',{{'s'}},'out',{{'s'}}),struct('in',[0,2],'out',[0,2]), ...
               struct('in',2,'out',true),struct('in',[1;1],'out',0),struct('otype',2));
Cops = {u.C.ops{1,1}, u3.C.ops{1,1}, fix_empty_params(Zop.ops{1,1})};
ncps = [2,3,2];         ntrials = [15,15,6];
for c=1:numel(Cops)
    Cop = Cops{c};      d = ncps(c);
    X = Xs(1:d);        Y = apply_op(Cop,X);
    p = size(Cop,1);
    for trial=1:ntrials(c)
        % the region of the copies, and of t (row d+1)
        ord = randperm(d);
        dom = cell(d+1,2);      dom(ord(1),:) = {0,2};
        for k=2:d
            q = cps(ord(randi(k-1)));
            lims_k = {0,q; q,2; q,q};
            dom(ord(k),:) = lims_k(randi(3),:);
        end
        q = cps(randi(d));      cc = sort(round(4*rand(1,2)-1,1));
        lims_t = {cc(1),cc(2); cc(1),q; q,cc(2); q,q};
        dom(d+1,:) = lims_t(randi(4),:);
        % a random weight, different for each row of Cop
        coef = randi([-3,3],1,5);
        W = (coef(1)+coef(2)*s_1+coef(3)*s_2^2+coef(4)*t*cps(d)+coef(5)*t^2)*ones(1,p) + (s_1*t)*(1:p)/p;
        KCop = int(Cop,W,[cps(1:d);t],dom);
        % the iterated integral: t, then the copies from the last in ord
        F = W*Y;
        for v=[d+1,flip(ord)]
            if v==d+1,  var = t;    else,   var = cps(v);   end
            [L,U] = dom{v,:};
            if ~isnumeric(L) && isequal(L,U)
                F = subs(F,var,L);
            else
                F = int(F,var,L,U);
            end
        end
        label = sprintf('random region, operator %d, trial %d: dom = {%s}',c,trial, ...
                        strjoin(cellfun(@lim2str,dom','UniformOutput',false),','));
        check(apply_val(KCop,X),double(F),label);
    end
end
end

function str = lim2str(L)
% A limit as text, for the labels of case_var_limits_random
if isnumeric(L),    str = num2str(L);
else,               str = L.varname{1};
end
end

function case_var_limits_errors()
pvar s s_1 s_2 s_3 t
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
x3 = polyopvar('x3',s,[0,2]);
u = otimes(x1,x2,2);                Cop2 = u.C.ops{1,1};
u3 = otimes(u,x3,2);                Cop3 = u3.C.ops{1,1};
Cop1 = x1.C.ops{1,1};               Cop1p = otimes(Cop1,Cop1,[1,2]);
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;0,s_2]),'own limits');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;s_2,s_2]),'own limits');
assert_err(@() int(Cop2,[],[s_1;s_2],[s_2,2;s_1,2]),'circular');
assert_err(@() int(Cop2,[],[s_1;s_2],[s_2,s_2;s_1,s_1]),'circular');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;0,s_3]),'named as in');
assert_err(@() int(Cop1p,[],{'s_1','t'},{0,2;0,'s'}),'named as in');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;0,2*s_1]),'pvar');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;1,s_1]),'full domain');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;1,1]),'full domain');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;2,s_1]),'full domain');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;s_1,0]),'full domain');
assert_err(@() int(Cop2,[],[s_1;s_2;t],[0,2;0,t;0,2]),'full domain');
assert_err(@() int(Cop2,[],[s_1;s_2],[0,2;0,s_1;0,2]),'N x 2');
% s_3 between s_1 and s_2, which are independent: not ordered
assert_err(@() int(Cop3,[],[s_1;s_2;s_3],[0,2;0,2;s_1,s_2]),'not ordered');
% other variables in a cycle
assert_err(@() int(Cop1p,[],{'s','t','r'},{0,2;0,'r';0,'t'}),'circular');
end
% DJ (CC), 11/19/2026: end.

% DJ (CC), 11/22/2026: start. The form int(Cop,W,intmat).
function case_intmat_input()
% The region as an intmat, int(Cop,W,{row}), against the iterated
% polynomial integral, for random regions of D>=d copies: the d copies of
% Cop first, in a random order (the first over [a,b], each other below,
% above or at an earlier copy), then the other copies (also over [a,b]),
% over which W is integrated. For factors without integrals (2 and 3
% copies), a pointwise output (1 copy, as for the pointwise blocks of
% piesos_quadvar), and a basis operator with integral terms (2 copies);
% then the errors of malformed regions
rng(20261122);
pvar s s_1 s_2 s_3 s_4
cps = [s_1;s_2;s_3;s_4];
x1 = polyopvar('x1',s,[0,2]);       x2 = polyopvar('x2',s,[0,2]);
x3 = polyopvar('x3',s,[0,2]);
u = otimes(x1,x2,2);                u3 = otimes(u,x3,2);
Cop2 = u.C.ops{1,1};                Cop3 = u3.C.ops{1,1};
Cop1 = x1.C.ops{1,1};               Cop1p = otimes(Cop1,Cop1,[1,2]);    % x1(s)^2
Xs = {1+2*s-s^2, s+3*s^3, 2-s+s^2};
Zop = monopmat(1,struct('in',{{'s'}},'out',{{'s'}}),struct('in',[0,2],'out',[0,2]), ...
               struct('in',2,'out',true),struct('in',[1;1],'out',0),struct('otype',2));
Cops = {Cop2, Cop3, Cop1p, fix_empty_params(Zop.ops{1,1})};
Xc = {Xs(1:2), Xs, Xs([1,1]), Xs(1:2)};
ds = [2,3,1,2];         ntrials = [12,12,12,4];
for c=1:numel(Cops)
    Cop = Cops{c};      d = ds(c);      X = Xc{c};
    Y = subs(apply_op(Cop,X),s,s_1);    % a pointwise output is in s, copy s_1
    p = size(Cop,1);
    for trial=1:ntrials(c)
        D = min(d+randi([0,2]),4);
        ord = [randperm(d),d+randperm(D-d)];
        row = zeros(1,2*D);     row(2*ord(1)+[-1,0]) = [0,D+1];
        for k=2:D
            q = ord(randi(k-1));
            lims_k = [0,q; q,D+1; q,q; 0,D+1];
            row(2*ord(k)+[-1,0]) = lims_k(randi(3+(ord(k)>d)),:);
        end
        coef = randi([-3,3],1,4);
        W = (coef(1)+coef(2)*s_1+coef(3)*cps(D)^2+coef(4)*s_1*cps(D))*ones(1,p) + cps(D)*(1:p)/p;
        KCop = int(Cop,W,{row});
        % the iterated integral, from the last copy in ord
        F = W*Y;
        lim_vals = [{0};num2cell(cps(1:D));{2}];
        for j=flip(ord)
            L = row(2*j-1);     U = row(2*j);
            if L==U
                F = subs(F,cps(j),cps(L));
            else
                F = int(F,cps(j),lim_vals{L+1},lim_vals{U+1});
            end
        end
        check(apply_val(KCop,X),double(F),sprintf('operator %d, trial %d: region %s',c,trial,mat2str(row)));
    end
end
% errors
assert_err(@() int(Cop2,[],{[0,3,0,1],[0,1]}),'1 x M cell');
assert_err(@() int(Cop2,[],{[0,2]}),'D>=2');
assert_err(@() int(Cop2,[],{[0,3,0,1;0,3,1,3]}),'1 x 2*D');
assert_err(@() int(Cop2,[],{[0,3,0,1.5]}),'integers');
assert_err(@() int(Cop2,[],{[0,4,0,3,0,1]}),'full domain');     % copy 2 of Cop up to copy 3
assert_err(@() int(Cop2,[],{[0,3,0,0]}),'full domain');         % at a
assert_err(@() int(Cop2,[],{[0,3,3,3]}),'full domain');         % at b
assert_err(@() int(Cop2,[],{[0,3,2,3]}),'own limits');
assert_err(@() int(Cop2,[],{[2,3,1,3]}),'circular');
assert_err(@() int(Cop1p,[],{[0,4,3,4,2,4]}),'circular');       % other copies 2, 3
assert_err(@() int(Cop3,[],{[0,4,0,4,1,2]}),'not ordered');
end
% DJ (CC), 11/22/2026: end.


%% ---- helpers ----------------------------------------------------------

function val = eval_at(V,X)
% V at the states X (struct of polynomials by state name), as a number
names = V.varname;
xv = cell(numel(names),1);
for k=1:numel(names)
    xv{k} = X.(names{k});
end
val = double(eval_polyopvar(V,xv));
end

function val = apply_val(KCop,X)                                            % DJ (CC), 11/19/2026
% The functional KCop applied to the factors X (cell of polynomials in the
% input variables), as a number
val = double(apply_op(KCop,X));
end

function T = fix_empty_params(T)                                            % DJ (CC), 11/19/2026
% Zeros for the empty parameters params{k}=[] that monopmat leaves for the
% excluded terms of the 'sopvar' factors of the 'tensopvar' T, as apply_op
% and otimes cannot handle them (as in piesos_quadvar)
if isa(T,'tensopvar')
    for k=1:numel(T.ops)
        T.ops{k} = fix_empty_params(T.ops{k});
    end
elseif isa(T,'sopvar')
    params = T.params;
    nL = prod(cellfun(@numel,T.ZL));    nR = prod(cellfun(@numel,T.ZR));
    for k=1:numel(params)
        if isempty(params{k}),  params{k} = sparse(T.dims(1)*nL,T.dims(2)*nR);  end
    end
    T = sopvar(params,T.vars,T.ZL,T.ZR,T.dom,T.dims);
end
end

function check(val,orc,label)
relerr = max(abs(val(:)-orc(:)))/max(max(abs(orc(:))),1);                   % DJ (CC), 11/19/2026
assert(relerr<1e-8,sprintf('%s: V=%.10g, oracle=%.10g (relative error %.2e)',label,val,orc,relerr));
end

function assert_err(fn,pattern)
try
    fn();
catch err
    assert(contains(err.message,pattern),sprintf('Unexpected error message: "%s" (expected to contain "%s")',err.message,pattern));
    return
end
error('No error thrown (expected "%s").',pattern);
end

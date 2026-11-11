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

% NOTES
% - The test needs eval_polyopvar (Testfolder/polyopvar);
% - The copies are named 'si_j'; the only copy of a pointwise output (degree
%   1, otype(1)=1) is also named 'si', as before. A single-factor tensopvar
%   has otype [2,2] by default, so its copy is named 'si_1' only.

tests = {@case_values_1D, @case_values_2D, @case_names, @case_no_weight_unchanged};
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

function check(val,orc,label)
relerr = abs(val-orc)/max(abs(orc),1);
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

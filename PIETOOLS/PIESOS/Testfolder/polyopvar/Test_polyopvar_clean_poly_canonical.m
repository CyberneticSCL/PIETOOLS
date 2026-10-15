% Test_polyopvar_clean_poly_canonical.m -- verify that
% @polyopvar/clean_poly.m gives the same functional for equal
% functionals, with coefficients of type 'intop': the same monomials,
% intmat, and kernels, and that the difference of equal functionals is
% recognized as zero (is_zero).
%
% Cases: (1) x o x with the full square [0,3,0,3] and kernel k, with
% kernel k(s2,s1), the chain [0,3,0,1] with k(s1,s2)+k(s2,s1), and the
% chain [0,2,0,3] with this kernel; (2) x o y for two states, where the
% square is the sum of the chains [0,3,0,1] and [0,2,0,3] with the same
% kernel, which cannot be relabeled; (3) a quartic functional in x with
% the full cube, and the same functional with the copies of x relabeled;
% (4) the same in two spatial variables, relabeling the copies in both;
% (5) a functional that is not equal is not zero; (6) x o x with a
% multiplier [0,3,1,1], whose kernel k(s1,s2) is the same functional as
% k(s1,s1), and as the same kernel on [2,2,0,3]; (7) a vector-valued
% state, int x'*M*x with M not symmetric: the kernels on the elements
% (1,2) and (2,1) are interchangeable, only their sum matters.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - Equal functionals must have identical clean_poly output because
%   piesos_eq equates the kernels cell by cell; same_poly compares
%   monomials, intmat, kernels.
% - Values are evaluated by eval_functional below, not apply_poly, which
%   does not substitute states in functionals of degree >1.

tol = 1e-8;

%% Case 1: x o x, one spatial variable
pvar s s_dum1 s_dum2
x = polyopvar('x1',s,[0,1]);
k = s_dum1^2+3*s_dum1*s_dum2^3;
ksw = subs(k,[s_dum1;s_dum2],[s_dum2;s_dum1]);
mk = @(K,im) set_coeff(x*x,intop({K},{im},{'s'},[0,1]));
forms = {mk(k,[0,3,0,3]), mk(ksw,[0,3,0,3]), mk(k+ksw,[0,3,0,1]), mk(k+ksw,[0,2,0,3])};
xv = rand_poly([1,1],{'s'},2,4);
vals = zeros(1,4);
for j=1:4
    vals(j) = double(cleanpoly(eval_functional(forms{j},{xv}),1e-12));
end
if max(vals)-min(vals)>tol*max(1,max(abs(vals)))
    error("Test_polyopvar_clean_poly_canonical failed: Case 1, the four forms are not equal functionals.")
end
for j=2:4
    if ~same_poly(clean_poly(forms{1}),clean_poly(forms{j}))
        error("Test_polyopvar_clean_poly_canonical failed: Case 1, form %d gives another clean_poly than form 1.",j)
    end
    [~,is_zero] = clean_poly(forms{1}-forms{j});
    if ~is_zero
        error("Test_polyopvar_clean_poly_canonical failed: Case 1, form 1 - form %d is not zero.",j)
    end
end

%% Case 2: x o y, two states: the square is the sum of two chains with the same kernel
y = polyopvar('y1',s,[0,1]);
fsq = set_coeff(x*y,intop({k},{[0,3,0,3]},{'s'},[0,1]));
fch = set_coeff(x*y,intop({k;k},{[0,3,0,1;0,2,0,3]},{'s'},[0,1]));
fo = set_coeff(x*y,intop({k},{[0,3,0,1]},{'s'},[0,1]));
[~,is_zero] = clean_poly(fsq-fch);
if ~is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 2, the square is not the sum of the chains.")
end
[~,is_zero] = clean_poly(fsq-fo);
if is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 2, the square is not one chain for two states.")
end
if ~same_poly(clean_poly(fsq),clean_poly(fch))
    error("Test_polyopvar_clean_poly_canonical failed: Case 2, other form after clean_poly.")
end

%% Case 3: quartic functional in x: the cube, with the copies of x relabeled
pvar s_dum3 s_dum4
vv = [s_dum1;s_dum2;s_dum3;s_dum4];
rng(1);
k4 = rand_poly([1,1],{'s_dum1','s_dum2','s_dum3','s_dum4'},2,6);
perm = [3,1,4,2];
k4p = subs(k4,vv,vv(perm));
mk4 = @(K) set_coeff(x*x*x*x,intop({K},{[0,5,0,5,0,5,0,5]},{'s'},[0,1]));
f4 = mk4(k4);       g4 = mk4(k4p);
[~,is_zero] = clean_poly(f4-g4);
if ~is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 3, relabeled quartic functional is not equal.")
end
if ~same_poly(clean_poly(f4),clean_poly(g4))
    error("Test_polyopvar_clean_poly_canonical failed: Case 3, relabeled quartic functional gives another clean_poly.")
end
xv = rand_poly([1,1],{'s'},1,3);
err = cleanpoly(eval_functional(f4,{xv})-eval_functional(clean_poly(f4),{xv}),1e-10);
if ~isdouble(err) || abs(double(err))>tol*max(1,abs(double(cleanpoly(eval_functional(f4,{xv}),1e-10))))
    error("Test_polyopvar_clean_poly_canonical failed: Case 3, clean_poly changed the value.")
end

%% Case 4: two spatial variables, copies relabeled in both
pvar s1 s2
x2 = polyopvar('x2',[s1,s2],[0,1]);
dn = {'s1_dum1','s1_dum2','s2_dum1','s2_dum2'};
w = polynomial(dn(:));
kk2 = rand_poly([1,1],dn,2,6);
kk2p = subs(kk2,w,w([2,1,4,3]));
mk2 = @(K,I1,I2) set_coeff(x2*x2,intop({K},{I1,I2},{'s1','s2'},[0,1;0,1]));
f2 = mk2(kk2,[0,3,0,3],[0,3,0,3]);
g2 = mk2(kk2p,[0,3,0,3],[0,3,0,3]);
[~,is_zero] = clean_poly(f2-g2);
if ~is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 4, relabeled functional is not equal.")
end
if ~same_poly(clean_poly(f2),clean_poly(g2))
    error("Test_polyopvar_clean_poly_canonical failed: Case 4, relabeled functional gives another clean_poly.")
end
% [2,3,0,3] along s1 is only one of the two chains of the square: another functional
h2 = mk2(kk2p,[2,3,0,3],[0,3,0,3]);
[~,is_zero] = clean_poly(f2-h2);
if is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 4, different functionals are equal.")
end

%% Case 5: a functional that is not equal is not zero
[~,is_zero] = clean_poly(forms{1}-mk(k,[0,3,0,1]));
if is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 5, a different functional was recognized as zero.")
end

%% Case 6: multiplier: the kernel k(s1,s2) on [0,3,1,1] (s2=s1) is the functional with kernel k(s1,s1)
kd = subs(k,s_dum2,s_dum1);
f6 = mk(k,[0,3,1,1]);
g6 = mk(kd,[0,3,1,1]);
h6 = mk(k,[2,2,0,3]);                   % s1=s2, with s2 the integrated variable
[~,is_zero] = clean_poly(f6-g6);
if ~is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 6, k(s1,s2) and k(s1,s1) on [0,3,1,1] are not equal.")
end
[~,is_zero] = clean_poly(f6-h6);
if ~is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 6, [0,3,1,1] and [2,2,0,3] are not equal.")
end
if ~same_poly(clean_poly(f6),clean_poly(g6)) || ~same_poly(clean_poly(f6),clean_poly(h6))
    error("Test_polyopvar_clean_poly_canonical failed: Case 6, other form after clean_poly.")
end
[~,is_zero] = clean_poly(f6-mk(kd+1,[0,3,1,1]));
if is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 6, a different functional was recognized as zero.")
end

%% Case 7: vector-valued state with two elements, int x'*M*x: the weights M(1,2) and M(2,1) are interchangeable
xv2 = polyopvar('y1',2,s,[0,1]);
M = [s+1, 2*s^2; -3*s, 5];                          % not symmetric
Msym = [s+1, (2*s^2-3*s)/2; (2*s^2-3*s)/2, 5];      % the same functional
Mdif = [s+1, 2*s^2; 2*s^2, 5];                      % another functional
f7 = innerprod(xv2,M,xv2);
g7 = innerprod(xv2,Msym,xv2);
h7 = innerprod(xv2,Mdif,xv2);
[~,is_zero] = clean_poly(f7-g7);
if ~is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 7, asymmetric and symmetric weights of x'*M*x are not equal.")
end
if ~same_poly(clean_poly(f7),clean_poly(g7))
    error("Test_polyopvar_clean_poly_canonical failed: Case 7, other form after clean_poly.")
end
[~,is_zero] = clean_poly(f7-h7);
if is_zero
    error("Test_polyopvar_clean_poly_canonical failed: Case 7, a different functional was recognized as zero.")
end
disp('Test_polyopvar_clean_poly_canonical passed.')


%% ---- local functions ----------------------------------------------------
function V = set_coeff(V,K)
% Set the coefficient of the monomial of V to the functional K (as in piesos_sosvar)
V.C = tensopmat({K});
V.dims.out = V.C.dims.out;      V.vars.out = V.C.vars.out;
V.dom.out = V.C.dom.out;        V.depmat.out = V.C.depmat.out;
V.chkval = true;
end

function Y = eval_functional(f,xvals)
% Value of f = sum_j K_j (x o ... o x) at xvals: apply K_j to the product of the states, with
% in each factor the spatial variable v renamed to v_k, k = number of factors so far depending on v
Y = 0;
svars = f.vars.in;
for j=1:size(f.Z,1)
    cnt = zeros(1,numel(svars));
    xval = 1;
    for st = nonzeros(f.Z(j,:))'
        dep = f.depmat.in(st,:);
        cnt(dep) = cnt(dep)+1;
        newnames = arrayfun(@(i) [svars{i},'_',num2str(cnt(i))],find(dep),'UniformOutput',false);
        xval = kron(xval,subs(xvals{st},polynomial(svars(dep)'),polynomial(newnames(:))));
    end
    Y = Y + apply_op(f.C.ops{1,j},xval);
end
end

function tf = same_poly(A,B)
% Same monomials, and for each block the same intmat and kernels
tf = isequal(full(A.Z),full(B.Z)) && isequal(size(A.C.ops),size(B.C.ops));
if ~tf
    return
end
for j=1:numel(A.C.ops)
    a = A.C.ops{j};     b = B.C.ops{j};
    if isempty(a) && isempty(b)
        continue
    end
    if ~isa(a,'intop') || ~isa(b,'intop') || ~isequal(a.intmat,b.intmat) || ~isequal(size(a.params),size(b.params))
        tf = false;
        return
    end
    ka = kernels(a);    kb = kernels(b);
    for r=1:numel(ka)
        e = cleanpoly(ka{r}-kb{r},1e-10);
        if ~isdouble(e) || any(any(abs(double(e))>1e-9))
            tf = false;
            return
        end
    end
end
end

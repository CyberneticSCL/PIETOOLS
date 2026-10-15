% Test_polyopvar_clean_poly_idempotent.m -- verify that
% @polyopvar/clean_poly.m is idempotent on functionals of repeated
% states, clean_poly(clean_poly(f)) = clean_poly(f), and does not change
% the functional.
%
% For each of the functionals f below, with coefficients of type 'intop':
% (1) clean_poly is applied four times; the results of the second, third
%     and fourth call must be identical to that of the first: the same
%     monomials Z, the same intmat, monomials Z.in, and parameters of
%     every functional (so also the same kernels() and intmat);
% (2) the cleaned functionals must have the same value as f at random
%     polynomial states (apply_op of the coefficients on copies of the
%     states; apply_poly does not substitute states in functionals of
%     degree >1).
% Cases: x o x with the tied full square [0,3,0,3] and with the two
% orderings (chains) of the variables, a quartic functional in x (4
% factors, with ties), x o x in two spatial variables (x(s1,s2)), a
% quartic functional in two spatial variables, and a functional of the
% monomial x o y o x o y of two states.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - The test checks that clean_poly does not swap the dummy variables
%   tied in the order of integration (full square [0,3,0,3]) on every
%   call (@intop/combine_terms.m); piesos_eq and piesos_ineq compare
%   kernels without symmetrizing.

tol = 1e-8;

%% Case 1: x o x in one spatial variable, tied and ordered variables
pvar s s_dum1 s_dum2
x = polyopvar('x1',s,[0,1]);
xv = rand_poly([1,1],{'s'},2,5);
intmat_arr = {[0,3,0,3],[0,3,0,1],[0,2,0,3]};
for j=1:numel(intmat_arr)
    f = set_coeff(x*x,intop({s_dum1*(1+s_dum2^2)},intmat_arr(j),{'s'},[0,1]));
    check_idempotent(f,{xv},tol,sprintf('Case 1, intmat %s',mat2str(intmat_arr{j})))
end

%% Case 2: quartic functional in x, with ties, (x o x o x o x), 3 terms
rng(1);
dumvars = {'s_dum1','s_dum2','s_dum3','s_dum4'};
intmat = [0,5,0,5,0,5,0,5; 0,5,0,5,0,1,0,2; 0,2,0,3,0,4,0,5];
kern = cell(3,1);
for k=1:3
    kern{k} = rand_poly([1,1],dumvars,1,6);
end
f = set_coeff(x*x*x*x,intop(kern,{intmat},{'s'},[0,1]));
xv = rand_poly([1,1],{'s'},1,3);
check_idempotent(f,{xv},tol,'Case 2')

%% Case 3: x o x in two spatial variables, ordered differently along s1 and s2
pvar s1 s2
x2 = polyopvar('x2',[s1,s2],[0,1]);
dumvars = {'s1_dum1','s1_dum2','s2_dum1','s2_dum2'};
xv2 = rand_poly([1,1],{'s1','s2'},2,5);
% intmat along s1, along s2, and the number of kernels
cases = {{[0,2,0,3;0,3,0,1],[0,3,0,1]}, ...                % two orders along s1, twisted along s2
         {[0,2,0,3;0,3,0,1;0,3,0,3],[0,3,0,1]}, ...        % plus a tied term along s1
         {[0,3,0,3],[0,3,0,1;0,2,0,3]}, ...                % tied along s1, orders along s2
         {[0,3,0,3],[0,3,0,3]}};                           % tied along both
for j=1:numel(cases)
    nt = [size(cases{j}{1},1),size(cases{j}{2},1)];
    kern = cell(nt);
    for k=1:numel(kern)
        kern{k} = rand_poly([1,1],dumvars,2,6);
    end
    f = set_coeff(x2*x2,intop(kern,cases{j},{'s1','s2'},[0,1;0,1]));
    check_idempotent(f,{xv2},tol,sprintf('Case 3.%d',j))
end

%% Case 4: quartic functional in two spatial variables
dumvars = {'s1_dum1','s1_dum2','s1_dum3','s1_dum4','s2_dum1','s2_dum2','s2_dum3','s2_dum4'};
I1 = [0,5,0,5,0,5,0,5; 0,5,0,5,0,1,0,2];
I2 = [0,5,0,5,0,1,0,5; 0,2,0,3,0,4,0,5];
kern = cell(2,2);
for k=1:numel(kern)
    kern{k} = rand_poly([1,1],dumvars,1,4);
end
f = set_coeff(x2*x2*x2*x2,intop(kern,{I1,I2},{'s1','s2'},[0,1;0,1]));
xv2 = rand_poly([1,1],{'s1','s2'},1,3);
check_idempotent(f,{xv2},tol,'Case 4')

%% Case 5: x o y o x o y, two states, in one spatial variable
y = polyopvar('y1',s,[0,1]);
dumvars = {'s_dum1','s_dum2','s_dum3','s_dum4'};
intmat = [0,5,0,5,0,5,0,5; 0,2,0,3,0,4,0,5; 0,5,0,5,0,1,0,2];
kern = cell(3,1);
for k=1:3
    kern{k} = rand_poly([1,1],dumvars,1,6);
end
f = x*y*x*y;
f = set_coeff(f,intop(kern,{intmat},{'s'},[0,1]));
xv = rand_poly([1,1],{'s'},1,3);
yv = rand_poly([1,1],{'s'},1,3);
check_idempotent(f,{xv,yv},tol,'Case 5')

disp('Test_polyopvar_clean_poly_idempotent passed.')


%% ---- local functions ----------------------------------------------------
function V = set_coeff(V,K)
% Set the coefficient of the monomial of V to the functional K (as in piesos_sosvar)
V.C = tensopmat({K});
V.dims.out = V.C.dims.out;      V.vars.out = V.C.vars.out;
V.dom.out = V.C.dom.out;        V.depmat.out = V.C.depmat.out;
V.chkval = true;
end

function check_idempotent(f,xvals,tol,label)
% clean_poly applied four times: calls 2-4 leave the result of call 1 unchanged,
% and the functional keeps its value
W = cell(1,4);
W{1} = clean_poly(f);
for k=2:4
    W{k} = clean_poly(W{k-1});
    if ~same_poly(W{k},W{1})
        error("Test_polyopvar_clean_poly_idempotent failed: %s, call %d of clean_poly changed the functional.",label,k)
    end
end
Y0 = eval_functional(f,xvals);
Y1 = eval_functional(W{1},xvals);
err = cleanpoly(Y0-Y1,1e-10);
nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
    error("Test_polyopvar_clean_poly_idempotent failed: %s, clean_poly changed the value of the functional.",label)
end
end

function Y = eval_functional(f,xvals)
% Value of f = sum_j K_j (x o ... o x), K_j = f.C.ops{1,j} of type 'intop', at the states xvals,
% by applying K_j to the product of the states, with in each factor the spatial variable v renamed to v_k,
% k = number of factors so far depending on v (as in Test_intop_combine_terms_ND_random)
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
% Same monomials, and for each block the same intmat, monomials, parameters, and kernels
tf = isequal(full(A.Z),full(B.Z)) && isequal(size(A.C.ops),size(B.C.ops));
if ~tf
    return
end
for j=1:numel(A.C.ops)
    a = A.C.ops{j};     b = B.C.ops{j};
    if isempty(a) && isempty(b)
        continue
    end
    if ~isa(a,'intop') || ~isa(b,'intop') || ~same_op(a,b)
        tf = false;
        return
    end
    ka = kernels(a);    kb = kernels(b);
    for k=1:numel(ka)
        err = cleanpoly(ka{k}-kb{k},1e-12);
        if ~isdouble(err) || any(any(double(err)))
            tf = false;
            return
        end
    end
end
end

function tf = same_op(A,B)
tf = isequal(A.intmat,B.intmat) && isequal(A.Z.in,B.Z.in) && isequal(A.Zd,B.Zd) ...
    && isequal(size(A.params),size(B.params)) && all(cellfun(@isequal,A.params(:),B.params(:)));
end

% Test_polyopvar_subs_nested.m -- substitution of some, or all, of the
% states in a monomial acted on by a nested operator (a factor of the
% 'tensopvar' is itself a 'tensopvar' acting on several of the factors of
% the monomial): the sum of monomials with the same Z and the same number
% of factors but multipliers m_k from R depending on different variables,
% x*m1*y + x*m2*y. subs_factor composes the coefficient operator with an
% operator that has a factor for each leaf factor of the monomial
% (@tensopvar/mtimes, nested A).
%
% Verification is semantic and independent of subs: the polynomial is
% evaluated (eval_polyopvar, test-only, one input for each leaf factor of
% the operator) after the substitution at the remaining states, and
% compared with the product of the polynomials x*(m1+m2)*y with the
% substituted values, as polynomials.
% Parts:
%   (1) substitution of x, y, both, and x then y, by constants
%       (double/polynomial);
%   (2) substitution by position, by a constant;
%   (3) substitution by polynomial functions of new states (a sum of
%       degree-one monomials in the new state, so that the operator of
%       the substituted factor has several terms), of x, y, and x then y;
%   (4) a constant for one state, a function for the other.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;
ncases = 0;
pvar s1 s2
m1 = 1-s1^2;        m2 = 2+s2+s1*s2;        m3 = 3-s2;      m4 = 1+s1+s2^2;
polyopvar x(s1,s2) y(s1,s2)
Id = struct('x',x.C.ops{1,1},'y',y.C.ops{1,1});
% Each row: Z of the monomials, factors of each monomial, coefficients
specs = {   [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [2,-1]
            [0,1,2],    {{m1,'x','y'},{m3,'x','y'}},                  [1,1]
            [1,2,0],    {{'x','y',m1},{'x','y',m3}},                  [1,1]
            [1,0,0,2],  {{'x',m1,m3,'y'},{'x',m3,m1,'y'}},            [1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'},{'x',m4,'y'}},     [1,1,1]   };
rng(1);

%% (1) constants, by name
for k=1:size(specs,1)
    [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
    cx = rand_const([s1,s2]);       cy = rand_const([s1,s2]);
    lbl = sprintf('(1.%d)',k);
    g = subs(f,{'x'},cx);
    check(g,struct('y',rand_const([s1,s2])),Zrow,facs,cf,struct('x',cx),tol,[lbl ' x']);
    g = subs(f,{'y'},cy);
    check(g,struct('x',rand_const([s1,s2])),Zrow,facs,cf,struct('y',cy),tol,[lbl ' y']);
    % both, at once, and one after the other
    g = subs(f,{'x';'y'},[cx;cy]);
    check(g,struct(),Zrow,facs,cf,struct('x',cx,'y',cy),tol,[lbl ' x,y']);
    g = subs(subs(f,{'x'},cx),{'y'},cy);
    check(g,struct(),Zrow,facs,cf,struct('x',cx,'y',cy),tol,[lbl ' x then y']);
    g = subs(subs(f,{'y'},cy),{'x'},cx);
    check(g,struct(),Zrow,facs,cf,struct('x',cx,'y',cy),tol,[lbl ' y then x']);
    ncases = ncases+5;
end

%% (2) by position in Z
for k=1:size(specs,1)
    [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
    cx = rand_const([s1,s2]);       cy = rand_const([s1,s2]);
    lbl = sprintf('(2.%d)',k);
    g = subs(f,find(Zrow==1),cx);
    check(g,struct('y',rand_const([s1,s2])),Zrow,facs,cf,struct('x',cx),tol,[lbl ' x']);
    g = subs(f,find(Zrow==2),cy);
    check(g,struct('x',rand_const([s1,s2])),Zrow,facs,cf,struct('y',cy),tol,[lbl ' y']);
    ncases = ncases+2;
end

%% (3) polynomial functions of new states: x = val1(w1) (a sum of dval operators acting on w1), y = val2(w2)
for k=1:size(specs,1)
    for dval=1:2
        rng(100*k+dval);
        [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
        valx = build_val('w1',dval);        valy = build_val('w2',dval);
        w1 = rand_const([s1,s2]);           w2 = rand_const([s1,s2]);
        vx = eval_polyopvar(valx,{w1});     vy = eval_polyopvar(valy,{w2});
        lbl = sprintf('(3.%d,%d)',k,dval);
        g = subs(f,{'x'},valx);
        check(g,struct('w1',w1,'y',rand_const([s1,s2])),Zrow,facs,cf,struct('x',vx),tol,[lbl ' x']);
        g = subs(f,{'y'},valy);
        check(g,struct('w2',w2,'x',rand_const([s1,s2])),Zrow,facs,cf,struct('y',vy),tol,[lbl ' y']);
        g = subs(subs(f,{'x'},valx),{'y'},valy);
        check(g,struct('w1',w1,'w2',w2),Zrow,facs,cf,struct('x',vx,'y',vy),tol,[lbl ' x then y']);
        % by position
        g = subs(f,find(Zrow==1),valx);
        check(g,struct('w1',w1,'y',rand_const([s1,s2])),Zrow,facs,cf,struct('x',vx),tol,[lbl ' x by position']);
        ncases = ncases+4;
    end
end

%% (4) a constant for one state, a function for the other
for k=1:size(specs,1)
    rng(300+k);
    [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
    valy = build_val('w2',2);
    cx = rand_const([s1,s2]);           w2 = rand_const([s1,s2]);
    vy = eval_polyopvar(valy,{w2});
    lbl = sprintf('(4.%d)',k);
    g = subs(subs(f,{'x'},cx),{'y'},valy);
    check(g,struct('w2',w2),Zrow,facs,cf,struct('x',cx,'y',vy),tol,[lbl ' x const, y function']);
    valx = build_val('w1',2);           cy = rand_const([s1,s2]);
    w1 = rand_const([s1,s2]);           vx = eval_polyopvar(valx,{w1});
    g = subs(subs(f,{'y'},cy),{'x'},valx);
    check(g,struct('w1',w1),Zrow,facs,cf,struct('x',vx,'y',cy),tol,[lbl ' y const, x function']);
    ncases = ncases+2;
end

fprintf('Test_polyopvar_subs_nested passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [f,facs,cf,Zrow] = build_nested(spec,Id)
% Sum of the monomials, with Z=spec{1} and factors spec{2}{l}, coefficients spec{3}, as a nested operator
Zrow = spec{1};     facs = spec{2};     cf = spec{3};
f = cf(1)*build_f(Zrow,facs{1},Id);
for l=2:numel(facs)
    f = f + cf(l)*build_f(Zrow,facs{l},Id);
end
f = combine_terms(f);
Top = f.C.ops{1,1};
if ~isa(Top,'tensopvar') || ~any(cellfun('isclass',Top.ops(:),'tensopvar'))
    error("Test_polyopvar_subs_nested failed: setup: the operator should be nested.")
end
end

function check(g,vals_g,Zrow,facs,cf,subvals,tol,label) %#ok
% Evaluate g at the values of its remaining states (vals_g: struct by name; substituted-out states that remain
% declared get the value 0, and must not appear in Z) and compare with the product of the factors, with x and y
% replaced by the substituted values subvals.x, subvals.y, if given, else by their value in vals_g
Yref = 0;
for l=1:numel(facs)
    Yref = Yref + cf(l)*prod_facs(facs{l},subvals,vals_g);
end
xv = cell(numel(g.varname),1);
for kk=1:numel(xv)
    nm = g.varname{kk};
    if isfield(vals_g,nm)
        xv{kk} = vals_g.(nm);
    else
        if any(g.Z(:)==kk)
            error("Test_polyopvar_subs_nested failed (%s): the substituted state '%s' still appears in Z.",label,nm)
        end
        xv{kk} = polynomial(0);
    end
end
Yg = eval_polyopvar(g,xv);
assert_close(Yg,Yref,tol,label);
end

function Y = prod_facs(facs,subvals,vals_g)
% Product of the values of the factors, as polynomials: x -> its value, y -> its value, polynomial -> itself
Y = 1;
for j=1:numel(facs)
    if ischar(facs{j})
        if isfield(subvals,facs{j}) && ~isempty(subvals.(facs{j}))
            Y = Y*subvals.(facs{j});
        else
            Y = Y*vals_g.(facs{j});
        end
    else
        Y = Y*facs{j};
    end
end
end

function assert_close(Y,Yref,tol,label)
% Output variables of the nested operator are named s_1,...: renamed back
Y = polynomial(Y);      Yref = polynomial(Yref);
Y.varname = regexprep(Y.varname,'_1$','');
D = Y-Yref;
if ~isequal(size(D),size(Yref))
    error("Test_polyopvar_subs_nested failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol*max(1,full(max(abs(nonzeros(Yref.coefficient)))))
    error("Test_polyopvar_subs_nested failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient)) || max(abs(nonzeros(Yref.coefficient)))<1e-6
    error("Test_polyopvar_subs_nested %s: reference value is zero; test would be vacuous.",label)
end
end

function p = rand_const(svs)
% Random polynomial in the spatial variables svs (cell or array of pvars), degree 2 in each
p = polynomial(randn);
for i=1:numel(svs)
    si = svs(i);
    if iscell(si),  si = si{1};     end
    p = p + randn*si + randn*si^2;
end
end

function val = build_val(name,dval)
% polyopvar val(w)(s) = sum_k B_k w(s), dval random degree-one monomials in the new state w, in s1,s2
svars = {'s1','s2'};
dom = struct('in',[0,1;0,1],'out',[0,1;0,1]);
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,2),'out',true(1,2));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {name};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,2);  val.depmat.out = ones(1,2);
val = initialize(val);
end

function f = build_f(Zrow,facs,Id)
% polyopvar with monomial Z=Zrow in the states x,y, acted on by the product of the factors facs: the
% identity operators of the states for 'x','y', multipliers from R for polynomials
K = [];
for j=1:numel(facs)
    if ischar(facs{j})
        Fj = Id.(facs{j});
    else
        vars = cell(1,0);
        if ismember('s1',facs{j}.varname),  vars = [vars,{'s1'}];   end
        if ismember('s2',facs{j}.varname),  vars = [vars,{'s2'}];   end
        Fj = mat2tensopvar(facs{j},struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vars),1)));
    end
    if isempty(K),  K = Fj;
    else,           K = otimes(K,Fj,[1,2]);     end
end
f = polyopvar();
f.varname = {'x';'y'};      f.Z = sparse(Zrow);     f.C = tensopmat(K);
f.dims.in = [1;1];
f.vars.in = {'s1','s2'};    f.dom.in = [0,1;0,1];   f.depmat.in = true(2,2);
f = initialize(f);
end

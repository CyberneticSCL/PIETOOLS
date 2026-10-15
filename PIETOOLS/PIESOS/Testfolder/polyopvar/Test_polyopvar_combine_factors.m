% Test_polyopvar_combine_factors.m -- @polyopvar/combine_factors.m,
% combining consecutive factors that are multipliers from R (consecutive
% zeros in Z) in tensor-PI operators with pointwise output into a single
% factor: Z=[1,0,0,2] becomes Z=[1,0,2], the multiplier the elementwise
% product of those of the factors, on the union of their variables.
% Parts:
%   (1) explicit multipliers on different variables: Z, the multiplier
%       of the new factor, its variables, and the value, compared with
%       the product of the polynomials directly;
%   (2) random operators (several terms, vector-valued states,
%       multipliers on random subsets of the variables): the value
%       before and after, and Z;
%   (3) several monomials at once, and monomials left alone;
%   (4) idempotence, and operations on the result;
%   (5) operators for which factors cannot be combined (Kronecker
%       output): left as they are, also if other monomials, or other
%       operators acting on the same monomial, can be;
%   (6) monomials of which all factors act on the constant state:
%       converted to a 'polynomial' block;
%   (7) clean_poly calls combine_factors: monomials that coincide
%       afterwards are merged.
% Oracle: the product of the multiplier polynomials and state values,
% formed directly (prod_facs), compared with the value by apply_poly.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - Random cases redraw (up to 20 attempts) when the random operator is
%   zero, and assert_close rejects a zero reference: a zero value would
%   make the comparison vacuous.
% - In (7), the merged factors must depend on both variables (m4 does): a
%   sum of operators whose factors depend on different variables is
%   nested, which apply_poly cannot evaluate.

tol = 1e-8;
pvar s1 s2
ncases = 0;

%% (1) Explicit multipliers
m1 = 1-s1^2;        m2 = 2+s2+s1*s2;        m3 = 3-s2;      m4 = 1+s1+s2^2;
polyopvar x(s1,s2) y(s1,s2)
Id = struct('x',x.C.ops{1,1},'y',y.C.ops{1,1});
mults = {m1,m2,m3,m4};
vals = struct('x',1+2*s1+s2,'y',3-s1*s2+s2^2);
specs = {   [1,0,0,2],          {'x',m1,m2,'y'},            [1,0,2]
            [1,0,0,0,2],        {'x',m1,m2,m3,'y'},         [1,0,2]
            [0,0,1],            {m1,m2,'x'},                [0,1]
            [1,0,0],            {'x',m2,m3},                [1,0]
            [0,0,1,0,0,2,0,0],  {m1,m2,'x',m3,m4,'y',m1,m3}, [0,1,0,2,0]
            [1,0,2],            {'x',m2,'y'},               [1,0,2]
            [0,0],              {m1,m4},                    [0]
            [1,2],              {'x','y'},                  [1,2]      };
for k=1:size(specs,1)
    Zrow = specs{k,1};      facs = specs{k,2};      Zexp = specs{k,3};
    f = build_f(Zrow,facs,Id);
    g = combine_factors(f);
    Zg = full(g.Z);     Zg(end+1:numel(Zexp)) = 0;
    if ~isequal(Zg(1:numel(Zexp)),Zexp) || any(Zg(numel(Zexp)+1:end)) || size(g.Z,1)~=1
        error("Test_polyopvar_combine_factors failed (1.%d): Z = [%s], expected [%s].",k,num2str(full(g.Z)),num2str(Zexp))
    end
    if num_factors(g)~=numel(Zexp)
        error("Test_polyopvar_combine_factors failed (1.%d): %d factors, expected %d.",k,num_factors(g),numel(Zexp))
    end
    assert_close(eval_f(g,vals),prod_facs(Zrow,facs,vals),tol,sprintf('(1.%d)',k));
    ncases = ncases+1;
end
% The multiplier of the new factor depends on the variables of all factors in the run
f = build_f([1,0,0,2],{'x',m1,m3,'y'},Id);          % m1(s1), m3(s2)
g = combine_factors(f);
op = g.C.ops{1,1};
fctr = op.ops{1,2};
if ~isequal(fctr.vars.out,{'s1','s2'}) || ~isequal(op.depmat.out(2,:),[1,1]) || ~isequal(size(op.ops,2),3)
    error("Test_polyopvar_combine_factors failed (1): the new factor does not depend on the variables of both factors.")
end
f = build_f([1,0,0],{'x',m1,m3},Id);
g = combine_factors(f);
assert_close(apply_sopvar(g.C.ops{1,1}.ops{1,2},1),m1*m3,tol,'(1) multiplier of the new factor');
ncases = ncases+2;

%% (2) Random operators
svar_all = {s1,s2};
zrows = {[1,0,0,2],[0,0,1],[2,0,0],[1,0,0,0,2],[0,0,1,0,0,2,0,0],[1,0,2],[0,0]};
for N = 1:2
    svs = svar_all(1:N);
    for nterms = 1:2
        for iz = 1:numel(zrows)
            sz = [2;1];
            for attempt=1:20        % random operators may be zero: draw again
                rng(100*N+10*nterms+iz+1000*attempt);
                f = make_rand_f(N,zrows{iz},sz,nterms);
                xv = {rand_pol(svs,2);rand_pol(svs,1)};
                Yref = apply_poly(f,xv);
                if is_nonzero(Yref),    break,      end
            end
            g = combine_factors(f);
            assert_close(apply_poly(g,xv),Yref,tol,sprintf('(2) N=%d,nterms=%d,Z=[%s]',N,nterms,num2str(zrows{iz})));
            Zexp = compress_zeros(zrows{iz});
            Zg = full(g.Z);     Zg(end+1:numel(Zexp)) = 0;
            if ~isequal(Zg(1:numel(Zexp)),Zexp) || any(Zg(numel(Zexp)+1:end))
                error("Test_polyopvar_combine_factors failed (2): Z=[%s] became [%s], expected [%s].",num2str(zrows{iz}),num2str(full(g.Z)),num2str(Zexp))
            elseif any(zrows{iz}) && size(g.C.ops{1,1}.ops,1)~=nterms
                error("Test_polyopvar_combine_factors failed (2): the number of terms changed.")
            elseif ~any(zrows{iz}) && ~isa(g.C.ops{1,1},'polynomial')
                % all factors act on the constant state: a block of type polynomial
                error("Test_polyopvar_combine_factors failed (2): Z=[%s] not converted to a polynomial.",num2str(zrows{iz}))
            end
            ncases = ncases+1;
        end
    end
end

%% (3) Several monomials at once, with monomials that are left alone
N = 1;      svs = svar_all(1);      sz = [2;1];
for attempt=1:20
    rng(300+attempt);
    fs = make_rand_f(N,[1,0,0,2],sz,1) + make_rand_f(N,[0,0,1],sz,1) + make_rand_f(N,[1,2],sz,1) + make_rand_f(N,[2,0,0],sz,1);
    xv = {rand_pol(svs,2);rand_pol(svs,1)};
    Yref = apply_poly(fs,xv);
    if is_nonzero(Yref),    break,      end
end
gs = combine_factors(fs);
assert_close(apply_poly(gs,xv),Yref,tol,'(3) several monomials');
[~,nops_f] = num_factors(fs);       [~,nops_g] = num_factors(gs);
if size(gs.Z,1)~=size(fs.Z,1) || sum(nops_g)~=sum(nops_f)-3
    error("Test_polyopvar_combine_factors failed (3): number of factors %d, expected %d.",sum(nops_g),sum(nops_f)-3)
end
% Nothing to combine: returned as is
g0 = combine_factors(gs);
if ~isequal(g0.Z,gs.Z)
    error("Test_polyopvar_combine_factors failed (3): not idempotent.")
end
ncases = ncases+2;

%% (4) Operations on the result
f = build_f([1,0,0,2],{'x',m1,m2,'y'},Id);
fx = build_f(1,{'x'},Id);
g = combine_factors(f);
assert_close(eval_f(g.*fx,vals),prod_facs([1,0,0,2],{'x',m1,m2,'y'},vals)*vals.x,tol,'(4) times');
assert_close(eval_f(g+f,vals),2*prod_facs([1,0,0,2],{'x',m1,m2,'y'},vals),tol,'(4) plus, original and combined');
assert_close(eval_f(subs(g,{'x'},known_fun(1+s2,{'s1','s2'})),vals),(1+s2)*m1*m2*vals.y,tol,'(4) subs');
ncases = ncases+3;

%% (5) Operators for which the factors cannot be combined are left as they are
% Kronecker output (otype(1)=2): multipliers cannot be multiplied; f is returned unchanged, no error
sz = [2;1];
zrows5 = {[1,0,0,2],[0,0,1],[0,0,1,0,0,2],[1,0,2]};
for iz = 1:numel(zrows5)
    for attempt=1:20
        rng(400+10*iz+attempt);
        f = make_rand_f(1,zrows5{iz},sz,2,[2,2],[2;ones(numel(zrows5{iz})-1,1)]);   % small output
        xv = {rand_pol(svar_all(1),2);rand_pol(svar_all(1),1)};
        Yref = apply_poly(f,xv);
        if is_nonzero(Yref),    break,      end
    end
    g = combine_factors(f);
    assert_unchanged(f,g,sprintf('(5) Kronecker, Z=[%s]',num2str(zrows5{iz})));
    assert_close(apply_poly(g,xv),Yref,tol,sprintf('(5) Kronecker, Z=[%s]',num2str(zrows5{iz})));
    ncases = ncases+1;
end

% Several monomials: only the monomial with the Kronecker output is left, the other one is reduced
for attempt=1:20
    rng(500+attempt);
    fpw = make_rand_f(1,[1,0,0,2],sz,1);
    fkr = make_rand_f(1,[2,0,0,1],sz,1,[2,2],[2;1;1;1]);
    fm = fpw;                   fm.Z = sparse([1,0,0,2;2,0,0,1]);
    Cm = fpw.C;                 Cm.ops = [fpw.C.ops,fkr.C.ops];
    Cm.dims.in = [fpw.C.dims.in;fkr.C.dims.in];     Cm.depmat.in = [fpw.C.depmat.in;fkr.C.depmat.in];
    Cm.chkval = false;
    fm.C = Cm;                  fm.depmat.out = [];     fm.chkval = false;
    fm = initialize(fm);        % output function space of the Kronecker term differs from that of the pointwise one
    xv = {rand_pol(svar_all(1),2);rand_pol(svar_all(1),1)};
    Yref = apply_poly(fm,xv);
    if is_nonzero(Yref),    break,      end
end
gm = combine_factors(fm);
[~,nops_g] = num_factors(gm);       Zg = full(gm.Z);
if ~isequal(nops_g(:)',[3,4]) || ~isequal(Zg(1,1:3),[1,0,2]) || ~isequal(Zg(2,:),[2,0,0,1]) || any(Zg(1,4:end))
    error("Test_polyopvar_combine_factors failed (5): mixed monomials: Z = %s, nops = %s.",mat2str(Zg),mat2str(nops_g(:)'))
end
assert_close(apply_poly(gm,xv),Yref,tol,'(5) mixed monomials');
ncases = ncases+2;

% Several operators acting on the same monomial (rows of C), one of which cannot be combined: the
% monomial is left as it is for all of them, as they share Z
rng(600);
fkr1 = make_rand_f(1,[1,0,0,2],sz,1,[2,2],[2;1;1;1]);      % same Z as fpw, output of size 2
for order = 1:2
    fs = fpw;
    if order==1,    fs.C = [fpw.C;fkr1.C];
    else,           fs.C = [fkr1.C;fpw.C];      end
    fs.dims.out = [2;2];
    fs.depmat.out = [];     fs.chkval = false;
    fs = initialize(fs);
    xv = {rand_pol(svar_all(1),2);rand_pol(svar_all(1),1)};
    Yref = apply_poly(fs,xv);
    gs = combine_factors(fs);
    assert_unchanged(fs,gs,sprintf('(5) stacked operators, order %d',order));
    assert_close(apply_poly(gs,xv),Yref,tol,sprintf('(5) stacked operators, order %d',order));
    ncases = ncases+1;
end

%% (6) Monomials of which all factors act on the constant state: a block of type 'polynomial'
% (tensopvar2polynomial): the value, the number of factors, and the type of the block
facs6 = {   {m1,m4},        [0,0]
            {m1,m2,m3},     [0,0,0]
            {m4,m3},        [0,0]       };
for k=1:size(facs6,1)
    facs = facs6{k,1};      Zrow = facs6{k,2};
    f = build_f(Zrow,facs,Id);
    g = combine_factors(f);
    if ~isa(g.C.ops{1,1},'polynomial') || ~isequal(size(g.C.ops{1,1}),[1,1]) || ~isequal(full(g.Z),0) || num_factors(g)~=1
        error("Test_polyopvar_combine_factors failed (6.%d): not a polynomial block with a single factor.",k)
    end
    if ~isvalid(g.C) || ~g.chkval
        error("Test_polyopvar_combine_factors failed (6.%d): not a valid object.",k)
    end
    assert_close(eval_f(g,vals),prod_facs(Zrow,facs,vals),tol,sprintf('(6.%d)',k));
    ncases = ncases+1;
end
% a Kronecker product of the factors in different variables: a polynomial too
mk = @(m,v) mat2tensopvar(m,struct('in',{cell(1,0)},'out',{{v}}),struct('in',zeros(0,2),'out',[0,1]));
Mv1 = [1-s1^2;2];       Mv2 = [3-s2;s2];
fkr = build_kron_f([0,0],mk(Mv1,'s1'),mk(Mv2,'s2'));
gkr = combine_factors(fkr);
if ~isa(gkr.C.ops{1,1},'polynomial') || ~isequal(full(gkr.Z),0) || ~isvalid(gkr.C)
    error("Test_polyopvar_combine_factors failed (6.4): Kronecker product in different variables not converted.")
end
assert_close(gkr.C.ops{1,1},kron(Mv1,Mv2),tol,'(6.4) matrix of the Kronecker product');
assert_close(apply_poly(gkr,{[];[]}),kron(Mv1,Mv2),tol,'(6.4) value');
ncases = ncases+2;
% ... and in the same variable: not a polynomial (a function of two copies of s1), left as it is
fsh = build_kron_f([0,0],mk(Mv1,'s1'),mk([2-s1;s1],'s1'));
gsh = combine_factors(fsh);
assert_unchanged(fsh,gsh,'(6.5) Kronecker product in the same variable');
if isa(gsh.C.ops{1,1},'polynomial')
    error("Test_polyopvar_combine_factors failed (6.5): converted to a polynomial.")
end
ncases = ncases+1;
% several operators acting on the monomial (rows of C): converted for all, or for none (as Z is shared)
fpw1 = build_f([0,0],{m1,m4},Id);
fkr1 = build_kron_f([0,0],mk(1-s1^2,'s1'),mk(2-s1^2,'s1'));        % scalar, the same variable: not convertible
for order = 1:2
    fs = fpw1;
    if order==1,    fs.C = [fpw1.C;fkr1.C];
    else,           fs.C = [fkr1.C;fpw1.C];     end
    fs.dims.out = [1;1];    fs.depmat.out = [];     fs.chkval = false;
    fs = initialize(fs);
    gs = combine_factors(fs);
    assert_unchanged(fs,gs,sprintf('(6.6) stacked operators, order %d',order));
    ncases = ncases+1;
end
% several operators, all convertible: all are converted
fpw2 = build_f([0,0],{m2,m3},Id);
fs = fpw1;      fs.C = [fpw1.C;fpw2.C];     fs.dims.out = [1;1];    fs.depmat.out = [];     fs.chkval = false;
fs = initialize(fs);
gs = combine_factors(fs);
if ~isa(gs.C.ops{1,1},'polynomial') || ~isa(gs.C.ops{2,1},'polynomial') || ~isvalid(gs.C)
    error("Test_polyopvar_combine_factors failed (6.7): not all operators converted.")
end
assert_close(gs.C.ops{2,1},m2*m3,tol,'(6.7) second row');
ncases = ncases+1;
% with another monomial that has a state, and idempotence
fm6 = build_f([0,0],{m1,m4},Id) + build_f([1,0,0,2],{'x',m2,m3,'y'},Id);
gm6 = combine_factors(fm6);
is_k = ~any(gm6.Z,2);
if nnz(is_k)~=1 || ~isa(gm6.C.ops{1,find(is_k)},'polynomial') || isa(gm6.C.ops{1,find(~is_k)},'polynomial')
    error("Test_polyopvar_combine_factors failed (6.8): constant monomial not converted, or the other one is.")
end
assert_close(eval_f(gm6,vals),prod_facs([0,0],{m1,m4},vals)+prod_facs([1,0,0,2],{'x',m2,m3,'y'},vals),tol,'(6.8)');
gg6 = combine_factors(gm6);
if ~isequal(gg6.Z,gm6.Z) || ~isequal(class(gg6.C.ops{1,find(is_k)}),'polynomial')
    error("Test_polyopvar_combine_factors failed (6.8): not idempotent.")
end
ncases = ncases+2;

%% (7) clean_poly calls combine_factors: monomials that coincide once the factors are combined are merged
fa = build_f([1,0,0,2],{'x',m1,m2,'y'},Id);        fb = build_f([1,0,2],{'x',m4,'y'},Id);
fab = fa+fb;
if size(fab.Z,1)~=2
    error("Test_polyopvar_combine_factors failed (7): setup: monomials with different numbers of factors must differ.")
end
gab = clean_poly(fab);
if size(gab.Z,1)~=1 || ~isequal(full(gab.Z(1,1:3)),[1,0,2]) || num_factors(gab)~=3
    error("Test_polyopvar_combine_factors failed (7): monomials not merged: Z = %s.",mat2str(full(gab.Z)))
end
assert_close(eval_f(gab,vals),prod_facs([1,0,0,2],{'x',m1,m2,'y'},vals)+prod_facs([1,0,2],{'x',m4,'y'},vals),tol,'(7) value');
% a constant monomial of several factors is converted by clean_poly, and its block cleaned
gk = clean_poly(build_f([0,0],{m1,m4},Id));
if ~isa(gk.C.ops{1,1},'polynomial')
    error("Test_polyopvar_combine_factors failed (7): constant monomial not converted by clean_poly.")
end
assert_close(eval_f(gk,vals),m1*m4,tol,'(7) constant');
gz = clean_poly(build_f([0,0],{1e-14*m1,m4},Id));       % the block is zero up to the tolerance: removed
if ~is_zero_poly(gz)
    error("Test_polyopvar_combine_factors failed (7): a polynomial block below the tolerance was not removed.")
end
ncases = ncases+4;

fprintf('Test_polyopvar_combine_factors passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function f = build_kron_f(Zrow,F1,F2)
% f = Kop*Z with Kop the Kronecker product (otype(1)=2) of the 'tensopvar' objects F1, F2, multipliers from R
K = otimes(F1,F2,[2,2]);
f = polyopvar();
f.varname = {'x';'y'};      f.Z = sparse(Zrow);     f.C = tensopmat(K);
f.dims.in = [1;1];
f.vars.in = {'s1','s2'};    f.dom.in = [0,1;0,1];   f.depmat.in = true(2,2);
f = initialize(f);
end

function tf = is_zero_poly(f)
% true if f has no monomials left that contribute, or a zero block
tf = size(f.Z,1)<=1 && (isempty(f.C.ops) || all(cellfun(@(c) isempty(c) || ~any(nonzeros(polynomial(c).coefficient)),f.C.ops(:))));
end

function f = build_f(Zrow,facs,Id)
% f = Kop*Z with Kop the pointwise product of the factors facs{j}: the identity on state 'x' or
% 'y' (names given), or a multiplier m(s) (polynomial given) from R, in the variables of m
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
% Z indices 1,2 refer to the states 'x','y'
f = polyopvar();
f.varname = {'x';'y'};      f.Z = sparse(Zrow);     f.C = tensopmat(K);
f.dims.in = [1;1];
f.vars.in = {'s1','s2'};    f.dom.in = [0,1;0,1];   f.depmat.in = true(2,2);
f = initialize(f);
end

function Y = prod_facs(Zrow,facs,vals)
names = {'x','y'};
Y = polynomial(1);
for j=1:numel(facs)
    if ischar(facs{j}),     Y = Y*vals.(facs{j});
    else,                   Y = Y*facs{j};      end
end
end

function Y = eval_f(f,vals)
xv = cell(numel(f.varname),1);
for k=1:numel(xv)
    xv{k} = vals.(f.varname{k});
end
Y = apply_poly(f,xv);
end

function f = make_rand_f(N,Zrow,sz,nterms,ot,dout)
% f = Kop*Z in states a (size sz(1)) and b (size sz(2)), depending on all N variables, with Kop random
% with nterms terms and pointwise (default) output of size 2; the multipliers from R at the zeros of
% Zrow depend on random subsets of the variables. For Kronecker output (ot(1)=2), the output sizes of
% the factors are dout (default 1+(1:d))
if nargin<5,    ot = [1,2];     end
if nargin<6,    dout = [];      end
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
d = numel(Zrow);
is_R = Zrow==0;
din = ones(d,1);        din(~is_R) = sz(Zrow(~is_R));
depin = true(d,N);      depin(is_R,:) = false;
depout = true(d,N);     depout(is_R,:) = rand(nnz(is_R),N)>0.4;
if ot(1)==1
    dims = [2*ones(d,1),din];       mout = 2;
else
    if isempty(dout),   dout = 1+(1:d)';    end
    dims = [dout(:),din];           mout = prod(dims(:,1));
end
Kop = rand_tensopvar(dims,vars,dom,struct('in',depin,'out',depout),nterms,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = Zrow;                 f.varname = {'a';'b'};
f.dims.in = sz(:);          f.dims.out = mout;
f.vars.in = svars;          f.dom.in = dom.in;          f.depmat.in = true(2,N);
f.C = tensopmat(Kop);
f.vars.out = svars;         f.dom.out = dom.out;        f.depmat.out = Kop.deps.out;
f = initialize(f);
end

function z = compress_zeros(z)
% Replace each run of consecutive zeros by a single zero
keep = true(size(z));
keep(2:end) = ~(z(2:end)==0 & z(1:end-1)==0);
z = z(keep);
end

function p = rand_pol(svs,m)
p = polynomial(randn(m,1));
for i=1:numel(svs)
    p = p + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
end
end

function g = known_fun(m,vars)
M = mat2tensopvar(m,struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vars),1)));
g = polyopvar();
g.varname = cell(0,1);      g.Z = sparse(1,0);      g.C = tensopmat(M);
g = initialize(g);
end

function tf = is_nonzero(Y)
Yc = polynomial(Y).coeff;
tf = ~isempty(Yc) && max(abs(full(Yc(:))))>=1e-6;
end

function assert_close(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_combine_factors (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("Test_polyopvar_combine_factors failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function assert_unchanged(f,g,label)
% g must have the same monomials, numbers of factors and terms as f
[~,nops_f] = num_factors(f);        [~,nops_g] = num_factors(g);
if ~isequal(g.Z,f.Z) || ~isequal(nops_f,nops_g) || ~isequal(size(g.C.ops),size(f.C.ops))
    error("Test_polyopvar_combine_factors failed (%s): the polynomial changed.",label)
end
for k=1:numel(f.C.ops)
    if ~isequal(size(g.C.ops{k}.ops),size(f.C.ops{k}.ops))
        error("Test_polyopvar_combine_factors failed (%s): the operator changed.",label)
    end
end
end

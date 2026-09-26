%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR_REGISTRY checks 'plus' and 'mtimes' of 'copvar' and 'cdopvar'
% on operands built over DIFFERENT variable registries. Until 09/26/2026
% both refused them (plus/mtimes:registryMismatch), so an R^n -> R^n
% operator - empty registry: -gam*Iw, Dzu, Dyw - could not meet an L_2 one,
% which blocked Hinf_gain with Tw ~= 0, Dzu*Z, D12*Z and Z*Dyw in the
% executives. They now restate both operands over the sorted union
% registry ('merge_copvar_registry'), blocks untouched.
%
% Every check is against SEMANTICS, never a round trip (CLAUDE.md S4):
%
%   plus     kernels(A+B) = kernels(A) + kernels(B), blockwise, against the
%            operands' ORIGINAL blocks, via 'pi_blk_kernels';
%   mtimes   (A*B)x = A(Bx) for random polynomial test functions x, both
%            sides applied from block kernels by the kernel form in the
%            'sopvar.m' header (apply_blk below). The right side uses the
%            ORIGINAL blocks of A and B, so no container code enters it.
%            The oracle is calibrated on same-registry products, which
%            'mtimes' already handled, and shown to reject a wrong result.
%
% Decision blocks are evaluated at random values drawn once per NAME.
% Metadata is checked by variable NAME: the registry is the sorted union,
% each name keeps its domain, and every space of the result is the same set
% of names, with the same component count, as in the operand it comes from.
% 'verify' runs on every result.
%
% Registry pairs, each over copvar/cdopvar combinations (a product of two
% cdopvar is refused by design):
%   (a) R^n only (empty registry) with L2[s]: A + X'*Y (the -gam*Iw + PB'*Tw
%       pattern) and Dzu*Z, Z'*Dzu';
%   (b) {s} with {t};
%   (c) 1-D {s1} with 2-D {s1,s2};
% plus: a variable given two domains raises cls:domConflict from all four
% methods, as concatenation does; spaces that differ as name sets still
% raise spaceMismatch on the union; a fixed x decision product leaves every
% block on the container's decision variable list, including decision
% blocks on an EMPTY list.
%
% Initial coding MMP, 09/26/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20260926);
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
% Distinct domains per variable, so a dom row landing on the wrong name shows.
DOM = struct('s',[0,1],'t',[-1,2],'s1',[0,1],'s2',[-1,1]);
E0 = {{}};                  % one R^n space
Ra = {{},{}};               % two R^n spaces
nchk = 0;

% ---------------------------------------------------------- (a) R^n, L2[s]
% Sum: A maps [R^1,R^2] -> [R^1,R^2] over an EMPTY registry; X'*Y maps the
% same spaces but keeps X's registry {s}, as PB'*Tw does in Hinf_gain.
A  = gen(DOM,Ra,Ra,[1;2],[1;2],0,[true false; true true]);
Ad = gen(DOM,Ra,Ra,[1;2],[1;2],3,[true false; true true]);
X  = gen(DOM,{{},{'s'}},Ra,[1;2],[1;2],0);
Y  = gen(DOM,{{},{'s'}},Ra,[1;2],[1;2],0);
Yd = gen(DOM,{{},{'s'}},Ra,[1;2],[1;2],4);
B  = X'*Y;          Bd = X'*Yd;
assert(isempty(A.vars) && isequal(B.vars,{'s'}) && isequal(Bd.vars,{'s'}),...
    'setup (a): expected registries {} and {s}');
% Product: Dzu is R^2 <- R^1 over an empty registry, Z is R^1 <- [R^1, L2^2[s]].
Dzu  = gen(DOM,E0,E0,2,1,0);            Dzud = gen(DOM,E0,E0,2,1,2);
Z    = gen(DOM,E0,{{},{'s'}},1,[1;2],0);
Zd   = gen(DOM,E0,{{},{'s'}},1,[1;2],3);
dmap = draw_dvals({Ad,Bd,Yd,Dzud,Zd});

% Oracle calibration on same-registry products (worked before 09/26/2026),
% and a wrong result it must reject.
assert(sem_product(B,X',Y,dmap) && sem_product(Bd,X',Yd,dmap),'oracle: X''*Y');
assert(~sem_product(2*B,X',Y,dmap),'oracle: accepted 2*(X''*Y) as X''*Y');
assert(~sem_lincomb(B,{A,B},[1,1],dmap),'oracle: accepted B as A + B');
nchk = nchk+3;

nchk = nchk + check_sum('(a) A + X''*Y',A,B,dmap);
nchk = nchk + check_sum('(a) A + X''*Yd',A,Bd,dmap);
nchk = nchk + check_sum('(a) Ad + X''*Y',Ad,B,dmap);
nchk = nchk + check_sum('(a) Ad + X''*Yd',Ad,Bd,dmap);
nchk = nchk + check_prod('(a) Dzu*Z',Dzu,Z,dmap);
nchk = nchk + check_prod('(a) Dzu*Zd',Dzu,Zd,dmap);
nchk = nchk + check_prod('(a) Dzud*Z',Dzud,Z,dmap);
nchk = nchk + check_prod('(a) Z''*Dzu''',Z',Dzu',dmap);
nchk = nchk + check_prod('(a) Zd''*Dzu''',Zd',Dzu',dmap);
nchk = nchk + check_prod('(a) Z''*Dzud''',Z',Dzud',dmap);
fprintf('  passed: (a) R^n with L2[s]      (%d checks so far)\n',nchk);

% ---------------------------------------------------------- (b) {s}, {t}
As  = gen(DOM,{{},{'s'}},E0,[1;2],1,0);     Asd = gen(DOM,{{},{'s'}},E0,[1;2],1,3);
Bt  = gen(DOM,E0,{{},{'t'}},1,[1;2],0);     Btd = gen(DOM,E0,{{},{'t'}},1,[1;2],4);
% Sums over R^n spaces, registries {s} and {t}; and over L2[s], registries
% {s} and {s,t} (F*G passes through L2[t]).
Xs = gen(DOM,{{},{'s'}},Ra,[1;2],[1;2],0);  Ys = gen(DOM,{{},{'s'}},Ra,[1;2],[1;2],0);
Ysd = gen(DOM,{{},{'s'}},Ra,[1;2],[1;2],3);
Xt = gen(DOM,{{},{'t'}},Ra,[1;2],[1;2],0);  Yt = gen(DOM,{{},{'t'}},Ra,[1;2],[1;2],0);
Ks = Xs'*Ys;    Ksd = Xs'*Ysd;      Kt = Xt'*Yt;
Ls  = gen(DOM,{{'s'}},{{'s'}},2,2,0);       Lsd = gen(DOM,{{'s'}},{{'s'}},2,2,5);
F   = gen(DOM,{{'s'}},{{'t'}},2,1,0);       G   = gen(DOM,{{'t'}},{{'s'}},1,2,0);
Lst = F*G;
assert(isequal(Ks.vars,{'s'}) && isequal(Kt.vars,{'t'}) && isequal(Lst.vars,{'s','t'}),...
    'setup (b): expected registries {s}, {t}, {s,t}');
dmap = draw_dvals({Asd,Btd,Ysd,Ksd,Lsd},dmap);
assert(sem_product(Lst,F,G,dmap),'oracle: F*G through L2[t]');
nchk = nchk+1;

nchk = nchk + check_sum('(b) Ks + Kt',Ks,Kt,dmap);
nchk = nchk + check_sum('(b) Ksd + Kt',Ksd,Kt,dmap);
nchk = nchk + check_sum('(b) Ls + F*G',Ls,Lst,dmap);
nchk = nchk + check_sum('(b) Lsd + F*G',Lsd,Lst,dmap);
nchk = nchk + check_prod('(b) As*Bt',As,Bt,dmap);
nchk = nchk + check_prod('(b) Asd*Bt',Asd,Bt,dmap);
nchk = nchk + check_prod('(b) As*Btd',As,Btd,dmap);
nchk = nchk + check_prod('(b) Bt''*As''',Bt',As',dmap);
nchk = nchk + check_prod('(b) Btd''*As''',Btd',As',dmap);
fprintf('  passed: (b) {s} with {t}          (%d checks so far)\n',nchk);

% ---------------------------------------------------------- (c) 1-D, 2-D
A1  = gen(DOM,{{},{'s1'}},{{},{'s1'}},[1;1],[1;1],0);
A1d = gen(DOM,{{},{'s1'}},{{},{'s1'}},[1;1],[1;1],3);
F2  = gen(DOM,{{},{'s1'}},{{'s1','s2'}},[1;1],1,0);
F2d = gen(DOM,{{},{'s1'}},{{'s1','s2'}},[1;1],1,2);
G2  = gen(DOM,{{'s1','s2'}},{{},{'s1'}},1,[1;1],0);
FG  = F2*G2;        FGd = F2d*G2;
assert(isequal(A1.vars,{'s1'}) && isequal(FG.vars,{'s1','s2'}),...
    'setup (c): expected registries {s1} and {s1,s2}');
dmap = draw_dvals({A1d,F2d,FGd},dmap);
assert(sem_product(FG,F2,G2,dmap) && sem_product(A1*A1,A1,A1,dmap),'oracle: 1-D and 2-D');
nchk = nchk+1;

nchk = nchk + check_sum('(c) A1 + F2*G2',A1,FG,dmap);
nchk = nchk + check_sum('(c) A1d + F2*G2',A1d,FG,dmap);
nchk = nchk + check_sum('(c) A1 + F2d*G2',A1,FGd,dmap);
nchk = nchk + check_sum('(c) A1d + F2d*G2',A1d,FGd,dmap);
nchk = nchk + check_prod('(c) A1*F2',A1,F2,dmap);
nchk = nchk + check_prod('(c) A1d*F2',A1d,F2,dmap);
nchk = nchk + check_prod('(c) A1*F2d',A1,F2d,dmap);
nchk = nchk + check_prod('(c) F2''*A1',F2',A1,dmap);
nchk = nchk + check_prod('(c) F2d''*A1',F2d',A1,dmap);
nchk = nchk + check_prod('(c) F2''*A1d',F2',A1d,dmap);
fprintf('  passed: (c) {s1} with {s1,s2}    (%d checks so far)\n',nchk);

% ---------------------------------------------------------- decision lists
% Fixed x decision: every block of the product on the container's list,
% which is the decision factor's own list, name for name and in order.
for PQ = {{Dzud,Z},{Dzu,Zd},{As,Btd},{A1d,F2}}
    [P,Q] = deal(PQ{1}{:});
    R = P*Q;
    D = P;      if isa(Q,'cdopvar'),    D = Q;      end
    assert(isa(R,'cdopvar') && all_blocks_on_list(R) && isequal(R.Zd(:),D.Zd(:)),...
        'fixed x decision product: blocks not on the decision factor''s list');
    nchk = nchk+1;
end
% A decision factor whose sdopvar blocks are on an EMPTY list: the product
% keeps an empty list and its decision blocks stay on it.
C0 = Z.C;
for ii = 1:numel(C0)
    if ~isempty(C0{ii}),    C0{ii} = sopvar2sdopvar(C0{ii});    end
end
Z0 = cdopvar(C0);
assert(isempty(Z0.Zd) && isa(Z0.C{1,2},'sdopvar'),'setup: Z0 should be decision blocks on no list');
R = Dzu*Z0;
assert(isempty(R.Zd) && all_blocks_on_list(R) && sem_product(R,Dzu,Z0,dmap),...
    'fixed x empty-list decision product');
nchk = nchk+2;
fprintf('  passed: decision variable lists  (%d checks so far)\n',nchk);

% ---------------------------------------------------------- errors
% A variable on two domains, from each of the four methods.
DOM2 = DOM;     DOM2.s = [0,2];
L01  = gen(DOM ,{{'s'}},{{'s'}},1,1,0);
L02  = gen(DOM2,{{'s'}},{{'s'}},1,1,0);
L02d = gen(DOM2,{{'s'}},{{'s'}},1,1,3);
expect_error(@() L01 + L02,'copvar:domConflict',"Variable 's'");     % @copvar/plus
expect_error(@() L01 * L02,'copvar:domConflict',"Variable 's'");     % @copvar/mtimes
expect_error(@() L01 + L02d,'cdopvar:domConflict',"Variable 's'");   % @cdopvar/plus
expect_error(@() L02d * L01,'cdopvar:domConflict',"Variable 's'");   % @cdopvar/mtimes
% ... also when the conflict sits in a registry entry no space of the
% operand uses: Ks keeps {s} on [0,1] although it maps R^n -> R^n.
K2 = gen(DOM2,{{},{'s'}},Ra,[1;2],[1;2],0);     K2 = K2'*K2;
expect_error(@() Ks + K2,'copvar:domConflict',"Variable 's'");
% Spaces that differ as name sets are still refused once on the union.
P1  = gen(DOM,E0,{{'s'}},1,1,0);
P2  = gen(DOM,E0,{{'t'}},1,1,0);        P2d = gen(DOM,E0,{{'t'}},1,1,2);
Q2  = gen(DOM,{{'t'}},E0,1,1,0);        Q2d = gen(DOM,{{'t'}},E0,1,1,2);
expect_error(@() P1 + P2,'plus:spaceMismatch','');
expect_error(@() P1 + P2d,'plus:spaceMismatch','');
expect_error(@() P1 * Q2,'mtimes:spaceMismatch','');
expect_error(@() P1 * Q2d,'mtimes:spaceMismatch','');
nchk = nchk+9;

fprintf('test_copvar_registry passed (%d checks).\n',nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = gen(DOM,out,in,dmo,dmi,ndec,occ)
% A random container between the named spaces, each variable on its DOM
% domain: 'copvar' for ndec = 0, else 'cdopvar' with ndec decision
% variables. Degree 2 in 1-D, 1 when a space has two variables.
if nargin<7 || isempty(occ),    occ = true(numel(out),numel(in));   end
v = reshape(unique([out{:}, in{:}]),1,[]);
deg = 2 - (max([0, cellfun(@numel,[out, in])])>=2);
if isempty(v)
    dom = [0,1];
else
    dom = cell2mat(cellfun(@(n) DOM.(n),v(:),'UniformOutput',false));
end
sp = struct('out',{out},'in',{in},'vars',{v});
dm = struct('out',dmo,'in',dmi);
if ndec==0,     P = rand_copvar(sp,dm,dom,deg,0.7,occ);
else,           P = rand_cdopvar(sp,dm,dom,deg,ndec,0.7,occ);
end
end


function n = check_sum(lbl,A,B,dmap)
% A + B and B + A across registries: class, metadata by name, kernels
% against the operands' original blocks, verify, the decision list; and
% A - B, which goes through 'plus'.
for pr = {{A,B},{B,A}}
    [P,Q] = deal(pr{1}{:});
    R = P + Q;
    want = 'copvar';
    if isa(P,'cdopvar') || isa(Q,'cdopvar'),    want = 'cdopvar';   end
    assert(strcmp(class(R),want),'%s: plus returned class %s',lbl,class(R));
    assert(meta_by_name(R,P,Q,'plus'),'%s: plus metadata by name',lbl);
    assert(sem_lincomb(R,{P,Q},[1,1],dmap),'%s: kernels(A+B) ~= kernels(A)+kernels(B)',lbl);
    v = verify(R);
    assert(v.true,'%s: plus result fails verify',lbl);
    assert(~isa(R,'cdopvar') || all_blocks_on_list(R),'%s: plus left a block off the list',lbl);
    Rm = P - Q;
    assert(sem_lincomb(Rm,{P,Q},[1,-1],dmap),'%s: kernels(A-B) ~= kernels(A)-kernels(B)',lbl);
end
n = 12;
end


function n = check_prod(lbl,A,B,dmap)
% A*B across registries: class, metadata by name, (A*B)x = A(Bx) from the
% original blocks, verify, and the decision list of a fixed x decision
% product is the decision factor's.
R = A*B;
want = 'copvar';
if isa(A,'cdopvar') || isa(B,'cdopvar'),    want = 'cdopvar';   end
assert(strcmp(class(R),want),'%s: mtimes returned class %s',lbl,class(R));
assert(meta_by_name(R,A,B,'mtimes'),'%s: mtimes metadata by name',lbl);
assert(sem_product(R,A,B,dmap),'%s: (A*B)x ~= A(Bx)',lbl);
v = verify(R);
assert(v.true,'%s: mtimes result fails verify',lbl);
if isa(R,'cdopvar')
    D = A;      if isa(B,'cdopvar'),    D = B;      end
    assert(all_blocks_on_list(R) && isequal(R.Zd(:),D.Zd(:)),...
        '%s: product blocks not on the decision factor''s list',lbl);
end
n = 5;
end


function tf = meta_by_name(R,A,B,kind)
% Registry = sorted union of the operands' names, each name on its operand
% domain; spaces as name sets: plus keeps both sides of A (= B), mtimes
% takes the output spaces of A and the input spaces of B.
u = reshape(unique([reshape(A.vars,1,[]), reshape(B.vars,1,[])]),1,[]);
tf = isequal(reshape(R.vars,1,[]),u) && size(R.dom,1)==numel(u);
if ~tf,     return,     end
for k = 1:numel(u)
    for Xc = {A,B}
        r = find(strcmp(Xc{1}.vars,u{k}));
        if ~isempty(r),     tf = tf && isequal(R.dom(k,:),Xc{1}.dom(r,:));     end
    end
end
if strcmp(kind,'plus')
    tf = tf && same_spaces(R,'out',A,'out') && same_spaces(R,'in',A,'in') ...
            && same_spaces(R,'out',B,'out') && same_spaces(R,'in',B,'in');
else
    tf = tf && same_spaces(R,'out',A,'out') && same_spaces(R,'in',B,'in');
end
end


function tf = same_spaces(R,sr,X,sx)
% Every space of R on side sr equals that of X on side sx as a SET of
% names, with the same component count.
[mr,dr] = side(R,sr);       [mx,dx] = side(X,sx);
tf = size(mr,1)==size(mx,1) && isequal(dr(:),dx(:));
for k = 1:size(mr,1)
    if ~tf,     return,     end
    tf = isequal(names(R,mr(k,:)),names(X,mx(k,:)));
end
end


function [m,d] = side(P,s)
if strcmp(s,'out'),     m = P.space_out;    d = P.dim_out;
else,                   m = P.space_in;     d = P.dim_in;
end
end


function c = names(P,mask)
c = reshape(sort(P.vars(mask)),1,[]);
end


function tf = sem_product(R,A,B,dmap)
% (A*B)x = A(Bx) for two random polynomial inputs, R's blocks on the left
% and A's and B's original blocks on the right, each applied by apply_blk.
M = size(R.C,1);    N = size(R.C,2);    K = size(A.C,2);
tf = size(B.C,1)==K && size(B.C,2)==N && size(A.C,1)==M;
for trial = 1:2
    if ~tf,     return,     end
    x = cell(1,N);
    for j = 1:N,    x{j} = rand_testfun(names(B,B.space_in(j,:)),B.dim_in(j));     end
    z = cell(1,K);
    for k = 1:K
        z{k} = polynomial(zeros(B.dim_out(k),1));
        for j = 1:N
            if ~isempty(B.C{k,j}),  z{k} = z{k} + apply_blk(B.C{k,j},x{j},dmap);   end
        end
    end
    for i = 1:M
        lhs = polynomial(zeros(R.dim_out(i),1));    rhs = lhs;
        for j = 1:N
            if ~isempty(R.C{i,j}),  lhs = lhs + apply_blk(R.C{i,j},x{j},dmap);     end
        end
        for k = 1:K
            if ~isempty(A.C{i,k}),  rhs = rhs + apply_blk(A.C{i,k},z{k},dmap);     end
        end
        sc = max([1, maxcoef(lhs), maxcoef(rhs)]);
        if maxcoef(lhs-rhs) > 1e-9*sc,  tf = false;     return,     end
    end
end
end


function y = apply_blk(B,x,dmap)
% (Bx) for one block, from its kernels, by the form in the 'sopvar.m' header:
%
%   (Bx)(S2,S3) = sum_gam int_S1 int_S3dum I_gam(S3-S3dum)
%                           K_gam(S2,S3; S3dum,S1) x(S3dum,S1)
%
% index 1 a multiplier delta(s-s'), 2 the lower integral s'<=s, 3 the upper,
% as in 'pi_kernel_pairing'. The parameter cell is indexed over S3 sorted.
O = B.vars.out;     I = B.vars.in;
S3 = intersect(O,I);        S1 = setdiff(I,S3);
Kc = kern(B,dmap);
if numel(Kc)~=3^numel(S3)
    error('test_copvar_registry: %d parameter cells for %d shared variables.',numel(Kc),numel(S3));
end
xd = x;
for t = 1:numel(I)
    xd = subs(xd,polynomial(I(t)),polynomial({[I{t} '_dum']}));
end
y = polynomial(zeros(B.dims(1),1));
for q = 1:numel(Kc)
    gam = pi_gamma_index(q,numel(S3));
    e = Kc{q}*xd;
    for t = 1:numel(S3)
        d = S3{t};      lim = B.dom.in(strcmp(I,d),:);
        sp = polynomial({[d '_dum']});      s = polynomial({d});
        switch gam(t)
            case 1,     e = subs(e,sp,s);
            case 2,     e = int(e,sp,lim(1),s);
            case 3,     e = int(e,sp,s,lim(2));
        end
    end
    for t = 1:numel(S1)
        d = S1{t};      lim = B.dom.in(strcmp(I,d),:);
        e = int(e,polynomial({[d '_dum']}),lim(1),lim(2));
    end
    y = y + e;
end
end


function x = rand_testfun(v,m)
% m x 1 random polynomial in the variables v: constant, linear and square
% terms per variable, and a cross term in two variables.
x = polynomial(randn(m,1));
for k = 1:numel(v)
    s = polynomial(v(k));
    x = x + randn(m,1)*s + randn(m,1)*s^2;
end
if numel(v)>=2,     x = x + randn(m,1)*(polynomial(v(1))*polynomial(v(2)));   end
end


function dmap = draw_dvals(Ps,dmap)
% One random value per decision variable NAME, kept across calls.
if nargin<2,    dmap = containers.Map('KeyType','char','ValueType','double');  end
for k = 1:numel(Ps)
    z = Ps{k}.Zd;
    for i = 1:numel(z)
        if ~isKey(dmap,z{i}),   dmap(z{i}) = randn;     end
    end
end
end


function K = kern(B,dmap)
% Kernels of one block, a decision block at the drawn values.
if isa(B,'sdopvar')
    z = B.Zd;
    dval = zeros(numel(z),1);
    for i = 1:numel(z),     dval(i) = dmap(z{i});   end
    K = pi_blk_kernels(B,dval);
else
    K = pi_blk_kernels(B,[]);
end
end


function m = maxcoef(X)
if isnumeric(X)
    m = full(max(abs(X(:))));
else
    X = polynomial(X);
    cf = X.coefficient;
    m = full(max(abs(cf(:))));
end
if isempty(m),  m = 0;  end
end


function tf = sem_lincomb(R,Xs,cs,dmap)
% Every block of R equals sum_k cs(k)*Xs{k} blockwise, as kernels.
tf = isequal(size(R.C),size(Xs{1}.C));
if ~tf,     return,     end
for ii = 1:numel(R.C)
    Kr = {};
    if ~isempty(R.C{ii}),   Kr = kern(R.C{ii},dmap);    end
    Ks = cell(1,numel(Xs));
    for k = 1:numel(Xs)
        if ~isempty(Xs{k}.C{ii}),   Ks{k} = kern(Xs{k}.C{ii},dmap);    end
    end
    ng = max([numel(Kr),cellfun(@numel,Ks)]);
    sc = 1;
    for g = 1:ng
        acc = 0;
        for k = 1:numel(Xs)
            if ~isempty(Ks{k}),     acc = acc + cs(k)*Ks{k}{g};     sc = max(sc,maxcoef(Ks{k}{g}));   end
        end
        if ~isempty(Kr),    acc = Kr{g} - acc;      else,   acc = -acc;     end
        if maxcoef(acc) > 1e-9*sc,  tf = false;     return,     end
    end
end
end


function tf = all_blocks_on_list(P)
tf = true;
for ii = 1:numel(P.C)
    if isa(P.C{ii},'sdopvar') && ~isequal(P.C{ii}.Zd(:),P.Zd(:))
        tf = false;     return
    end
end
end


function expect_error(f,id,msg)
% f must raise error 'id', and its message must contain 'msg'.
try
    f();
catch ME
    if ~strcmp(ME.identifier,id)
        error('test_copvar_registry: expected error ''%s'', got ''%s'': %s',id,ME.identifier,ME.message);
    end
    if ~contains(ME.message,msg)
        error('test_copvar_registry: error ''%s'' does not mention "%s": %s',id,msg,ME.message);
    end
    return
end
error('test_copvar_registry: expected error ''%s'', but none was raised.',id);
end

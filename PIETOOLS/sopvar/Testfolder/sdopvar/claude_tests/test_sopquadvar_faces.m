function R = test_sopquadvar_faces()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_SOPQUADVAR_FACES() asserting tests of options.psatz in
% sopquadvar / possopvar (MMP 09/28/2026): the value check and the face
% codes 2k+1 / 2k+2, k over SORTED S3, mirroring copquadvar's.
%
%  (i)   value checks, n3 = 0..3: 2 (poslpivar_2d's ball), 2n3+3, 1.5,
%        2n3+2.5, [3 4], -3 and NaN are rejected by the psatz check, as is
%        type 'sym' with any face code; 0, 1, true, false, int8(1) and every
%        code 3..2n3+2 are accepted; int32/uint8/single codes build exactly
%        the operator of the double code (two basis operators, degree 0);
%  (ii)  semantics at degree 0 by independent quadrature: for random
%        Q >= 0 put into the declared variables (by Qcell),
%            <x, Pop x>  =  int g(th) (Z x)(th)' Q (Z x)(th) dth,
%        the left side EXACT from the stored kernels (pi_sdopvar_kernels,
%        pi_kernel_pairing), the right side by Gauss quadrature straight
%        from the definition (Z_alpha x)(th) = int I_alpha(th - s) x(s) ds,
%        with g built HERE from the header: code 2k+1 (th_k-a_k)/L_k, 2k+2
%        (b_k-th_k)/L_k, th_k the k-th entry of sort(vars), 1 the product,
%        0 none. n3 = 1, 2, 3 on a non-unit box, names given UNSORTED so
%        that sorted and caller order differ. Controls that must NOT match:
%        the other face of the same direction, and (n3 >= 2) the face read
%        in caller order, vars{k};
%  (iii) semantics at degree >= 1, closed form: the bilinear identity
%            <y, Pop x> = sum_ij int g (Z_i y)' Q_ij (Z_j x) dth
%        for a random (indefinite) Q, both sides exact, the basis rebuilt
%        from its documented description as test_possopvar does; n3 = 1
%        (degree 2), 2 (degree 1, unsorted names; and with sep + include),
%        3 (degree 1 in theta, 0 in s, unsorted names);
%  (iv)  CROSS-CHECK against copquadvar: for a single space possopvar(c) and
%        poscopvar(c) must define the same operator. Their Gram bases may
%        be ordered differently, so they are compared by semantics: the
%        blocks are matched by multi-index, both are given the same Q
%        (fully random at degree 0; at degree >= 1 from the family
%        Q_ij = M_ij 1 1' + delta_ij D_i I, M >= 0, D > 0, which is invariant
%        under any reordering within a block), and the resulting kernels
%        are compared cell by cell (multiplier directions evaluated on the
%        diagonal, where delta puts them). Controls: poscopvar with the other
%        face, and with the caller-order code, must differ.
%
% Initial coding MMP, 09/28/2026 ((ii) adapted from test_copquadvar_faces,
%   (iii) from test_possopvar).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rs = rng;   rng(4243);
w0 = warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
R = struct('check',{},'err',{},'tol',{});
Dnu = [0 1; -1 2; 0.5 3];
VU = {{'s1'}, {'y','x'}, {'z','x','y'}};        % unsorted from n3 = 2 on
T0 = tic;

% % % (i) value checks.
for n3 = 0:3
    vars = arrayfun(@(i) sprintf('s%d',i),1:n3,'UniformOutput',false);
    dom = Dnu(1:n3,:);
    % Two basis operators (multiplier, all-lower) keep the builds cheap; the
    % face still shows, in the kernel of both.
    o = struct();   if n3>0,    o.include = [ones(1,n3); 2*ones(1,n3)];     end
    for bad = {2, 2*n3+3, 1.5, 2*n3+2.5, [3 4], -3, NaN}
        ok = errs_psatz(@() possopvar(mkprog(vars,dom),1,vars,dom,0,setf(o,'psatz',bad{1})));
        R(end+1) = rec(sprintf('(i) n3=%d psatz=%s rejected',n3,mat2str(bad{1})),~ok,0); %#ok<AGROW>
    end
    for c = 3:2*n3+2
        ok = errs_psatz(@() sopquadvar(mkprog(vars,dom),1,vars,dom,0,setf(setf(o,'psatz',c),'type','sym')));
        R(end+1) = rec(sprintf('(i) n3=%d psatz=%d, type sym rejected',n3,c),~ok,0); %#ok<AGROW>
    end
    for c = [{0,1,true,false,int8(1)}, num2cell(3:2*n3+2)]
        ok = ~errs(@() possopvar(mkprog(vars,dom),1,vars,dom,0,setf(o,'psatz',c{1})));
        R(end+1) = rec(sprintf('(i) n3=%d psatz=%s(%g) accepted',n3,class(c{1}),double(c{1})),~ok,0); %#ok<AGROW>
    end
    % Integer classes pass the value check, so they must build the operator
    % of the double code; without the conversion to double, floor((c-1)/2)
    % rounds in integer arithmetic and picks the next direction.
    for c = 3:2*n3+2
        [~,P0,Q0] = possopvar(mkprog(vars,dom),1,vars,dom,0,setf(o,'psatz',c));
        for cls = {@int32,@uint8,@single}
            [~,P1,Q1] = possopvar(mkprog(vars,dom),1,vars,dom,0,setf(o,'psatz',cls{1}(c)));
            R(end+1) = rec(sprintf('(i) n3=%d psatz=%s(%d) equals double',n3,func2str(cls{1}),c), ...
                           ~(isequal(P0,P1) && isequal(Q0,Q1)),0); %#ok<AGROW>
        end
    end
end
fprintf('[test_sopquadvar_faces] (i) done, %.1f s\n',toc(T0));

% % % (ii) semantics at degree 0, quadrature, random Q >= 0.
for n3 = 1:3
    vars = VU{n3};      dom = Dnu(1:n3,:);      % dom rows follow 'vars'
    nq = 6;             % Gauss exact to degree 11 per direction; integrand <= 8
    [Xq,wq] = grid(dom,nq);
    for code = [0 1 3:2*n3+2]
        [g,gw,gc] = weight(code,vars,dom);
        [~,Pop,Qcell,al] = possopvar(mkprog(vars,dom),1,vars,dom,0,struct('psatz',code));
        nb = size(al,1);
        G = randn(nb);  Q = G*G';                                           % random PSD
        dval = assign_q(Qcell,Q,Pop.Zd);
        [xp,xh] = randpoly(vars,4);
        lhs = pi_kernel_pairing(pi_sdopvar_kernels(Pop,dval),xp,xp,Pop.vars.in,Pop.dom.in);
        Zx = zeros(size(Xq,1),nb);                                          % (Z_alpha x)(th)
        for c = 1:nb,   Zx(:,c) = basis_apply(al(c,:),xh,Xq,dom,nq);    end
        q = sum((Zx*Q).*Zx,2);
        rhs = sum(wq.*g(Xq).*q);
        R(end+1) = rec(sprintf('(ii) n3=%d psatz=%d: <x,Px> = int g (Zx)''Q(Zx)',n3,code), ...
                       abs(lhs-rhs)/max(abs(rhs),eps),1e-10); %#ok<AGROW>
        R(end+1) = rec(sprintf('(ii) n3=%d psatz=%d: <x,Px> >= 0',n3,code), ...
                       max(0,-lhs)/max(abs(rhs),eps),1e-12); %#ok<AGROW>
        if ~isempty(gw)
            rw = sum(wq.*gw(Xq).*q);    d = abs(lhs-rw)/max(abs(rw),eps);
            R(end+1) = rec(sprintf('(ii) n3=%d psatz=%d: control, other face differs (%.2g)',n3,code,d), ...
                           d<1e-3,0); %#ok<AGROW>
        end
        if ~isempty(gc)
            rc = sum(wq.*gc(Xq).*q);    d = abs(lhs-rc)/max(abs(rc),eps);
            R(end+1) = rec(sprintf('(ii) n3=%d psatz=%d: control, caller-order face differs (%.2g)',n3,code,d), ...
                           d<1e-3,0); %#ok<AGROW>
        end
    end
end
fprintf('[test_sopquadvar_faces] (ii) done, %.1f s\n',toc(T0));

% % % (iii) degree >= 1, closed form, random indefinite Q.
% n3 = 3: theta-degree 1 only (s-degree 0); the product weight there costs
% ~65 s at joint degree 1 with s-degree 1 and is unchanged code.
J1 = struct('int',1,'mult',0,'joint',1);
cs = { {1, {'s1'},      Dnu(1,:),   2,  struct(),                                [0 1 3 4]}
       {2, {'y','x'},   Dnu(1:2,:), 1,  struct(),                                [0 1 3:6]}
       {2, {'y','x'},   Dnu(1:2,:), 1,  struct('sep',[true false],'include',[1 1;4 2;4 3;1 3]), [3 6]}
       {3, {'z','x','y'}, Dnu,      J1, struct(),                                [1 4 7]} };
for ic = 1:numel(cs)
    [n3,vars,dom,deg,op,codes] = deal(cs{ic}{:});
    for code = codes
        op.psatz = code;
        [~,Pop,Qcell,al] = possopvar(mkprog(vars,dom),1,vars,dom,deg,op);
        dvars = cellstr(string(Pop.Zd(:)));
        dval = randn(numel(dvars),1);
        xp = randpoly(vars,3);      yp = randpoly(vars,3);
        lhs = pi_kernel_pairing(pi_sdopvar_kernels(Pop,dval),xp,yp,Pop.vars.in,Pop.dom.in);
        rhs = pair_factored(Qcell,al,xp,yp,vars,dom,deg,code,dval,dvars);
        R(end+1) = rec(sprintf('(iii) n3=%d case %d psatz=%d: <y,Px> identity',n3,ic,code), ...
                       abs(lhs-rhs)/max(1,abs(rhs)),1e-9); %#ok<AGROW>
    end
end
fprintf('[test_sopquadvar_faces] (iii) done, %.1f s\n',toc(T0));

% % % (iv) cross-check against poscopvar, single space.
cs = { {{'s1'},        [-1 2],     1}
       {{'y','x'},     Dnu(1:2,:), 0}
       {{'y','x'},     Dnu(1:2,:), 1}
       {{'s1','s2'},   Dnu(1:2,:), 1}
       {{'z','x','y'}, Dnu,        0} };
for ic = 1:numel(cs)
    [vars,dom,deg] = deal(cs{ic}{:});
    n3 = numel(vars);   sv = sort(vars);
    pos = cellfun(@(v) find(strcmp(vars,v)),sv);        % caller column of sorted var r
    doms = struct('vars',{vars},'dom',dom);             % copquadvar: name-keyed domain
    K2 = containers.Map('KeyType','double','ValueType','any');
    K1s = containers.Map('KeyType','double','ValueType','any');
    for code = [0 1 3:2*n3+2]
        [~,P1,Q1,al] = possopvar(mkprog(vars,dom),1,vars,dom,deg,struct('psatz',code));
        K1s(code) = {P1,Q1};
        [~,P2,Q2,bl] = poscopvar(mkprog(sv,dom(pos,:)),1,vars,doms,deg,struct('psatz',code));
        B2 = P2.C{1};
        ok = isequal(P1.vars.in(:),B2.vars.in(:)) && isequal(P1.vars.out(:),B2.vars.out(:)) ...
             && isequal(P1.dom.in,B2.dom.in);
        R(end+1) = rec(sprintf('(iv) case %d psatz=%d: same vars and dom',ic,code),~ok,0); %#ok<AGROW>
        % Blocks by multi-index: al is over 'vars', bl(:,2:end) over sorted.
        [tf,p] = ismember(al(:,pos),bl(:,2:end),'rows');
        nb = size(al,1);
        ok = all(tf) && numel(unique(p))==nb && size(bl,1)==nb;
        R(end+1) = rec(sprintf('(iv) case %d psatz=%d: same basis operators',ic,code),~ok,0); %#ok<AGROW>
        if ~ok,     continue,   end
        [QA,QB] = shared_q(Q1,Q2,p,deg);
        d1 = assign_q(Q1,QA,P1.Zd);
        d2 = assign_q(Q2,QB,B2.Zd);
        K1 = pi_sdopvar_kernels(P1,d1);
        Kc = pi_sdopvar_kernels(B2,d2);
        K2(code) = {Kc,QA,QB};
        R(end+1) = rec(sprintf('(iv) case %d deg %d psatz=%d: possopvar = poscopvar',ic,deg_of(deg),code), ...
                       kdiff(K1,Kc,P1.vars.in),1e-10); %#ok<AGROW>
    end
    % Controls: poscopvar's other face, and its caller-order code, differ.
    for code = 3:2*n3+2
        E = K1s(code);  P1 = E{1};      Q1 = E{2};
        k = floor((code-1)/2);  par = 2-mod(code,2);    % par 1 lower, 2 upper
        % Other face of the same direction, and the code poscopvar gives
        % the face of vars{k} (sorted index r), i.e. a caller-order reading.
        r = find(pos==k);
        alt = [2*k+3-par, 2*r+par];
        alt = unique(alt(alt~=code));
        for ca = alt
            E = K2(ca);     Kc = E{1};      QA = E{2};
            d1 = assign_q(Q1,QA,P1.Zd);
            K1 = pi_sdopvar_kernels(P1,d1);     d = kdiff(K1,Kc,P1.vars.in);
            R(end+1) = rec(sprintf('(iv) case %d psatz=%d vs poscopvar %d: control differs (%.2g)',ic,code,ca,d), ...
                           d<1e-3,0); %#ok<AGROW>
        end
    end
end
fprintf('[test_sopquadvar_faces] (iv) done, %.1f s\n',toc(T0));

warning(w0);    rng(rs);
ok = true;
for k = 1:numel(R)
    pass = R(k).err<=R(k).tol;     ok = ok && pass;
    fprintf('%-66s %10.2e  %s\n',R(k).check,R(k).err,char(string(pass)));
end
assert(ok,'test_sopquadvar_faces: a check failed.');
fprintf('test_sopquadvar_faces: all %d checks pass (%.1f s).\n',numel(R),toc(T0));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = rec(check,err,tol)
r = struct('check',check,'err',double(err),'tol',tol);
end

function s = setf(s,f,v)
s.(f) = v;
end

function tf = errs(f)
% True if F() errors.
tf = false;
try,    f();    catch,  tf = true;  end
end

function tf = errs_psatz(f)
% True if F() errors in the psatz check itself (its message names 'psatz'),
% not somewhere else.
tf = false;
try,    f();    catch ME,   tf = contains(ME.message,"'psatz'");    end
end

function d = deg_of(deg)
if isstruct(deg),   d = deg.joint;  else,   d = deg;    end
end

function [g,gw,gc] = weight(code,vars,dom)
% The weight of CODE from sopquadvar's header, on points T whose columns
% follow 'vars': 0 none, 1 the product, 2k+1 (th_k-a_k)/L_k, 2k+2
% (b_k-th_k)/L_k with th_k the k-th entry of sort(vars). GW: the other face
% of that direction. GC: the same face of vars{k} (caller order), [] where
% that is the same variable.
gw = [];    gc = [];
switch code
    case 0,     g = @(T) ones(size(T,1),1);                         return
    case 1,     g = @(T) prod((T-dom(:,1)').*(dom(:,2)'-T),2);      return
end
k = floor((code-1)/2);  sv = sort(vars);    j = find(strcmp(vars,sv{k}));
[g,gw] = face(j,mod(code,2)==1,dom);
if j~=k,    gc = face(k,mod(code,2)==1,dom);    end
end

function [g,gw] = face(j,lower,dom)
L = dom(j,2)-dom(j,1);
lo = @(T) (T(:,j)-dom(j,1))/L;      hi = @(T) (dom(j,2)-T(:,j))/L;
if lower,   g = lo;     gw = hi;
else,       g = hi;     gw = lo;
end
end

function prog = mkprog(vars,dom)
% lpiprogram refuses n3 > 2 (and n3 = 0); for those what it would return.
if numel(vars)>=1 && numel(vars)<=2
    prog = lpiprogram(polynomial(vars(:)),[],dom);  return
end
prog = sosprogram(polynomial([]),dpvar(zeros(0,1)));
if ~isempty(vars)
    prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(strcat(vars(:),'_dum'))];
end
prog.dom = dom;
end

function [xp,xh] = randpoly(vars,nt)
% Random scalar polynomial of degree <= 2 per variable, as a 'polynomial'
% in VARS and as a handle on points whose columns follow VARS.
n3 = numel(vars);
E = unique(randi([0 2],nt,n3),'rows');  c = randn(size(E,1),1);
xp = polynomial(c,E,vars(:),[1 1]);
xh = @(X) sum(c'.*prod(reshape(X,[],1,n3).^reshape(E,1,[],n3),3),2);
end

function dval = assign_q(Qcell,Q,Zd)
% Values of the decision variables Zd that put the symmetric matrix Q into
% the Gram matrix named by Qcell (blocks in Qcell order).
nb = size(Qcell,1);     sz = zeros(1,nb);
for i = 1:nb,   sz(i) = size(Qcell{i,1},1);     end
off = [0,cumsum(sz)];
names = cell(off(end));
for i = 1:nb
    for j = 1:nb
        names(off(i)+1:off(i+1),off(j)+1:off(j+1)) = cellstr(string(Qcell{i,j}));
    end
end
zd = cellstr(string(Zd(:)));
[tf,loc] = ismember(names,zd);
assert(all(tf(:)),'test_sopquadvar_faces: a Qcell name is not in Zd.');
dval = zeros(numel(zd),1);  dval(loc(:)) = Q(:);
assert(norm(reshape(dval(loc),size(Q))-Q,'fro')<=1e-12*max(1,norm(Q,'fro')), ...
       'test_sopquadvar_faces: Qcell naming is not symmetric.');
end

function [QA,QB] = shared_q(QcA,QcB,p,deg)
% One Gram matrix for two parameterizations whose block i of A is block p(i)
% of B. At degree 0 each block is 1x1 and Q is fully random; otherwise
% Q_ij = M_ij 1 1' + delta_ij D_i I, PSD for M >= 0, D > 0, and invariant
% under any reordering of the monomials inside a block.
nb = numel(p);  sz = zeros(1,nb);   szB = zeros(1,nb);
for i = 1:nb,   sz(i) = size(QcA{i,1},1);   szB(p(i)) = size(QcB{p(i),1},1);    end
assert(isequal(sz,szB(p)),'test_sopquadvar_faces: block sizes differ.');
G = randn(nb);  M = G*G'/nb;    D = 0.5+rand(1,nb);
if ~isstruct(deg) && deg==0,    M = M + diag(D);    D = zeros(1,nb);    end
QA = blocks(M,D,sz,1:nb);
QB = zeros(size(QA));   offB = [0,cumsum(szB)];     offA = [0,cumsum(sz)];
for i = 1:nb
    for j = 1:nb
        QB(offB(p(i))+1:offB(p(i)+1),offB(p(j))+1:offB(p(j)+1)) = ...
            QA(offA(i)+1:offA(i+1),offA(j)+1:offA(j+1));
    end
end
end

function Q = blocks(M,D,sz,~)
off = [0,cumsum(sz)];   Q = zeros(off(end));
for i = 1:numel(sz)
    for j = 1:numel(sz)
        B = M(i,j)*ones(sz(i),sz(j));
        if i==j,    B = B + D(i)*eye(sz(i));    end
        Q(off(i)+1:off(i+1),off(j)+1:off(j+1)) = B;
    end
end
end

function e = kdiff(K1,K2,vars)
% Largest coefficient of K1-K2 over all cells, relative to the largest of
% K1. In a multiplier direction delta(s_k-s_k') puts the kernel on the
% diagonal, so there s_k' := s_k first: the comparison is of the operator,
% not of how its kernel happens to be stored.
n3 = numel(vars);   num = 0;    den = 0;
for k = 1:numel(K1)
    gam = pi_gamma_index(k,n3);
    A = K1{k};  B = K2{k};
    for d = find(gam==1)
        A = subs(A,polynomial({[vars{d} '_dum']}),polynomial(vars(d)));
        B = subs(B,polynomial({[vars{d} '_dum']}),polynomial(vars(d)));
    end
    cD = coefs(A-B);    cA = coefs(A);          % subs may return a double
    if ~isempty(cD),    num = max(num,full(max(abs(cD(:)))));   end
    if ~isempty(cA),    den = max(den,full(max(abs(cA(:)))));   end
end
e = num/max(den,eps);
end

function c = coefs(p)
if isa(p,'polynomial'),     c = p.coefficient;  else,   c = p;  end
end

function y = basis_apply(alpha,x,Th,dom,nq)
% (Z_alpha x)(th) = int I_alpha(th - s) x(s) ds, degree-0 basis (kernel 1);
% alpha, dom and the columns of Th follow 'vars'.
[xg,wg] = gl01(nq);     N = numel(alpha);   y = zeros(size(Th,1),1);
for k = 1:size(Th,1)
    th = Th(k,:);   pts = cell(1,N);    wts = cell(1,N);
    for d = 1:N
        switch alpha(d)
            case 1,     pts{d} = th(d);     wts{d} = 1;
            case 2,     pts{d} = dom(d,1)+(th(d)-dom(d,1))*xg;  wts{d} = (th(d)-dom(d,1))*wg;
            case 3,     pts{d} = th(d)+(dom(d,2)-th(d))*xg;     wts{d} = (dom(d,2)-th(d))*wg;
            otherwise,  error('basis_apply:alpha','Unexpected multi-index %d.',alpha(d));
        end
    end
    [P,Wt] = tens(pts,wts);
    y(k) = sum(Wt.*x(P));
end
end

function [X,w] = grid(dom,n)
[xg,wg] = gl01(n);  N = size(dom,1);    pts = cell(1,N);    wts = cell(1,N);
for d = 1:N,    L = dom(d,2)-dom(d,1);  pts{d} = dom(d,1)+L*xg;  wts{d} = L*wg;  end
[X,w] = tens(pts,wts);
end

function [X,w] = tens(pts,wts)
% Tensor grid of per-direction nodes/weights; columns in direction order.
N = numel(pts);     X = zeros(1,0);     w = 1;
for d = 1:N
    n = numel(pts{d});  m = size(X,1);
    X = [repmat(X,n,1), repelem(pts{d}(:),m,1)];
    w = repmat(w,n,1).*repelem(wts{d}(:),m,1);
end
end

function [x,w] = gl01(n)
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*V(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end

function val = pair_factored(Qcell,alpha_list,xfun,yfun,vars,dom,deg,code,dval,dvars)
% sum_ij int_theta g(theta) (Z_alpha_i y)(theta)' Q_ij (Z_alpha_j x)(theta),
% in closed form, g from the header definition (see WEIGHT), the basis
% rebuilt from its documented description over [theta; s] in 'vars' order.
n3 = numel(vars);
vars_int = strcat(vars,'_int');
nblk = size(alpha_list,1);
Zy = cell(1,nblk);  Zx = cell(1,nblk);
for i = 1:nblk
    Zi = basis_monomials(alpha_list(i,:),deg,vars_int,vars);
    Zy{i} = apply_basis_op(Zi,yfun,alpha_list(i,:),vars,vars_int,dom);
    Zx{i} = apply_basis_op(Zi,xfun,alpha_list(i,:),vars,vars_int,dom);
end
gfun = polynomial(1);
if code==1
    for k = 1:n3
        thk = polynomial(vars_int(k));
        gfun = gfun*(thk-dom(k,1))*(dom(k,2)-thk);
    end
elseif code>=3
    k = floor((code-1)/2);  sv = sort(vars);    j = find(strcmp(vars,sv{k}));
    thj = polynomial(vars_int(j));      L = dom(j,2)-dom(j,1);
    if mod(code,2)==1,  gfun = (thj-dom(j,1))/L;
    else,               gfun = (dom(j,2)-thj)/L;
    end
end
expr = polynomial(0);
for i = 1:nblk
    for j = 1:nblk
        names = cellstr(string(Qcell{i,j}));
        [tf,loc] = ismember(names,dvars);   assert(all(tf(:)));
        expr = expr + Zy{i}.'*reshape(dval(loc),size(names))*Zx{j};
    end
end
expr = gfun*expr;
for k = 1:n3
    expr = int(expr,polynomial(vars_int(k)),dom(k,1),dom(k,2));
end
val = double(expr);
end

function Zi = basis_monomials(alpha,deg,vars_int,vars_out)
% Z^alpha(theta,s) as documented in sopquadvar: per-variable caps 'int' in
% theta and 'mult' in s (0 where alpha_k = 1), total degree <= 'joint',
% exponents over [theta; s] in 'vars' order, rows sorted.
n3 = numel(vars_out);
if isnumeric(deg),  deg = struct('int',deg);    end
if ~isfield(deg,'mult') || isempty(deg.mult),   deg.mult = deg.int;     end
ci = deg.int.*ones(1,n3);   cm = deg.mult.*ones(1,n3);  cm(alpha==1) = 0;
caps = [ci,cm];
if ~isfield(deg,'joint') || isempty(deg.joint),     jc = sum(caps);
else,                                               jc = deg.joint;
end
grids = arrayfun(@(c) (0:c)',caps,'UniformOutput',false);
sub = cell(1,numel(caps));  [sub{:}] = ndgrid(grids{:});
E = cell2mat(cellfun(@(s) s(:),sub,'UniformOutput',false));
E = sortrows(E(sum(E,2)<=jc,:));
T = size(E,1);
Zi = polynomial(speye(T),E,[vars_int(:);vars_out(:)],[T,1]);
end

function Zw = apply_basis_op(Zi,w,alpha,vars,vars_int,dom)
% (Z_alpha w)(theta) = int_s I_alpha(theta-s) Z^alpha(theta,s) w(s) ds in
% closed form; alpha 4 is the full-domain integral of a separable direction.
T = size(Zi,1);     Zw = polynomial(zeros(T,1));
for u = 1:T
    expr = Zi(u)*w;
    for k = 1:numel(vars)
        sk = polynomial(vars(k));   thk = polynomial(vars_int(k));
        switch alpha(k)
            case 1,     expr = subs(expr,sk,thk);
            case 2,     expr = int(expr,sk,dom(k,1),thk);
            case 3,     expr = int(expr,sk,thk,dom(k,2));
            case 4,     expr = int(expr,sk,dom(k,1),dom(k,2));
        end
    end
    Zw(u) = expr;
end
end

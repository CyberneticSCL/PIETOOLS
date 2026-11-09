function R = test_heatNd_pie(dostock)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_HEATND_PIE(DOSTOCK) asserting semantic tests of the PIE operators
% T, A built by HEATND_PIE, for N = 1, 2, 3, each object tested against its
% OWN definition, never against an inverse routine (CLAUDE.md sec. 4):
%
%  (a) T v for a random non-separable polynomial v equals the unique
%      solution u of D^delta u = v with the BCs, built from the 1-D ODE
%      u'' = s^p (particular s^(p+2)/((p+1)(p+2)) plus the affine part
%      fixed by the BC) - no kernel formula involved. Implies D^delta T v =
%      v and the BCs;
%  (b) T D^delta phi = phi for a polynomial phi IN the domain (each 1-D
%      factor satisfies its BC), i.e. T is also a left inverse;
%  (c) BC traces of T v evaluated directly: Dirichlet values at the face,
%      Neumann by a one-sided 7-point stencil exact to degree 6;
%  (d) A v equals r T v + sum_i prod_{j~=i} T_j v, built from (a)'s exact
%      1-D solutions with the identity in direction i; (d0) A0 v (the r = 0
%      part HEATND_LPI builds from) equals the sum alone;
%  (e) eigenfunction: T phi1 = prod(-1/lambda_i) phi1 and A phi1 =
%      (r - lambda1) T phi1 for the lowest mode (quadrature, nq = 20);
%  (f) N = 2 only: T against the paper's printed Ex. 3 kernel (four terms);
%  (g) DOSTOCK (default true), N = 1, 2: STOCK convert(PDE,'pie') T and A
%      (opvar / opvar2d), each evaluated from ITS OWN class definition by
%      the local quadrature below (no converter), applied to the same
%      polynomial at the same points, against the container T and A.
% Operators are applied with HEATND_APPLY (sopvar definition by quadrature).
% A non-default BC set and a non-unit domain are included at N = 1, 2.
%
% OUTPUT R: table-like struct array, one row per (case, check), with the
% relative max error; the function errors on the first failed assertion
% after printing every row.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): check (d0) of pie.A0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(dostock),    dostock = true;     end
info = heatNd_path();
fprintf('test_heatNd_pie: maxNumCompThreads = %d, MATLAB %s\n',info.threads,info.version);
rng(20260927);
tol = 1e-10;
cases = { {1,{'DD'},[0 1]}, {1,{'DN'},[0 1]}, {1,{'ND'},[0 1]}, {1,{'DD'},[-1 2]}, ...
          {2,{'DD','DN'},[0 1;0 1]}, {2,{'ND','DD'},[0 2;-1 1]}, ...
          {3,{'DD','DN','DN'},[0 1;0 1;0 1]}, {3,{'DN','ND','DD'},[0 1;1 2;0 3]} };
R = struct('case',{},'check',{},'err',{},'tol',{});
r = 3.7;
for c = 1:numel(cases)
    [N,bc,dom] = cases{c}{:};
    pie = heatNd_pie(N,r,bc,dom);
    tag = sprintf('N=%d %s dom=%s',N,strjoin(bc,'x'),mat2str(dom));
    S = dom(:,1)' + rand(6,N).*(dom(:,2)-dom(:,1))';       % interior points
    % Random non-separable polynomial, degree <= 3 per variable.
    nt = 5;     Pw = randi([0 3],nt,N);     cf = randn(nt,1);
    v = @(X) polyval_terms(X,Pw,cf);
    uex = @(X) exact_inv(X,Pw,cf,bc,dom,1:N);
    % (a)
    Tv = heatNd_apply(pie.T,v,S);
    R(end+1) = row(tag,'(a) T v = D^-delta v',relerr(Tv,uex(S)),tol);   %#ok<AGROW>
    % (b) phi = prod_i phi_i(s_i), phi_i in the 1-D domain, times (1+s_i/3).
    [ph,d2ph] = domain_poly(bc,dom);
    Tph = heatNd_apply(pie.T,d2ph,S);
    R(end+1) = row(tag,'(b) T D^delta phi = phi',relerr(Tph,ph(S)),tol);    %#ok<AGROW>
    % (c) traces.
    R(end+1) = row(tag,'(c) BC traces of T v',bc_trace(pie.T,v,bc,dom,S),1e-9); %#ok<AGROW>
    % (d)
    Av = heatNd_apply(pie.A,v,S);
    Aex = r*uex(S);
    for i = 1:N,    Aex = Aex + exact_inv(S,Pw,cf,bc,dom,setdiff(1:N,i));   end
    R(end+1) = row(tag,'(d) A v = rTv + sum_i prod_{j~=i}T_j v',relerr(Av,Aex),tol); %#ok<AGROW>
    R(end+1) = row(tag,'(d0) A0 v = sum_i prod_{j~=i}T_j v', ...
                   relerr(heatNd_apply(pie.A0,v,S),Aex - r*uex(S)),tol);             %#ok<AGROW>
    % (e)
    phi1 = pie.exact.eigfun;
    Tp = heatNd_apply(pie.T,phi1,S,20);     Ap = heatNd_apply(pie.A,phi1,S,20);
    R(end+1) = row(tag,'(e) T phi1 = prod(-1/lam) phi1',relerr(Tp,pie.exact.Tfac*phi1(S)),tol); %#ok<AGROW>
    R(end+1) = row(tag,'(e) A phi1 = (r-lam1) T phi1',relerr(Ap,(r-pie.exact.lambda1)*Tp),tol); %#ok<AGROW>
    % (f) paper Ex. 3 (N = 2, DD x DN on [0,1]^2).
    if N==2 && isequal(bc,{'DD','DN'}) && isequal(dom,[0 1;0 1])
        R(end+1) = row(tag,'(f) T vs paper Ex. 3 kernel',relerr(Tv,paper_ex3(v,S)),tol); %#ok<AGROW>
    end
    % (g) stock.
    if dostock && N<=2
        [Ts,As] = stock_pie(N,r,bc,dom);
        R(end+1) = row(tag,'(g) T vs stock convert(PDE,pie)',relerr(Tv,stock_apply(Ts,v,S)),tol); %#ok<AGROW>
        R(end+1) = row(tag,'(g) A vs stock convert(PDE,pie)',relerr(Av,stock_apply(As,v,S)),tol); %#ok<AGROW>
    end
end
% Negative control: the evaluator must see a lower/upper swap in one direction.
pie = heatNd_pie(2,r);   B = pie.T.C{1,1};  prm = B.params;
prm = prm([1 3 2],:);                                       % swap cells in s1
Tsw = copvar({sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims)});
S = rand(6,2);  Pw = randi([0 3],5,2);  cf = randn(5,1);
e = relerr(heatNd_apply(Tsw,@(X) polyval_terms(X,Pw,cf),S), ...
           exact_inv(S,Pw,cf,pie.bc,pie.dom,1:2));
R(end+1) = row('N=2 control','lower/upper swapped in s1 (must FAIL)',e,-0.1);

fprintf('\n%-34s %-42s %10s %8s\n','case','check','rel err','result');
ok = true;
for k = 1:numel(R)
    if R(k).tol>0,  pass = R(k).err<=R(k).tol;  else,   pass = R(k).err>=-R(k).tol;  end
    ok = ok && pass;
    fprintf('%-34s %-42s %10.2e %8s\n',R(k).case,R(k).check,R(k).err,tern(pass,'pass','FAIL'));
end
assert(ok,'test_heatNd_pie: at least one check failed (see table).');
fprintf('\ntest_heatNd_pie: all %d checks pass.\n',numel(R));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = row(c,k,e,t),  s = struct('case',c,'check',k,'err',e,'tol',t);  end
function s = tern(b,x,y),   if b, s = x; else, s = y; end,  end
function e = relerr(a,b),   e = max(abs(a(:)-b(:)))/max(max(abs(b(:))),eps);   end

function y = polyval_terms(X,Pw,cf)
% sum_t cf(t) prod_i X(:,i).^Pw(t,i).
y = zeros(size(X,1),1);
for t = 1:numel(cf),    y = y + cf(t)*prod(X.^Pw(t,:),2);   end
end

function y = exact_inv(X,Pw,cf,bc,dom,dirs)
% sum_t cf(t) prod_i f_i(X_i): f_i = U_i (1-D inverse of d^2/ds^2 under
% bc{i}, applied to s^p) for i in DIRS, f_i = s^p otherwise.
y = zeros(size(X,1),1);
for t = 1:numel(cf)
    m = cf(t)*ones(size(X,1),1);
    for i = 1:size(X,2)
        p = Pw(t,i);    x = X(:,i);
        if ismember(i,dirs)
            [al,be] = affine_part(p,bc{i},dom(i,:));
            m = m.*(x.^(p+2)/((p+1)*(p+2)) + al*x + be);
        else
            m = m.*x.^p;
        end
    end
    y = y + m;
end
end

function [al,be] = affine_part(p,bc,ab)
% u = s^(p+2)/((p+1)(p+2)) + al s + be with u'' = s^p and the 1-D BC.
P  = @(s) s.^(p+2)/((p+1)*(p+2));   dP = @(s) s.^(p+1)/(p+1);
M = zeros(2);   rhs = zeros(2,1);
for e = 1:2
    c = ab(e);
    if bc(e)=='D',  M(e,:) = [c 1];     rhs(e) = -P(c);
    else,           M(e,:) = [1 0];     rhs(e) = -dP(c);
    end
end
z = M\rhs;  al = z(1);  be = z(2);
end

function [ph,d2ph] = domain_poly(bc,dom)
% phi(s) = prod_i q_i(s_i), q_i satisfying bc{i} exactly; d2ph = D^delta phi.
N = numel(bc);  q = cell(1,N);  q2 = cell(1,N);
for i = 1:N
    a = dom(i,1);   b = dom(i,2);
    switch bc{i}        % cubic q with q'' known; coefficients high to low
        case 'DD',  c = conv(conv([1 -a],[-1 b]),[1/3 1]);          % (s-a)(b-s)(1+s/3)
        case 'DN',  c = conv([1 -a],[-1 2*b-a]);                    % q(a)=0, q'(b)=0
        case 'ND',  c = conv([-1 b],[1 b-2*a]);                     % q'(a)=0, q(b)=0
    end
    q{i} = c;   q2{i} = polyder(polyder(c));
end
ph   = @(X) prodpoly(X,q);
d2ph = @(X) prodpoly(X,q2);
end

function y = prodpoly(X,q)
y = ones(size(X,1),1);
for i = 1:numel(q),     y = y.*polyval(q{i},X(:,i));    end
end

function e = bc_trace(T,v,bc,dom,S)
% Largest BC residual of T v on each face, relative to max |T v| inside.
% Dirichlet: value at the face; Neumann: 7-point one-sided derivative,
% exact for degree <= 6 (T v has degree <= 5 in each variable here).
N = numel(bc);  sc = max(abs(heatNd_apply(T,v,S)));  e = 0;
k = (0:6)';     V = (-k').^((0:6)');    cD = V\[0;1;0;0;0;0;0];     % f'(x0) ~ sum c f(x0-k h)/h
for i = 1:N
    for side = 1:2
        x0 = dom(i,side);   h = 0.02*(dom(i,2)-dom(i,1));   sg = 3-2*side;  % +1 at a (step inward)
        P = S(1:3,:);   P(:,i) = x0;
        if bc{i}(side)=='D'
            res = heatNd_apply(T,v,P);
        else
            F = zeros(3,7);
            for t = 1:7
                Pt = P;     Pt(:,i) = x0 + sg*(t-1)*h;
                F(:,t) = heatNd_apply(T,v,Pt);
            end
            res = (F*cD)/h;        % derivative along -sg; zero either way
        end
        e = max(e,max(abs(res))/sc);
    end
end
end

function y = paper_ex3(v,S)
% Paper Ex. 3, (T v)(x,y) = int_0^x int_0^y th(1-x)eta v + int_0^x int_y^1
% th(1-x) y v + int_x^1 int_0^y x(1-th) eta v + int_x^1 int_y^1 x(1-th) y v.
[xg,wg] = gl01(12);     y = zeros(size(S,1),1);
for k = 1:size(S,1)
    x = S(k,1);     yy = S(k,2);
    for I = 1:2
        if I==1,    t = x*xg;   wt = x*wg;      Kx = @(th) th*(1-x);
        else,       t = x+(1-x)*xg;  wt = (1-x)*wg;  Kx = @(th) x*(1-th);
        end
        for J = 1:2
            if J==1,    e = yy*xg;  we = yy*wg;     Ky = @(et) et;
            else,       e = yy+(1-yy)*xg;  we = (1-yy)*wg;  Ky = @(et) yy+0*et;
            end
            [TT,EE] = ndgrid(t,e);  [WT,WE] = ndgrid(wt,we);
            y(k) = y(k) + sum(WT(:).*WE(:).*Kx(TT(:)).*Ky(EE(:)).*v([TT(:),EE(:)]));
        end
    end
end
end

function [Ts,As] = stock_pie(N,r,bc,dom)
% STOCK PIE from pde_var / convert(PDE,'pie'), as in Sec. 7.1.
pvar t s1 s2
sv = [s1;s2];   X = sv(1:N);
u = pde_var('state',1,X,dom);
rhs = r*u;
for i = 1:N,    rhs = rhs + diff(u,X(i),2);     end
PDE = diff(u,t)==rhs;
for i = 1:N
    s = X(i);
    for side = 1:2
        if bc{i}(side)=='D',    PDE = [PDE; subs(u,s,dom(i,side))==0];              %#ok<AGROW>
        else,                   PDE = [PDE; subs(diff(u,s),s,dom(i,side))==0];      %#ok<AGROW>
        end
    end
end
evalc('PIE = convert(PDE,''pie'');');
Ts = PIE.T;     As = PIE.A;
end

function y = stock_apply(P,f,S)
% Stock opvar (1-D) / opvar2d (2-D) PI operator on L2 applied by
% quadrature, from each class's own definition (opvar.m / opvar2d.m):
% 1-D  (P v)(s) = R0 v(s) + int_a^s R1(s,t) v(t) dt + int_s^b R2(s,t) v(t) dt,
% 2-D  (P v)(x,y) = sum_{i,j} R22{i,j} applied with row i the x-direction
%      and column j the y-direction, 1 multiplier, 2 int_a^x, 3 int_x^b.
[xg,wg] = gl01(12);     y = zeros(size(S,1),1);
if isa(P,'opvar')
    cells = {P.R.R0, P.R.R1, P.R.R2};   nm = [elnames(P.var1); elnames(P.var2)];
    I = P.I;
    for k = 1:size(S,1)
        s = S(k,1);
        y(k) = pev(cells{1},nm,[s s])*f(s);
        for g = 2:3
            if g==2,    t = I(1)+(s-I(1))*xg;   w = (s-I(1))*wg;
            else,       t = s+(I(2)-s)*xg;      w = (I(2)-s)*wg;
            end
            y(k) = y(k) + sum(w.*pev(cells{g},nm,[s+0*t,t]).*f(t));
        end
    end
    return
end
nm = [elnames(P.var1); elnames(P.var2)];   I = P.I;     % [x y x_dum y_dum]
for k = 1:size(S,1)
    x = S(k,:);
    for i = 1:3
        for j = 1:3
            Rij = P.R22{i,j};
            if isempty(Rij) || (isa(Rij,'polynomial') && all(Rij.coefficient(:)==0)),  continue,  end
            [t1,w1] = rule1(i,x(1),I(1,:),xg,wg);
            [t2,w2] = rule1(j,x(2),I(2,:),xg,wg);
            [T1,T2] = ndgrid(t1,t2);    [W1,W2] = ndgrid(w1,w2);
            Pt = [T1(:),T2(:)];
            y(k) = y(k) + sum(W1(:).*W2(:).*pev(Rij,nm,[repmat(x,numel(T1),1),Pt]).*f(Pt));
        end
    end
end
end

function nm = elnames(pv)
% Names of the elements of a polynomial vector, in ELEMENT order (varname
% of the vector is the sorted union, which need not be element order).
nm = cell(numel(pv),1);
for i = 1:numel(pv),    e = pv(i);  nm{i} = e.varname{1};  end
end

function [t,w] = rule1(g,x,ab,xg,wg)
if g==1,        t = x;  w = 1;
elseif g==2,    t = ab(1)+(x-ab(1))*xg;     w = (x-ab(1))*wg;
else,           t = x+(ab(2)-x)*xg;         w = (ab(2)-x)*wg;
end
end

function v = pev(p,names,X)
% Scalar polynomial p evaluated at the rows of X, columns named NAMES
% (variables of p missing from NAMES are an error; unused NAMES ignored).
if isnumeric(p),    v = p*ones(size(X,1),1);    return,     end
[tf,loc] = ismember(p.varname,names);
if ~all(tf),    error('pev:var','Polynomial variable not supplied.'),   end
M = ones(size(X,1),size(p.degmat,1));
for j = 1:numel(loc)
    M = M.*X(:,loc(j)).^(full(p.degmat(:,j))');     % ' and .^ share precedence
end
v = M*full(p.coefficient(:,1));
end

function [x,w] = gl01(n)
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*V(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end

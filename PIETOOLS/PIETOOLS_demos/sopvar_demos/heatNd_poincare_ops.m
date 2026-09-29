function O = heatNd_poincare_ops(pie,doI0,O)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% O = HEATND_POINCARE_OPS(PIE,DOI0) the fixed operators of the PIE-form
% Poincare diagnostic (HEATND_POINCARE) for a HEATND_PIE instance, and, if
% DOI0, the SDP-free identity checks I0 - the first diagnostic of the T and
% A that the heat LPI (HEATND_LPI) is built from.
%
% With u = T v (v = D^delta u, HEATND_PIE; T = T_N...T_1, T_i the 1-D
% Green operator of d^2/ds_i^2 under bc{i}):
%   A_i  = d^2/ds_i^2 o T = prod_{j~=i} T_j       (identity in s_i),
%   A0   = sum_i A_i                               (A = r T + A0),
%   D_iT = d/ds_i o T = (dT_i/ds_i) prod_{j~=i} T_j,
% where dT_i/ds_i has the kernel dG/ds of the 1-D Green function on [a,b]:
%   DD  G = -(min(s,th)-a)(b-max(s,th))/L:  (th-a)/L for th <= s (lower
%       cell), (th-b)/L for th > s (upper);   [0,1]: th and th-1;
%   DN  G = a - min(s,th):                    0 lower, -1 upper;
%   ND  G = max(s,th) - b:                    1 lower,  0 upper.
% Identities (reference study, integration by parts; the boundary term
% [T w d_i T v] vanishes for DD, DN and ND, one factor being zero at each
% end):
%   T*A_i = A_i*T = -(D_iT)*(D_iT) <= 0,   T*A0 = -sum_i (D_iT)*(D_iT),
% so <v,-T*A0 v> = ||grad u||^2, <v,T*T v> = ||u||^2, <v,A0*A0 v> =
% ||Lap u||^2, and per direction T*A_i = T_i (x) T_ihat^2, A_i*A_i = I_i (x)
% T_ihat^2 (T_ihat = prod_{j~=i} T_j).
%
% I0 (DOI0 true; no SDP, seconds at N <= 3). Each object against its OWN
% definition (CLAUDE.md sec. 4); a failure implicates the T/A construction
% (BC-to-direction map, kernel coefficients, sorted-variable cell indexing,
% ZL/ZR) and the LPIs should not be run:
%  (k) kernels: every cell of T, A_i and D_iT, evaluated from the stored
%      coefficients (sopvar definition: K_g = ZL(s)' C_g ZR(th)) at random
%      (s,th), against the product of the textbook 1-D Green kernels above
%      (G, dG/ds, or the identity in direction i);
%  (c) coefficient identities through the class algebra (relative max
%      coefficient; zero iff the operator is zero, since the monomials are
%      independent per cell and the multiplier form is canonical, sopvar.m):
%      A - r T - A0, A0 - pie.A0, T*A_i - A_i*T, T*A_i + (D_iT)*(D_iT), and
%      T*A0 + sum_i (D_iT)*(D_iT);
%  (q) semantics by quadrature (HEATND_APPLY), v a random sum of monomials
%      of degree <= 2 per variable, u = T v the exact solution of D^delta u
%      = v with the BCs (1-D polynomial inverses, no kernel involved):
%      D_iT v = d_i u and A_i v at random points; <v,T*A_i v> = -||d_i u||^2,
%      <v,(D_iT)*(D_iT) v> = ||d_i u||^2, <v,T*T v> = ||u||^2 and
%      <v,A_i*A_i v> = ||A_i v||^2 on a tensor Gauss grid (all integrands
%      are polynomials of degree <= 8 per variable: exact at nq = 6).
% Tolerances: (k), (c) 1e-12 (the reference's; stock measured 2.2e-16 to
% 3.3e-16); (q) 1e-10.
%
% INPUT
% - pie:  struct from HEATND_PIE (any N, BC set, box and r);
% - doI0: (optional, default false) run I0 and store it in O.I0;
% - O:    (optional) an O to run I0 ON instead of building one (negative
%         controls in TEST_HEATND_POINCARE: a deliberately broken A_i/D_iT).
% OUTPUT struct O
% - N, vars, dom, bc, r, mu (1 x N: exact per-direction Poincare constants
%   pie.exact.lambda: pi^2/L^2 DD, pi^2/(4L^2) DN/ND), lambda1 (= sum mu);
% - T (pie.T), A0 (sum of the A_i built here), A{i}, DT{i} (1x1 'copvar');
% - dK (1 x N struct: ZL, ZR, C = {mult, lower, upper} of dT_i/ds_i);
% - I0 (DOI0): struct array (check, err, tol, ok), and O.I0ok.
% Compositions (T*T, (D_iT)*(D_iT), ...) are formed by HEATND_POINCARE for
% the selected target only.
%
% Cost: 2N+1 tensor-product 'copvar' of 3^N cells, 2^N nonzero, of size
% <= 2^N x 2^N; no decision variables. I0 adds O(N) compositions and a
% quadrature of ns * 3^N * nq^N kernel evaluations per application.
%
% Initial coding MMP, 09/27/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(doI0),   doI0 = false;   end
if nargin>=3 && ~isempty(O)                     % I0 on a given (possibly broken) O
    O.I0 = i0_checks(pie,O);    O.I0ok = all([O.I0.ok]);
    return
end
N = pie.N;  vars = pie.vars;    dom = pie.dom;  bc = pie.bc;
Id = struct('ZL',0,'ZR',0,'C',{{1,0,0}});      % identity in one direction
dK = struct('ZL',cell(1,N),'ZR',[],'C',[]);
for i = 1:N
    a = dom(i,1);   b = dom(i,2);   L = b-a;
    switch bc{i}            % dG/ds, coefficients of th^0, th^1 (ZR = [0;1])
        case 'DD',  dK(i).ZL = 0;   dK(i).ZR = [0;1];   dK(i).C = {[0 0], [-a 1]/L, [-b 1]/L};
        case 'DN',  dK(i).ZL = 0;   dK(i).ZR = 0;       dK(i).C = {0, 0, -1};
        case 'ND',  dK(i).ZL = 0;   dK(i).ZR = 0;       dK(i).C = {0, 1, 0};
        otherwise,  error('heatNd_poincare_ops:bc','BC ''%s'' is not DD, DN or ND.',bc{i})
    end
end
A = cell(1,N);  DT = cell(1,N);     A0 = [];
for i = 1:N
    fac = num2cell(pie.op1);    fac{i} = Id;    % prod_{j~=i} T_j (x) I_i
    A{i} = copvar({kronop(fac,vars,dom)});
    fac{i} = dK(i);                             % (dT_i/ds_i) (x) prod_{j~=i} T_j
    DT{i} = copvar({kronop(fac,vars,dom)});
    if isempty(A0),     A0 = A{i};  else,   A0 = A0 + A{i};     end
end
O = struct('N',N,'vars',{vars},'dom',dom,'bc',{bc},'r',pie.r,'mu',pie.exact.lambda, ...
           'lambda1',pie.exact.lambda1,'T',pie.T,'A0',A0,'A',{A},'DT',{DT},'dK',dK, ...
           'I0',[],'I0ok',[]);
if doI0
    O.I0 = i0_checks(pie,O);
    O.I0ok = all([O.I0.ok]);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = i0_checks(pie,O)
% The I0 checks of the header; one row per check.
N = O.N;    dom = O.dom;    T = O.T;
R = struct('check',{},'err',{},'tol',{},'ok',{});
tk = 1e-12;     tq = 1e-10;
% % % (k) kernels against the textbook Green products.
for g = 1:3^N
    gam = gamma_of(g,N);
    S = dom(:,1)' + rand(5,N).*(dom(:,2)-dom(:,1))';
    Th = dom(:,1)' + rand(5,N).*(dom(:,2)-dom(:,1))';
    ek(1) = kerr(T,g,S,Th,@(s,t,k,gk) green(O.bc{k},dom(k,:),s,t,gk));
    for i = 1:N
        ek(1+i)   = kerr(O.A{i},g,S,Th,@(s,t,k,gk) tern(k==i,double(gk==1)*ones(size(s)), ...
                                                  green(O.bc{k},dom(k,:),s,t,gk)));
        ek(1+N+i) = kerr(O.DT{i},g,S,Th,@(s,t,k,gk) tern(k==i,dgreen(O.bc{k},dom(k,:),s,t,gk), ...
                                                  green(O.bc{k},dom(k,:),s,t,gk)));
    end
    if g==1,    EK = ek;    else,   EK = max(EK,ek);    end
    clear ek
end
R(end+1) = row('(k) T kernels vs Green product',EK(1),tk);
for i = 1:N
    R(end+1) = row(sprintf('(k) A_%d kernels vs identity (x) Green',i),EK(1+i),tk);     %#ok<AGROW>
    R(end+1) = row(sprintf('(k) D_%dT kernels vs dG/ds (x) Green',i),EK(1+N+i),tk);     %#ok<AGROW>
end
% % % (c) coefficient identities through the class algebra.
R(end+1) = row('(c) A - rT - A0 (A0 = sum A_i)',cn(pie.A - O.r*T - O.A0)/cn(pie.A),tk);
R(end+1) = row('(c) A0 - pie.A0',cn(O.A0 - pie.A0)/cn(pie.A0),tk);
GG = [];
for i = 1:N
    TA = T'*O.A{i};     G = O.DT{i}'*O.DT{i};
    R(end+1) = row(sprintf('(c) T*A_%d - A_%d*T',i,i),cn(TA - O.A{i}'*T)/cn(TA),tk);     %#ok<AGROW>
    R(end+1) = row(sprintf('(c) T*A_%d + (D_%dT)*(D_%dT)',i,i,i),cn(TA + G)/cn(TA),tk);  %#ok<AGROW>
    if isempty(GG),     GG = G;     else,   GG = GG + G;    end
end
TA0 = T'*O.A0;
R(end+1) = row('(c) T*A0 + sum_i (D_iT)*(D_iT)',cn(TA0 + GG)/cn(TA0),tk);
% % % (q) semantics by quadrature.
nt = 4;     Pw = randi([0 2],nt,N);     cf = randn(nt,1);
v = @(X) polyterms(X,Pw,cf);
S = dom(:,1)' + rand(6,N).*(dom(:,2)-dom(:,1))';
nq = 6;     [Xq,wq] = qgrid(dom,nq);
uq = exactu(Xq,Pw,cf,O.bc,dom,1:N,[]);      vq = v(Xq);
R(end+1) = row('(q) <v,T*T v> = ||u||^2', ...
               rel(sum(wq.*vq.*heatNd_apply(T'*T,v,Xq,nq)),sum(wq.*uq.^2)),tq);
for i = 1:N
    di = exactu(S,Pw,cf,O.bc,dom,1:N,i);                    % d_i u at S
    R(end+1) = row(sprintf('(q) D_%dT v = d_%d u',i,i),rel(heatNd_apply(O.DT{i},v,S,nq),di),tq); %#ok<AGROW>
    ai = exactu(S,Pw,cf,O.bc,dom,setdiff(1:N,i),[]);        % prod_{j~=i} T_j v
    R(end+1) = row(sprintf('(q) A_%d v',i),rel(heatNd_apply(O.A{i},v,S,nq),ai),tq);     %#ok<AGROW>
    dq = exactu(Xq,Pw,cf,O.bc,dom,1:N,i);   gi = sum(wq.*dq.^2);
    R(end+1) = row(sprintf('(q) <v,T*A_%d v> = -||d_%d u||^2',i,i), ...
                   rel(sum(wq.*vq.*heatNd_apply(T'*O.A{i},v,Xq,nq)),-gi),tq);            %#ok<AGROW>
    R(end+1) = row(sprintf('(q) <v,(D_%dT)*(D_%dT) v> = ||d_%d u||^2',i,i,i), ...
                   rel(sum(wq.*vq.*heatNd_apply(O.DT{i}'*O.DT{i},v,Xq,nq)),gi),tq);      %#ok<AGROW>
    aq = exactu(Xq,Pw,cf,O.bc,dom,setdiff(1:N,i),[]);
    R(end+1) = row(sprintf('(q) <v,A_%d*A_%d v> = ||A_%d v||^2',i,i,i), ...
                   rel(sum(wq.*vq.*heatNd_apply(O.A{i}'*O.A{i},v,Xq,nq)),sum(wq.*aq.^2)),tq); %#ok<AGROW>
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = kronop(ops,vars,dom)
% 1x1 'sopvar' on L2[vars], the tensor product of 1-D operators ops{i}
% (fields ZL, ZR, C{1:3} = mult, lower, upper), op i acting in vars{i}.
% Same construction as HEATND_PIE's local KRONOP (params = kron of the 1-D
% coefficients, first variable slowest, as ZL(s) = kron_i ZL_i(s_i)); I0
% (c) checks the A_i built here against pie.A0.
N = numel(ops);
params = cell([3*ones(1,N),1]);
for g = 1:3^N
    gam = gamma_of(g,N);
    M = 1;
    for i = 1:N,    M = kron(M,ops{i}.C{gam(i)});   end
    params{g} = sparse(M);
end
ZL = cellfun(@(o) o.ZL(:),ops,'UniformOutput',false);
ZR = cellfun(@(o) o.ZR(:),ops,'UniformOutput',false);
P = sopvar(params,struct('in',{vars},'out',{vars}),ZL,ZR, ...
           struct('in',dom,'out',dom),[1 1]);
end


function gam = gamma_of(g,N)
% Multi-index of linear cell g, direction 1 fastest (sopvar.m).
c = cell(1,N);  [c{:}] = ind2sub([3*ones(1,N),1],g);    gam = [c{1:N}];
end


function e = kerr(X,g,S,Th,kfun)
% Max |K_g(s,th) - prod_k kfun(s_k,th_k,k,gamma_k)| over the rows of S, Th,
% relative to max(1, |analytic|); K_g from the stored coefficients of the
% 1x1 block (K_g = ZL(s)' C_g ZR(th), ZL over vars.out, ZR over vars.in).
B = X.C{1,1};   N = size(S,2);  gam = gamma_of(g,N);
[~,io] = ismember(B.vars.out,sort(B.vars.out));     [~,ii] = ismember(B.vars.in,sort(B.vars.in));  % S, Th: sorted
C = full(B.params{g});
ref = ones(size(S,1),1);
for k = 1:N,    ref = ref.*kfun(S(:,k),Th(:,k),k,gam(k));  end
if isempty(C) || ~any(C(:))
    Kc = zeros(size(ref));
else
    ZL = monvec(S(:,io),B.ZL);  ZR = monvec(Th(:,ii),B.ZR);     % columns in block order
    Kc = sum((ZL*C).*ZR,2);
end
e = max(abs(Kc-ref))/max(1,max(abs(ref)));
end


function k = green(bc,ab,s,t,gk)
% Textbook 1-D Green kernel of d^2/ds^2 under bc on [a,b], G = -(min(s,th)
% -a)(b-max(s,th))/L (DD), a-min (DN), max-b (ND), as the POLYNOMIAL of cell
% gk (1 multiplier: none; 2 lower, min = th; 3 upper, min = s), so random
% (s,th) off the cell's region test the stored polynomial itself.
a = ab(1);  b = ab(2);  L = b-a;
if gk==1,   k = zeros(size(s));     return,     end
if gk==2,   lo = t;     hi = s;     else,   lo = s;     hi = t;     end
switch bc
    case 'DD',  k = -(lo-a).*(b-hi)/L;
    case 'DN',  k = a - lo;
    case 'ND',  k = hi - b;
end
end


function k = dgreen(bc,ab,s,t,gk)
% d/ds of the cell polynomial of GREEN.
a = ab(1);  b = ab(2);  L = b-a;
if gk==1,   k = zeros(size(s));     return,     end
switch bc
    case 'DD',  if gk==2,   k = (t-a)/L;    else,   k = (t-b)/L;    end
    case 'DN',  k = -double(gk==3)*ones(size(s));
    case 'ND',  k = double(gk==2)*ones(size(s));
end
end


function Z = monvec(X,Zc)
% Rows: kron over variables of X(:,i).^Zc{i}, first variable slowest
% (as HEATND_APPLY).
n = size(X,1);  Z = ones(n,1);
for i = 1:numel(Zc)
    Zi = X(:,i).^reshape(Zc{i},1,[]);
    Z = reshape(reshape(Zi,n,[],1).*reshape(Z,n,1,[]),n,[]);
end
end


function n = cn(X)
% Largest |coefficient| of a 1x1 'copvar' over all cells (0 if zero).
B = X.C{1,1};   n = 0;
if isempty(B),  return,     end
for g = 1:numel(B.params)
    c = B.params{g};
    if ~isempty(c) && nnz(c),   n = max(n,full(max(abs(c(:)))));   end
end
end


function y = polyterms(X,Pw,cf)
% sum_t cf(t) prod_i X(:,i).^Pw(t,i).
y = zeros(size(X,1),1);
for t = 1:numel(cf),    y = y + cf(t)*prod(X.^Pw(t,:),2);   end
end


function y = exactu(X,Pw,cf,bc,dom,dirs,dd)
% sum_t cf(t) prod_i f_i(X_i) with f_i = U_i(s^p) (the 1-D solution of u''
% = s^p under bc{i}) for i in DIRS, s^p otherwise; in direction DD (if not
% empty, and in DIRS) the derivative U_i'. U = s^(p+2)/((p+1)(p+2)) + al s
% + be, the affine part fixed by the two BCs (no Green kernel involved).
y = zeros(size(X,1),1);
for t = 1:numel(cf)
    m = cf(t)*ones(size(X,1),1);
    for i = 1:size(X,2)
        p = Pw(t,i);    x = X(:,i);
        if ismember(i,dirs)
            [al,be] = affine_part(p,bc{i},dom(i,:));
            if isequal(i,dd),   m = m.*(x.^(p+1)/(p+1) + al);
            else,               m = m.*(x.^(p+2)/((p+1)*(p+2)) + al*x + be);
            end
        else
            m = m.*x.^p;
        end
    end
    y = y + m;
end
end


function [al,be] = affine_part(p,bc,ab)
% u = s^(p+2)/((p+1)(p+2)) + al s + be with u'' = s^p and the 1-D BC.
P = @(s) s.^(p+2)/((p+1)*(p+2));    dP = @(s) s.^(p+1)/(p+1);
M = zeros(2);   rhs = zeros(2,1);
for e = 1:2
    c = ab(e);
    if bc(e)=='D',  M(e,:) = [c 1];     rhs(e) = -P(c);
    else,           M(e,:) = [1 0];     rhs(e) = -dP(c);
    end
end
z = M\rhs;  al = z(1);  be = z(2);
end


function [X,w] = qgrid(dom,n)
% Tensor Gauss-Legendre grid on the box, n nodes per direction.
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     wg = 2*V(1,ix)'.^2;     x = (x+1)/2;    wg = wg/2;
N = size(dom,1);    X = zeros(1,0);     w = 1;
for d = 1:N
    L = dom(d,2)-dom(d,1);  p = dom(d,1)+L*x;   q = L*wg;   m = size(X,1);
    X = [repmat(X,n,1), repelem(p,m,1)];    w = repmat(w,n,1).*repelem(q,m,1);
end
end


function s = row(c,e,t),    s = struct('check',c,'err',e,'tol',t,'ok',e<=t);  end
function e = rel(a,b),      e = max(abs(a(:)-b(:)))/max(max(abs(b(:))),eps);   end
function x = tern(c,a,b),   if c, x = a; else, x = b; end,  end

function out = mincap_dual(P,D,d,M,mult)
% OUT = MINCAP_DUAL(P,D,d,M,MULT) the minimum-cap dual of the proof program
% (minimum_cap_primal_dual.md (5), algebraic_hierarchy_duality.md (6)) at
% one finite level, in its POINTWISE form:
%
%   beta = sup tr(Pq*T)  s.t.  S_k >= F(s_k) T F(s_k)',  S_k >= 0,
%                              sum_k w_k tr(S_k) <= 1,
%
% i.e. the integral of the positive-part trace of the covariance field
% Z_T(s) = F(s) T F(s)' is replaced by a 2-point Gauss rule on M uniform
% cells (2M nodes s_k, weights w_k). F(s) is the canonical lift of degree D
% on [a,b] (pure: running moments z_j and full moments m_j, j = 0..D;
% mixed, MULT = 1: the input itself as a further component) applied to an
% orthonormal (shifted Legendre) input basis of degree <= d, q = d+1.
% beta approximates the level value of the hierarchy, which increases to
% the minimum cap m_D(P) of a bounded positive weight at lift degree D as
% d and the mesh grow, and is +inf when no bounded weight exists.
%
% P: struct with a, b, R0c (coefficients of the multiplier in s, index
%    alpha+1 <-> s^alpha; [] for none), R1c (matrix, (alpha+1,beta+1) <->
%    s^alpha theta^beta, the kernel on theta < s), R2c (same, theta > s).
% OUT: beta, status, T (in the Legendre basis), Pq, q, r, M, nodes.
% Scalar inputs (n = 1).
%
% Initial coding MMP, 10/08/2026 (cell form); pointwise form and Legendre
% basis the same day: the cell Gramians B_E are numerically rank deficient
% on fine meshes and the solver stalls.
a = P.a;    b = P.b;    L = b-a;    q = d+1;
nz = D+1;   r = (mult~=0) + 2*nz;
% shifted Legendre basis, orthonormal on [a,b]: coefficient rows in t^0..t^d
C = zeros(q,q);                             % C(j+1,:) = coefficients of p_j in powers of t
u = [-(a+b)/L, 2/L];                        % u = (2t-a-b)/L as a polynomial in t
Lg = cell(1,q);     Lg{1} = 1;
if q>1,     Lg{2} = u;  end
for j = 2:d                                 % (j+1) P_{j+1} = (2j+1) u P_j - j P_{j-1}
    Lg{j+1} = ((2*j-1)*conv(u,Lg{j}));
    tmp = zeros(size(Lg{j+1}));     tmp(1:numel(Lg{j-1})) = (j-1)*Lg{j-1};
    Lg{j+1} = (Lg{j+1} - tmp)/j;
end
for j = 0:d
    C(j+1,1:numel(Lg{j+1})) = Lg{j+1}*sqrt((2*j+1)/L);
end
% monomial pairings <t^i, P t^j> and Pq = C Mm C'
Mm = zeros(q);
for i = 0:d
    for j = 0:d,    Mm(i+1,j+1) = pair(P,i,j);  end
end
Pq = C*Mm*C';   Pq = (Pq+Pq')/2;
% nodes
edges = linspace(a,b,M+1);  g = 1/sqrt(3);
nodes = zeros(1,2*M);   wts = zeros(1,2*M);
for E = 1:M
    c0 = (edges(E)+edges(E+1))/2;   h = edges(E+1)-edges(E);
    nodes(2*E-1:2*E) = c0 + h/2*[-g,g];     wts(2*E-1:2*E) = h/2;
end
K = numel(nodes);
% lift at the nodes: F_k (r x q), entries from the monomial lift (s^e - a^e)/e etc.
Fk = cell(1,K);
for k = 1:K
    s = nodes(k);
    Fm = zeros(r,q);                        % monomial inputs t^0..t^d
    row = 0;
    if mult,    Fm(1,:) = s.^(0:d);     row = 1;    end
    for kk = 0:D
        for j = 0:d
            e = kk+j+1;
            Fm(row+kk+1,j+1) = (s^e - a^e)/e;
            Fm(row+nz+kk+1,j+1) = (b^e - a^e)/e;
        end
    end
    Fk{k} = Fm*C';                          % Legendre inputs
end
% SeDuMi primal form: x = [tvec (free, q(q+1)/2); tau (>=0); S_1 (r^2); E_1 (r^2); ...]
nt = q*(q+1)/2;
[iu,ju] = find(triu(ones(q)));
nx = nt + 1 + 2*K*r^2;
[ir,jr] = find(triu(ones(r)));  nr = numel(ir);
neq = K*nr + 1;
Ai = [];    Aj = [];    Av = [];
rowc = 0;
for k = 1:K
    offS = nt+1+(k-1)*2*r^2;    offE = offS + r^2;    F = Fk{k};
    for e = 1:nr
        i = ir(e);  j = jr(e);  rowc = rowc+1;
        % E_k(i,j) - S_k(i,j) + (F T F')(i,j) = 0 on the symmetric parts
        if i==j
            Ai(end+1) = rowc;   Aj(end+1) = offE+(j-1)*r+i;    Av(end+1) = 1;      %#ok<AGROW>
            Ai(end+1) = rowc;   Aj(end+1) = offS+(j-1)*r+i;    Av(end+1) = -1;     %#ok<AGROW>
        else
            Ai(end+1:end+2) = rowc;   Aj(end+1:end+2) = offE+[(j-1)*r+i,(i-1)*r+j];  Av(end+1:end+2) = 0.5;  %#ok<AGROW>
            Ai(end+1:end+2) = rowc;   Aj(end+1:end+2) = offS+[(j-1)*r+i,(i-1)*r+j];  Av(end+1:end+2) = -0.5; %#ok<AGROW>
        end
        for t = 1:nt                        % (F T F')(i,j) = sum_pq F(i,p) T(p,q) F(j,q)
            p = iu(t);  qq = ju(t);
            if p==qq,   cf = F(i,p)*F(j,p);
            else,       cf = F(i,p)*F(j,qq) + F(i,qq)*F(j,p);
            end
            if cf~=0,   Ai(end+1) = rowc;   Aj(end+1) = t;  Av(end+1) = cf;    end   %#ok<AGROW>
        end
    end
end
rowc = rowc+1;                              % quadrature trace row
for k = 1:K
    offS = nt+1+(k-1)*2*r^2;
    Ai(end+1:end+r) = rowc;     Aj(end+1:end+r) = offS+((1:r)-1)*r+(1:r);   Av(end+1:end+r) = wts(k);   %#ok<AGROW>
end
Ai(end+1) = rowc;   Aj(end+1) = nt+1;   Av(end+1) = 1;
A = sparse(Ai,Aj,Av,neq,nx);
bb = zeros(neq,1);  bb(end) = 1;
pP = max(norm(Pq),eps);
c = zeros(nx,1);
for t = 1:nt
    i = iu(t);  j = ju(t);
    c(t) = -(2-(i==j))*Pq(i,j)/pP;
end
Kc = struct('f',nt,'l',1,'s',repmat(r,1,2*K));
T = zeros(q);   status = '';
if isempty(which('mosekopt'))
    pars = struct('fid',0,'eps',1e-9);
    [x,~,info] = sedumi(A,bb,c,Kc,pars);
    beta = -c'*x*pP;
    if info.dinf,   beta = Inf;     status = 'sedumi dinf';    else,   status = 'sedumi';   end
    for t = 1:nt,   T(iu(t),ju(t)) = x(t);  T(ju(t),iu(t)) = x(t);   end
else
    prob = Sedumi2Mosek(A,bb,c,Kc);
    [~,res] = mosekopt('minimize echo(0)',prob);
    info = res.sol.itr;
    status = [info.prosta ' / ' info.solsta];
    if contains(info.prosta,'DUAL_INFEASIBLE'),     beta = Inf;
    elseif contains(info.solsta,'OPTIMAL'),         beta = -info.pobjval*pP;
    else,                                           beta = NaN;
    end
    if isfield(info,'xx') && numel(info.xx)>=nt
        for t = 1:nt,   T(iu(t),ju(t)) = info.xx(t);  T(ju(t),iu(t)) = info.xx(t);   end
    end
end
out = struct('beta',beta,'status',status,'T',T,'Pq',Pq,'q',q,'r',r,'M',M,'nodes',nodes,'D',D,'d',d,'mult',mult);
end


function v = pair(P,i,j)
% <t^i, P t^j> on [a,b] from the coefficient tables
a = P.a;    b = P.b;
I = @(n) (b^(n+1) - a^(n+1))/(n+1);         % int_a^b s^n ds
v = 0;
if ~isempty(P.R0c)
    for al = 0:numel(P.R0c)-1
        v = v + P.R0c(al+1)*I(i+al+j);
    end
end
if ~isempty(P.R1c)                          % theta < s
    [nA,nB] = size(P.R1c);
    for al = 0:nA-1
        for be = 0:nB-1
            if P.R1c(al+1,be+1)==0, continue, end
            v = v + P.R1c(al+1,be+1)*(I(i+al+j+be+1) - a^(j+be+1)*I(i+al))/(j+be+1);
        end
    end
end
if ~isempty(P.R2c)                          % theta > s
    [nA,nB] = size(P.R2c);
    for al = 0:nA-1
        for be = 0:nB-1
            if P.R2c(al+1,be+1)==0, continue, end
            v = v + P.R2c(al+1,be+1)*(b^(j+be+1)*I(i+al) - I(i+al+j+be+1))/(j+be+1);
        end
    end
end
end

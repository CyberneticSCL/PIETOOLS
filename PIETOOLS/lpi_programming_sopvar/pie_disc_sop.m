function [M,info] = pie_disc_sop(Pop,N,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [M,INFO] = PIE_DISC_SOP(POP,N,DOM) the Galerkin matrix of a 1-D 4-PI
% operator (an 'opvar', or a double) in the orthonormal basis
%
%     R^{m0} (+) span{ p_j(s) e_k : j = 0..N, k = 1..m1 }   (p_j shifted Legendre on DOM),
%
% component-outer, coefficient-inner for the L2 part. The entries are the
% pairings <basis_i, Pop basis_j>, computed by Gauss-Legendre quadrature
% with the basis evaluated by its three-term recurrence; the multiplier R0
% and the kernels R1 (theta < s), R2 (theta > s) are polynomial, so with
% N + dmax + 3 nodes (dmax the largest total degree of an entry) every
% quadrature is exact. Because the basis is orthonormal the matrix of a
% composition is the product of the matrices up to the truncation, and the
% L2 norms of inputs, outputs and states are the Euclidean norms of their
% coefficient vectors: no mass matrix is needed.
%
% INPUT
% - Pop:  'opvar' (fields P, Q1, Q2, R.R0, R.R1, R.R2, var1, var2, I), or a
%         double matrix (a finite-dimensional operator);
% - N:    polynomial degree of the basis (N+1 functions per component);
% - dom:  [a b] (default Pop.I).
%
% OUTPUT
% - M:    (m0 + m1 (N+1)) x (n0 + n1 (N+1)) matrix;
% - info: evalbasis (handle: values of p_0..p_N at a column of points, one
%         row per point), dims, a, b, q = N+1, nq (quadrature nodes).
%
% Initial coding MMP, 10/08/2026 (PIESIM's matrices are built on a rescaled
% PIE inside its own pipeline and are not the plain coefficient maps; this
% routine is the executives' own discretization for the numerical witness).
% MMP, 10/08/2026: Pairings by quadrature on the recurrence-evaluated basis
%                  instead of closed-form monomial moments of the basis'
%                  monomial coefficients: those coefficients grow as 4^N
%                  and the moment sums lost every digit beyond N = 12
%                  (spurious eigenvalues at N = 16, singular T at N >= 24
%                  on the io1 plant), while the executives' default is 24.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(Pop,'double')
    M = Pop;    info = struct('evalbasis',[],'dims',[size(Pop,1) 0; 0 0; size(Pop,2) 0],'a',NaN,'b',NaN,'q',NaN,'nq',NaN);
    return
end
if ~isa(Pop,'opvar'),   error('pie_disc_sop:class','POP should be an opvar or a double.');  end
if nargin<3 || isempty(dom),    dom = Pop.I;    end
a = dom(1);     b = dom(2);     q = N+1;
d = Pop.dim;    m0 = d(1,1);    m1 = d(2,1);    n0 = d(1,2);    n1 = d(2,2);
sname = Pop.var1.varname{1};    tname = Pop.var2.varname{1};
% coefficient tables of every entry, and the largest total degree
Q1c = cell(m0,n1);  Q2c = cell(m1,n0);  R0c = cell(m1,n1);  R1c = cell(m1,n1);  R2c = cell(m1,n1);
dmax = 0;
for i = 1:m0,   for k = 1:n1,   Q1c{i,k} = pcoef1(entry(Pop.Q1,i,k),sname);    dmax = max(dmax,numel(Q1c{i,k})-1);    end,    end
for k = 1:m1,   for i = 1:n0,   Q2c{k,i} = pcoef1(entry(Pop.Q2,k,i),sname);    dmax = max(dmax,numel(Q2c{k,i})-1);    end,    end
for k = 1:m1
    for l = 1:n1
        R0c{k,l} = pcoef1(entry(Pop.R.R0,k,l),sname);           dmax = max(dmax,numel(R0c{k,l})-1);
        R1c{k,l} = pcoef2(entry(Pop.R.R1,k,l),sname,tname);     dmax = max(dmax,sum(size(R1c{k,l}))-2);
        R2c{k,l} = pcoef2(entry(Pop.R.R2,k,l),sname,tname);     dmax = max(dmax,sum(size(R2c{k,l}))-2);
    end
end
nq = N + max(dmax,0) + 3;                 % exact for integrands of degree <= 2nq-1 >= 2N+2dmax+5
[sq,wq] = gauss_legendre(nq,a,b);
Pq = legvals(sq,N,a,b);                   % nq x q, p_j(s_m)
M = zeros(m0+m1*q, n0+n1*q);
% finite -> finite
if m0>0 && n0>0,    M(1:m0,1:n0) = double(Pop.P);    end
% L2 -> finite: z = int Q1(s) x(s) ds
for i = 1:m0
    for k = 1:n1
        c = Q1c{i,k};   if isempty(c),  continue,   end
        f = monoval(sq,c);                                  % Q1_ik(s_m)
        M(i,n0+(k-1)*q+(1:q)) = (wq(:).*f(:))'*Pq;
    end
end
% finite -> L2: x(s) = Q2(s) u
for k = 1:m1
    for i = 1:n0
        c = Q2c{k,i};   if isempty(c),  continue,   end
        f = monoval(sq,c);
        M(m0+(k-1)*q+(1:q),i) = Pq'*(wq(:).*f(:));
    end
end
% L2 -> L2: multiplier and kernels
if m1>0 && n1>0
    % inner quadrature per outer node, on [a,s_m] (R1) and [s_m,b] (R2)
    [xr,wr] = gauss_legendre(nq,-1,1);                      % reference rule
    for k = 1:m1
        for l = 1:n1
            if isempty(R0c{k,l}) && isempty(R1c{k,l}) && isempty(R2c{k,l}),   continue,   end
            Mm = zeros(q);
            if ~isempty(R0c{k,l})
                f = monoval(sq,R0c{k,l});
                Mm = Mm + Pq'*((wq(:).*f(:)).*Pq);
            end
            if ~isempty(R1c{k,l}) || ~isempty(R2c{k,l})
                F = zeros(nq,q);                            % F(m,j) = int_a^b K(s_m,th) p_j(th) dth
                for m = 1:nq
                    s = sq(m);
                    if ~isempty(R1c{k,l}) && s>a
                        th = (s-a)/2*xr + (s+a)/2;      w = (s-a)/2*wr;
                        F(m,:) = F(m,:) + (w(:).*kerval(s,th,R1c{k,l}))'*legvals(th,N,a,b);
                    end
                    if ~isempty(R2c{k,l}) && s<b
                        th = (b-s)/2*xr + (b+s)/2;      w = (b-s)/2*wr;
                        F(m,:) = F(m,:) + (w(:).*kerval(s,th,R2c{k,l}))'*legvals(th,N,a,b);
                    end
                end
                Mm = Mm + Pq'*(wq(:).*F);
            end
            M(m0+(k-1)*q+(1:q), n0+(l-1)*q+(1:q)) = Mm;
        end
    end
end
info = struct('evalbasis',@(x) legvals(x,N,a,b),'dims',d,'a',a,'b',b,'q',q,'nq',nq);
end


function V = legvals(x,N,a,b)
% values of the orthonormal shifted Legendre polynomials p_0..p_N on [a,b] at
% the points x (column), one row per point, by the three-term recurrence
x = x(:);   L = b-a;    u = (2*x-a-b)/L;
V = zeros(numel(x),N+1);
V(:,1) = 1;
if N>=1,    V(:,2) = u;     end
for j = 1:N-1
    V(:,j+2) = ((2*j+1)*u.*V(:,j+1) - j*V(:,j))/(j+1);
end
V = V.*sqrt((2*(0:N)+1)/L);
end


function [x,w] = gauss_legendre(n,a,b)
% Golub-Welsch nodes and weights on [a,b]
beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
J = diag(beta,1)+diag(beta,-1);
[V,Dg] = eig(J);
[x,i] = sort(diag(Dg));     w = 2*V(1,i).^2;
x = (b-a)/2*x + (a+b)/2;    w = (b-a)/2*w(:);
end


function f = monoval(x,c)
% sum_al c(al+1) x^al at the points x
x = x(:);   f = zeros(size(x));
for al = 0:numel(c)-1
    if c(al+1)~=0,  f = f + c(al+1)*x.^al;  end
end
end


function f = kerval(s,th,Cm)
% sum_{al,be} Cm(al+1,be+1) s^al th^be at the points th (column), s scalar
[nA,nB] = size(Cm);
sv = s.^(0:nA-1);                           % 1 x nA
tv = th(:).^(0:nB-1);                       % nq x nB
f = tv*(Cm'*sv');                           % nq x 1
end


function e = entry(X,i,k)
% entry (i,k) of a (possibly empty / double / polynomial) matrix
if isempty(X),  e = [];  return,   end
if isa(X,'double'),     e = X(i,k);     return,     end
e = X(i,k);
end


function c = pcoef1(p,s)
% coefficient vector of a polynomial in s (index alpha+1 <-> s^alpha); [] if zero
if isempty(p),  c = [];  return,   end
if isa(p,'double')
    if all(p(:)==0),    c = [];     else,   c = p;  end
    return
end
p = polynomial(p);
if isempty(p.coefficient) || ~any(p.coefficient(:)),  c = [];  return,     end
is = find(strcmp(p.varname,s));
dg = zeros(size(p.degmat,1),1);     if ~isempty(is),    dg = p.degmat(:,is);    end
others = setdiff(1:numel(p.varname),is);
if ~isempty(others) && any(any(p.degmat(:,others)>0))
    error('pie_disc_sop:vars','A multiplier or finite-coupling entry depends on a variable other than %s.',s);
end
c = zeros(max(dg)+1,1);
for k = 1:numel(dg),    c(dg(k)+1) = c(dg(k)+1) + p.coefficient(k);    end
end


function Cm = pcoef2(p,s,th)
% coefficient matrix of a polynomial in (s,theta): (alpha+1,beta+1) <-> s^alpha theta^beta; [] if zero
if isempty(p),  Cm = [];  return,   end
if isa(p,'double')
    if all(p(:)==0),    Cm = [];    else,   Cm = p;     end
    return
end
p = polynomial(p);
if isempty(p.coefficient) || ~any(p.coefficient(:)),  Cm = [];  return,    end
is = find(strcmp(p.varname,s));  it = find(strcmp(p.varname,th));
nterm = size(p.degmat,1);
da = zeros(nterm,1);    db = zeros(nterm,1);
if ~isempty(is),    da = p.degmat(:,is);    end
if ~isempty(it),    db = p.degmat(:,it);    end
Cm = zeros(max(da)+1,max(db)+1);
for k = 1:nterm,    Cm(da(k)+1,db(k)+1) = Cm(da(k)+1,db(k)+1) + p.coefficient(k);    end
end

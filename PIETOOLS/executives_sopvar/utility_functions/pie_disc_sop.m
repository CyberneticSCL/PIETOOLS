function [M,info] = pie_disc_sop(Pop,N,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [M,INFO] = PIE_DISC_SOP(POP,N,DOM) returns the Galerkin matrix M of the
% 1-D PI operator POP: R^n0 x L2^n1[a,b] -> R^m0 x L2^m1[a,b],
%
%   (Pop x)_0    = P x_0 + int_a^b Q1(s) x_1(s) ds,
%   (Pop x)_1(s) = Q2(s) x_0 + R0(s) x_1(s) + int_a^s R1(s,th) x_1(th) dth
%                  + int_s^b R2(s,th) x_1(th) dth,
%
% in the orthonormal basis of each space made of the unit vectors of its
% finite part and the functions p_j(s) e_k, j = 0..N, of its L2 part, p_j
% the Legendre polynomial of degree j normalized on [a,b] and e_k the k-th
% unit vector, ordered finite part first, then k outer and j inner: M(r,c)
% is the pairing <basis_r, Pop basis_c>. The pairings are computed by
% Gauss-Legendre quadrature on N + dmax + 3 nodes, dmax a bound on the
% degree of the parameters of POP, with the basis evaluated by its
% three-term recurrence; for polynomial parameters the quadrature is
% exact. The Euclidean norm of a coefficient vector is the norm of the
% function it represents.
%
% INPUT
% - Pop:  'opvar' object with Pop.dim = [m0 n0; m1 n1], parameters P
%         (double), Q1, Q2, R.R0 (polynomials in Pop.var1, or double) and
%         R.R1, R.R2 (polynomials in Pop.var1 and Pop.var2, or double); or
%         a 'double' matrix, returned as M;
% - N:    double, nonnegative integer, the largest degree of the basis;
% - dom:  (optional) 1x2 double [a b], the domain of the basis (default
%         Pop.I).
%
% OUTPUT
% - M:    (m0 + m1*(N+1)) x (n0 + n1*(N+1)) double;
% - info: 'struct' with fields
%           evalbasis  'function_handle', V = evalbasis(x): the
%                      numel(x) x (N+1) values of p_0, ..., p_N at the
%                      points x;
%           dims       Pop.dim;
%           a, b       double, the domain;
%           q          double, N+1;
%           nq         double, the number of quadrature nodes;
%         for a 'double' POP, evalbasis is [], a, b, q and nq are NaN, and
%         dims is [size(Pop,1) 0; 0 0; size(Pop,2) 0].
%
% NOTES
% POP of another class is an error, as is an entry of Q1, Q2 or R.R0 with a
% positive degree in a variable other than Pop.var1.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie_disc_sop
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), the licence block added, step comments
%                  in the body and headers of the local functions; code
%                  unchanged. The deleted header text was a remark on
%                  compositions and on the mass matrix.

% A matrix is its own Galerkin matrix
if isa(Pop,'double')
    M = Pop;    info = struct('evalbasis',[],'dims',[size(Pop,1) 0; 0 0; size(Pop,2) 0],'a',NaN,'b',NaN,'q',NaN,'nq',NaN);
    return
end
% Otherwise an opvar: domain, basis size, dimensions and variable names
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
% Quadrature nodes and weights on [a,b], and the basis at the nodes
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
            % Multiplier: sum over the nodes of wq * p_i R0 p_j
            if ~isempty(R0c{k,l})
                f = monoval(sq,R0c{k,l});
                Mm = Mm + Pq'*((wq(:).*f(:)).*Pq);
            end
            % Kernels: at each node s_m the inner integrals of R1 and R2
            % against p_j, then the outer quadrature against p_i
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
% Basis evaluator and sizes
info = struct('evalbasis',@(x) legvals(x,N,a,b),'dims',d,'a',a,'b',b,'q',q,'nq',nq);
end


function V = legvals(x,N,a,b)
% V = LEGVALS(X,N,A,B) returns the values of the orthonormal Legendre
% polynomials p_0, ..., p_N on [A,B] at the points X, by the three-term
% recurrence.
% INPUT  x: double, the points; N: double, the largest degree;
%        a, b: double, the interval.
% OUTPUT V: numel(x) x (N+1) double, V(r,j+1) = p_j(x(r)).
% Points mapped to u in [-1,1]; Legendre P_j(u) by
% (j+1) P_{j+1} = (2j+1) u P_j - j P_{j-1}; p_j = sqrt((2j+1)/(b-a)) P_j
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
% [X,W] = GAUSS_LEGENDRE(N,A,B) returns the N-point Gauss-Legendre nodes and
% weights on [A,B], by the Golub-Welsch method.
% INPUT  n: double, the number of nodes; a, b: double, the interval.
% OUTPUT x: n x 1 double, the nodes in increasing order;
%        w: n x 1 double, the weights.
% Nodes: the eigenvalues of the Jacobi matrix of the Legendre recurrence;
% weights: twice the squared first entries of its eigenvectors; then the
% map from [-1,1] to [a,b]
beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
J = diag(beta,1)+diag(beta,-1);
[V,Dg] = eig(J);
[x,i] = sort(diag(Dg));     w = 2*V(1,i).^2;
x = (b-a)/2*x + (a+b)/2;    w = (b-a)/2*w(:);
end


function f = monoval(x,c)
% F = MONOVAL(X,C) returns sum_al c(al+1)*x.^al at the points X.
% INPUT  x: double, the points; c: double vector, the coefficients.
% OUTPUT f: numel(x) x 1 double.
x = x(:);   f = zeros(size(x));
for al = 0:numel(c)-1
    if c(al+1)~=0,  f = f + c(al+1)*x.^al;  end
end
end


function f = kerval(s,th,Cm)
% F = KERVAL(S,TH,CM) returns sum_{al,be} Cm(al+1,be+1)*s^al*th.^be at the
% points TH, for a scalar S.
% INPUT  s: double scalar; th: double, the points; Cm: double matrix, the
%        coefficients.
% OUTPUT f: numel(th) x 1 double.
[nA,nB] = size(Cm);
sv = s.^(0:nA-1);                           % 1 x nA
tv = th(:).^(0:nB-1);                       % nq x nB
f = tv*(Cm'*sv');                           % nq x 1
end


function e = entry(X,i,k)
% E = ENTRY(X,I,K) returns the entry (I,K) of the matrix X, or [] when X is
% empty.
% INPUT  X: double or 'polynomial' matrix, possibly empty; i, k: double.
% OUTPUT e: double or 'polynomial' scalar, or [].
if isempty(X),  e = [];  return,   end
if isa(X,'double'),     e = X(i,k);     return,     end
e = X(i,k);
end


function c = pcoef1(p,s)
% C = PCOEF1(P,S) returns the coefficient vector of the polynomial P in the
% variable S, c(al+1) the coefficient of s^al, or [] when P is zero; a
% positive degree of P in another variable is an error.
% INPUT  p: double or 'polynomial' scalar, possibly empty; s: char, the
%           variable name.
% OUTPUT c: double column, the double P itself, or [].
% A double: itself, or [] when zero
if isempty(p),  c = [];  return,   end
if isa(p,'double')
    if all(p(:)==0),    c = [];     else,   c = p;  end
    return
end
% A polynomial: [] when zero; the exponents of S, no other variable allowed
p = polynomial(p);
if isempty(p.coefficient) || ~any(p.coefficient(:)),  c = [];  return,     end
is = find(strcmp(p.varname,s));
dg = zeros(size(p.degmat,1),1);     if ~isempty(is),    dg = p.degmat(:,is);    end
others = setdiff(1:numel(p.varname),is);
if ~isempty(others) && any(any(p.degmat(:,others)>0))
    error('pie_disc_sop:vars','A multiplier or finite-coupling entry depends on a variable other than %s.',s);
end
% Sum the coefficients of each exponent
c = zeros(max(dg)+1,1);
for k = 1:numel(dg),    c(dg(k)+1) = c(dg(k)+1) + p.coefficient(k);    end
end


function Cm = pcoef2(p,s,th)
% CM = PCOEF2(P,S,TH) returns the coefficient matrix of the polynomial P in
% the variables S and TH, Cm(al+1,be+1) the coefficient of s^al*th^be, or
% [] when P is zero; the exponents of other variables are not read.
% INPUT  p:  double or 'polynomial' scalar, possibly empty; s, th: char,
%            the variable names.
% OUTPUT Cm: double matrix, the double P itself, or [].
% A double: itself, or [] when zero
if isempty(p),  Cm = [];  return,   end
if isa(p,'double')
    if all(p(:)==0),    Cm = [];    else,   Cm = p;     end
    return
end
% A polynomial: [] when zero; the exponents of S and TH of each term
p = polynomial(p);
if isempty(p.coefficient) || ~any(p.coefficient(:)),  Cm = [];  return,    end
is = find(strcmp(p.varname,s));  it = find(strcmp(p.varname,th));
nterm = size(p.degmat,1);
da = zeros(nterm,1);    db = zeros(nterm,1);
if ~isempty(is),    da = p.degmat(:,is);    end
if ~isempty(it),    db = p.degmat(:,it);    end
% Sum the coefficients of each pair of exponents
Cm = zeros(max(da)+1,max(db)+1);
for k = 1:nterm,    Cm(da(k)+1,db(k)+1) = Cm(da(k)+1,db(k)+1) + p.coefficient(k);    end
end

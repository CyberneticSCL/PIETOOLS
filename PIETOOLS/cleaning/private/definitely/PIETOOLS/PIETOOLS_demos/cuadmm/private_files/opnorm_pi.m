function [M,info] = opnorm_pi(Op,D,nq)
% opnorm_pi -- Galerkin matrix of a full 4-PI opvar on R^n0 x L2[a,b]^n1:
%
%     [ y_o ]   [ P        Q1<.>  ] [ x_o ]
%     [ y(s)] = [ Q2(s)    R      ] [ u   ]
%
% in the basis {e_k} (+) {phi_i}, with phi orthonormal in L2[a,b].  Because the
% basis is orthonormal, norm(M) IS the induced L2->L2 operator norm restricted
% to that subspace -- a LOWER bound on ||Op||, nondecreasing in D.  Report it
% with a D-refinement so the reader can see it plateau; a single D is not a
% norm, it is a truncation.
%
% Ordering: the ODE block first, then L2 component-major
% (index = n0 + (p-1)*(D+1) + i).

if nargin<2 || isempty(D),  D  = 12;     end
if nargin<3 || isempty(nq), nq = D + 12; end

dom = Op.I;  a = dom(1);  b = dom(2);
s   = Op.var1;   th = Op.var2;

P = Op.P;  Q1 = Op.Q1;  Q2 = Op.Q2;
R0 = Op.R.R0;  R1 = Op.R.R1;  R2 = Op.R.R2;

n0 = size(P,1);
if isempty(R0), n1 = 0; else, n1 = size(R0,1); end
nb = D+1;
N  = n0 + n1*nb;
M  = zeros(N);

[sq,wq] = gaussleg(nq,a,b);
Phi     = legval(sq,a,b,D);                 % nq x nb

% --- ODE-ODE
if n0>0, M(1:n0,1:n0) = double(P); end

% --- ODE <- L2   : int Q1(p,q)(s) phi_j(s) ds
if n0>0 && n1>0 && ~isempty(Q1)
    for p = 1:n0
        for q = 1:n1
            g = pgridval(Q1(p,q),{s},{sq});
            M(p, n0+(q-1)*nb+(1:nb)) = (wq(:).*g(:))' * Phi;
        end
    end
end

% --- L2 <- ODE   : int phi_i(s) Q2(p,q)(s) ds
if n0>0 && n1>0 && ~isempty(Q2)
    for p = 1:n1
        for q = 1:n0
            g = pgridval(Q2(p,q),{s},{sq});
            M(n0+(p-1)*nb+(1:nb), q) = Phi' * (wq(:).*g(:));
        end
    end
end

% --- L2 <- L2
for p = 1:n1
    for q = 1:n1
        r0 = sub_or_empty(R0,p,q);
        r1 = sub_or_empty(R1,p,q);
        r2 = sub_or_empty(R2,p,q);
        if isempty(r0) && isempty(r1) && isempty(r2), continue; end
        M(n0+(p-1)*nb+(1:nb), n0+(q-1)*nb+(1:nb)) = ...
            opgram(r0,r1,r2,s,th,dom,D,nq);
    end
end

info = struct('n0',n0,'n1',n1,'D',D,'nq',nq,'N',N,'dom',dom);
end


function e = sub_or_empty(A,p,q)
% Return A(p,q), or [] when the parameter is absent or identically zero, so
% opgram can skip the term rather than integrate a zero kernel.
if isempty(A), e = []; return; end
e = A(p,q);
if isa(e,'double')
    if e==0, e = []; end
    return
end
c = e.coefficient;
if isempty(c) || ~any(any(c~=0)), e = []; end
end

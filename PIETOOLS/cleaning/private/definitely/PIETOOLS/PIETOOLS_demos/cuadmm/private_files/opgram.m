function G = opgram(R0,R1,R2,s,th,dom,D,nq)
% opgram -- Galerkin matrix of the 3-PI operator
%     (R u)(s) = R0(s) u(s) + int_a^s R1(s,th) u(th) dth
%                           + int_s^b R2(s,th) u(th) dth
% in the L2-ORTHONORMAL shifted-Legendre basis of degree <= D on dom=[a,b].
% With an orthonormal basis the coefficient Euclidean norm IS the L2 norm, so
% norm(G) is the induced L2->L2 norm of R restricted to that subspace -- a
% LOWER bound on ||R||, increasing in D.
%
% QUADRATURE, NOT SYMBOLIC INTEGRATION.  The first version built the Legendre
% basis as polynomial objects and called int() twice.  Exact in principle, but
% Legendre expressed in the monomial basis is badly conditioned: the identity
% operator already showed ||G-I|| = 4e-6 at D=8.  Gauss-Legendre with the basis
% EVALUATED by recurrence at nodes is exact for polynomials and perfectly
% conditioned.  The inner integrals have a variable limit, so each is mapped to
% a fixed interval:  int_a^s f dth  =  (s-a) int_0^1 f(a+(s-a)t) dt,
% which is polynomial in (s,t) and therefore exact under a tensor rule.

if nargin<8 || isempty(nq), nq = D + 12; end      % generous: exact for our degrees
a = dom(1);  b = dom(2);

[sq,wq] = gaussleg(nq,a,b);        % outer rule in s
[tq,vq] = gaussleg(nq,0,1);        % inner rule on the mapped variable
Phi = legval(sq,a,b,D);            % nq x (D+1), orthonormal
G   = zeros(D+1);

% --- multiplier term
if ~isempty(R0)
    r0 = pgridval(R0,{s},{sq});
    G  = G + Phi' * (wq(:).*r0(:) .* Phi);
end

% --- int_a^s  R1(s,th) u(th) dth
if ~isempty(R1)
    G = G + tri_term(R1,s,th,sq,wq,tq,vq,Phi,a,b,D,'lower');
end

% --- int_s^b  R2(s,th) u(th) dth
if ~isempty(R2)
    G = G + tri_term(R2,s,th,sq,wq,tq,vq,Phi,a,b,D,'upper');
end
end


function Gt = tri_term(R,s,th,sq,wq,tq,vq,Phi,a,b,D,which)
% Accumulate int_a^b phi_i(s) [ int R(s,th) phi_j(th) dth ] ds over the
% triangular domain, via the affine map of the inner variable.
nq = numel(sq);
Gt = zeros(D+1);
for q = 1:nq
    sv = sq(q);
    if strcmp(which,'lower')
        thq = a + (sv-a)*tq;   jac = (sv-a);
    else
        thq = sv + (b-sv)*tq;  jac = (b-sv);
    end
    if jac == 0, continue; end
    rv  = pgridval(R,{s,th},{repmat(sv,size(thq)),thq});     % nq x 1
    Psi = legval(thq,a,b,D);                              % nq x (D+1)
    inner = jac * ( (vq(:).*rv(:))' * Psi );              % 1 x (D+1)
    Gt = Gt + wq(q) * (Phi(q,:)' * inner);
end
end



function [M,info] = op2gal(Pop,N)
% PROVENANCE.  scratchpad/reach1d/op2gal.m verbatim.  Used by T1 as a
% GRAM-FREE second opinion on an accepted certificate: eig(op2gal(Pop)) > 0
% cannot be satisfied by a Gram parametrisation artefact.  Read the README
% before asserting anything about -Dop through this route -- the t<s indicator
% costs ~8% relative quadrature error even on the analytic P = I certificate
% (measured, reach1d/cal2.log), so only Pop's strict positivity is assertable.
%
% Gauss-Legendre discretisation of the QUADRATIC FORM of a square opvar, so
% that sign-definiteness can be tested WITHOUT reference to any Gram matrix.
%   <(z,u),Pop(z,u)> = z'Pz + int z'Q1 u + int u'Q2 z + int u'R0 u
%                       + int int_{t<s} u(s)'R1(s,t)u(t) + int int_{t>s} ...R2
% With v_i = sqrt(w_i) u(s_i) the form is v'Mv with M as built below, so
% eig((M+M')/2) has the same signs as the operator's spectrum (up to the
% quadrature error introduced by the t<s / t>s indicator, which is why
% e2e_verify reports M at two values of N).
if nargin<2||isempty(N), N=24; end
d = Pop.dim;  m = d(1,1); n = d(2,1);
assert(d(1,1)==d(1,2) && d(2,1)==d(2,2),'op2gal: square operators only');
a = Pop.I(1); b = Pop.I(2);
[x,w] = gleg(N,a,b);
v1 = Pop.var1.varname{1};  v2 = Pop.var2.varname{1};
sw = sqrt(w(:));
M = zeros(m+n*N);
if m>0, M(1:m,1:m) = full(Pop.P); end
if m>0 && n>0
    Q1v = polyeval2(Pop.Q1,v1,v2,[x(:) x(:)]);
    Q2v = polyeval2(Pop.Q2,v1,v2,[x(:) x(:)]);
    for i=1:N
        ri = m+(i-1)*n+(1:n);
        M(1:m,ri) = sw(i)*Q1v(:,:,i);
        M(ri,1:m) = sw(i)*Q2v(:,:,i);
    end
end
if n>0
    R0v = polyeval2(Pop.R.R0,v1,v2,[x(:) x(:)]);
    [SI,TJ] = ndgrid(1:N,1:N);
    pts = [x(SI(:)).' ; x(TJ(:)).'].';
    R1v = polyeval2(Pop.R.R1,v1,v2,pts);
    R2v = polyeval2(Pop.R.R2,v1,v2,pts);
    for i=1:N
        ri = m+(i-1)*n+(1:n);
        M(ri,ri) = M(ri,ri) + R0v(:,:,i);
        for j=1:N
            if j==i, continue; end
            cj = m+(j-1)*n+(1:n);
            k = i+(j-1)*N;            % ndgrid index for (s_i,t_j)
            if j<i, Kk = R1v(:,:,k); else, Kk = R2v(:,:,k); end
            M(ri,cj) = M(ri,cj) + sw(i)*sw(j)*Kk;
        end
    end
end
M = (M+M')/2;
info.N=N; info.m=m; info.n=n; info.x=x; info.w=w;
end

function [x,w] = gleg(N,a,b)
% Golub-Welsch nodes/weights on [a,b].
k = 1:N-1;  be = k./sqrt(4*k.^2-1);
J = diag(be,1)+diag(be,-1);
[V,D] = eig(J);
[x,ord] = sort(diag(D));  V = V(:,ord);
w = 2*(V(1,:).^2).';
x = 0.5*(b-a)*x + 0.5*(a+b);
w = 0.5*(b-a)*w;
end

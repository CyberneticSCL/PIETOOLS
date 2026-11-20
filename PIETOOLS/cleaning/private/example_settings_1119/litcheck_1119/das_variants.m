% Spectral abscissa of the Das (2018, arXiv 1812.07186, Sec. 6.2) ODE-PDE
% system as coded in PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m, and variants.
% Chebyshev collocation on [0,1], Dirichlet BCs z(0) = z(1) = 0.
A = [-1.2142,  1.9649,  0.2232,  0.5616;
     -1.8042, -0.7260, -0.3479,  5.4355;
     -0.2898,  0.7381, -1.7606,  0.8294;
     -0.9417, -5.3399, -1.0704, -0.7590];
B = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];
C = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
fprintf('eig(A): %s\n', mat2str(eig(A).',5));
fprintf('max Re eig(A) = %.4f\n', max(real(eig(A))));
lams = [0, 2, 5, 8, pi^2-1, 9.8];
for N = [40 56]
  [D,x] = cheb(N);  s = (x+1)/2;  Ds = 2*D;  Dss = Ds^2;   % s(1)=1, s(N+1)=0
  in = 2:N;  n = N-1;  i0 = N+1;  i1 = 1;                  % interior; s=0; s=1
  for v = 1:5
    ab = zeros(size(lams));
    for k = 1:numel(lams)
      lam = lams(k);
      L = lam*eye(n) + Dss(in,in);
      % flux operators on interior values (boundary values are zero)
      F0 = Ds(i0,in);   F1 = Ds(i1,in);
      switch v
        case 1, sB = 1;  sC = 1;  F = F0;   % as coded: B*z_s(0), +C*xo
        case 2, sB = -1; sC = 1;  F = F0;   % ODE input -z_s(0)
        case 3, sB = 1;  sC = 1;  F = F1;   % ODE input z_s(1)
        case 4, sB = -1; sC = 1;  F = F1;   % ODE input -z_s(1)
        case 5, sB = 0;  sC = 1;  F = F0;   % no feedback into the ODE (cascade)
      end
      M = [A, sB*B(:,1)*F, sB*B(:,2)*F;
           sC*ones(n,1)*C(1,:), L, zeros(n);
           sC*ones(n,1)*C(2,:), zeros(n), L];
      ab(k) = max(real(eig(M)));
    end
    fprintf('N=%d variant %d: lam = %s -> max Re = %s\n', N, v, mat2str(lams,4), mat2str(ab,5));
  end
end
% threshold in lam of the coded system (variant 1), N = 56
[D,x] = cheb(56); Ds = 2*D; Dss = Ds^2; in = 2:56; n = 55; F0 = Ds(57,in);
f = @(lam) max(real(eig([A, B(:,1)*F0, B(:,2)*F0; ones(n,1)*C(1,:), lam*eye(n)+Dss(in,in), zeros(n); ones(n,1)*C(2,:), zeros(n), lam*eye(n)+Dss(in,in)])));
g = arrayfun(f, linspace(-20,9.8,150));
fprintf('coded: max Re over lam in [-20, 9.8]: min %.4f at lam = %.3f\n', min(g), -20+(9.8+20)*(find(g==min(g),1)-1)/149);
if any(g<0), fprintf('coded: stable for lam <= %.4f (grid)\n', max(linspace(-20,9.8,150).*(g<0)-1e9*(g>=0))); else, fprintf('coded: unstable on the whole grid\n'); end

function [D,x] = cheb(N)
% Trefethen, Spectral Methods in MATLAB, p. 54
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X - X';
D = (c*(1./c)')./(dX + eye(N+1));
D = D - diag(sum(D,2));
end

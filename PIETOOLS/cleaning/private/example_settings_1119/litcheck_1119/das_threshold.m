% Stability threshold in lam of the Das ODE-PDE system with the flux sign
% flipped (ODE input -B*z_s(0)) and of the coded system; Chebyshev, N = 56, 72.
A = [-1.2142,  1.9649,  0.2232,  0.5616;
     -1.8042, -0.7260, -0.3479,  5.4355;
     -0.2898,  0.7381, -1.7606,  0.8294;
     -0.9417, -5.3399, -1.0704, -0.7590];
B = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];
C = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
for N = [56 72]
  [D,~] = cheb(N);  Ds = 2*D;  Dss = Ds^2;  in = 2:N;  n = N-1;  F0 = Ds(N+1,in);
  for sB = [1 -1]
    f = @(lam) max(real(eig([A, sB*B(:,1)*F0, sB*B(:,2)*F0; ...
          ones(n,1)*C(1,:), lam*eye(n)+Dss(in,in), zeros(n); ...
          ones(n,1)*C(2,:), zeros(n), lam*eye(n)+Dss(in,in)])));
    lo = -30;  hi = 9.86;
    if f(hi) < 0, fprintf('N=%d sB=%+d: stable at lam = 9.86\n', N, sB); continue; end
    if f(lo) > 0, fprintf('N=%d sB=%+d: unstable at lam = -30\n', N, sB); continue; end
    for it = 1:50, mid = (lo+hi)/2; if f(mid) < 0, lo = mid; else, hi = mid; end, end
    fprintf('N=%d sB=%+d: threshold lam* = %.6f (max Re at pi^2-1: %.5f)\n', N, sB, lo, f(pi^2-1));
  end
end
function [D,x] = cheb(N)
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X - X';
D = (c*(1./c)')./(dX + eye(N+1));
D = D - diag(sum(D,2));
end

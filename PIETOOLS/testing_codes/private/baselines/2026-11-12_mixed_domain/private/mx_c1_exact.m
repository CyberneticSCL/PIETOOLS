function [gi,h2] = mx_c1_exact()
% [GI,H2] = MX_C1_EXACT() returns the exact H-infinity gain and H2 norm of
% the control C1 (MX_SYS2). The modes 2 sin(m pi s1) sin(n pi s2), m, n odd,
% carry the residue a_m a_n, a_m = 8/(m^2 pi^2) > 0, so the impulse
% response is f(t)^2 with f(t) = sum_{m odd} a_m exp(-m^2 pi^2 t). Every
% residue is positive, so the gain is G(0) = int f^2 dt; H2^2 = int f^4 dt.
% Prints the gain, its double-series check and the H2 norm.
% INPUT
% none.
% OUTPUT
% - gi:     double, the H-infinity gain;
% - h2:     double, the H2 norm.

% Initial coding MMP, 11/12/2026.
m = (1:2:40001)';
a = 8./(m.^2*pi^2);  l = m.^2*pi^2;
f = @(t) sum(a.*exp(-l.*t),1);
gi = integral(@(t) f(t).^2,0,Inf,'AbsTol',1e-15,'RelTol',1e-12);
h2 = sqrt(integral(@(t) f(t).^4,0,Inf,'AbsTol',1e-15,'RelTol',1e-12));
% Cross-check of the gain: sum over m, n odd <= 2001 of a_m a_n/((m^2+n^2) pi^2)
mm = (1:2:2001)';  A = 8./(mm.^2*pi^2);
gs = sum(sum((A*A')./((mm.^2+mm'.^2)*pi^2)));
fprintf('C1 exact: gain %.10f (double series to 2001: %.10f), H2 norm %.10f\n',gi,gs,h2);
end

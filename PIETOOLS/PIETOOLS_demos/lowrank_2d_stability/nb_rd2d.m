function PIE = nb_rd2d(n,lam)
% nb_rd2d(n,lam) -- n-state 2-D reaction-diffusion (heat) equation on the unit
% square, Dirichlet on all four edges:
%       x_t = x_{s1s1} + x_{s2s2} + lam*x,   x in R^n.
% Spelling follows PIETOOLS' PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq.
%
% STABILITY BOUNDARY (analytic reference).  The Dirichlet Laplacian on (0,1)^2
% has eigenfunctions sin(j*pi*s1)*sin(k*pi*s2), j,k >= 1, with eigenvalues
% -(j^2+k^2)*pi^2.  The generator A = Delta + lam*I is diagonal in that basis
% (no advection term to couple the modes) with spectrum lam - (j^2+k^2)*pi^2,
% whose supremum is attained at j=k=1.  Hence
%       exponentially stable  <=>  lam < lam* = 2*pi^2 = 19.7392.
% Each state decouples, so lam* is independent of n.
%
% lam defaults to 2, i.e. 10.1% of lam* -- a slack operating point.
if nargin<2, lam = 2; end
pvar s1 s2
clear stateNameGenerator
x = pde_var('state',n,[s1;s2],[0,1;0,1]);
sys = [diff(x,'t')==diff(x,s1,2)+diff(x,s2,2)+lam*x;
       subs(x,s1,0)==zeros(n,1);   subs(x,s1,1)==zeros(n,1);
       subs(x,s2,0)==zeros(n,1);   subs(x,s2,1)==zeros(n,1)];
PIE = convert(sys,'pie');
end

function [sol,M] = bl_b_fisher(R,use_bnd)
% Nonlinear local stability on the L2 ball of radius R (arXiv 2604.01115).
[sol,I0] = fisher_prog(R,use_bnd,true,'mosek');
M = struct('info',I0);
end

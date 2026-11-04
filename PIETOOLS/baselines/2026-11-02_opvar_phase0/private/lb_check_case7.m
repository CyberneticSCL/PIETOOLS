function lb_check_case7()
% LB_CHECK_CASE7() measures whether the documented parameter override
% 'lam=...' changes PDE library example 7 (Parabolic_Eq_Gahlawat), and
% whether the catalog setter, which rebuilds c from lam, does.
%
% OUTPUT (printed)
%   fingerprints of PIE.A (lb_fp) for the default, for {'lam=5;'}, and for
%   the catalog setter at lam=5, with the verdict.

M = lb_meta();
m = M(strcmp({M.lib},'PDE') & strcmp({M.id},'7'));
P0 = lb_build_system(m, [], true);
m2 = m;   m2.setter = 'lam=%.10g;';
P1 = lb_build_system(m2, 5, true);
P2 = lb_build_system(m, 5, true);
f0 = lb_fp(P0.A);  f1 = lb_fp(P1.A);  f2 = lb_fp(P2.A);
fprintf('\nCASE7 |A_lam=5 alone - A_default|_1 = %.3g ; |A_setter - A_default|_1 = %.3g\n', ...
    sum(abs(f1-f0)), sum(abs(f2-f0)));
end

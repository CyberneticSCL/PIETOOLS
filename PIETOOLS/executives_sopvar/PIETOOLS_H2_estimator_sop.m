function [prog,Lop,gam,P,Z,W,info] = PIETOOLS_H2_estimator_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,LOP,GAM,P,Z,W,INFO] = PIETOOLS_H2_ESTIMATOR_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_H2_estimator: the H2-optimal observer gain
% L = P^{-1} Z from the LPI P >= 0 (+ eppos), Z and W free,
% [-gam I_z, C; C', (T'P)A + (A'P)T + (T'Z)C2 + (C2'Z')T] <= 0,
% [W, -(B'P + D21' Z'); (.)', P] >= 0, gam >= trace(W), on the container
% path; L by the stock getObserver on the solved P and Z. 1-D only.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
[prog,Lop,gam,info] = h2_build_sop(PIE,settings,'est','H2_estimator_sop');
P = safe_op(info.P);    Z = safe_op(info.Z);    W = safe_op(info.W);
end

function X = safe_op(C)
X = C;
if isa(C,'copvar')
    try,    X = cop2opvar_sop(C);   catch,  end
end
end

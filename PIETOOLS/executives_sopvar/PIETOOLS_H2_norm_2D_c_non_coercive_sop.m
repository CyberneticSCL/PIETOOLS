function [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm_2D_c_non_coercive_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,R,Q,INFO] = PIETOOLS_H2_NORM_2D_C_NON_COERCIVE_SOP(PIE,SETTINGS)
% the container version of PIETOOLS_H2_norm_2D_c_non_coercive: R >= 0,
% TQ = R, W >= 0 on R^nz, [-gam I_w, B'; B, Q'A' + AQ] <= 0,
% [W, CQ; Q'C', R] >= 0, gam >= trace(W) (H2_BUILD_SOP 'c'). OPTIONS is
% accepted for signature parity and not used.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
[prog,W,gam,info] = h2_build_sop(PIE,settings,'c','H2_norm_2D_c_non_coercive_sop');
R = info.R;     Q = info.Q;
end

function [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm_2D_o_non_coercive_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,R,Q,INFO] = PIETOOLS_H2_NORM_2D_O_NON_COERCIVE_SOP(PIE,SETTINGS)
% the container version of PIETOOLS_H2_norm_2D_o_non_coercive: R >= 0,
% T'Q = R, W >= 0 on R^nw, [-gam I_z, C; C', Q'A + A'Q] <= 0,
% [W, B'Q; Q'B, R] >= 0, gam >= trace(W) (H2_BUILD_SOP 'o'). OPTIONS is
% accepted for signature parity and not used.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
[prog,W,gam,info] = h2_build_sop(PIE,settings,'o','H2_norm_2D_o_non_coercive_sop');
R = info.R;     Q = info.Q;
end

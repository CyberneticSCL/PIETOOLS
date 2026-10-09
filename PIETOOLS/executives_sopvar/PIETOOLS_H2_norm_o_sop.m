function [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm_o_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,R,Q,INFO] = PIETOOLS_H2_NORM_O_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_H2_norm_o: the observability-side LPI, R >= 0, T'Q = R,
% W >= 0 on the input space, [-gam I_z, C; C', Q'A + A'Q] <= 0,
% [W, B'Q; Q'B, R] >= 0, gam >= trace(W), on the container path. GAM is the
% stock's return value. A 2-D PIE is routed to
% PIETOOLS_H2_norm_2D_o_non_coercive_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm_2D_o_non_coercive_sop(PIE,settings);
    return
end
[prog,W,gam,info] = h2_build_sop(PIE,settings,'o','H2_norm_o_sop');
R = info.R;     Q = info.Q;
end

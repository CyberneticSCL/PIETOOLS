function [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm_c_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,R,Q,INFO] = PIETOOLS_H2_NORM_C_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_H2_norm_c: the H2 norm bound through the
% controllability-side LPI, R >= 0, TQ = R, W >= 0 on the output space,
% [-gam I_w, B'; B, Q'A' + AQ] <= 0, [W, CQ; Q'C', R] >= 0, gam >= trace(W),
% on the container path (executives_sopvar/README.md; the LPI as the stock
% file and its cx_exec transcription write it). GAM is the stock's return
% value (the LPI objective). OPTIONS is accepted for signature parity and
% not used (the stock's 'h2' option selects lpi_ineq). A 2-D PIE is routed to
% PIETOOLS_H2_norm_2D_c_non_coercive_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm_2D_c_non_coercive_sop(PIE,settings);
    return
end
[prog,W,gam,info] = h2_build_sop(PIE,settings,'c','H2_norm_c_sop');
R = info.R;     Q = info.Q;
end

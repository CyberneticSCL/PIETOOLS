function [prog,R,gam,info] = PIETOOLS_Hinf_gain_dual_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,R,GAM,INFO] = PIETOOLS_HINF_GAIN_DUAL_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_Hinf_gain_dual: the L2 gain bound by the dual KYP LPI
% in its non-coercive form, TQ = R >= 0, K33 = Q'A' + AQ, and the negativity
% constraint on the container path (executives_sopvar/README.md). Measured
% on io1 (10/08/2026): this form reaches the closed-form gain, where the
% primal non-coercive form stalls 2.5e-6 above it.
%
% INPUT / OUTPUT as PIETOOLS_Hinf_gain_sop; a boundary disturbance is
% refused, as stock; a 2-D PIE is routed to PIETOOLS_Hinf_gain_dual_2D_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,R,gam,info] = PIETOOLS_Hinf_gain_dual_2D_sop(PIE,settings);
    return
end
[prog,R,gam,info] = hinf_build_sop(PIE,settings,'Qd','Hinf_gain_dual_sop');
end

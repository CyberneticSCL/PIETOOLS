function [prog,R,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,R,GAM,INFO] = PIETOOLS_HINF_GAIN_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_Hinf_gain: the L2 gain bound of the PIE
%
%   T x_t = A x + B1 w,   z = C1 x + D11 w
%
% by the primal KYP LPI in its non-coercive form, R = T'Q >= 0 with Q a
% free operator, and the negativity constraint on the container path
% (executives_sopvar/README.md). With a boundary disturbance (Tw ~= 0) the
% coercive form is used, as the stock executive does.
%
% INPUT
% - PIE:      'pie_struct' (1-D; a 2-D PIE is routed to
%             PIETOOLS_Hinf_gain_2D_non_coercive_sop);
% - settings: lpisettings struct (default lpisettings('light')); the field
%             settings.sop (EXEC_SOP_SETTINGS) controls the slack sizing,
%             the degree loop, the dual read-back and the numerical witness.
%
% OUTPUT
% - prog:     the solved LPI program;
% - R:        the storage operator, as an 'opvar' (the container in INFO);
% - gam:      the certified gain bound (NaN when the solve did not certify);
% - info:     status, certified, gain_lb (the numerical lower bound of the
%             frequency response), omega, gap, dual (the dual kernel of the
%             negativity constraint), dual_alignment, degrees, terms, shape,
%             hist (the degree loop), storage, K, N.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,R,gam,info] = PIETOOLS_Hinf_gain_2D_non_coercive_sop(PIE,settings);
    return
end
[prog,R,gam,info] = hinf_build_sop(PIE,settings,'Q','Hinf_gain_sop');
end

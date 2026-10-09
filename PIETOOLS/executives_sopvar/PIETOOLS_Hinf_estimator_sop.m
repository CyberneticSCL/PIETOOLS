function [prog,L,gam,P,Z,info] = PIETOOLS_Hinf_estimator_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,L,GAM,P,Z,INFO] = PIETOOLS_HINF_ESTIMATOR_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_Hinf_estimator: the H-infinity optimal
% observer gain L = P^{-1} Z from the LPI P >= 0 (+ eppos), Z free,
% [-gam I_w, -Dzw', -(PB + Z Dyw)' T; (.)', -gam I_z, Cz; (.)', Cz',
%  (PA + Z Cy)' T + T' (PA + Z Cy)] <= 0, with the Tw terms of the stock when
% the disturbance enters at the boundary, on the container path
% (SYNTH_BUILD_SOP 'est'). L by the stock getObserver on the solved P and Z.
% A 2-D PIE is routed to PIETOOLS_Hinf_estimator_2D_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,L,gam,P,Z,info] = PIETOOLS_Hinf_estimator_2D_sop(PIE,settings);
    return
end
[prog,L,gam,info] = synth_build_sop(PIE,settings,'est','Hinf_estimator_sop');
P = info.P;     Z = info.Z;
end

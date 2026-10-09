function [prog,L,gam,P,Z,info] = PIETOOLS_Hinf_estimator_2D_sop(PIE,settings,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,L,GAM,P,Z,INFO] = PIETOOLS_HINF_ESTIMATOR_2D_SOP(PIE,SETTINGS,GAIN)
% the container version of PIETOOLS_Hinf_estimator_2D: the same LPI as the
% 1-D estimator (SYNTH_BUILD_SOP 'est') on a 2-D PIE, the storage and the
% stock slack from settings_2d (POSLPIVAR_SETTINGS_2D_SOP), Z with the
% largest entry of settings_2d.Zop_deg as its degree cap, L by
% getObserver_2D. A nonzero GAIN poses the feasibility test at that gamma
% (the stock's third argument; the stock bisection options are not
% reproduced: the degree loop of the container path replaces them).
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
if nargin>=3 && ~isempty(gain) && gain~=0,  settings.sop.gam_fixed = gain;  end
[prog,L,gam,info] = synth_build_sop(PIE,settings,'est','Hinf_estimator_2D_sop');
P = info.P;     Z = info.Z;
end

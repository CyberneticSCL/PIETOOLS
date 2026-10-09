function [prog,P,gam,info] = PIETOOLS_Hinf_gain_2D_sop(PIE,settings,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,GAM,INFO] = PIETOOLS_HINF_GAIN_2D_SOP(PIE,SETTINGS,GAIN) the
% container version of PIETOOLS_Hinf_gain_2D: the coercive H-infinity gain
% LPI, P >= 0 (+ eppos I), [-gam I_w, D', (PB)'T; D, -gam I_z, C; T'PB, C',
% (PA)'T + T'PA] <= 0, with the Tw terms when the disturbance enters at the
% boundary (HINF_BUILD_SOP 'P'; storage and stock slack from settings_2d).
% A nonzero GAIN poses the feasibility test at that gamma (the stock's
% third argument; the stock bisection options are not reproduced).
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
if nargin>=3 && ~isempty(gain) && gain~=0,  settings.sop.gam_fixed = gain;  end
[prog,P,gam,info] = hinf_build_sop(PIE,settings,'P','Hinf_gain_2D_sop');
end

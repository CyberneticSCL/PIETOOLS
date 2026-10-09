function [prog,R,gam,info] = PIETOOLS_Hinf_gain_2D_non_coercive_sop(PIE,settings,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,R,GAM,INFO] = PIETOOLS_HINF_GAIN_2D_NON_COERCIVE_SOP(PIE,SETTINGS,GAIN)
% the container version of PIETOOLS_Hinf_gain_2D_non_coercive: the Q form,
% R = T'Q >= 0 (Q an lpivar sized from R), [-gam I_w, D', B'Q; D, -gam I_z,
% C; Q'B, C', A'Q + Q'A] <= 0 (HINF_BUILD_SOP 'Q'); with a boundary
% disturbance the coercive form, as the stock. A nonzero GAIN poses the
% feasibility test at that gamma.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
if nargin>=3 && ~isempty(gain) && gain~=0,  settings.sop.gam_fixed = gain;  end
[prog,R,gam,info] = hinf_build_sop(PIE,settings,'Q','Hinf_gain_2D_non_coercive_sop');
end

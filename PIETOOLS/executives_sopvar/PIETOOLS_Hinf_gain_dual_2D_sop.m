function [prog,P,gam,info] = PIETOOLS_Hinf_gain_dual_2D_sop(PIE,settings,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,GAM,INFO] = PIETOOLS_HINF_GAIN_DUAL_2D_SOP(PIE,SETTINGS,GAIN) the
% container version of PIETOOLS_Hinf_gain_dual_2D: the coercive dual LPI,
% P >= 0 (+ eppos I; [1e-4;0;0;1e-6] when the settings carry no eppos, as
% the stock), [-gam I_z, D, CPT'; D', -gam I_w, B'; TPC', B, TPA' + APT']
% <= 0 (HINF_BUILD_SOP 'Pd'). Boundary disturbances are refused, as the
% stock. A nonzero GAIN poses the feasibility test at that gamma.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
if nargin>=3 && ~isempty(gain) && gain~=0,  settings.sop.gam_fixed = gain;  end
[prog,P,gam,info] = hinf_build_sop(PIE,settings,'Pd','Hinf_gain_dual_2D_sop');
end

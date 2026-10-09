function [prog,P,gam,info] = PIETOOLS_Hinf_gain_dual_coercive_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,GAM,INFO] = PIETOOLS_HINF_GAIN_DUAL_COERCIVE_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_Hinf_gain_dual_coercive: the dual KYP LPI
% with K33 = TPA' + APT', P >= 0, on the container path. MEASURED
% (10/08/2026, io1): the stock executive returns gamma = 8252 at light and
% 14905 at heavy, and every tailored slack is reported primal infeasible;
% this form is not usable on that plant as shipped, and the executive
% reports its status rather than a bound.
%
% INPUT / OUTPUT as PIETOOLS_Hinf_gain_sop; a boundary disturbance is
% refused; a 2-D PIE is routed to PIETOOLS_Hinf_gain_dual_2D_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,P,gam,info] = PIETOOLS_Hinf_gain_dual_2D_sop(PIE,settings);
    return
end
[prog,P,gam,info] = hinf_build_sop(PIE,settings,'Pd','Hinf_gain_dual_coercive_sop');
end

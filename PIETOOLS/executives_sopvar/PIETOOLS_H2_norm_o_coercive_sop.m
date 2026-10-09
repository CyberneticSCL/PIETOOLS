function [prog,Wo,gam,info] = PIETOOLS_H2_norm_o_coercive_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,WO,GAM,INFO] = PIETOOLS_H2_NORM_O_COERCIVE_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_H2_norm_o_coercive: W >= 0,
% (A'W)T + T'(WA) + C'C <= 0, gam >= trace(B' W B), GAM = sqrt(gam) as the
% stock returns it. A 2-D PIE is routed to PIETOOLS_H2_norm_2D_o_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,Wo,gam,info] = PIETOOLS_H2_norm_2D_o_sop(PIE,settings);
    return
end
[prog,Wo,gam,info] = h2_build_sop(PIE,settings,'oco','H2_norm_o_coercive_sop');
end

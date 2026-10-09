function [prog,Wc,gam,info] = PIETOOLS_H2_norm_c_coercive_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,WC,GAM,INFO] = PIETOOLS_H2_NORM_C_COERCIVE_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_H2_norm_c_coercive: W >= 0,
% (AW)T' + T(WA') + BB' <= 0, gam >= trace(C W C'), the norm bound
% GAM = sqrt(gam) as the stock returns it. A 2-D PIE is routed to
% PIETOOLS_H2_norm_2D_c_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,Wc,gam,info] = PIETOOLS_H2_norm_2D_c_sop(PIE,settings);
    return
end
[prog,Wc,gam,info] = h2_build_sop(PIE,settings,'cco','H2_norm_c_coercive_sop');
end

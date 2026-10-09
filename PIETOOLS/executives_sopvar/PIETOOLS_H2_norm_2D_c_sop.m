function [prog,W,gam,info] = PIETOOLS_H2_norm_2D_c_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,INFO] = PIETOOLS_H2_NORM_2D_C_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_H2_norm_2D_c: the controllability Gramian LPI,
% W >= 0 (+ eppos on R^n and L2[x,y], as the stock), (AW)T' + T(WA') + BB'
% <= 0, gam >= trace(C W C'), GAM = sqrt(gam) as the 1-D coercive executive
% returns it (H2_BUILD_SOP 'cco'). OPTIONS is accepted for signature parity
% and not used.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
[prog,W,gam,info] = h2_build_sop(PIE,settings,'cco','H2_norm_2D_c_sop');
end

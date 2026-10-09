function [prog,W,gam,info] = PIETOOLS_H2_norm_2D_o_sop(PIE,settings,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,INFO] = PIETOOLS_H2_NORM_2D_O_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_H2_norm_2D_o: the observability Gramian LPI,
% W >= 0 (+ eppos on R^n and L2[x,y]), (A'W)T + T'(WA) + C'C <= 0,
% gam >= trace(B' W B), GAM = sqrt(gam) (H2_BUILD_SOP 'oco'). OPTIONS is
% accepted for signature parity and not used.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
[prog,W,gam,info] = h2_build_sop(PIE,settings,'oco','H2_norm_2D_o_sop');
end

function [prog,Kop,gam,P,Z,W,info] = PIETOOLS_H2_control_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,KOP,GAM,P,Z,W,INFO] = PIETOOLS_H2_CONTROL_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_H2_control: the H2-optimal state feedback
% u = K x, K = Z P^{-1}, from the LPI P >= 0 (+ eppos), Z and W free,
% [-gam I_w, B'; B, (TP)A' + (AP)T' + (B2 Z)T' + (T Z')B2'] <= 0,
% [W, CP + D12 Z; (.)', P] >= 0, gam >= trace(W), on the container path.
% K is recovered by the stock getController on the solved P and Z (opvars).
% 1-D only, as the stock.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
[prog,Kop,gam,info] = h2_build_sop(PIE,settings,'ctrl','H2_control_sop');
P = safe_op(info.P);    Z = safe_op(info.Z);    W = safe_op(info.W);
end

function X = safe_op(C)
X = C;
if isa(C,'copvar')
    try,    X = cop2opvar_sop(C);   catch,  end
end
end

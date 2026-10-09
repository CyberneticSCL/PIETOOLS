function [prog,K,gam,P,Z,info] = PIETOOLS_Hinf_control_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,K,GAM,P,Z,INFO] = PIETOOLS_HINF_CONTROL_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_Hinf_control: the H-infinity optimal state
% feedback u = K x, K = Z P^{-1}, from the LPI P >= 0 (+ eppos), Z free,
% [-gam I_z, Dzw, (CP + Dzu Z) T'; (.)', -gam I_w, B'; (.)', B,
%  (AP + Bu Z) T' + T (AP + Bu Z)'] <= 0, on the container path
% (SYNTH_BUILD_SOP 'ctrl'; executives_sopvar/README.md). K is recovered by
% the stock getController on the solved P and Z (opvars). 1-D only, as the
% stock; boundary disturbances or inputs are refused, as the stock.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
[prog,K,gam,info] = synth_build_sop(PIE,settings,'ctrl','Hinf_control_sop');
P = info.P;     Z = info.Z;
end

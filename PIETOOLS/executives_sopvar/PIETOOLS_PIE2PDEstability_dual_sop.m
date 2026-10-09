function [prog,P,info] = PIETOOLS_PIE2PDEstability_dual_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,INFO] = PIETOOLS_PIE2PDESTABILITY_DUAL_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_PIE2PDEstability_dual: P >= 0 (+ eppos2 T T'),
% TQ = P, AQ + Q'A' + epneg P <= 0, on the container path. INPUT / OUTPUT as
% PIETOOLS_PDEstability_sop; a 2-D PIE is routed to
% PIETOOLS_stability_dual_2D_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings)
    settings = lpisettings('heavy');    settings.eppos = 1e-4;  settings.eppos2 = 1e-6;  settings.epneg = 0;
end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,P,info] = PIETOOLS_stability_dual_2D_sop(PIE,settings);
    return
end
[prog,P,info] = stability_build_sop(PIE,settings,'Qd','PIE2PDEstability_dual_sop');
end

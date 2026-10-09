function [prog,P,info] = PIETOOLS_stability_2D_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,INFO] = PIETOOLS_STABILITY_2D_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_stability_2D: P >= 0 (+ eppos I, the four-entry eppos
% of settings_2d), T'PA + A'PT + 2 epneg T'PT <= 0, on the container path
% (STABILITY_BUILD_SOP 'P'; storage and stock slack from settings_2d through
% POSLPIVAR_SETTINGS_2D_SOP). Without SETTINGS the stock defaults: light 2-D,
% eppos = 1e-2 on every space, epneg = 0.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings)
    settings = lpisettings('light');
    settings.sos_opts.simplify = 1;
    settings.settings_2d.eppos = 1e-2*ones(4,1);    settings.settings_2d.epneg = 0;
end
[prog,P,info] = stability_build_sop(PIE,settings,'P','stability_2D_sop');
end

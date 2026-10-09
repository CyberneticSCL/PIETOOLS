function [prog,P,info] = PIETOOLS_PDEstability_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,INFO] = PIETOOLS_PDESTABILITY_SOP(PIE,SETTINGS) the container
% version of PIETOOLS_PDEstability: a Lyapunov certificate of stability of
% the PIE T x_t = A x, V = <Tx,PTx> with P >= 0, T'PA + A'PT + epneg T'PT
% <= 0, on the container path (executives_sopvar/README.md).
%
% INPUT
% - PIE:      'pie_struct' (a 2-D PIE is routed to PIETOOLS_stability_2D_sop);
% - settings: lpisettings struct (default: 'heavy' with eppos 1e-4, eppos2
%             1e-6, epneg 0, as the stock); settings.sop as EXEC_SOP_SETTINGS.
% OUTPUT
% - prog:  the solved program;  P: the storage operator (opvar);
% - info:  status, certified, epneg, maxre (the numerical spectrum of (A,T)),
%          dual, degrees, terms, shape, hist.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings)
    settings = lpisettings('heavy');    settings.eppos = 1e-4;  settings.eppos2 = 1e-6;  settings.epneg = 0;
end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,P,info] = PIETOOLS_stability_2D_sop(PIE,settings);
    return
end
[prog,P,info] = stability_build_sop(PIE,settings,'P','PDEstability_sop');
end

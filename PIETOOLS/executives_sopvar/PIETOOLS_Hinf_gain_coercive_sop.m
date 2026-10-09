function [prog,P,gam,info] = PIETOOLS_Hinf_gain_coercive_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,GAM,INFO] = PIETOOLS_HINF_GAIN_COERCIVE_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_Hinf_gain_coercive: the L2 gain bound by
% the primal KYP LPI with the storage V = <Tx,PTx>, P >= 0, and the
% negativity constraint on the container path (executives_sopvar/README.md).
% The slack is sized from P with the weight raised by settings.sop.dw_P
% (default 2; measured in 1-D and 2-D, lpi_programming_sopvar/README.md).
%
% INPUT / OUTPUT as PIETOOLS_Hinf_gain_sop, with P the storage operator; a
% 2-D PIE is routed to PIETOOLS_Hinf_gain_2D_sop.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
PIE = initialize(PIE);
if PIE.dim==2
    [prog,P,gam,info] = PIETOOLS_Hinf_gain_2D_sop(PIE,settings);
    return
end
[prog,P,gam,info] = hinf_build_sop(PIE,settings,'P','Hinf_gain_coercive_sop');
end

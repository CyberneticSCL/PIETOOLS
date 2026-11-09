function [names, execs] = lb_toggles()
% [NAMES, EXECS] = LB_TOGGLES() lists the executive toggles of the stock
% flow in the order PIETOOLS_auto_execute tests them, with the executive
% each one selects.
%
% OUTPUT
%   names:  cellstr of toggle names ('stability', 'stability_dual', ...)
%   execs:  cellstr of the executive functions, same order
%
% NOTES
%   - Mapping of executives/PIETOOLS_auto_execute.m. PIETOOLS_DDE.m maps
%     H2_control to PIETOOLS_H2_controller and PIETOOLS_DDF.m maps
%     stability to PIETOOLS_stability; neither function exists, so the
%     auto_execute names are used for all three libraries.

T = {'stability',      'PIETOOLS_PIE2PDEstability'
     'stability_dual', 'PIETOOLS_PIE2PDEstability_dual'
     'Hinf_gain',      'PIETOOLS_Hinf_gain'
     'Hinf_gain_dual', 'PIETOOLS_Hinf_gain_dual'
     'Hinf_estimator', 'PIETOOLS_Hinf_estimator'
     'Hinf_control',   'PIETOOLS_Hinf_control'
     'H2_norm',        'PIETOOLS_H2_norm_o'
     'H2_norm_dual',   'PIETOOLS_H2_norm_c'
     'H2_estimator',   'PIETOOLS_H2_estimator'
     'H2_control',     'PIETOOLS_H2_control'
     'well_posedness', 'PIETOOLS_well_posedness'};
names = T(:,1)';
execs = T(:,2)';
end

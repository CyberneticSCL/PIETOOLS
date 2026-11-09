function [fn,fa] = bl_sop(exec)                                             % CC, 11/07/2026
% [FN,FA] = BL_SOP(EXEC) names the container executive FN, and its form
% argument FA, that pose the LPI of the stock executive PIETOOLS_<EXEC>; a
% builder runs the case as feval(FN,PIE,st,FA{:}).
%
% INPUT
%   exec:  char, a stock executive name without 'PIETOOLS_', as bl_cases
%          lists them ('PIE2PDEstability', 'Hinf_gain_2D', ...).
% OUTPUT
%   fn:    char, the executive in executives_sopvar;
%   fa:    cell, its form argument ({} for an executive with one form).
%
% Why: convert returns a pie_struct of 'copvar' operators since 3df33dcc
% (10/29/2026), which every stock executive rejects (Top.dim), and the opvar
% removal retires executives/. The names in bl_cases stay the case labels.
% Forms as the consolidation of 10/28/2026 maps the stock executives (headers
% of PIETOOLS_stability_sop, _Hinf_gain_sop, _H2_norm_sop). A stock 2-D
% executive maps to the form of its LPI: 'stability_2D' was the P form (K.f =
% 0, blocks [8 424]); the stock 2-D gain and H2 executives are coercive, except
% _non_coercive. The container executives pose other SDPs than the stock ones
% (block partition, slack sizing), so values banked on the stock path
% (bl_expect.tsv, results/2026-09-24) do not carry over.

M = {'PIE2PDEstability',          'PIETOOLS_stability_sop',      {'Q'};
     'PIE2PDEstability_dual',     'PIETOOLS_stability_sop',      {'Qd'};
     'PDEstability',              'PIETOOLS_stability_sop',      {'P'};
     'PDEstability_dual',         'PIETOOLS_stability_sop',      {'Pd'};
     'stability_2D',              'PIETOOLS_stability_sop',      {'P'};
     'stability_dual_2D',         'PIETOOLS_stability_sop',      {'Pd'};
     'Hinf_gain',                 'PIETOOLS_Hinf_gain_sop',      {'Q'};
     'Hinf_gain_coercive',        'PIETOOLS_Hinf_gain_sop',      {'P'};
     'Hinf_gain_dual',            'PIETOOLS_Hinf_gain_sop',      {'Qd'};
     'Hinf_gain_dual_coercive',   'PIETOOLS_Hinf_gain_sop',      {'Pd'};
     'Hinf_gain_2D',              'PIETOOLS_Hinf_gain_sop',      {'P'};
     'Hinf_gain_2D_non_coercive', 'PIETOOLS_Hinf_gain_sop',      {'Q'};
     'Hinf_gain_dual_2D',         'PIETOOLS_Hinf_gain_sop',      {'Pd'};
     'H2_norm_c',                 'PIETOOLS_H2_norm_sop',        {'c'};
     'H2_norm_o',                 'PIETOOLS_H2_norm_sop',        {'o'};
     'H2_norm_c_coercive',        'PIETOOLS_H2_norm_sop',        {'cco'};
     'H2_norm_o_coercive',        'PIETOOLS_H2_norm_sop',        {'oco'};
     'H2_norm_2D_c',              'PIETOOLS_H2_norm_sop',        {'cco'};
     'Hinf_estimator',            'PIETOOLS_Hinf_estimator_sop', {};
     'Hinf_control',              'PIETOOLS_Hinf_control_sop',   {};
     'H2_estimator',              'PIETOOLS_H2_estimator_sop',   {};
     'H2_control',                'PIETOOLS_H2_control_sop',     {};
     'well_posedness',            'PIETOOLS_well_posedness_sop', {}};
j = find(strcmp(M(:,1),exec),1);
if isempty(j), error('bl_sop:exec','no container executive for PIETOOLS_%s',exec); end
fn = M{j,2};    fa = M{j,3};
end

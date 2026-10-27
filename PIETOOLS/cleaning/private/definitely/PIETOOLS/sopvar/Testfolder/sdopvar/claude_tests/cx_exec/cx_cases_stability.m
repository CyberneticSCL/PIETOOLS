function C = cx_cases_stability(which)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = CX_CASES_STABILITY(WHICH) returns the cx_compare case structs of the
% stability family, mirroring PIETOOLS_demos/cuadmm/bl_cases.m (ids, plants,
% settings; Mosek in 1-D as bl_settings uses).
%
% WHICH: '1d' (default) the eight 1-D feasibility cases, solved by
%        cx_compare; '2d' the two 2-D cases, STRUCTURE ONLY - never pass them
%        to cx_compare, which solves (use cx_stability_check(C,'nosolve')).
%
% thresh: cx_compare bisects plant{2} assuming feasibility at LARGE values
% (cx_bisect is written for gamma). For 'rd' a larger frac is LESS stable,
% so that bisection cannot bracket; cx_stability_thresh does it with the
% orientation reversed. Kept on for stab_rd1/stabpde_rd1 as specified.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(which),  which = '1d';   end
mk = @(id,exec,setname,solver,plant,thresh) struct('id',id,'exec',exec, ...
        'setname',setname,'solver',solver,'kind','feas','plant',{plant},'thresh',thresh);
switch lower(which)
    case '1d'
        C = [mk('stab_rd1',     'PIE2PDEstability',     'light','mosek',{'rd',0.5},  true)
             mk('stabdual_rd1', 'PIE2PDEstability_dual','light','mosek',{'rd',0.5},  false)
             mk('stabpde_rd1',  'PDEstability',         'light','mosek',{'rd',0.5},  true)
             mk('stabpded_rd1', 'PDEstability_dual',    'light','mosek',{'rd',0.5},  false)
             mk('stab_rd1_hv',  'PIE2PDEstability',     'heavy','mosek',{'rd',0.5},  false)
             mk('stab_tr1',     'PIE2PDEstability',     'light','mosek',{'tr',0.5},  false)
             mk('stab_wave1',   'PIE2PDEstability',     'light','mosek',{'wave',0.5},false)
             mk('wellposed_rd1','well_posedness',       'light','mosek',{'syn1'},    false)];
    case '2d'
        C = [mk('stab2_rd',     'stability_2D',         'light','sedumi',{'rd2',0.5},false)
             mk('stab2dual_rd', 'stability_dual_2D',    'light','sedumi',{'rd2',0.5},false)];
    otherwise
        error('cx_cases_stability:which','WHICH is ''1d'' or ''2d''.')
end
end

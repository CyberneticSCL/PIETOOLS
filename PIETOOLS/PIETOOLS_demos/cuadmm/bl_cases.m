function C = bl_cases()
% bl_cases -- the baseline case registry.
%
% The axis is BREADTH OF APPLICATION CLASS, not spatial dimension: stability
% (four variants), Hinf gain (four), H2 norm (four), estimator and controller
% synthesis, well-posedness, nonlinear local stability, then 2-D. 1-D dominates
% deliberately -- it is where a method bug is cheapest to find, and every bug
% found so far has been dimension-independent.
%
% Every plant is written INLINE from the pde_var API. No Examples_Library call:
% those files run evalin('base',...) side effects and at least one has an
% unguarded input() prompt, neither of which survives -batch.
%
% Each builder returns the SOLVED program, because every executive calls
% lpisolve internally. Timing therefore comes from Mosek's own optimizer clock
% (solinfo.info.cpusec, set from MSK_DINF_OPTIMIZER_TIME at sossolve.m:367) and
% from a separate direct re-solve of the dumped bytes -- not from wall time
% around the executive, which would include the PDE->PIE conversion.
%
% ORDERED BY PRIORITY: if the overnight run is cut short, the cases that carry
% the most information have already been measured.

L = {};
add = @(id,cls,dim,kind,bld,args) struct('id',id,'cls',cls,'dim',dim, ...
                                         'kind',kind,'builder',bld,'args',{args});

% ---- 1-D stability: all four executives on one plant, so the four differ only
%      in the analysis, plus two other physics to show the class is not tuned
%      to the Dirichlet Laplacian.
L{end+1} = add('stab_rd1',      'stability','1D','feas', 'bl_b_stab1', {'rd',  0.5,'light','PIE2PDEstability'});
L{end+1} = add('stabdual_rd1',  'stability','1D','feas', 'bl_b_stab1', {'rd',  0.5,'light','PIE2PDEstability_dual'});
L{end+1} = add('stabpde_rd1',   'stability','1D','feas', 'bl_b_stab1', {'rd',  0.5,'light','PDEstability'});
L{end+1} = add('stabpded_rd1',  'stability','1D','feas', 'bl_b_stab1', {'rd',  0.5,'light','PDEstability_dual'});
L{end+1} = add('stab_rd1_hv',   'stability','1D','feas', 'bl_b_stab1', {'rd',  0.5,'heavy','PIE2PDEstability'});
L{end+1} = add('stab_rd1_tight','stability','1D','feas', 'bl_b_stab1', {'rd',  0.95,'heavy','PIE2PDEstability'});
L{end+1} = add('stab_tr1',      'stability','1D','feas', 'bl_b_stab1', {'tr',  0.5,'light','PIE2PDEstability'});
L{end+1} = add('stab_wave1',    'stability','1D','feas', 'bl_b_stab1', {'wave',0.5,'light','PIE2PDEstability'});

% ---- 1-D Hinf gain: primal/dual x coercive/non-coercive is the full 2x2.
L{end+1} = add('hinf_rd1',      'hinf','1D','obj', 'bl_b_io1', {'Hinf_gain',              'light'});
L{end+1} = add('hinfco_rd1',    'hinf','1D','obj', 'bl_b_io1', {'Hinf_gain_coercive',     'light'});
L{end+1} = add('hinfdu_rd1',    'hinf','1D','obj', 'bl_b_io1', {'Hinf_gain_dual',         'light'});
L{end+1} = add('hinfduco_rd1',  'hinf','1D','obj', 'bl_b_io1', {'Hinf_gain_dual_coercive','light'});
L{end+1} = add('hinf_rd1_hv',   'hinf','1D','obj', 'bl_b_io1', {'Hinf_gain',              'heavy'});

% ---- 1-D H2 norm: controllability and observability gramians, both forms.
L{end+1} = add('h2c_rd1',       'h2','1D','obj', 'bl_b_io1', {'H2_norm_c',         'light'});
L{end+1} = add('h2cco_rd1',     'h2','1D','obj', 'bl_b_io1', {'H2_norm_c_coercive','light'});
L{end+1} = add('h2o_rd1',       'h2','1D','obj', 'bl_b_io1', {'H2_norm_o',         'light'});
L{end+1} = add('h2oco_rd1',     'h2','1D','obj', 'bl_b_io1', {'H2_norm_o_coercive','light'});

% ---- 1-D synthesis: these need channels the gain cases do not have, so the
%      plant gains a control input and a sensed output.
L{end+1} = add('est_rd1',       'estimator', '1D','obj', 'bl_b_syn1', {'Hinf_estimator','light'});
L{end+1} = add('h2est_rd1',     'estimator', '1D','obj', 'bl_b_syn1', {'H2_estimator',  'light'});
L{end+1} = add('ctrl_rd1',      'controller','1D','obj', 'bl_b_syn1', {'Hinf_control',  'light'});
L{end+1} = add('h2ctrl_rd1',    'controller','1D','obj', 'bl_b_syn1', {'H2_control',    'light'});
L{end+1} = add('wellposed_rd1', 'wellposed', '1D','feas','bl_b_syn1', {'well_posedness','light'});

% ---- nonlinear local stability, three problem forms at a fixed radius.
L{end+1} = add('nl_fisher_opt',   'nonlinear','1D','obj', 'bl_b_fisher',{4.0,true});
L{end+1} = add('nl_fisher_gamfix','nonlinear','1D','feas','bl_b_fisher',{4.0,2.0});
L{end+1} = add('nl_fisher_nobnd', 'nonlinear','1D','feas','bl_b_fisher',{4.0,false});

% ---- 2-D. Present, not dominant. The Hinf gain plant has a CLOSED-FORM L2
%      gain (Ex_2D_ReactionDiffusion_DDDD), so that row has a ground truth.
L{end+1} = add('stab2_rd',      'stability','2D','feas','bl_b_stab2',{0.5,'light','stability_2D'});
L{end+1} = add('stab2dual_rd',  'stability','2D','feas','bl_b_stab2',{0.5,'light','stability_dual_2D'});
% psatz generators [3;4;5;6] (poslpivar_2d, 0cf1add1): with them Mosek certifies at
% 1e-08 where psatz-off fails, so these two are the regression for that patch.
L{end+1} = add('stab2_rd_psz',  'stability','2D','feas','bl_b_stab2',{0.5,'light','stability_2D',[3;4;5;6]});
L{end+1} = add('stab2_rd_psz90','stability','2D','feas','bl_b_stab2',{0.9,'light','stability_2D',[3;4;5;6]});
L{end+1} = add('hinf2_rd',      'hinf','2D','obj', 'bl_b_io2',  {'Hinf_gain_2D',             'light'});
L{end+1} = add('hinf2nc_rd',    'hinf','2D','obj', 'bl_b_io2',  {'Hinf_gain_2D_non_coercive','light'});
L{end+1} = add('hinf2du_rd',    'hinf','2D','obj', 'bl_b_io2',  {'Hinf_gain_dual_2D',        'light'});
L{end+1} = add('h2_2dc_rd',     'h2','2D','obj',   'bl_b_io2',  {'H2_norm_2D_c',             'light'});

% ---- SCALING ARM. n decoupled states multiplies the block sizes without
%      changing the physics, so m grows with the analysis held fixed. This is
%      the arm that tests the only claim that matters for cuADMM: that it keeps
%      going where an interior-point method runs out of memory.
for n = [2 4 8 16 24 32]
    L{end+1} = add(sprintf('scale_stab_n%02d',n),'scaling','1D','feas', ...
                      'bl_b_stab1',{'rd',0.5,'heavy','PIE2PDEstability',n});
end
for n = [2 4 8]
    L{end+1} = add(sprintf('scale_hinf_n%02d',n),'scaling','1D','obj', ...
                      'bl_b_io1',{'Hinf_gain','heavy',n});
end
C = [L{:}];
end

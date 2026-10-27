function [cases,cases2d] = cx_cases_h2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CASES,CASES2D] = CX_CASES_H2() returns the H2 family's comparison cases.
%
% CASES   1-D, for cx_compare: the six H2 executives at light settings with
%         Mosek, on the baseline plants they were assigned (io1 for the four
%         norm executives, syn1 for synthesis). All are min-gamma ('obj'):
%         the container is posed at fixed gamma and bisected, the stock
%         program is pinned at the same decision variable gam and bisected.
%         For the coercive forms that variable is the SQUARED bound (stock
%         prints sqrt(gam)), so the printed stock gam* is sqrt of the scale
%         the pinned and cx brackets are on; for the others it is the bound.
% CASES2D 2-D, STRUCTURE ONLY (never solved; see cx_h2_struct2d): the two
%         coercive 2-D norm executives on io2 at light, SeDuMi, at a fixed
%         gamma. gam is the decision variable, the squared bound (stock
%         options.h2 branch: gam = options.h2^2, 2D_c L136, 2D_o L138).
%
% Plants carry distributed disturbances, so several cases are degenerate
% by construction (map_h2-norm.json family_notes): c_coercive is
% infeasible at every gamma (measured to gam = 4151 on both paths); o,
% o_coercive and estimator take the trace over an empty finite block, so
% gam there is not an H2 bound. They are kept because they are the assigned
% baseline rows, and the comparison is of SDPs, not of H2 values.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mk = @(id,ex,pl) struct('id',id,'exec',ex,'setname','light','solver','mosek', ...
                        'kind','obj','plant',{pl},'thresh',[]);
cases = [ mk('h2c_io1',    'H2_norm_c',          {'io1'})
          mk('h2cco_io1',  'H2_norm_c_coercive', {'io1'})
          mk('h2o_io1',    'H2_norm_o',          {'io1'})
          mk('h2oco_io1',  'H2_norm_o_coercive', {'io1'})
          mk('h2ctrl_syn1','H2_control',         {'syn1'})
          mk('h2est_syn1', 'H2_estimator',       {'syn1'}) ];

mk2 = @(id,ex) struct('id',id,'exec',ex,'setname','light','solver','sedumi', ...
                      'plant',{{'io2'}},'gam',1);
cases2d = [ mk2('h2_2dc_io2','H2_norm_2D_c')
            mk2('h2_2do_io2','H2_norm_2D_o') ];
end

function lb_validate(stage)
% LB_VALIDATE(STAGE) runs the validation subset of the library baseline
% into results.csv (rows are reused by the full run, which resumes).
%
% INPUT
%   stage:  'A' - six examples on 'light' and 'heavy': two bisections
%                 (PDE 6 on lam, DDE05 on tau), two gains (DDE10, PDE 20),
%                 one NDS stability solve (NDS10) and one 2-D stability
%                 example (PDE 35), with the 2-D skip threshold lowered to
%                 180 s to bound the validation time
%           'B' - tier calibration: single solves of two 1-D rows of
%                 PIE size [1 1] and [2 4] on all six tiers
%
% NOTES
%   - Stage A also runs lb_check_case7 (no solve).

switch stage
    case 'A'
        lb_setup();
        lb_check_case7();
        lb_run({'PDE:6','DDE:DDE05','DDE:DDE10','PDE:20','NDSDDF:NDS10','PDE:35'}, ...
               {'light','heavy'}, struct('max2d',180));
    case 'B'
        lb_run({'DDE:DDE05/stability_dual','NDSDDF:NDS04/stability_dual'}, 'all');
end
end

function R = test_endtoend_sop(which_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_ENDTOEND_SOP(WHICH) runs the Tier 1 end-to-end checks of
% lpi_programming_sopvar/tests (until 10/21/2026 the examples of this
% folder), each of which asserts its own outcomes (see its header):
%
%   'E1'  e2e_hinf_gain_1d      1-D H-infinity gain, gamma a decision
%                               variable in one SDP, against stock
%                               PIETOOLS_Hinf_gain              ~15 s
%   'E2'  e2e_stability_1d      1-D stability, extracted certificate  ~3 s
%   'E3'  e2e_stability_2d      2-D stability, SeDuMi, extraction only
%                               (not a certificate)              ~5 min
%   'E4'  e2e_heat3d_stability  3-D heat at kappa = 14.0, certified;
%                               computes its own row set (sparse QR,
%                               ~6 min) unless given one   ~15-20 min, ~11-21 GB
%
% WHICH: cellstr or comma list (default {'E1','E2'}; 'all' is E1-E4).
% Times measured on the development workstation (Tier 1 README).
%
% Initial coding MMP, 09/29/2026
% MMP, 10/21/2026: E1-E4 moved from examples/ (hinf_gain_1d_sop,
%                stability_1d_sop, stability_2d_sop, heat3d_stability_sop)
%                to this folder as e2e_*: they check, through test helpers;
%                the examples now show standard usage only. Same checks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(which_),     which_ = {'E1','E2'};   end
if ischar(which_),  which_ = strsplit(which_,',');  end
if isequal(which_,{'all'}),     which_ = {'E1','E2','E3','E4'};     end
R = struct();
for k = 1:numel(which_)
    t = tic;
    switch which_{k}
%       case 'E1',  R.E1 = hinf_gain_1d_sop();                              % MMP, 10/21/2026 (was)
%       case 'E2',  R.E2 = stability_1d_sop();                              % MMP, 10/21/2026 (was)
%       case 'E3',  R.E3 = stability_2d_sop();                              % MMP, 10/21/2026 (was)
%       case 'E4',  R.E4 = heat3d_stability_sop();                          % MMP, 10/21/2026 (was)
        case 'E1',  R.E1 = e2e_hinf_gain_1d();                              % MMP, 10/21/2026
        case 'E2',  R.E2 = e2e_stability_1d();                              % MMP, 10/21/2026
        case 'E3',  R.E3 = e2e_stability_2d();                              % MMP, 10/21/2026
        case 'E4',  R.E4 = e2e_heat3d_stability();                          % MMP, 10/21/2026
        otherwise,  error('test_endtoend_sop:which','Unknown example ''%s''.',which_{k})
    end
    fprintf('test_endtoend_sop: %s passed (%.1f s)\n',which_{k},toc(t));
end
end

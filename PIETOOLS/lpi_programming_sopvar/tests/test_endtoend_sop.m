function R = test_endtoend_sop(which_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_ENDTOEND_SOP(WHICH_) runs, each with its default arguments, the
% end-to-end checks of the container LPI path that WHICH_ names, and errors
% at the first failed check:
%   'E1'  e2e_hinf_gain_1d: the H-infinity gain LPI of a 1-D PIE with the
%         gain gamma a decision variable, solved once with MOSEK; gamma
%         against stock PIETOOLS_Hinf_gain, the extracted certificate
%         checked by its action;
%   'E2'  e2e_stability_1d: the stability LPI of a 1-D PIE, solved with
%         MOSEK, the extracted certificate checked by its action;
%   'E3'  e2e_stability_2d: the stability LPI of the 2-D heat equation,
%         solved with SeDuMi, the extracted operators checked against the
%         residual of the solver (not a certificate);
%   'E4'  e2e_heat3d_stability: exponential stability at rate 14.0 of the
%         3-D heat equation, certified with MOSEK on a maximal independent
%         subset of the equality rows (found by one sparse QR), the
%         extracted certificate checked by its action.
%
% INPUT
% - which_: (optional) 'cellstr' subset of {'E1','E2','E3','E4'}, a 'char'
%           comma-separated list of the same names, or 'all' for the four;
%           default {'E1','E2'}.
%
% OUTPUT
% - R:      'struct' with one field per check run (E1, E2, E3, E4), the
%           output 'struct' of that check; one printed line per check
%           passed, with its time.
%
% NOTES
% A name other than 'E1', 'E2', 'E3', 'E4' is an error
% ('test_endtoend_sop:which').
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_endtoend_sop
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/29/2026
% MMP, 10/21/2026: E1-E4 moved from examples/ (hinf_gain_1d_sop,
%                stability_1d_sop, stability_2d_sop, heat3d_stability_sop)
%                to this folder as e2e_*: they check, through test helpers;
%                the examples now show standard usage only. Same checks.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. The deleted header text was
%                the measured run time and memory of each check and a
%                history note.

% Default selection; a comma list or 'all' becomes a cellstr of names
if nargin<1 || isempty(which_),     which_ = {'E1','E2'};   end
if ischar(which_),  which_ = strsplit(which_,',');  end
if isequal(which_,{'all'}),     which_ = {'E1','E2','E3','E4'};     end
R = struct();
% Run each check with its defaults (it errors at its first failure), then
% print its time
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

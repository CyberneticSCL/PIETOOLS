function prog_sol = lpisolve(prog,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG_SOL = LPISOLVE(PROG,OPTS) takes an LPI program structure 'prog', and
% solves the SDP it poses with options 'opts', returning the solved LPI
% program structure.
%
% INPUT
% - prog:       An LPI program structure (based on SOSTOOLS structure)
%               specifying an LPI to solve;
% - opts:       (optional) 'struct' of options to be passed to SOSTOOLS for
%               solving the optimization program specified by lpi; or an
%               LPI settings structure, whose field sos_opts is used.
%               Default: prog.settings.sos_opts, else the solver defaults.
%               An empty field solver is set to 'mosek' when mosekopt is
%               on the path, else to 'sedumi';
%
% OUTPUT
% - prog_sol:   'struct', the solved LPI program;
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% DJ, 10/18/2024 - Initial coding
% MMP, 11/04/2026 - Moved from lpi_programming/ to lpi_programming_sopvar/,
%   code unchanged: it references no operator class and serves both
%   operator families, while lpi_programming/ is retired with the opvar
%   classes.
% MMP, 11/08/2026 - Comment pass (maintainer, 11/08/2026): OUTPUT names
%   prog_sol; code unchanged.
% MMP, 11/12/2026 - Without OPTS, the solver options of the program
%   settings (LPIPROGRAM, LPISETSETTINGS); OPTS may also be a settings
%   structure (maintainer, 11/12/2026: "Settings on the program"). Every
%   tier sets sos_opts.simplify = 0 and no solver, as SOSSOLVE does without
%   options, so a call without OPTS solves the same program as before.
% MMP, 11/18/2026 - An empty sos_opts.solver, which every tier carries since the
%   layout of 11/17/2026 (section 5: '' for MOSEK if found, else SeDuMi), is
%   resolved before SOSSOLVE as EXEC_SOP_SETTINGS does; SOSSOLVE refuses ''.
%   Found by DEMO3 (lpisolve(prob) with the program's tier settings) in the
%   demo pass before the commit of 11/18/2026. The 11/12 entry above is
%   stale on this point: the tiers set a solver field, empty.
% MMP, 11/19/2026: Release comment pass (comment-only).

% % Solve using 'sossolve' with desired settings, if specified.
% if nargin>=2                                                              % MMP, 11/12/2026 (was)
%     prog_sol = sossolve(prog,opts);                                       % MMP, 11/12/2026 (was)
% else                                                                      % MMP, 11/12/2026 (was)
%     prog_sol = sossolve(prog);                                            % MMP, 11/12/2026 (was)
% end                                                                       % MMP, 11/12/2026 (was)
% The options: OPTS; sos_opts of a settings structure OPTS; without OPTS,   % MMP, 11/12/2026
% sos_opts of prog.settings; else none (the solver defaults)                % MMP, 11/12/2026
if nargin<2 || isempty(opts)                                                % MMP, 11/12/2026
    opts = [];                                                              % MMP, 11/12/2026
    if isfield(prog,'settings') && isfield(prog.settings,'sos_opts')        % MMP, 11/12/2026
        opts = prog.settings.sos_opts;                                      % MMP, 11/12/2026
    end                                                                     % MMP, 11/12/2026
elseif isfield(opts,'sos_opts')                                             % MMP, 11/12/2026
    opts = opts.sos_opts;                                                   % MMP, 11/12/2026
end                                                                         % MMP, 11/12/2026
% An empty solver resolved: MOSEK when mosekopt is on the path, else        % MMP, 11/18/2026
% SeDuMi; the solver call refuses ''                                        % MMP, 11/18/2026
if isstruct(opts) && isfield(opts,'solver') && isempty(opts.solver)         % MMP, 11/18/2026
    if ~isempty(which('mosekopt')),     opts.solver = 'mosek';              % MMP, 11/18/2026
    else,                               opts.solver = 'sedumi';             % MMP, 11/18/2026
    end                                                                     % MMP, 11/18/2026
end                                                                         % MMP, 11/18/2026
% Solve the program, with the options when there are any
if isempty(opts),   prog_sol = sossolve(prog);                              % MMP, 11/12/2026
else,               prog_sol = sossolve(prog,opts);                         % MMP, 11/12/2026
end                                                                         % MMP, 11/12/2026

end
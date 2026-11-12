function prog_sol = lpisolve(prog,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG_SOL = LPISOLVE(PROG,OPTS) takes an LPI program structure 'prog', and
% passes this to sossolve to solve with options 'opts', returning the
% solved LPI program structure.
%
% INPUT
% - prog:       An LPI program structure (based on SOSTOOLS structure)
%               specifying an LPI to solve;
% - opts:       (optional) 'struct' of options to be passed to SOSTOOLS for
%               solving the optimization program specified by lpi; or a
%               settings structure (LPISETTINGS_SOP), whose field sos_opts
%               is passed. Default: prog.settings.sos_opts (LPIPROGRAM,
%               LPISETSETTINGS), else the defaults of SOSSOLVE;
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

% % Solve using 'sossolve' with desired settings, if specified.
% if nargin>=2                                                              % MMP, 11/12/2026 (was)
%     prog_sol = sossolve(prog,opts);                                       % MMP, 11/12/2026 (was)
% else                                                                      % MMP, 11/12/2026 (was)
%     prog_sol = sossolve(prog);                                            % MMP, 11/12/2026 (was)
% end                                                                       % MMP, 11/12/2026 (was)
% The options: OPTS; sos_opts of a settings structure OPTS; without OPTS,   % MMP, 11/12/2026
% sos_opts of prog.settings; else none (the SOSSOLVE defaults)              % MMP, 11/12/2026
if nargin<2 || isempty(opts)                                                % MMP, 11/12/2026
    opts = [];                                                              % MMP, 11/12/2026
    if isfield(prog,'settings') && isfield(prog.settings,'sos_opts')        % MMP, 11/12/2026
        opts = prog.settings.sos_opts;                                      % MMP, 11/12/2026
    end                                                                     % MMP, 11/12/2026
elseif isfield(opts,'sos_opts')                                             % MMP, 11/12/2026
    opts = opts.sos_opts;                                                   % MMP, 11/12/2026
end                                                                         % MMP, 11/12/2026
if isempty(opts),   prog_sol = sossolve(prog);                              % MMP, 11/12/2026
else,               prog_sol = sossolve(prog,opts);                         % MMP, 11/12/2026
end                                                                         % MMP, 11/12/2026

end
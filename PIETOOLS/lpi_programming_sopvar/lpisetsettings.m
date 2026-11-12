function prog = lpisetsettings(prog,settings)
% PROG = LPISETSETTINGS(PROG,SETTINGS) sets the settings of the LPI program
% PROG: LPI_INEQ_SOP takes from them the degree options a call does not
% give (settings.lpi.ineq), and LPISOLVE the solver options when a call
% gives none (settings.sos_opts).
%
% INPUT
% - prog:      'struct', the LPI program (LPIPROGRAM);
% - settings:  'struct' of LPI settings (LPISETTINGS_SOP); or a 'char' or
%              'string' naming a tier: 'stripped', 'light', 'heavy',
%              'veryheavy', 'custom' or 'extreme'; or [] (or omitted) for
%              'light';
%
% OUTPUT
% - prog:      'struct', PROG with field settings, the settings structure.
%
% NOTES
% One field can also be set directly, as in prog.settings.lpi.ineq.dw = 2.
% Without a field lpi, LPI_INEQ_SOP reads that of 'light'; without a field
% sos_opts, LPISOLVE leaves the solver options to SOSSOLVE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpisetsettings
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding MMP - 11/12/2026: the settings of a program, so that a
%   hand-written LPI chooses its tier once instead of passing options to
%   every call (maintainer, 11/12/2026: "there should also be an option to
%   change the settings associated with a prog structure"). LPIPROGRAM
%   calls it with its optional last input.

if ~isa(prog,'struct')
    error("LPI program structure should be specified as object of type 'struct'.")
end
if nargin<2 || isempty(settings)
    settings = lpisettings_sop();
elseif ischar(settings) || isstring(settings)
    settings = lpisettings_sop(char(settings));
elseif ~isa(settings,'struct')
    error("Settings should be a 'struct' (LPISETTINGS_SOP) or the name of a tier.")
end
prog.settings = settings;

end

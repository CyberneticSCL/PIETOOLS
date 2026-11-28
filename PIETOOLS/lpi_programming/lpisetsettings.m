function prog = lpisetsettings(prog,settings)
% PROG = LPISETSETTINGS(PROG,SETTINGS) sets the field settings of the LPI
% program PROG: the degree options (settings.lpi.ineq) and the solver
% options (settings.sos_opts) that the LPI routines take when a call gives
% none.
%
% INPUT
% - prog:      'struct', the LPI program;
% - settings:  'struct' of LPI settings; or a 'char' or 'string' naming a
%              tier: 'stripped', 'light', 'heavy', 'veryheavy', 'custom'
%              or 'extreme'; or [] (or omitted) for 'light';
%
% OUTPUT
% - prog:      'struct', PROG with field settings, the settings structure.
%
% NOTES
% One field can also be set directly, as in prog.settings.lpi.ineq.dw = 2.
% Without a field lpi, the degree options of 'light' apply; without a field
% sos_opts, the default solver options apply.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

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
% Initial coding MMP - 11/12/2026

% Check PROG
if ~isa(prog,'struct')
    error("LPI program structure should be specified as object of type 'struct'.")
end
% The settings: the 'light' tier when omitted or empty, a tier by name, or
% a 'struct' as given
if nargin<2 || isempty(settings)
    settings = lpisettings();
elseif ischar(settings) || isstring(settings)
    settings = lpisettings(char(settings));
elseif ~isa(settings,'struct')
    error("Settings should be a 'struct' (LPISETTINGS) or the name of a tier.")
end
% Store them in PROG
prog.settings = settings;

end

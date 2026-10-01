function prog = lpi_eq_sop(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_EQ_SOP(PROG,P,OPTS) takes an LPI optimization program
% structure 'prog' and a (PI operator) decision variable P of either
% operator family, and adds equality constraints enforcing P==0. It is
% 'lpi_eq' with the 'sopvar' family added, so that one LPI body serves
% both:
%
%   class of P                              handled by
%   'cdopvar', 'copvar'                     lpi_eq_cdopvar
%   'sdopvar', 'sopvar'                     lpi_eq_sdopvar
%   anything else ('dpvar', 'dopvar',       lpi_eq, unchanged
%   'dopvar2d', 'opvar', 'opvar2d',
%   'polynomial', 'double', ...)
%
% Each routine keeps its own checks and errors: a fixed 'copvar' or
% 'sopvar' is refused by the container routines (no decision variables),
% and 'polynomial'/'double' by 'lpi_eq', as each does today.
%
% INPUT
% - prog:   'struct' specifying the LPI program structure to modify (see
%           'lpiprogram_sop' or 'lpiprogram');
% - P:      the operator or decision variable to set to zero;
% - opts:   (optional) 'symmetric' if P is self-adjoint, to constrain one
%           parameter (container: one block) per adjoint pair; passed to
%           the routine unchanged, and only when given.
%
% OUTPUT
% - prog:   the program with the constraints P==0 added: the same program,
%           field for field, that the routine in the table above builds
%           when called directly (tests/test_lpi_eq_sop).
%
% NOTES
% Cost: two to four class tests, then the routine's own cost. The
% container routines impose the rows in joined 'soseq' calls, about one
% per numel(prog.decvartable) nonzeros (see 'lpi_eq_sdopvar').
%
% See also LPI_EQ, LPI_EQ_CDOPVAR, LPI_EQ_SDOPVAR, LPIPROGRAM_SOP,
%          LPIGETSOL_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_eq_sop
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
% Initial coding MMP, 09/29/2026. Dispatcher over both operator families,
%                in a folder parallel to lpi_programming and suffixed _sop
%                so that it shadows nothing (every folder is on the path).
%                'lpi_eq' refuses the new classes ("Input must be of type
%                'dopvar', 'dopvar2d', or 'dpvar'"), so an LPI body had to
%                call lpi_eq_cdopvar or lpi_eq_sdopvar by class.
% MMP, 09/30/2026: Comment only: lpi_eq_sdopvar has one output and three
%                inputs since its modes were split into collect_eq_rows and
%                impose_eq_rows; the note on a second output is obsolete.

% Only opts the caller gave is passed on: lpi_eq, lpi_eq_cdopvar and
% lpi_eq_sdopvar each test nargin for it.
if isa(P,'cdopvar') || isa(P,'copvar')
    if nargin>=3,   prog = lpi_eq_cdopvar(prog,P,opts);
    else,           prog = lpi_eq_cdopvar(prog,P);
    end
elseif isa(P,'sdopvar') || isa(P,'sopvar')
%   % One output: lpi_eq_sdopvar imposes the rows (a second output would    % MMP, 09/30/2026 (was)
%   % return them unimposed).                                               % MMP, 09/30/2026 (was)
    % lpi_eq_sdopvar collects and imposes the rows.                         % MMP, 09/30/2026
    if nargin>=3,   prog = lpi_eq_sdopvar(prog,P,opts);
    else,           prog = lpi_eq_sdopvar(prog,P);
    end
else
    if nargin>=3,   prog = lpi_eq(prog,P,opts);
    else,           prog = lpi_eq(prog,P);
    end
end

end

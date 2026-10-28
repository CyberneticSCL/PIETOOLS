function [prog,tag] = lpi_eq_sop(prog,P,opts)                               % MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,TAG] = LPI_EQ_SOP(PROG,P,OPTS) adds to the LPI program PROG the
% equality constraints P == 0 on the decision variables of P: for a
% 'cdopvar' P, that every block is zero (LPI_EQ_CDOPVAR); for an 'sdopvar'
% P, that every stored coefficient of every parameter is zero
% (LPI_EQ_SDOPVAR); for any other P, the constraints of LPI_EQ.
%
% INPUT
% - prog:   'struct', the LPI program to modify (from LPIPROGRAM_SOP or
%           LPIPROGRAM), whose field decvartable holds every decision
%           variable of P;
% - P:      'cdopvar' or 'sdopvar' object, the decision operator to set to
%           zero; or a 'dpvar', 'opvar', 'dopvar', 'opvar2d' or 'dopvar2d'
%           object;
% - opts:   (optional) char 'symmetric', for a self-adjoint P: one
%           parameter of each adjoint pair is constrained (for a 'cdopvar',
%           one block of each pair of blocks (i,j) and (j,i), and one
%           parameter of each pair within a diagonal block). For a P that
%           is not self-adjoint the constraints are then weaker than
%           P == 0.
%
% OUTPUT
% - prog:   'struct', PROG with the constraints P == 0 added;
% - tag:    'double', for a 'cdopvar' P and a PROG with the field sopeq (a
%           cell): the index in prog.sopeq of the record of the rows this
%           call added (the TAG input of LPIGETDUAL_SOP); 0 otherwise.
%
% NOTES
% A 'copvar', 'sopvar', 'polynomial' or 'double' P, which has no decision
% variables, is an error.
%
% See also LPI_EQ, LPI_EQ_CDOPVAR, LPI_EQ_SDOPVAR, LPIPROGRAM_SOP,
%          LPIGETSOL_SOP, LPIGETDUAL_SOP.
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
% MMP, 10/02/2026: Comment only: the container rows are written by
%                'lpi_soseq', not soseq.
% MMP, 10/08/2026: second output TAG, the index in prog.sopeq of the row
%                  record lpi_eq_cdopvar keeps for this equality, so that
%                  'lpigetdual_sop' can assemble its dual kernel; 0 on the
%                  single-block and legacy paths, which keep no record.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body; code
%                unchanged. The 10/08/2026 entry, which stood in the body,
%                moved to this log. The deleted header text was the
%                comparison with lpi_eq, the per-routine remarks on errors,
%                a test reference and a cost note with its 10/02/2026 edit.

% Only opts the caller gave is passed on: lpi_eq, lpi_eq_cdopvar and
% lpi_eq_sdopvar each test nargin for it.
% TAG is 0 unless lpi_eq_cdopvar records the rows of this equality
tag = 0;                                                                    % MMP, 10/08/2026
if isa(P,'cdopvar') || isa(P,'copvar')
    % Container: every block, by lpi_eq_cdopvar; TAG the index of the row
    % record it appends to prog.sopeq when the program has that field
    if nargin>=3,   prog = lpi_eq_cdopvar(prog,P,opts);
    else,           prog = lpi_eq_cdopvar(prog,P);
    end
    if isfield(prog,'sopeq'),   tag = numel(prog.sopeq);    end              % MMP, 10/08/2026
elseif isa(P,'sdopvar') || isa(P,'sopvar')
%   % One output: lpi_eq_sdopvar imposes the rows (a second output would    % MMP, 09/30/2026 (was)
%   % return them unimposed).                                               % MMP, 09/30/2026 (was)
    % lpi_eq_sdopvar collects and imposes the rows.                         % MMP, 09/30/2026
    if nargin>=3,   prog = lpi_eq_sdopvar(prog,P,opts);
    else,           prog = lpi_eq_sdopvar(prog,P);
    end
else
    % Every other class: lpi_eq
    if nargin>=3,   prog = lpi_eq(prog,P,opts);
    else,           prog = lpi_eq(prog,P);
    end
end

end

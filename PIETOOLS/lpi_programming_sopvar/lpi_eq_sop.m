function [prog,tag] = lpi_eq_sop(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,TAG] = LPI_EQ_SOP(PROG,P,OPTS) adds to the LPI program PROG the
% equality constraints P == 0 on the decision variables of P: for a
% 'cdopvar' P, that every block is zero; for an 'sdopvar' P, that every
% stored coefficient of every parameter is zero; for a 'dpvar' P, P == 0 as
% an equality of the SOS program; any other P is an error.
%
% INPUT
% - prog:   'struct', the LPI program to modify, whose field decvartable
%           holds every decision variable of P;
% - P:      'cdopvar' or 'sdopvar' object, the decision operator to set to
%           zero; or a 'dpvar' object;
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
%           call added; 0 otherwise.
%
% NOTES
% A 'copvar', 'sopvar', 'polynomial' or 'double' P, which has no decision
% variables, is an error, as is a P of any other class
% ('lpi_eq_sop:class').
%
% See also SOSEQ, LPI_EQ_CDOPVAR, LPI_EQ_SDOPVAR, LPIPROGRAM,
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
% Initial coding MMP, 09/29/2026

% OPTS is passed on only when given; TAG is 0 unless the rows of this
% equality are recorded in prog.sopeq
tag = 0;
if isa(P,'cdopvar') || isa(P,'copvar')
    % Container: every block zero; TAG the index of the row record appended
    % to prog.sopeq when the program has that field
    if nargin>=3,   prog = lpi_eq_cdopvar(prog,P,opts);
    else,           prog = lpi_eq_cdopvar(prog,P);
    end
    if isfield(prog,'sopeq'),   tag = numel(prog.sopeq);    end
elseif isa(P,'sdopvar') || isa(P,'sopvar')
    % 'sdopvar': every stored coefficient of every parameter zero
    if nargin>=3,   prog = lpi_eq_sdopvar(prog,P,opts);
    else,           prog = lpi_eq_sdopvar(prog,P);
    end
elseif isa(P,'dpvar')
    % 'dpvar': P == 0 as an equality of the SOS program
    prog = soseq(prog,P);
else
    % Every other class: an error
    error('lpi_eq_sop:class','P should be a cdopvar, sdopvar or dpvar object; got %s.',class(P))
end

end

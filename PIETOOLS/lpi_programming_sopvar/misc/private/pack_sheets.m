function CG = pack_sheets(V,mrow,ncol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CG = PACK_SHEETS(V,MROW,NCOL) reshapes each row of V into an MROW x NCOL
% matrix and lays the results out side by side: sheet sg,
% reshape(V(sg,:),MROW,NCOL), occupies columns (sg-1)*NCOL+1 : sg*NCOL of
% CG.
%
% INPUT
% - V:      'double', nsheet x (mrow*ncol), sparse or full, row sg the
%           column-major vec of sheet sg;
% - mrow:   'double', the number of rows of a sheet;
% - ncol:   'double', the number of columns of a sheet.
%
% OUTPUT
% - CG:     'double', sparse, mrow x (nsheet*ncol).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pack_sheets
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. No help text deleted. The
%                file had no log: initial coding MP, 08/22/2026, as a local
%                function of 'possopvar' (the file since renamed
%                'sopquadvar'), moved to private/ by MMP, 09/21/2026.

% Number of sheets; the nonzeros of V by sheet and position in the sheet
nsheet = size(V,1);
[sg,lin,val] = find(V);
% Row and column of each nonzero within its sheet (column-major)
irow = mod(lin-1,mrow)+1;
icol = floor((lin-1)/mrow)+1;
% Sheet sg shifted right by (sg-1)*ncol columns
CG = sparse(irow,(sg-1)*ncol+icol,val,mrow,nsheet*ncol);

end



function Cup = liftIndex(C, m, n, mapS, mapT, dsU, dtU)
%LIFTINDEX Lift coefficient matrix by reindexing nnz entries using maps.
% Lifts C such that 
%       (I_m\otimes Zs_old')*C*(I_n\otimes Zt_old) 
%               = (I_m \otimes Zs_union')*Cup*(I_n \otimes Zt_union)
% mapS: dsOld×1 old->union map (values in 1..dsU)
% i.e., Zs_old = Zs_union(mapS)
% mapT: dtOld×1 old->union map (values in 1..dtU)
% i.e., Zt_old = Zt_union(mapT)
%
% INPUT
% - C:          m*dsOld x n*dtOld sparse coefficient matrix;
% - m, n:       integers, the matrix dimensions;
% - mapS, mapT: the old-to-union maps above;
% - dsU, dtU:   integers, the sizes of the union bases;
%
% OUTPUT
% - Cup:        m*dsU x n*dtU sparse coefficient matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - liftIndex
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding SS - 01/29/2026

C = sparse(C);

dsOld = numel(mapS);
dtOld = numel(mapT);

[i, j, v] = find(C);

a  = mod(i(:)-1, dsOld) + 1;
ib = floor((i(:)-1) / dsOld);

b  = mod(j(:)-1, dtOld) + 1;
jb = floor((j(:)-1) / dtOld);

i2 = ib*dsU + mapS(a);
j2 = jb*dtU + mapT(b);

Cup = sparse(i2, j2, v(:), m*dsU, n*dtU);
end

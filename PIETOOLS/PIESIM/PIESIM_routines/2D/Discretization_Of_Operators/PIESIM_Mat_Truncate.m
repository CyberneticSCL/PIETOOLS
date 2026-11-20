%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Mat_Truncate.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncates matrix to the correct dimension by skipping the columns of rows
% corresponding to the unresolved coefficients of differentiated PDE states comparative to full 
% PDE states that appear in the fundamental state
%
% Index k_x + (N_x+1)*k_y + 1 of a 2D coefficient vector (x index fastest)
% is kept if k_x <= N_x - p(1) and k_y <= N_y - p(2).
%
% Inputs:
% A - matrix to be truncated 
% N   - 1x2 vector [N_x, N_y], polynomial orders of the Chebyshev
%       discretization in x and y
% p - 2-vector [p_x; p_y], degree of differentiability corresponding to the
% column states (if dir='col') or row states (if dir='row') % 
% dir ('col' or 'row') - direction of truncation (columns or rows)
%
% Outputs:
% Atrun - truncated matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_Mat_Truncate
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 1/31/2026

function Atrun = PIESIM_Mat_Truncate(A,N,p,dir)

% Keep the first N_x+1-p(1) of each block of N_x+1 x-coefficients, in the
% first N_y+1-p(2) blocks
MatrixSize=prod(N+1);
nSkip=p(1);
nKeep=N(1)+1-p(1);
nBlocksKeep = N(2)+1-p(2);
pattern = [true(nKeep,1); false(nSkip,1)];
keep = repmat(pattern, nBlocksKeep, 1);
% automatically drop the last p blocks, where p is degree of
% differentiability of column states
keep = [keep; false(MatrixSize - length(keep), 1)];

% Select the kept columns or rows
if (dir=='col')
Atrun = A(:, keep);
else
Atrun = A(keep,:);
end

end
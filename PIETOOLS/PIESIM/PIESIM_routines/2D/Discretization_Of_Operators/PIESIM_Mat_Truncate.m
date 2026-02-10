%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Mat_Truncate.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncates matrix to the correct dimension by skipping the columns of rows
% corresponding to the unresolved coefficients of differentiated PDE states comparative to full 
% PDE states that appear in the fundamental state
%
% Inputs:
% A - matrix to be truncated 
% N   - polynomial order of Chebyshev discretization polynomial
% p - degree of differentiability corresponding to the
% column states (if dir='col') or row states (if dir='row') % 
% dir ('col' or 'row') - direction of truncation (columns or rows)

% Outputs:
% Atrun - truncated matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 1/31/2026

function Atrun = PIESIM_Mat_Truncate(A,N,p,dir)

MatrixSize=prod(N+1);
nSkip=p(1);
nKeep=N(1)+1-p(1);
nBlocksKeep = N(2)+1-p(2);
pattern = [true(nKeep,1); false(nSkip,1)];
keep = repmat(pattern, nBlocksKeep, 1);
% automatically drop the last p blocks, where p is degree of
% differentiability of column states
keep = [keep; false(MatrixSize - length(keep), 1)];

if (dir=='col')
Atrun = A(:, keep);
else
Atrun = A(keep,:);
end

end
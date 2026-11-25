function P_reshaped = reshape(P1, M, N)
%RESHAPE  Reshape a quadPoly object.
%
%   P_RESHAPED = RESHAPE(P1,M,N) returns the M-by-N quadPoly whose
%   elements are taken columnwise from P1. P1 must have M*N elements.
%
% INPUTS
% - P1:         m x n 'quadPoly' object (or 'double' array);
% - M, N:       positive integers with M*N = m*n; the form
%               RESHAPE(P1,[M,N]) is not supported;
%
% OUTPUTS
% - P_reshaped: M x N 'quadPoly' object, the entries of P1 taken
%               columnwise, with the same monomial bases ('double' for a
%               'double' P1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - reshape
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
% AT, 02/18/2026: Initial coding;
 

% A 'double' array is reshaped by the builtin reshape
if isa(P1, 'double')
    P_reshaped = reshape(P1, M, N);
    return
end

if ~isa(P1, 'quadPoly')
    error('quadPoly:reshape:badType', 'reshape supports quadPoly or numeric inputs.');
end


% The number of elements is kept
sz = P1.dim;

if sz(1)*sz(2) ~= M*N
    error('quadPoly:reshape:dimMismatch', 'Number of elements must not change.');
end

% extract monomials and coefficients.
Zt = P1.Zt;
Zs = P1.Zs;
nt = P1.nt;
ns = P1.ns;
C = P1.C;


n_Zt = 1;% length of right monomials
for idx = 1:size(Zt,2)
    n_Zt = n_Zt*size(Zt{idx}, 1);
end

n_Zs = 1;% length of left monomials
for idx = 1:size(Zs,2)
    n_Zs = n_Zs*size(Zs{idx}, 1);
end

[i, j, v] = find(C);   % find nonzero entries

% Within-block indices and block indices
a  = mod(i-1, n_Zs) + 1;        % 1..n_Zs
ib = floor((i-1) / n_Zs);       % 0..m-1

b  = mod(j-1, n_Zt) + 1;        % 1..n_Zt
jb = floor((j-1) / n_Zt);       % 0..n-1

full_index = ib + jb*sz(1); % numbering of blocks

iblc2 = mod(full_index, M); %new row indx of blocks
jblc2 = floor(full_index/M); %new col indx of blocks

i2 = iblc2*n_Zs + a;             % rows: M blocks of size n_Zs
j2 = jblc2*n_Zt + b;             % cols: N blocks of size n_Zt

Cnew = sparse(i2, j2, v, M*n_Zs, N*n_Zt);

% define reshaped quadpoly 
P_reshaped = quadPoly(Cnew, Zs, Zt, [M, N], ns, nt);

end


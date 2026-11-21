%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fcgltran4d.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - fcgltran4d
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
% Initial coding YP - 07/05/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.

function B=fcgltran4d(A)
% 4D Fast Chebyshev Transform from nodal to spectral
%
% Performs the fast transform of data sampled at the
% Chebyshev-Gauss-Lobatto Nodes x=cos(k*pi/N);
%
% A - original data in columns
% B - transformed data in columns
%
% Written by Yulia Peet 07/05/22  
% Contact: ypeet@asu.edu
[N,M,P,Q]=size(A);
    F=ifft([A(1:N,:,:,:);A(N-1:-1:2,:,:,:)]);
    B=real([F(1,:,:,:); 2*F(2:N,:,:,:)]);
    G=permute(B,[2 1 3 4]);
    F=ifft([G(1:M,:,:,:);G(M-1:-1:2,:,:,:)]);
    B=permute(real([F(1,:,:,:); 2*F(2:M,:,:,:)]),[2 1 3 4]);
    G=permute(B,[3 2 1 4]);
    F=ifft([G(1:P,:,:,:);G(P-1:-1:2,:,:,:)]);
    B=permute(real([F(1,:,:,:); 2*F(2:P,:,:,:)]),[3 2 1 4]);
    G=permute(B,[4 2 3 1]);
    F=ifft([G(1:Q,:,:,:);G(Q-1:-1:2,:,:,:)]);
    B=permute(real([F(1,:,:,:); 2*F(2:Q,:,:,:)]),[4 2 3 1]);
end

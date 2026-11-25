function P_tr = trace(P1)
%TRACE  Sum of diagonal elements.
%   P_TR = TRACE(P1) is the sum of the diagonal elements of P1.
%
% INPUTS
% - P1:     n x n 'quadPoly' object, or a 'double' matrix;
%
% OUTPUTS
% - P_tr:   1 x 1 'quadPoly' object sum_i P1(i,i), in the monomial bases
%           of P1; for a 'double' P1, the 'double' trace(P1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - trace
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
 

if isa(P1, 'double')
    P_tr = trace(P1);
    return
end

if ~isa(P1, 'quadPoly')
    error('quadPoly:trace:badType', 'trace supports quadPoly and numeric inputs.');
end


sz = P1.dim;
if sz(1) ~= sz(2)
    error('quadPoly:trace:dimMismatch', 'quadPoly must be square.');
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


% compute C matrix 
C_tr = C(1:n_Zs, 1:n_Zt);
for idx = 2:sz(1)
    left_idx  = n_Zs*(idx - 1)+1;
    right_idx = n_Zt*(idx - 1)+1;
    C_tr = C_tr + C(left_idx:(left_idx + n_Zs-1), right_idx:(right_idx+n_Zt-1));
end

% define trace quadpoly 
P_tr = quadPoly(C_tr, Zs, Zt, [1, 1], ns, nt);

end


function P_repmat = repmat(P1, N, M)
% P_REPMAT = REPMAT(P1,N,M) tiles the quadPoly P1 N times vertically and
% M times horizontally, keeping its monomial bases.
%
% INPUTS
% - P1:         m x n 'quadPoly' object;
% - N, M:       (optional) positive integers, M = N if omitted; without
%               N and M, P1 is returned. The form REPMAT(P1,[N,M]) is
%               not supported;
%
% OUTPUTS
% - P_repmat:   m*N x n*M 'quadPoly' object, the N x M tiling of P1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - repmat
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

% Only 'quadPoly' objects; P1 itself without N, M; M = N by default
if ~isa(P1, 'quadPoly')
    error('quadPoly:repmat:badType', 'repmat supports quadPoly only.');
end

if nargin==1
    P_repmat=P1;
    return 
end



if nargin==2
    M = N;
end
 


% extract monomials and coefficients.
Zt = P1.Zt;
Zs = P1.Zs;
nt = P1.nt;
ns = P1.ns;
C = P1.C;

% Tile the coefficients (blocks of rows m*ds, columns n*dt), same bases
C_repmat = repmat(C, N, M); 

P_repmat = quadPoly(C_repmat, Zs, Zt, P1.dim.*[N, M], ns, nt); 

end
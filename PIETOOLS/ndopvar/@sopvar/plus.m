function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = plus(A,B) adds two sopvar operators A,B: L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% B: sopvar class object L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% Note: A,B must have identical dimensions and domains
%
% OUTPUT
% C = A+B:  L_2^q[Sj,...,Sm] to L_2^p[Si,...Sn]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - plus
%
% Copyright (C)2026  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 1_16_2026
% Allowing different monomials AT - 05/26/26
% Switch order ZL and ZR, DJ 06/08/2026

% initialize added output class

% Error handling: Checks to ensure A and B are compatible
if any(A.dims~=B.dims)
    error('Summands A and B have different dimensions');
end

if any(~strcmp(A.vars.in,B.vars.in)) ||any(~strcmp(A.vars.out,B.vars.out))
    error('Summands A and B map between different spaces');
end

if any(any(A.dom.in~=B.dom.in)) || any(any(A.dom.out~=B.dom.out))
    error('input or output variables in summands A and B have different domains');
end

if numel(A.params)~=numel(B.params)
    error('number of terms in summands is not equal -- one of them is probably malformed');
end
% dimensions of params in A and B, vars_in and vars_out


ZL1 = A.ZL;
ZL2 = B.ZL;
ZR1 = A.ZR;
ZR2 = B.ZR;

%find common basis
[ZR, C1R, C2R] = UnionBasisMonomials(ZR1, ZR2);
[ZL, C1L, C2L] = UnionBasisMonomials(ZL1, ZL2);


C1L = kron(eye(A.dims(1)), C1L);
C2L = kron(eye(B.dims(1)), C2L);
C1R = kron(eye(A.dims(2)), C1R);
C2R = kron(eye(B.dims(2)), C2R);



params = A.params;
for i=1:numel(A.params)  % linear indexing of multi-dimensional cell array
    params{i} = C1L'*A.params{i}*C1R +C2L'*B.params{i}*C2R;  % adding quadpoly objects. % Note: chatgpt says this approach is faster than cellfun
end


C = sopvar(params, A.vars, ZR, ZL, A.dom, A.dims);

end

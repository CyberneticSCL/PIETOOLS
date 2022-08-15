function [P] = ctranspose(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P] = ctranspose(P) transposes an operator 
% P: R^m0 x L2^mx x L2^my x L2^m2 to R^n0 x L2^nx x L2^ny x L2^n2
% Date: 02/04/21
% Version: 1.0
% 
% INPUT
% P: opvar2d class object
% 
% OUTPUT
% P: transpose of the input opvar2d with the same matlab structure as P
% 
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - ctranspose
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 02_04_2021  
%   ^ Based heavily on "@opvar"-ctranspose code by SS ^

if ~isa(P,'opvar2d')
    error('Input must be an opvar2d variable.');
end

ds = P.var1;
dt = P.var2;

R0x = P.R0x;    R0y = P.R0y;    R02 = P.R02;
Rxx = P.Rxx;    Rxy = P.Rxy;    Rx2 = P.Rx2;
                Ryy = P.Ryy;    Ry2 = P.Ry2;
R2x = P.R2x;    R2y = P.R2y;    R22 = P.R22;

P.R00 = P.R00';
P.R0x = P.Rx0';
P.R0y = P.Ry0';
P.R02 = P.R20';

P.Rx0 = R0x';
P.Rxx{1,1} = Rxx{1,1}';
P.Rxy = P.Ryx';
P.Rx2{1,1} = R2x{1,1}';

P.Ry0 = R0y';
P.Ryx = Rxy';
P.Ryy{1,1} = Ryy{1,1}';
P.Ry2{1,1} = R2y{1,1}';

P.R20 = R02';
P.R2x{1,1} = Rx2{1,1}';
P.R2y{1,1} = Ry2{1,1}';
P.R22{1,1} = R22{1,1}';

for i=2:3
    P.Rxx{i,1} = var_swap(Rxx{5-i,1}',ds(1),dt(1));
    P.Rx2{i,1} = var_swap(R2x{5-i,1}',ds(1),dt(1));
    P.R2x{i,1} = var_swap(Rx2{5-i,1}',ds(1),dt(1));
    
    P.Ryy{1,i} = var_swap(Ryy{1,5-i}',ds(2),dt(2));
    P.Ry2{1,i} = var_swap(R2y{1,5-i}',ds(2),dt(2));
    P.R2y{1,i} = var_swap(Ry2{1,5-i}',ds(2),dt(2));
    
    P.R22{i,1} = var_swap(R22{5-i,1}',ds(1),dt(1));
    P.R22{1,i} = var_swap(R22{1,5-i}',ds(2),dt(2));
    for j=2:3
        P.R22{i,j} = var_swap(var_swap(R22{5-i,5-j}',ds(1),dt(1)),ds(2),dt(2));
    end
end

end
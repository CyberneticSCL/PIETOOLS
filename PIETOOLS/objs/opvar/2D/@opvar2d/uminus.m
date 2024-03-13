function [P] = uminus(P) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P] = uminus(P) performs unitary minus of an operator 
% P: R^m0 x L2^mx x L2^my x L2^m2 to R^n0 x L2^nx x L2^ny x L2^n2
% Date: 02/01/21
% Version: 1.0
% 
% INPUTS
% P: opvar2d class object
% 
% OUTPUTS
% P: minus of P
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - uminus
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 02_01_2021  
%   ^ Based heavily on "@opvar"-uminus code by SS ^
%

if ~isa(P,'opvar2d')
    error('Input must be an opvar2d variable.');
end

P.R00 = -P.R00; P.R0x = -P.R0x; P.R0y = -P.R0y; P.R02 = -P.R02;
P.Rx0 = -P.Rx0;                 P.Rxy = -P.Rxy;
P.Ry0 = -P.Ry0; P.Ryx = -P.Ryx;
P.R20 = -P.R20;

for i=1:3
    P.Rxx{i,1} = -P.Rxx{i,1};
    P.Rx2{i,1} = -P.Rx2{i,1};
    P.R2x{i,1} = -P.R2x{i,1};
    
    P.Ryy{1,i} = -P.Ryy{1,i};
    P.Ry2{1,i} = -P.Ry2{1,i};
    P.R2y{1,i} = -P.R2y{1,i};
    
    for j=1:3
        P.R22{i,j} = -P.R22{i,j};
    end
end

end
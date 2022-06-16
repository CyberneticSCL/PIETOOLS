function PIE = rescalePIE(PIE,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = rescalePIE(PIE,I) takes in a PIE structure and spatial domain 
% I = [c,d], and outputs an equivalent PIE structure on I. 
%
% INPUT:
%
% T: a 4-PI operator
% [c,d]: interval, default [-1,1]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - rescalePIE
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
% Initial coding YP  - 10_06_2021
% YP 6/16/2022 - Separated PDE LHS operator (T) from the state map operator
% (T0)

if nargin==1
    I = [-1,1];
end

PIE.T = transl_mod(PIE.T,I);
if isfield(PIE,'T0')
PIE.T0 = transl_mod(PIE.T0,I);
end
PIE.Tw = transl_mod(PIE.Tw,I);
PIE.Tu = transl_mod(PIE.Tu,I);
PIE.A = transl_mod(PIE.A,I);
PIE.B1 = transl_mod(PIE.B1,I);
PIE.B2 = transl_mod(PIE.B2,I);
PIE.C1 = transl_mod(PIE.C1,I);
PIE.C2 = transl_mod(PIE.C2,I);
PIE.D11 = transl_mod(PIE.D11,I);
PIE.D12 = transl_mod(PIE.D12,I);
PIE.D21 = transl_mod(PIE.D21,I);
PIE.D22 = transl_mod(PIE.D22,I);

if isfield(PIE,'L')
PIE.L = transl_mod(PIE.L,I);
end
if isfield(PIE,'K')
PIE.K = transl_mod(PIE.K,I);
end
end
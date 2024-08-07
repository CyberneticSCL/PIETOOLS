function PIE = rescalePIE_2D(PIE,dom,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = rescalePIE_2D(PIE,I) takes in a PIE structure and spatial domain 
% I = [a,b;c,d], and outputs an equivalent PIE structure on I domain. 
%
% INPUT:
%
% T: a 2D PI operator
% dom: original domain
% I: new domain
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% Y. Peet at ypeet@asu.edu, S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - rescalePIE_2D
%
% Copyright (C)2024  M. Peet, Y. Peet, S. Shivakumar, D. Jagt
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
% Initial coding YP  - 06_23_2022

if nargin<3
    I = [-1,1;-1,1];
elseif nargin>3
    error("Incorrect number of inputs; only 3 inputs are allowed");
end

PIE.T = transl_2D(PIE.T,dom,I);
PIE.A = transl_2D(PIE.A,dom,I);

if isfield(PIE,'T0')
PIE.T0 = transl_2D(PIE.T0,dom,I);
end
PIE.Tw = transl_2D(PIE.Tw,dom,I);
PIE.Tu = transl_2D(PIE.Tu,dom,I);
PIE.B1 = transl_2D(PIE.B1,dom,I);
PIE.B2 = transl_2D(PIE.B2,dom,I);
PIE.C1 = transl_2D(PIE.C1,dom,I);
PIE.C2 = transl_2D(PIE.C2,dom,I);
PIE.D11 = transl_2D(PIE.D11,dom,I);
PIE.D12 = transl_2D(PIE.D12,dom,I);
PIE.D21 = transl_2D(PIE.D21,dom,I);
PIE.D22 = transl_2D(PIE.D22,dom,I);

if isfield(PIE,'L')
PIE.L = transl_2D(PIE.L,dom,I);
end
if isfield(PIE,'K')
PIE.K =transl_2D(PIE.K,dom,I);
end

end
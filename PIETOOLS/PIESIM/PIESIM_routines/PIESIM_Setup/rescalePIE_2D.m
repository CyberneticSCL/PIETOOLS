%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rescalePIE_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP - 4_16_2024

function PIE = rescalePIE_2D(PIE,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = rescalePIE_2D(PIE,I) takes in a PIE structure and spatial domain 
% I = [a,b;c,d], and outputs an equivalent PIE structure on I domain. 
%
% INPUT:
%
% PIE: a PIE structure, with 2D PI operators defined on an original domain
%      (T, Tw, Tu, A, B1, B2, C1, C2, D11, D12, D21, D22, and
%      PIE.misc.Tmap, Tumap, Twmap, PIE.L, PIE.K, if present)
% I: new domain (2x2), default [-1,1;-1,1]
%
% OUTPUT:
%
% PIE: the PIE structure with each operator translated to I
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
% YP 6/26/2025 - modified support for Tmap operator with new pie_struct
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: a given domain I is
%   used (was replaced by the default whenever nargin < 3); more than 2
%   inputs is an error; the maps PIE.misc.Tmap, Tumap, Twmap are translated
%   (was PIE.misc.T0); help and comments revised.


% Default domain [-1,1]^2
if nargin<2
    I = [-1,1;-1,1];
elseif nargin>2
    error("Incorrect number of inputs; only 2 inputs are allowed");
end

PIE.T = transl_2D(PIE.T,I);
PIE.A = transl_2D(PIE.A,I);

% Maps from the fundamental to the primary states, if present
if isfield(PIE.misc,'Tmap')
PIE.misc.Tmap = transl_2D(PIE.misc.Tmap,I);
end
if isfield(PIE.misc,'Tumap')
PIE.misc.Tumap = transl_2D(PIE.misc.Tumap,I);
end
if isfield(PIE.misc,'Twmap')
PIE.misc.Twmap = transl_2D(PIE.misc.Twmap,I);
end
% Remaining operators of the PIE
PIE.Tw = transl_2D(PIE.Tw,I);
PIE.Tu = transl_2D(PIE.Tu,I);
PIE.B1 = transl_2D(PIE.B1,I);
PIE.B2 = transl_2D(PIE.B2,I);
PIE.C1 = transl_2D(PIE.C1,I);
PIE.C2 = transl_2D(PIE.C2,I);
PIE.D11 = transl_2D(PIE.D11,I);
PIE.D12 = transl_2D(PIE.D12,I);
PIE.D21 = transl_2D(PIE.D21,I);
PIE.D22 = transl_2D(PIE.D22,I);

% Operators L and K, if present
if isfield(PIE,'L')
PIE.L = transl_2D(PIE.L,I);
end
if isfield(PIE,'K')
PIE.K =transl_2D(PIE.K,I);
end

end

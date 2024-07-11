function sos = lpi_eq_2d(sos,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sos = lpi_eq_2d(prog,P) sets up equality constraints for each component.
% P.R00 = 0     P.R0x = 0       P.R0y = 0       P.R02 = 0
% P.Rx0 = 0     P.Rxx{i} = 0    P.Rxy = 0       P.Rx2{i} = 0
% P.Ry0 = 0     P.Ryx = 0       P.Ryy{j} = 0    P.Ry2{j} = 0
% P.R20 = 0     P.R2x{i} = 0    P.R2y{j} = 0    P.R22{i,j} = 0
%
% INPUT
%   sos:    SOS program to modify.
%   P:      PI dopvar2d variable to set equal to 0.
%   opts:   Optional input. Set opts = 'symmetric' if P is known to be
%           symmetric, to enforce only parameters on and below the diagonal
%           to be zero (as the rest will follow by symmetry).
% OUTPUT 
%   sos:    Modified SOS program
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpi_eq_2d
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
% Initial coding MMP, SS, DJ  - 07_21_2021
% 07/24/2023 - DJ: Add option to exploit symmetry of operators;
%

% Check if symmetric option is specified.
if nargin>=3 &&  strcmp(opts,'symmetric')
    % Enforce only lower-triangular parameters to be zero
    % --> rest will follow by symmetry of the operator
    fset  = {'R00','Rx0','Ry0','Ryx','R20'};
    imax = 2;
    is_symmetric = true;
else
    % Enforce all parameters to be zero
    fset = {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'};
    imax = 3;
    is_symmetric = false;
end

% % Enforce parameters in operator to be zero.
% Start with parameters consisting of just one element
for f = fset
    if ~isempty(P.(f{:}))
        sos = soseq(sos, P.(f{:}));
    end
end
% Now for off-diagonal 3-PI operators
for i=1:3
    if ~isempty(P.R2x{i}) && any(any(P.R2x{i}.C))
        sos = soseq(sos, P.R2x{i});
    end
    if ~isempty(P.R2y{i}) && any(any(P.R2y{i}.C))
        sos = soseq(sos, P.R2y{i});
    end
if ~is_symmetric
    if ~isempty(P.Rx2{i}) && any(any(P.Rx2{i}.C))
        sos = soseq(sos, P.Rx2{i});
    end
    if ~isempty(P.Ry2{i}) && any(any(P.Ry2{i}.C))
        sos = soseq(sos, P.Ry2{i});
    end
end
end
% Finally, diagonal 3-PI (9-PI) operators
for i=1:imax
    if ~isempty(P.Rxx{i}) && any(any(P.Rxx{i}.C))
        sos = soseq(sos, P.Rxx{i});
    end    
    if ~isempty(P.Ryy{i}) && any(any(P.Ryy{i}.C))
        sos = soseq(sos, P.Ryy{i});
    end
    for j=1:imax
        if ~isempty(P.R22{i,j}) && any(any(P.R22{i,j}.C))
            sos = soseq(sos, P.R22{i,j});
        end
    end
end
if is_symmetric
    % correction for the 9-PI operator
    if ~isempty(P.R22{3,2}) && any(any(P.R22{3,2}.C))
        sos = soseq(sos, P.R22{3,2});
    end
end
end
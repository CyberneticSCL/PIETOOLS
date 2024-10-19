function sos = lpi_eq(sos,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sos = lpi_eq(prog,P) sets up equality constraints for each component.
% P.P=0
% P.Qi =0
% P.Ri = 0
% INPUT
%   sos:    SOS program to modify.
%   P:      PI dopvar variable to set equal to 0.
%   opts:   Optional input. Set opts = 'symmetric' if P is known to be
%           symmetric, to enforce only parameters on and below the diagonal
%           to be zero (as the rest will follow by symmetry).
% OUTPUT 
%   sos:    Modified SOS program
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpi_eq
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding MMP, SS, DJ  - 09/26/2021
% 07/24/2023 - DJ: Add option to exploit symmetry of operators;
%

% Check that the input is of appropriate type.
if isa(P,'polynomial') || isa(P,'double')
    error('Enforcing equality constraints on fixed values or polynomials is not supported.')
elseif isa(P,'dpvar')
    % If P is not an opvar, we can enforce equality using just soseq.
    sos = soseq(sos,P);
elseif isa(P,'opvar2d') || isa(P,'dopvar2d')
    % Pass 2D operator to associated lpi_eq function.
    if nargin>=3
        sos = lpi_eq_2d(sos,P,opts);
    else
        sos = lpi_eq_2d(sos,P);
    end
    return
elseif ~isa(P,'opvar') && ~isa(P,'dopvar')
    error('Input must be of type ''dopvar'' or ''dpvar''.')
end

% Check if symmetric option is specified.
if nargin>=3 && strcmpi(opts,'symmetric')
    % Enforce only lower-triangular parameters to be zero
    % --> rest will follow by symmetry of the operator
    i_set1 = {'P','Q2'};
    i_set2 = {'R0','R1'};
else
    % Enforce all parameters to be zero
    i_set1 = {'P','Q1','Q2'};
    i_set2 = {'R0','R1','R2'};
end

% Enforce parameters in operator to be zero.
for i = i_set1
    if ~isempty(P.(i{:}))
        sos = soseq(sos, P.(i{:}));
    end
end
for i = i_set2
    if ~isempty(P.R.(i{:}))
        C = P.R.(i{:}); 
        if ~all(all(C.C==0))
            sos = soseq(sos, C);
        end
    end
end
end
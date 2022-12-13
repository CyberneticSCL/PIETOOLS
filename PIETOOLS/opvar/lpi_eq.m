function sos = lpi_eq(sos,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sos = lpi_eq(prog,P) sets up equality constraints for each component.
% P.P=0
% P.Qi =0
% P.Ri = 0
% INPUT
%   prog: SOS program to modify.
%   P: PI dopvar variable
% OUTPUT 
%   sos: SOS program
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
%
if isa(P,'opvar2d') || isa(P,'dopvar2d')
    sos = lpi_eq_2d(sos,P);
    return
end
for i = {'P','Q1','Q2'}
    if ~isempty(P.(i{:}))
        sos = soseq(sos, P.(i{:}));
    end
end
for i = {'R0','R1','R2'}
    if ~isempty(P.R.(i{:}))
        C = P.R.(i{:}); 
        if ~all(all(C.C==0))
            sos = soseq(sos, C);
        end
    end
end
end
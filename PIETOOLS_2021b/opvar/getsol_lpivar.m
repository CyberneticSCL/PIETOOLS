function Psol = getsol_lpivar(sos,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = getsol_lpivar(sos,P) returns the solution for each component of P, 
% after solving the sos program prog.
% 
% INPUT 
%   prog: solved SOS program
%   P: PI dopvar variable
% 
% OUTPUT 
%   P: solved opvar class object
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - getsol_lpivar
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
% Initial coding MMP, SS, DJ  - 09/26/2021

opvar Psol; Psol.I = P.I; Psol.var1 = P.var1; Psol.var2 = P.var2;
Psol.dim = P.dim;

for i = {'P','Q1','Q2'}
    if ~isempty(P.(i{:}))
        Psol.(i{:}) = sosgetsol(sos, P.(i{:}));
    end
end
for i = {'R0','R1','R2'}
    if ~isempty(P.R.(i{:}))
        Psol.R.(i{:}) = sosgetsol(sos, P.R.(i{:}));
    end
end
Psol.P = double(Psol.P);
end
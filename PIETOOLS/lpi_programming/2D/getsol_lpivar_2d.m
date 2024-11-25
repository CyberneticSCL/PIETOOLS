function P = getsol_lpivar_2d(sos,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = getsol_lpivar(sos,P) returns the solution for each component of P, 
% after solving the sos program prog.
% 
% INPUT 
%   prog: solved SOS program
%   P: PI dopvar2d variable
% 
% OUTPUT 
%   P: solved opvar2d class object
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - getsol_lpivar_2d
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
% Initial coding MMP, SS, DJ - 07_21_2021 
%

for f = {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'}
    if ~isempty(P.(f{:}))
        P.(f{:}) = sosgetsol(sos, P.(f{:}));
    end
end
for i=1:3
    if ~isempty(P.Rxx{i,1})
        P.Rxx{i,1} = sosgetsol(sos, P.Rxx{i,1});
    end
    if ~isempty(P.Rx2{i,1})
        P.Rx2{i,1} = sosgetsol(sos, P.Rx2{i,1});
    end
    if ~isempty(P.R2x{i,1})
        P.R2x{i,1} = sosgetsol(sos, P.R2x{i,1});
    end
    
    if ~isempty(P.Ryy{1,i})
        P.Ryy{1,i} = sosgetsol(sos, P.Ryy{1,i});
    end
    if ~isempty(P.Ry2{1,i})
        P.Ry2{1,i} = sosgetsol(sos, P.Ry2{1,i});
    end
    if ~isempty(P.R2y{1,i})
        P.R2y{1,i} = sosgetsol(sos, P.R2y{1,i});
    end
    
    for j=1:3
        if ~isempty(P.R22{i,j})
            P.R22{i,j} = sosgetsol(sos, P.R22{i,j});
        end
    end
end

P = opvar2d(P);

end
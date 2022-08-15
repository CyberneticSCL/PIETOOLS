function P1 = clean_opvar(P1,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P1 = clean_opvar(P1,tol) removes terms from opvar2d P1 with coefficients
% with a value less than a tolerance tol
% Date: 08/27/21
% Version: 1.0
% 
% INPUT
% P1: opvar2d class objects
% tol: acceptable tolerance value, defaults to 1e-14. 
%       -Any term (of any component) of P1 with a coefficient value less 
%        than or equal to tol gets removed.
% 
% OUTPUT
% P1: opvar2d class object with only nontrivial terms retained
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - clean_opvar
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
% Initial coding DJ - 08_27_2021 

if nargin<2
    tol=1e-14;
end    

% For each parameter P1.fset{f}, set terms with values less than the
% tolerance to zero
fset = {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'};

for f=1:length(fset)
    PR = P1.(fset{f});
    if isa(PR,'double')
        % If object is double, get rid of elements close to zero
        PR(abs(PR)<tol) = 0;
        P1.(fset{f}) = PR;
    elseif isa(PR,'polynomial')
        % If object is polynomial, get rid of coefficients close to zero
        PRC = PR.coef;
        PRC(abs(PRC)<tol) = 0;
        PR.coef = PRC;
        P1.(fset{f}) = combine(PR); % Combine polynomial to get rid of unnecessary zeros
    else
        error('Parameters of an opvar2d class object must be of type "double" or "polynomial"')
    end
end

% Components in fsetc are comprised of cells, with each cell element a
% double or polynomial
fsetx = {'Rxx','Rx2','R2x'};
fsety = {'Ryy','Ry2','R2y'};
fset2 = {'R22'};
fsetc = [fsetx, fsety, fset2];
% Move through each cell, and set terms below threshold to zero
for f=1:length(fsetc)
    PR = P1.(fsetc{f});
    for indx=1:numel(PR)
        PRR = PR{indx};
        if isa(PRR,'double')
            PRR(abs(PRR)<tol) = 0;
            PR{indx} = PRR;
        elseif isa(PRR,'polynomial')
            PRRC = PRR.coef;
            PRRC(abs(PRRC)<tol) = 0;
            PRR.coef = PRRC;
            PR{indx} = combine(PRR);
        else
            error('Parameters of an opvar2d class object must be of type "double" or "polynomial"')
        end
    end
    P1.(fsetc{f}) = PR;
end

end
function P1 = clean_opvar(P1,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P1 = clean_opvar(P1,tol) removes terms from dopvar P1 with coefficients
% with a value less than a tolerance tol
% 
% INPUT
% P1:  dopvar class objects
% tol: acceptable tolerance value, defaults to 1e-14. 
%       -Any term (of any component) of P1 with a coefficient value less 
%        than or equal to tol gets removed.
% 
% OUTPUT
% P1: dopvar class object with only nontrivial terms retained
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - clean_opvar
%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 03_09_2022

if nargin<2
    tol=1e-14;
end    

% For each parameter P1.fset{f}, set terms with values less than the
% tolerance to zero
fset = {'P','Q1','Q2'};

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
    elseif isa(PR,'dpvar')
        % If object is dpvar, also get rid of coefficients close to zero
        PRC = PR.C;
        PRC(abs(PRC)<tol) = 0;
        PR.C = PRC;
        P1.(fset{f}) = compress(PR);    % Remove monomials that no longer contribute
    else
        error('Parameters of a dopvar class object must be of type "double", "polynomial", or "dpvar"')
    end
end

% Components in fsetc are comprised of cells, with each cell element a
% double or polynomial
fsetx = {'R0','R1','R2'};
% Move through each cell, and set terms below threshold to zero
for f=1:length(fsetx)
    PR = P1.R.(fsetx{f});
    if isa(PR,'double')
        PR(abs(PR)<tol) = 0;
        P1.R.(fsetx{f}) = PR;
    elseif isa(PR,'polynomial')
        PRC = PR.coef;
        PRC(abs(PRC)<tol) = 0;
        PR.coef = PRC;
        P1.R.(fsetx{f}) = combine(PR);
    elseif isa(PR,'dpvar')
        % If object is dpvar, also get rid of coefficients close to zero
        PRC = PR.C;
        PRC(abs(PRC)<tol) = 0;
        PR.C = PRC;
        P1.R.(fsetx{f}) = compress(PR);    % Remove monomials that no longer contribute
    else
        error('Parameters of a dopvar class object must be of type "double", "polynomial", or "dpvar"')
    end
end

end
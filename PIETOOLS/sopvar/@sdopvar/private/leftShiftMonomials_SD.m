function [varsC,ZC,CC] = leftShiftMonomials_SD( ...
    varsA,ZA,CA,varsB,ZB,CB)
% leftShiftMonomials_SD
%
% Decision-variable-dependent version of 'leftShiftMonomials_SS'. Computes
%
%   ((I_p \otimes ZA')*CA{I,J})*((I_q \otimes ZB')*CB{I})
%       = (I_p \otimes ZC')*CC{I,J}
%
% for each cell index I,J, where each CA{I,J} = A + B'*d is affine in the
% decision variables d.
%
% INPUTS
% - varsA, ZA:  variables and monomial bases of the left factor;
% - CA:         nI x nJ cell of structs, with 'A' the vectorized constant
%               part (mA x q, mA = p*prod(deg(ZA)+1)) and 'B' the rows
%               indexed by decision variables;
% - varsB, ZB:  variables and monomial bases of the right factor;
% - CB:         nI x 1 cell of numeric (q*NB) x r matrices;
%
% OUTPUTS
% - varsC, ZC:  union of the variables, and the resulting monomial basis;
% - CC:         nI x nJ cell of structs with fields 'A' and 'B', such that
%               the constant part is shifted as in 'leftShiftMonomials_SS'
%               and 'B' is mapped with the same linear map, Tmap{I}, i.e.
%               CC{I,J}.B = CA{I,J}.B*Tmap{I}.';
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - leftShiftMonomials_SD
%
% Copyright (C) 2026 PIETOOLS Team
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
% Initial coding MP, SS, DJ, AT, DB - 10/24/2026
%

[nI,nJ] = size(CA);

NB = prod(cellfun(@numel,ZB));
if isempty(ZB)
    NB = 1;
end

if numel(CB) ~= nI
    error('leftShiftMonomials_SD:cellDimensions', ...
        'CB must contain one cell for each row of CA.');
end

if mod(size(CB{1},1),NB) ~= 0
    error('leftShiftMonomials_SD:CBDimensions', ...
        'CB dimensions are inconsistent with ZB.');
end

q = size(CB{1},1)/NB;
nCoef = numel(CA{1,1}.A);

if mod(nCoef,q) ~= 0
    error('leftShiftMonomials_SD:CADimensions', ...
        'CA dimensions are inconsistent with CB.');
end

mA = nCoef/q;
CAconstant = cell(nI,nJ);

for i = 1:nI
    for j = 1:nJ
        param = CA{i,j};

        if numel(param.A) ~= nCoef || size(param.B,2) ~= nCoef
            error('leftShiftMonomials_SD:affineDimensions', ...
                'Affine coefficient dimensions are inconsistent.');
        end

        CAconstant{i,j} = reshape(param.A,mA,q);
    end
end

[varsC,ZC,CCconstant,Tmap] = leftShiftMonomials_SS( ...
    varsA,ZA,CAconstant,varsB,ZB,CB);

CC = cell(nI,nJ);

for i = 1:nI
    T = Tmap{i};

    for j = 1:nJ
        CC{i,j} = struct( ...
            'A',CCconstant{i,j}(:), ...
            'B',CA{i,j}.B*T.');
    end
end
end
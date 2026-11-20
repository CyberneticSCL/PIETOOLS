function [varsC,ZC,CC,Tmap] = leftShiftMonomials_SD( ...
    varsA,ZA,CA,varsB,ZB,CB)
% [VARSC,ZC,CC,TMAP] = LEFTSHIFTMONOMIALS_SD(VARSA,ZA,CA,VARSB,ZB,CB)
% computes, for a left factor affine in the decision variables,
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
%               the constant part is shifted as for a fixed factor, and
%               'B' is mapped with the same linear map, Tmap{I}, i.e.
%               CC{I,J}.B = CA{I,J}.B*Tmap{I}.';
% - Tmap:       (optional) nI x 1 cell of sparse matrices, the linear maps
%               of the shift, vec(CC) = Tmap{I}*vec(CA). If requested,
%               CC{I,J}.B is returned
%               UNMAPPED (= CA{I,J}.B, nd x size(Tmap{I},2)) and the caller
%               applies Tmap{I}: the mapped B has numel(CC{I,J}.A) columns,
%               57e6 for three variables at degree 4, and a sparse matrix
%               stores one column pointer per column, empty or not;
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding MP, SS, DJ, AT, DB - 10/24/2026

% Sizes: NB monomials of the right factor, q its row blocks, mA x q the
% constant part of each CA cell as a matrix
[nI,nJ] = size(CA);

NB = prod(cellfun(@numel,ZB));
if isempty(ZB)
    NB = 1;
end

% CB: one cell per row of CA, with a multiple of NB rows
if numel(CB) ~= nI
    error('leftShiftMonomials_SD:cellDimensions', ...
        'CB must contain one cell for each row of CA.');
end

if mod(size(CB{1},1),NB) ~= 0
    error('leftShiftMonomials_SD:CBDimensions', ...
        'CB dimensions are inconsistent with ZB.');
end

% q row blocks of CB; the nCoef = mA*q constant coefficients of each CA cell
q = size(CB{1},1)/NB;
nCoef = numel(CA{1,1}.A);

if mod(nCoef,q) ~= 0
    error('leftShiftMonomials_SD:CADimensions', ...
        'CA dimensions are inconsistent with CB.');
end

mA = nCoef/q;
% Constant parts as matrices, for the numeric routine
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

% Shift the constant parts and get the linear map Tmap{I} of the shift,
% vec(CC) = Tmap{I}*vec(CA)
[varsC,ZC,CCconstant,Tmap] = leftShiftMonomials_SS( ...
    varsA,ZA,CAconstant,varsB,ZB,CB);

% The decision coefficients take the same map; the rows of B are the
% decision variables and are only carried through
CC = cell(nI,nJ);

for i = 1:nI
    T = Tmap{i};

    for j = 1:nJ
        CC{i,j} = struct( ...
            'A',CCconstant{i,j}(:), ...
            'B',CA{i,j}.B);
        % Apply the map only if the caller does not take it (4th output)
        if nargout < 4
            CC{i,j}.B = CC{i,j}.B*T.';
        end
    end
end
end
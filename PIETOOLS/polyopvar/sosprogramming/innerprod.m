function V = innerprod(Z1,Z2,P)
% V = INNERPROD(Z1,Z2,P) computes the weighted innter product of two
% distributed polynomial functions, V = <Z1,P*Z2>
%
% INPUTS
% - Z1,Z2:  n x 1 'polyopvar' objects representing distributed polynomial 
%           functions. If Z2 = [], we assume Z2 = Z1;
% - P:      n x n object of type 'double', 'polynomial', or 'dpvar';
%
% OUTPUTS
% - V:      1 x 1 'polyopvar' object representing the inner product of the
%           distributed polynomial vectors with weight P.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - innerprod
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
% DJ, 03/05/2026: Initial coding

% Extract information from the left polynomial
if isa(Z1,'double') || isa(Z1,'polynomial') || isa(Z1,'dpvar')
    % Convert to polyopvar
    Z1_C = Z1;
    Z1 = polyopvar();
    Z1.C.ops{1} = Z1_C;
    Z1.degmat = zeros(1,0);
elseif ~isa(Z1,'polyopvar')
    error("Inputs must be of type 'polyopvar'.")
end
ZopL = Z1.C;
ZxL = Z1;
ZxL.C = tensopvar();

% Set the weight of the inner product
if nargin<=2
    % Assume weight 1;
    P = 1;
end
if all(size(P)==1)
    P = P*eye(size(Z1,1));
elseif size(P,1)~=size(Z1,1)
    error("Dimension of weight must match that of the vectors to multiply.")
end

% Extract information from the right polynomial
if isempty(Z2)
    % Take the symmetric inner product
    if size(P,2)~=size(Z1,1)
        error("Dimension of weight must match that of the vectors to multiply.")
    end
    if size(ZopL.ops,1)==size(ZxL.degmat,1)
        V = quad2lin(P,ZopL,ZxL);
    else
        V = 0;
        for i=1:size(ZopL.ops,2)
            ZopL_i = ZopL;
            ZopL_i.ops = ZopL_i.ops(1);
            ZxL_i = ZxL;
            ZxL_i.degmat = ZxL_i.degmat(1,:);
            for j=1:size(ZopL.ops,2)
                ZopR_j = ZopL;
                ZopR_j.ops = ZopR_j.ops(j);
                ZxR_j = ZxL;
                ZxR_j.degmat = ZxR_j.degmat(j,:);
                V = V + quad2lin(P,ZopL_i,ZxL_i,ZopR_j,ZxR_j);
            end
        end
    end
    return
elseif isa(Z2,'double') || isa(Z2,'polynomial') || isa(Z2,'dpvar')
    % Convert to polyopvar
    Z2_C = Z2;
    Z2 = polyopvar();
    Z2.C.ops{1} = Z2_C;
    Z2.degmat = zeros(1,0);
elseif ~isa(Z2,'polyopvar')
    error("Inputs must be of type 'polyopvar'.")
end
ZopR = Z2.C;
ZxR = Z2;
ZxR.C = tensopvar();

% Take the inner product
if size(P,2)~=size(Z2,1)
    if size(Z2,1)~=size(Z1,1)
        error("Dimensions of vectors must match for inner product.")
    else
        error("Dimension of weight must match that of the vectors to multiply.")
    end
end
if size(ZopL.ops,1)==size(ZxL.degmat,1) && size(ZopR.ops,1)==size(ZxR.degmat,1)
    V = quad2lin(P,ZopL,ZxL,ZopR,ZxR);
elseif size(ZopR.ops,1)==size(ZxR.degmat,1)
    V = 0;
    for j=1:size(ZopL.ops,2)
        ZopL_j = ZopL;
        ZopL_j.ops = ZopL_j.ops(j);
        ZxL_j = ZxL;
        ZxL_j.degmat = ZxL_j.degmat(j,:);
        V = V + quad2lin(P,ZopL_j,ZxL_j,ZopR,ZxR);
    end
elseif size(ZopL.ops,1)==size(ZxL.degmat,1)
    V = 0;
    for j=1:size(ZopR.ops,2)
        ZopR_j = ZopR;
        ZopR_j.ops = ZopR_j.ops(j);
        ZxR_j = ZxR;
        ZxR_j.degmat = ZxR_j.degmat(j,:);
        V = V + quad2lin(P,ZopL,ZxL,ZopR_j,ZxR_j);
    end
else
    V = 0;
    for i=1:size(ZopL.ops,2)
        ZopL_i = ZopL;
        ZopL_i.ops = ZopL_i.ops(j);
        ZxL_i = ZxL;
        ZxL_i.degmat = ZxL_i.degmat(j,:);
        for j=1:size(ZopR.ops,2)
            ZopR_j = ZopR;
            ZopR_j.ops = ZopR_j.ops(j);
            ZxR_j = ZxR;
            ZxR_j.degmat = ZxR_j.degmat(j,:);
            V = V + quad2lin(P,ZopL_i,ZxL_i,ZopR_j,ZxR_j);
        end
    end
end


end
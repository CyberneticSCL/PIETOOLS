function V = innerprod_v2(Z1,Z2,P)
    % V = innerprod_v2(Z1,Z2,P) computes the weighted inner product
    % V = <Z1,P*Z2>_{L2} in vector-valued form. It can be considered as
    % returning the individual terms from Lem. 7.
    %
    % For TDP factors from Def. 4, Z1 and Z2 remain vector-valued. The
    % coefficient matrix P is passed to quad2lin_v2 in one call; no
    % componentwise expansion is performed here.
    %
    % INPUTS
    % - Z1,Z2: n x 1 polyopvar distributed-polynomial vectors. If Z2 is
    %           empty, the symmetric inner product is used.
    % - P:     n x n double, polynomial, or dpvar weight.
    %
    % OUTPUT
    % - V:     scalar polyopvar representing <Z1,P*Z2>_{L2}.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - innerprod_v2
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 3 of the License, or
    % (at your option) any later version.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CR, 09/07/2026: Initial coding - Main diff to innerprod is that
    % quad2lin_v2 is called instead of quad2lin and the full P matrix is
    % passed in.

    narginchk(1,3);

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
    ZxL.C.ops = {};

    % Set the weight of the inner product. An empty weight retains the v2 
    % convenience of using the identity weight.
    if nargin<=2 || isempty(P)
        P = 1;
    end
    if all(size(P)==1)
        P = P*eye(size(Z1,1));
    elseif size(P,1)~=size(Z1,1)
        error("Dimension of weight must match that of the vectors to multiply.")
    end

    % Extract information from the right polynomial.
    if nargin<2 || isempty(Z2)
        % Take the symmetric inner product
        if size(P,2)~=size(Z1,1)
            error("Dimension of weight must match that of the vectors to multiply.")
        end

        % quad2lin_v2 recognizes the symmetric form and avoids generating
        % the two mirrored block terms separately.
        V = quad2lin_v2(P,ZopL,ZxL);                                      % CR, 09/07/2026
        return
    elseif isa(Z2,'double') || isa(Z2,'polynomial') || isa(Z2,'dpvar')
        % Convert to polyopvar.
        Z2_C = Z2;
        Z2 = polyopvar();
        Z2.C.ops{1} = Z2_C;
        Z2.degmat = zeros(1,0);
    elseif ~isa(Z2,'polyopvar')
        error("Inputs must be of type 'polyopvar'.")
    end

    ZopR = Z2.C;
    ZxR = Z2;
    ZxR.C.ops = {};

    % Check the right dimension after the right polynomial is available.
    if size(P,2)~=size(Z2,1)
        if size(Z2,1)~=size(Z1,1)
            error("Dimensions of vectors must match for inner product.")
        else
            error("Dimension of weight must match that of the vectors to multiply.")
        end
    end

    % The vector-valued backend replaces the original quad2lin branching:
    % one call passes the complete P block with both TDP vectors.
    V = quad2lin_v2(P,ZopL,ZxL,ZopR,ZxR); 
end

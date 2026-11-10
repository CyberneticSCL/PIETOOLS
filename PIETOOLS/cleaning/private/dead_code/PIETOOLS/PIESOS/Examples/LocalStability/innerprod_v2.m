function V = innerprod_v2(Z1,Z2,P,opts)
    % V = INNERPROD_V2(Z1,Z2,P,OPTS) returns the weighted inner product
    % V = <Z1,P*Z2>_{L2} in vector-valued form, which can be considered as
    % the individual terms from Lem. 7.
    %
    % For TDP factors from Def. 4, Z1 and Z2 remain vector-valued. The
    % coefficient matrix P is passed to quad2lin_v2 in one call; no
    % componentwise expansion is performed here.
    %
    % INPUTS
    % - Z1,Z2:  m x 1 and n x 1 'polyopvar' objects representing
    %           distributed-polynomial vectors (or 'double',
    %           'polynomial' or 'dpvar' constants). If Z2 is empty, the
    %           symmetric inner product <Z1,P*Z1> is used;
    % - P:      m x n 'double', 'polynomial' or 'dpvar' weight. A scalar
    %           weight is taken to be a multiple of the identity.
    %           Defaults to 1;
    % - opts:   (optional) char. Pass 'skip_combine' to skip the
    %           combine_terms folding of repeated states in quad2lin_v2.
    %           Defaults to '';
    %
    % OUTPUTS
    % - V:      scalar 'polyopvar' object representing <Z1,P*Z2>_{L2};
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
    % DJ, 10/15/2026: Updated to use new polyopvar functions;

    narginchk(1,4);
    if nargin<4 || isempty(opts)
        opts = '';                                                         % DJ, 10/05/2026
    end

    % Extract information from the left polynomial
    if isa(Z1,'double') || isa(Z1,'polynomial') || isa(Z1,'dpvar')
        % Convert to polyopvar
        Z1_C = Z1;
        Z1 = constant2polyopvar(Z1_C);  % degree-0 monomial, degmat not settable  % DJ, 10/14/2026
    elseif ~isa(Z1,'polyopvar')
        error("Inputs must be of type 'polyopvar'.")
    end

    ZopL = Z1.C;
    ZxL = Z1;
    % (ZxL is only used for its monomials and states, ZopL for the operators)

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
        % the two mirrored block terms separately. ZopR/ZxR passed as []
        % since opts is a trailing argument.                               % DJ, 10/05/2026
        V = quad2lin_v2(P,ZopL,ZxL,[],[],opts);                            % CR, 09/07/2026
        return
    elseif isa(Z2,'double') || isa(Z2,'polynomial') || isa(Z2,'dpvar')
        % Convert to polyopvar.
        Z2_C = Z2;
        Z2 = constant2polyopvar(Z2_C);                                      % DJ, 10/14/2026
    elseif ~isa(Z2,'polyopvar')
        error("Inputs must be of type 'polyopvar'.")
    end

    ZopR = Z2.C;
    ZxR = Z2;

    % Check the right dimension after the right polynomial is available.
    if size(P,2)~=size(Z2,1)
        if size(Z2,1)~=size(Z1,1)
            error("Dimensions of vectors must match for inner product.")
        else
            error("Dimension of weight must match that of the vectors to multiply.")
        end
    end

    % The vector-valued backend takes the complete P block with both TDP
    % vectors in one call.
    V = quad2lin_v2(P,ZopL,ZxL,ZopR,ZxR,opts);                             % DJ, 10/05/2026 
end

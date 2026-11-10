function C = DMB(A,B)
    % C = DMB(A,B) returns the 'polyopvar' object C representing the
    % product of the distributed monomials of A and B, as in the m=1 case
    % of Def. 4 in Paper1. The tensor-PI operators of A and B are combined
    % using the tensor-PI product otimes, so that C acts on the tensor
    % product of the monomials of A and B. This allows products such as
    % Zs{2} = Zx*Zx, where Zs{i} represents (Z*x)^(\otimes i).
    %
    % INPUTS
    % - A, B:   'polyopvar' objects representing a single monomial in the
    %           same scalar state variable, defined on the same spatial
    %           variables and domain. The 'tensopmat' A.C (B.C) must hold
    %           a single 'tensopvar' operator, in A.C.ops{1,1}, which may
    %           act on a tensor product of any degree;
    %
    % OUTPUTS
    % - C:      'polyopvar' object representing the tensor product of A
    %           and B, with kA*kB output components if the operators of A
    %           and B have kA and kB output components, respectively;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - DMB
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CR, 09/03/2026: Initial coding
    % CR, 09/07/2026: Use Def. 4 for vector-valued 3-PI factors. The output
    %                   component count is the product of the two counts.
    % DJ, 10/15/2026: Updated to use new polyopvar functions;

    narginchk(2,2);
    
    if ~isa(A,'polyopvar') || ~isa(B,'polyopvar')
        error('polyopvar_product:InvalidInput', ...
            'Both inputs must be polyopvar objects.');
    end

    if ~isequal(A.varname,B.varname)
    error('polyopvar_product:IncompatibleVariables', ...
        'Both operators must act on the same state variables.');
    end

    if ~isequal(A.dims.in,B.dims.in)
        error('polyopvar_product:Invalid_varsize', ...
            'Both polyopvar objects must contain the vector-dimension of the state variables.');
    end

    if ~isequal(size(A.degmat),size(B.degmat)) || ~isequal(size(A.degmat),[1,1])
        error('polyopvar_product:Invalid_degmat', ...
            'Both polyopvar objects must contain one tensopmat coefficient.');
    end

    if ~isequal(A.vars.in,B.vars.in)
        error('polyopvar_product:Invalid_pvarname', ...
            'Both polyopvar objects must contain the same spatial variables.');
    end


    if ~isequal(A.dom.in,B.dom.in)
        error('polyopvar_product:Incompatible_dom', ...
            'The polyopvar objects must have the same domain.');
    end

    if ~isequal(A.depmat.in,B.depmat.in) || ~isequal(A.dims.in,1)
        error('polyopvar_product:Incompatible_varmat', ...
            'Both input factors must describe one scalar state variable.');
    end

    if ~isequal(size(A.C.ops),[1,1]) || ~isequal(size(B.C.ops),[1,1]) || ...
            ~isa(A.C.ops{1},'tensopvar') || ~isa(B.C.ops{1},'tensopvar')
        error('polyopvar_product:InvalidCoefficient', ...
            'Each input must contain one tensopvar coefficient in C{1,1}.');
    end

    % Def. 4 maps the tensor-product input through R_A otimes R_B. Since
    % R_A and R_B have k_A and k_B output components, respectively, the new
    % TDP has k_A*k_B components while retaining one TDP coefficient.
    TA = A.C.ops{1};
    TB = B.C.ops{1};
    kA = size(TA,1);
    kB = size(TB,1);
    % Kronecker product of the components at the same point s, tensor product of the inputs. The new
    % otimes takes this as the pointwise product (type 1) of the rows repelem(TA,kB,1) and repmat(TB,kA,1)
    % (as `piesos_quadvar`), and type 2 would give a function of (s1,s2)
    TDP = otimes(repelem(TA,kB,1),repmat(TB,kA,1),[1,2]);                   % DJ, 10/14/2026
    if size(TDP,1) ~= kA*kB
        error('polyopvar_product:InvalidProductDimension', ...
            'Def. 4 product has an inconsistent output component count.');
    end

    % The monomial of the product has the factors of A, followed by those of B
    % (as the input of the tensor product otimes(TA,TB)), degmat is not settable
    C = A;
    C.Z = [A.Z,B.Z];                                                        % DJ, 10/14/2026

    % Rebuild the enclosing tensopmat from the TDP. Direct cell assignment
    % leaves vars, dom, and dependency arrays empty for later products.
    C.C = tensopmat(TDP);                                                % CR, 09/07/2026

    % The output space of the product is that of the new tensopmat (as in piesos_quadvar)
    C.dims.out = C.C.dims.out;
    C.vars.out = C.C.vars.out;
    C.dom.out = C.C.dom.out;
    C.depmat.out = C.C.depmat.out;
    C.chkval = true;

    % Stacking varmat would incorrectly create a second state variable.
    % C.varmat = [A.varmat; B.varmat];

end

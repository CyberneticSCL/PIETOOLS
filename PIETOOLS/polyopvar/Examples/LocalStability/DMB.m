function C = DMB(A,B)
    % C = DMB(A,B) Computes the m=1 case of Def. 4 in Paper1.
    % The product is formed between the TDP coefficient operators associated
    % with two distributed-monomial basis factors. This allows products such as
    % Zs{2} = Zx*Zx, where Zs{i} represents (Z*x)^(\otimes i). Hence, this function is restricted to
    % polyopvar objects having a single operator coefficient: a 1-by-1 C array
    % containing one tensopmat.  That tensopmat may act on a tensor product of
    % any degree.  Polyopvars with multiple C coefficients or separately stored
    % operator summands are not handled.  The tensopmat coefficients are
    % combined using the Tensor-PI product otimes.
    %
    % The intended inputs also use the same state-variable vector and state
    % dimensions, the same spatial variables and ordering, and the same
    % spatial domain. Each input is assumed to have one state-variable block
    % and one spatial-variable block (varmat == 1).
    %
    % INPUTS
    % - A, B  polyopvar objects with compatible spatial domains, each containing
    %          exactly one tensopmat operator coefficient in C{1,1}.
    %
    % OUTPUT
    % - C     polyopvar object representing the tensor product A otimes B.
    
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

    narginchk(2,2);
    
    if ~isa(A,'polyopvar') || ~isa(B,'polyopvar')
        error('polyopvar_product:InvalidInput', ...
            'Both inputs must be polyopvar objects.');
    end

    if ~isequal(A.varname,B.varname)
    error('polyopvar_product:IncompatibleVariables', ...
        'Both operators must act on the same state variables.');
    end

    if ~isequal(A.varsize,B.varsize)
        error('polyopvar_product:Invalid_varsize', ...
            'Both polyopvar objects must contain the vector-dimension of the state variables.');
    end

    if ~isequal(size(A.degmat),size(B.degmat)) || ~isequal(size(A.degmat),[1,1])
        error('polyopvar_product:Invalid_degmat', ...
            'Both polyopvar objects must contain one tensopmat coefficient.');
    end

    if ~isequal(A.pvarname,B.pvarname)
        error('polyopvar_product:Invalid_pvarname', ...
            'Both polyopvar objects must contain the same spatial variables.');
    end


    if ~isequal(A.dom,B.dom)
        error('polyopvar_product:Incompatible_dom', ...
            'The polyopvar objects must have the same domain.');
    end

    if ~isequal(A.varmat,B.varmat) || ~isequal(A.varmat,1)
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
    TDP = otimes(TA,TB,[true,true]);                                     % CR, 09/07/2026
    if size(TDP,1) ~= kA*kB
        error('polyopvar_product:InvalidProductDimension', ...
            'Def. 4 product has an inconsistent output component count.');
    end

    C = polyopvar();
    C.varname = A.varname;
    C.varsize = A.varsize(:);
    C.degmat = A.degmat + B.degmat;

    % Rebuild the enclosing tensopmat from the TDP. Direct cell assignment
    % leaves vars, dom, and dependency arrays empty for later products.
    C.C = tensopmat(TDP);                                                % CR, 09/07/2026
    
    C.pvarname = A.pvarname;
    C.dom = A.dom;
    
    % Stacking varmat would incorrectly create a second state variable.
    % C.varmat = [A.varmat; B.varmat];
    C.varmat = A.varmat;                                                 % CR, 09/07/2026

end

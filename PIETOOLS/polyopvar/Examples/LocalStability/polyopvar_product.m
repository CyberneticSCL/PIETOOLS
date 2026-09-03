function C = polyopvar_product(A,B)
    % C = polyopvar_product(A,B) Computes a special case of the tensor product between 
    % two polyopvars parameterized by tensopvar operators. This function allows products such as
    %
    %     Zs{i} = polyopvar_product(Zs{i-1},Zx);
    %
    % to be computed when Zs{i} represents (Z*x)^(\otimes i).  Hence, this function is restricted to
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
    % PIETOOLS - polyopvar_product
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CRR, 09/03/2026: Initial coding

    narginchk(2,2);
    
    % A.C.ops % 1 by 1 cell arrays
    % size(A.C) % {n by m tensopvar}
    % numel(A.C) % number of cells


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

    C = polyopvar();
    C.varname = unique([A.varname(:);B.varname(:)]);
    C.varsize = A.varsize(:);
    C.degmat = A.degmat + B.degmat;
    C.C.ops{1} = otimes(A.C.ops{1},B.C.ops{1});
    C.pvarname = A.pvarname;
    C.dom = A.dom;
    C.varmat = [A.varmat; B.varmat];

end

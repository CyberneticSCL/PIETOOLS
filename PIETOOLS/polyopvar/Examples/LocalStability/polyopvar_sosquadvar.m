function [prog,Pcell] = polyopvar_sosquadvar(prog,Z1,Z2,option)
    % Declares an appropriately sized block Gram matrix (represented as a cell) 
    % determined by the polyopvar bases. This is the polyopvar analogue of sosquadvar.
    % Z1 and Z2 are cell arrays containing polyopvar basis operators (as returned, for example, by Z*x).
    % A single Gram matrix is declared, so all returned blocks are coupled by the selected PSD/symmetry 
    % option. The output Pcell{i,j} is the (i,j)-th block of that Gram matrix.
    %
    % INPUTS
    % - prog    Current PIESOS/SOSTOOLS program.
    % - Z1      Cell array of left polyopvar basis operators.
    % - Z2      Cell array of right polyopvar basis operators.
    % - option  'pos', 'sym', or another option accepted by sosquadvar.
    %
    % OUTPUTS
    % - prog    Updated program containing the Gram decision variables.
    % - Pcell   Cell array of Gram-matrix blocks.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - ndopvar_sosquadvar
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CRR, 09/01/2026: Initial coding
    
    narginchk(4,4);
    
    if ~iscell(Z1) || ~iscell(Z2) || isempty(Z1) || isempty(Z2)
        error('polyopvar_sosquadvar:InvalidBasis', ...
            'Z1 and Z2 must be nonempty cell arrays.');
    end
    
    if numel(Z1) ~= numel(Z2)
        error('polyopvar_sosquadvar:InvalidBasis', ...
            'Z1 and Z2 must contain the same number of blocks.');
    end
    
    d = numel(Z1);  % degree of Z1 and Z2.
    m = zeros(d,1); % first dimension of each monomial basis.
    n = zeros(d,1); % second dimension of each monomial basis.

    for i = 1:d
        if ~isa(Z1{i},'polyopvar') || ~isa(Z2{i},'polyopvar')
            error('polyopvar_sosquadvar:InvalidBasis', ...
                'Each basis block must be a polyopvar object.');
        end

        m(i) = size(Z1{i},1);
        n(i) = size(Z2{i},1);
    end

    if strcmpi(option,'pos') && any(m ~= n)
        error('polyopvar_sosquadvar:InvalidDimension', ...
            'A positive Gram variable requires matching left and right dimensions.');
    end

    % Declare one matrix variable.  This is what couples all Pcell blocks
    % into one global PSD/symmetric Gram matrix.
    [prog,P] = sosquadvar(prog,{1},{1},sum(m),sum(n),option);
    P = P{1}; % P is no longer a cell.

    Pcell = cell(d,d);
    
    row0 = 0;
    for i = 1:d
        col0 = 0;
        for j = 1:d
            Pcell{i,j} = P( row0+(1:m(i)), col0+(1:n(j)) );
            col0 = col0 + n(j);
        end
        row0 = row0 + m(i);
    end
end

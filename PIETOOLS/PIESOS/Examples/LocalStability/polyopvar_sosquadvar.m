function [prog,Pcell] = polyopvar_sosquadvar(prog,Z1,Z2,option)
    % [PROG,PCELL] = POLYOPVAR_SOSQUADVAR(PROG,Z1,Z2,OPTION) declares a
    % block Gram matrix P, represented as a cell, with the block sizes
    % determined by the 'polyopvar' bases Z1 and Z2. This is the
    % 'polyopvar' analogue of sosquadvar. A single matrix is declared, so
    % all blocks are coupled by the selected PSD/symmetry option.
    %
    % INPUTS
    % - prog:   PIESOS/SOSTOOLS program structure;
    % - Z1:     cell array of d 'polyopvar' objects, representing the
    %           left basis operators applied to the state, as returned,
    %           for example, by Z*x;
    % - Z2:     cell array of d 'polyopvar' objects, representing the
    %           right basis operators;
    % - option: (optional) 'pos', 'sym' or 'free', not case-sensitive.
    %           For 'pos' (positive semidefinite) and 'sym' (symmetric),
    %           size(Z1{i},1)=size(Z2{i},1) is required for all i. The
    %           'free' option declares an unconstrained matrix variable.
    %           Defaults to 'free';
    %
    % OUTPUTS
    % - prog:   PIESOS/SOSTOOLS program structure with the Gram decision
    %           variables added;
    % - Pcell:  d x d cell, with Pcell{i,j} the (i,j)th block of P, of
    %           size size(Z1{i},1) x size(Z2{j},1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - polyopvar_sosquadvar
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CR, 09/01/2026: Initial coding
    % DJ, 10/05/2026: Require m=n for 'sym' case as well, and pass
    %                   "lower(option)" to sosquadvar;
    % DJ, 10/15/2026: Updated to use new polyopvar functions;

    narginchk(3,4);

    if nargin<4 || isempty(option)
        option = 'free';
    end
    
    if ~iscell(Z1) || ~iscell(Z2) || isempty(Z1) || isempty(Z2)
        error('polyopvar_sosquadvar_v2:InvalidBasis', ...
            'Z1 and Z2 must be nonempty cell arrays.');
    end
    
    if numel(Z1) ~= numel(Z2)
        error('polyopvar_sosquadvar_v2:InvalidBasis', ...
            'Z1 and Z2 must contain the same number of blocks.');
    end
    
    d = numel(Z1);  % degree of Z1 and Z2.
    m = zeros(d,1); % first dimension of each monomial basis.
    n = zeros(d,1); % second dimension of each monomial basis.

    for i = 1:d
        if ~isa(Z1{i},'polyopvar') || ~isa(Z2{i},'polyopvar')
            error('polyopvar_sosquadvar_v2:InvalidBasis', ...
                'Each basis block must be a polyopvar object.');
        end

        m(i) = size(Z1{i},1);
        n(i) = size(Z2{i},1);
    end

    if (strcmpi(option,'pos') || strcmpi(option,'sym')) && any(m ~= n)
        error('polyopvar_sosquadvar_v2:InvalidDimension', ...
            'A positive or symmetric Gram variable requires matching left and right dimensions.');
    end

    % Declare one matrix variable. This is what couples all Pcell blocks
    % into one global matrix P. For the free option, omitting the final
    % argument leaves P unconstrained; 'pos' and 'sym' impose the usual
    % SOS quadratic-form constraints.
    if strcmpi(option,'free')
        [prog,P] = sosquadvar(prog,{1},{1},sum(m),sum(n));
    else                                                                    % DJ, 10/05/2026
        % sosquadvar requires option to match case-sensitively
        [prog,P] = sosquadvar(prog,{1},{1},sum(m),sum(n),lower(option));    
    end
    P = P{1}; % Extract the matrix from the cell returned by sosquadvar.

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

function k = kernel_degree(F)
    % k = kernel_degree(F) returns the maximum total polynomial degree
    % appearing in the kernels of a polyopvar functional F.
    %
    % For an intop coefficient, the rows of Kop.params.degmat contain the
    % exponents of the kernel monomials in the independent spatial variables.
    % The total degree of a kernel monomial is the sum of the exponents in
    % its row. Therefore, k is the maximum row sum over all intop
    % coefficients in F.
    %
    % Constant and decision-variable coefficients do not contribute spatial
    % kernel degree and are assigned degree zero.
    %
    % INPUTS
    % - F    1 x 1 'polyopvar' functional whose kernel degree is required.
    %
    % OUTPUTS
    % - kg   Maximum total polynomial degree of the kernels in F.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - kernel_degree
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 3 of the License, or
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
    % CR, 09/28/2026: Initial coding.

    narginchk(1,1);

    if ~isa(F,'polyopvar') || ~all(size(F)==1)
        error('kernel_degree:InvalidInput', ...
            'F must be a 1-by-1 polyopvar functional.');
    end

    % Initialize to zero so that functionals containing only constant or
    % decision-variable terms return k = 0.
    k = 0;

    % A polyopvar may contain several coefficient operators. Only intop
    % coefficients carry explicit polynomial kernels in their params field.
    for k = 1:numel(F.C.ops)
        Kop = F.C.ops{k};

        if ~isa(Kop,'intop')
            continue
        end

        degmat = full(Kop.params.degmat);
        if isempty(degmat)
            continue
        end

        % Each row represents one kernel monomial. Sum across the spatial
        % variables to obtain its total polynomial degree.
        k = max(k,max(sum(degmat,2)));
    end

    k = double(k);
end

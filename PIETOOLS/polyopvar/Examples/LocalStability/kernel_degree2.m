function opdeg = kernel_degree2(F)
    % opdeg = kernel_degree2(F) returns the spatial monomial degree to use
    % as opdeg in SOS_DP when F must equate to an SOS-FDP.
    %
    % Constant and decision-variable coefficients do not contribute spatial
    % kernel degree and are assigned degree zero.
    %
    % INPUTS
    % - F    1 x 1 'polyopvar' functional whose SOS basis degree is required.
    %
    % OUTPUTS
    % - opdeg  Smallest spatial monomial degree suitable for SOS_DP.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - kernel_degree2
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
    % CR, 09/29/2026: Return the SOS_DP spatial basis degree rather than the
    %                   expanded FDP kernel degree; add optional distributed
    %                   degree and fix accumulator/index handling.
    % CR, 10/07/2026: Calculate the required spatial monomial degree for
    %                   each FDP coefficient operator before taking the
    %                   maximum; document each step of the calculation.


    % The current implementation supports FDPs in one distributed state;
    % each row of F.degmat then identifies one distributed monomial v^w.
    if ~isa(F,'polyopvar') || ~(size(F.degmat,2)==1) 
        error('kernel_degree2:InvalidInput', ...
            'F must be a 1-by-1 polyopvar functional.');
    end

    % There is one coefficient operator K_w for every row of F.degmat in
    % the linear FDP representation F(v)=sum_w K_w*v^(otimes w).
    nz = size(F.degmat,1);

    % Store the SOS_DP spatial monomial degree required by each K_w. The
    % entries remain zero for coefficients with no spatial kernel.
    opdeg_opts = zeros(1,nz);

    % A polyopvar may contain several coefficient operators. Only intop
    % coefficients carry explicit polynomial kernels in their params field.
    for idx = 1:nz
        % Extract the coefficient operator multiplying the distributed
        % monomial whose degree is recorded in row idx of F.degmat.
        Kop = F.C.ops{idx};

        % Numeric coefficients and empty integral kernels have spatial
        % degree zero, so their initialized candidate opdeg is retained.
        if ~isa(Kop,'intop') || isempty(Kop.params.degmat)
            continue
        end

        % Each row of Kop.params.degmat is the exponent vector of one
        % spatial monomial in the kernel. Summing across its columns gives
        % that monomial's total spatial degree.
        kernel_term_degrees = sum(full(Kop.params.degmat),2);

        % The raw kernel degree of K_w is the largest total degree among
        % all spatial monomials appearing in that coefficient operator.
        raw_degree = max(kernel_term_degrees);

        % The row F.degmat(idx,:) contains the exponents of the distributed
        % states in v^w. Their sum is the linear-form distributed degree w.
        distributed_degree = sum(F.degmat(idx,:));

        % SOS_DP represents a degree-2d SOS-FDP as a quadratic form. Thus,
        % the basis degree associated with a degree-w linear term is the
        % smallest integer d for which 2d can contain that term.
        quadratic_degree = ceil(distributed_degree/2);

        % Invert the kernel-degree expansion used by the quadratic-to-linear
        % SOS_DP map: raw_degree = 2*d*opdeg+1. The subtraction by one
        % removes the degree introduced by the outer integration, division
        % by 2*d recovers the spatial degree of one U-hat basis factor, and
        % ceil selects an integer degree large enough to contain the kernel.
        basis_degree = (raw_degree-1)/(2*quadratic_degree);

        % A spatial basis degree cannot be negative. Store the degree needed
        % for this coefficient operator before considering the other K_w.
        opdeg_opts(idx) = ceil(max(0,basis_degree));
    end

    % One SOS_DP basis is shared by every distributed-degree block, so its
    % monomial degree must be the maximum required by any coefficient K_w.
    opdeg = max(opdeg_opts);

end

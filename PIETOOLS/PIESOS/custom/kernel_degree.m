function opdeg = kernel_degree(F)
    % OPDEG = KERNEL_DEGREE(F) returns the smallest spatial monomial degree
    % opdeg of the basis operators U^i of a degree-2d SOS distributed
    % polynomial <Z_d(x),P*Z_d(x)> that can equate to F: for the coefficient
    % K_w of a monomial of distributed degree w, with maximal total kernel
    % degree D, opdeg = max over w of ceil(max(0,(D-1)/(2*ceil(w/2)))).
    % Coefficient
    % operators of F that are not 'intop' objects, such as constants or
    % decision variables, do not contribute spatial kernel degree and are
    % assigned degree zero.
    %
    % INPUTS
    % - F:      'polyopvar' object representing a functional in a single
    %           distributed state variable, whose SOS basis degree is
    %           required;
    %
    % OUTPUTS
    % - opdeg:  smallest nonnegative integer spatial monomial degree of the
    %           basis operators;
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
    % If you modify this code, document all changes carefully and
    % include date authorship, and a brief description of modifications
    %
    % For support, contact M. Peet, Arizona State University at mpeet@asu.edu
    % CR, 09/28/2026: Initial coding.


    % The current implementation supports FDPs in one distributed state;
    % each row of F.degmat then identifies one distributed monomial v^w.
    if ~isa(F,'polyopvar') || ~(size(F.degmat,2)==1) 
        error('kernel_degree:InvalidInput', ...
            'F must be a 1-by-1 polyopvar functional.');
    end

    % There is one coefficient operator K_w for every row of F.degmat in
    % the linear FDP representation F(v)=sum_w K_w*v^(otimes w).
    nz = size(F.degmat,1);

    % Store the spatial monomial degree required by each K_w. The
    % entries remain zero for coefficients with no spatial kernel.
    opdeg_opts = zeros(1,nz);

    % A polyopvar may contain several coefficient operators. Only intop
    % coefficients carry explicit polynomial kernels in their params field.
    for idx = 1:nz
        % Extract the coefficient operator multiplying the distributed
        % monomial whose degree is recorded in row idx of F.degmat.
        Kop = F.C.ops{1,idx};   % F.C is 1 x nz, column idx acts on monomial idx

        % Numeric coefficients and empty integral kernels have spatial
        % degree zero, so their initialized candidate opdeg is retained.
        if ~isa(Kop,'intop')
            continue
        end

        % Each monomial of each kernel of the operator has a total spatial
        % degree: the exponents of its monomials, summed over the variables.
        % The raw kernel degree of K_w is the largest total degree among
        % all spatial monomials appearing in that coefficient operator,
        % 0 for constant kernels.
        m = 0;
        Kc = kernels(Kop);
        for k=1:numel(Kc)
            Kk = Kc{k};
            if (isa(Kk,'polynomial') || isa(Kk,'dpvar')) && ~isempty(Kk.degmat)
                m = max(m,max(sum(full(Kk.degmat),2)));
            end
        end
        raw_degree = m;

        % The row F.degmat(idx,:) contains the exponents of the distributed
        % states in v^w. Their sum is the linear-form distributed degree w.
        distributed_degree = sum(F.degmat(idx,:));

        % A degree-2d SOS-FDP is a quadratic form. Thus,
        % the basis degree associated with a degree-w linear term is the
        % smallest integer d for which 2d can contain that term.
        quadratic_degree = ceil(distributed_degree/2);

        % Invert the kernel-degree expansion of the quadratic-to-linear
        % map: raw_degree = 2*d*opdeg+1. The subtraction by one
        % removes the degree introduced by the outer integration, division
        % by 2*d recovers the spatial degree of one U-hat basis factor, and
        % ceil selects an integer degree large enough to contain the kernel.
        basis_degree = (raw_degree-1)/(2*quadratic_degree);

        % A spatial basis degree cannot be negative. Store the degree needed
        % for this coefficient operator before considering the other K_w.
        opdeg_opts(idx) = ceil(max(0,basis_degree));
    end

    % One basis is shared by every distributed-degree block, so its
    % monomial degree must be the maximum required by any coefficient K_w.
    opdeg = max(opdeg_opts);

end



function opdeg = kernel_degree(F)
    % opdeg = kernel_degree(F) returns the spatial monomial degree to use
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
    % CR, 09/29/2026: Return the SOS_DP spatial basis degree rather than the
    %                   expanded FDP kernel degree; add optional distributed
    %                   degree and fix accumulator/index handling.


    if ~isa(F,'polyopvar') || ~(size(F.degmat,2)==1) 
        error('kernel_degree:InvalidInput', ...
            'F must be a 1-by-1 polyopvar functional.');
    end

    nz = size(F.degmat,1); % number of operators.
    deg_opts = zeros(1,nz); % max degrees of each operator.

    % A polyopvar may contain several coefficient operators. Only intop
    % coefficients carry explicit polynomial kernels in their params field.
    for idx = 1:nz
        Kop = F.C.ops{idx};

        if ~isa(Kop,'intop') || isempty(Kop.params.degmat)
            continue
        end
        
        % size = (number of unique terms in polynomial, number of variables in polynomial)
        % row entries = number of times each variable appears in each unique term
        degmat = full(Kop.params.degmat);

        if isempty(degmat)
            continue
        end
        
        sor = sum(degmat,2); % sum over rows = degree of each term.
        deg_opts(idx) = max(sor); % max. of sum over rows = degree of polynomial.
    end

    [max_deg, arg_deg] = max(deg_opts);
    d = F.degmat(arg_deg); % degree of FDP in linear form - should always be even for SOS, could be odd for anti-symmetric LFs.
    d = ceil(d/2); % degree of FDP in quadratic form - used in formula.

    opdeg = (max_deg-1) / (2*d); % this could be negative if deg_opts = 0.
    opdeg = max(0, opdeg); % maximum polynomial degree - could be fractional in case of non-symmetric LF, should be exact for SOS_DP.
    opdeg = ceil(opdeg);

end

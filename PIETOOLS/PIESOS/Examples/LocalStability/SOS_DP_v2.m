function [prog, DP, Pcell, Zop] = SOS_DP_v2(prog, d, opdeg, x, dom, opts)
    % [PROG,DP,PCELL,ZOP] = SOS_DP_V2(PROG,D,OPDEG,X,DOM,OPTS)
    % declares the degree-2*d SOS distributed polynomial
    %   DP = < Z_d(x), P Z_d(x) >_{L2}
    %      = \sum_i=1^d \sum_j=1^d <U^i x^i, Pmat U^j x^j>_{L_2}
    % (as in Def. 9 CDC paper) in the PIESOS program PROG, using
    % `piesos_sosvar`.
    %
    % INPUTS
    % - prog:   PIESOS program structure;
    % - d:      positive integer, maximal degree in Z_d. DP has degree
    %           2*d;
    % - opdeg:  nonnegative integer, degree of the spatial monomial
    %           basis in the operators U^i;
    % - x:      'polyopvar' object representing the fundamental state;
    % - dom:    1 x 2 array [a,b] specifying the spatial domain. Not
    %           used, the domain is that of x;
    % - opts:   (optional) struct with fields
    %               'psatz': order q of the psatz multipliers used to
    %                   enforce positivity only on the interval over
    %                   which the operator integrates, as in
    %                   `piesos_sosvar`. If psatz=0, no multiplier is
    %                   used. If psatz=q>0, P = sum_k g_k*P_k, with
    %                   g_k = (s-a)^k*(b-s)^(q-k) for k=0,...,q, and
    %                   P_k >= 0. Defaults to 0;
    %               'degs': degrees of the monomials in the operators
    %                   U^i, as in `piesos_sosvar`. Defaults to opdeg;
    %               'exclude': 3 x 1 boolean array specifying which
    %                   terms to exclude from the basis operators U^i,
    %                   as in `piesos_sosvar`. Defaults to [1;0;0], so
    %                   that U^i has no multiplier term.
    %               Other fields are passed on to `piesos_sosvar`;
    %
    % OUTPUTS
    % - prog:   PIESOS program structure with the decision variables
    %           defining the Gram operator P added;
    % - DP:     'polyopvar' object representing the inner product
    %           DP = < Z_d(x), P Z_d(x) >_{L2};
    % - Pcell:  q x q cell, with Pcell{i,j} the (i,j)th block of the
    %           Gram operator P, for q the number of monomials in Z_d;
    % - Zop:    'tensopmat' object representing the basis of operators
    %           U^i acting on the monomials x^i;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - SOS_DP_v2
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
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
    % CR, 09/01/2026: Initial coding.
    % CR, 09/07/2026: Re-routed innerprod via innerprod_v2.m.
    % DJ, 09/21/2026: Removed strict positivity constraint.
    % DJ, 10/06/2026: Reduce monomial degree of kernels for tensor product
    %                 of basis operators.
    % CR, 10/14/2026: Added opts structure as optional input so Pcell can
    %                 be constructed using a psatz multiplier.
    % DJ, 10/14/2026: Alternative formulation using `piesos_sosvar`
    % DJ, 10/15/2026: Updated to use new polyopvar functions;
    %% Declare the options for the SOS functional.
    if nargin<6 || isempty(opts)
        opts = struct();
    end
    if isfield(opts,'psatz') && (~isnumeric(opts.psatz) && ~islogical(opts.psatz))
        error("SOS_DP_v2: Invalid opts.psatz value!")
    end
    % Degrees of the monomials in the basis operators U^i; the multiplier term
    % is excluded, so that U^i only has the two integral terms (R0=0, R1
    % lower, R2 upper)
    if ~isfield(opts,'degs')
        opts.degs = opdeg;
    end
    if ~isfield(opts,'exclude')
        opts.exclude = [1;0;0];
    end

    %% Declare the SOS functional DP = <Zop*Z_d(x),P*Zop*Z_d(x)>, P>=0.

    % The basis of distributed monomials Z_d(x) = [x; x^2; ...; x^d]
    Z = dmonomials(x,(1:d)');                                               % DJ, 10/14/2026
    % The Gram operator P >= 0, and the SOS-FDP DP it defines (piesos_quadvar
    % exploits the symmetry of P: only the upper triangle is evaluated)
    [prog,DP,Pcell,Zop] = piesos_sosvar(prog,Z,opts);                       % DJ, 10/14/2026

end

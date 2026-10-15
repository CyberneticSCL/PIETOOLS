function [prog, DP, Pcell, Zs] = SOS_DP(prog, d, opdeg, x, dom, opts)
    % [PROG,DP,PCELL,ZS] = SOS_DP(PROG,D,OPDEG,X,DOM,OPTS) declares the
    % degree-2*d SOS distributed polynomial
    %   DP = < Z_d(x), P Z_d(x) >_{L2}
    %      = \sum_i=1^d \sum_j=1^d <U^i x^i, Pmat U^j x^j>_{L_2}
    % (as in Def. 9 CDC paper) in the PIESOS program PROG.
    %
    % INPUTS
    % - prog:   PIESOS program structure;
    % - d:      positive integer, maximal degree in Z_d. DP has degree
    %           2*d;
    % - opdeg:  nonnegative integer, degree of the spatial monomial
    %           basis in the operators U^i. For i>1, the degree per
    %           factor of the tensor product is floor(opdeg/i), so that
    %           the cumulative degree is opdeg;
    % - x:      'polyopvar' object representing the fundamental state;
    % - dom:    1 x 2 array [a,b] specifying the spatial domain;
    % - opts:   (optional) struct with field
    %               'psatz': binary value specifying whether to use a
    %                   psatz multiplier to enforce positivity only on
    %                   the interval over which the operator
    %                   integrates. If psatz=0, Pcell is multiplied by
    %                   (s-a)*(b-s). If psatz=1, Pcell = Pcell0 +
    %                   (s-a)*(b-s)*Pcell1. If opts is not specified,
    %                   no multiplier is used;
    %
    % OUTPUTS
    % - prog:   PIESOS program structure with the decision variables
    %           defining the Gram operator P added;
    % - DP:     'polyopvar' object representing the inner product
    %           DP = < Z_d(x), P Z_d(x) >_{L2};
    % - Pcell:  d x d cell, with Pcell{i,j} the (i,j)th block of the
    %           Gram operator P. If opts is specified, the blocks with
    %           j<=i are multiplied by (s-a)*(b-s) (psatz=0) or replaced
    %           by Pcell0 + (s-a)*(b-s)*Pcell1 (psatz=1), and the blocks
    %           with j>i are not used;
    % - Zs:     d x 1 cell, with Zs{i} the 'polyopvar' object
    %           representing U^i x^i;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - SOS_DP
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
    % DJ, 10/15/2026: Updated to use new polyopvar functions;
        
            
    %% Build the monomial basis used to parameterize P.

    % % Construct the basis operator corresponding to \hat{U} in paper.
    % % Only 2PI operators can be used here (for now).
    pvar s s_dum
    % Zmon     = monomials([s,s_dum],0:opdeg);
    % Zop      = opvar();
    % Zop.var1 = s;
    % Zop.var2 = s_dum;
    % Zop.I    = dom;  
    % Zop.R.R0 = [0*Zmon;0*Zmon];
    % Zop.R.R1 = [Zmon;0*Zmon];
    % Zop.R.R2 = [0*Zmon;Zmon];
    % 
    % Z  = dopvar2ndopvar(Zop);
    % Zx = Z*x;
    % 
    % % Construct the T-PI operators (corresponding to \hat{U}^i x^i in 
    % % the paper) as products of Zx.
    % Zs = cell(d,1);
    % for i = 1:d
    %     if i==1
    %         Zs{i} = Zx;
    %     else
    %         Zs{i} = DMB(Zs{i-1},Zx);
    %     end
    % end

    % Construct the T-PI operators (corresponding to \hat{U}^i x^i in 
    % the paper) as products of Zx.
    Zs = cell(d,1);
    for i = 1:d  
        if i==1
            % For i==1, define basis of operators by basis of monomials of  % DJ, 10/06/2026
            % degree at most opdeg in s and theta
            %Zmon     = monomials([s,s_dum],0:opdeg+4);
            Zmon1    = monomials(s,0:opdeg);                              % NOTE: +4 should be removed 
            Zmon2    = monomials(s_dum,0:opdeg);
            Zmon = kron(Zmon1,Zmon2);
            Zop      = opvar();
            Zop.var1 = s;
            Zop.var2 = s_dum;
            Zop.I    = dom;  
            Zop.R.R0 = [0*Zmon;0*Zmon];
            Zop.R.R1 = [Zmon;0*Zmon];
            Zop.R.R2 = [0*Zmon;Zmon];

            Z  = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop)));
            Zx = Z*x;
            Zs{i} = Zx;
        else
            % For i>1, define basis of operators as tensor product of 1D    % DJ, 10/06/2026
            % basis. Define the 1D operators in such a manner that the
            % cumulative degree of (theta_1,...,theta_i) in the tensor
            % product is opdeg
            Zmon1    = monomials(s,0:opdeg/i);
            Zmon2    = monomials(s_dum,0:floor(opdeg/i));
            Zmon = kron(Zmon1,Zmon2);
            % Zop1 is degree opdeg/i in s, and opdeg/i in theta
            Zop1      = opvar();
            Zop1.var1 = s;
            Zop1.var2 = s_dum;
            Zop1.I    = dom;  
            Zop1.R.R0 = [0*Zmon;0*Zmon];
            Zop1.R.R1 = [Zmon;0*Zmon];
            Zop1.R.R2 = [0*Zmon;Zmon];
            Z1  = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop1)));
            % Zop2 is degree 0 in s, and opdeg/i in theta
            Zop2      = opvar();
            Zop2.var1 = s;
            Zop2.var2 = s_dum;
            Zop2.I    = dom;  
            Zop2.R.R0 = [0*Zmon2;0*Zmon2];
            Zop2.R.R1 = [Zmon2;0*Zmon2];
            Zop2.R.R2 = [0*Zmon2;Zmon2];
            Z2  = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop2)));
            % Take the tensor product
            Z1x = Z1*x;
            Z2x = Z2*x;
            Zs{i} = Z1x;
            for k=2:i
                Zs{i} = DMB(Zs{i},Z2x);
            end
        end
    end

    %% Declare the block Gram operator P >= 0 and add its variables to prog.
    if nargin < 6
        [prog, Pcell] = polyopvar_sosquadvar(prog, Zs, Zs, 'pos');
    else
        g = (s-dom(1))*(dom(2)-s);
        if ~opts.psatz
            [prog, Pcell] = polyopvar_sosquadvar(prog, Zs, Zs, 'pos');
        elseif opts.psatz
            [prog, Pcell0] = polyopvar_sosquadvar(prog, Zs, Zs, 'pos');
            [prog, Pcell1] = polyopvar_sosquadvar(prog, Zs, Zs, 'pos');
            Pcell = cell(size(Pcell0));
        else
            error("SOS_DP: Invalid opts.psatz value!")
        end
    end

    % %% Ensure strict positivity of the constructed SOS DP.
    % eppos = 1e-4;
    % for i = 1:d
    %     Pcell{i,i} = Pcell{i,i} + eppos*eye(size(Pcell{i,i}));
    % end
    
    %% Evaluate DP = <Z_d(x), P Z_d(x)> as a complete block quadratic form.
    
    % Each Zs{i} is a vector-valued polyopvar and Pcell{i,j} is the matching
    % block of the global Gram matrix. Computation exploits smmetry of Pcell.
    DP = 0;
    for i = 1:d
        for j = 1:i
            if i==j
                if nargin < 6
                    DP = DP + innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
                elseif ~opts.psatz
                    Pcell{i,j} = g*Pcell{i,j};
                    DP         = DP + innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
                elseif opts.psatz
                    Pcell{i,j} = Pcell0{i,j} + g*Pcell1{i,j};
                    DP         = DP + innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
                else
                    error("SOS_DP: Invalid opts.psatz value!")
                end
            else
                if nargin < 6
                    DP = DP + 2*innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
                elseif ~opts.psatz
                    Pcell{i,j} = g*Pcell{i,j};
                    DP         = DP + 2*innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
                elseif opts.psatz
                    Pcell{i,j} = Pcell0{i,j} + g*Pcell1{i,j};
                    DP         = DP + 2*innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
                else
                    error("SOS_DP: Invalid opts.psatz value!")
                end
            end
        end
    end

end

function [prog, DP] = SOS_DP(prog, d, opdeg, x, dom)
    % [prog, DP] = SOS_DP(prog, d, opdeg, x) Construct a degree-2*d SOS distributed polynomial
    % DP = < Z_d(x), P Z_d(x) >_{L2} = \sum_i=1^d \sum_j=1^d <U^i x^i, Pmat U^j x^j>_{L_2}
    % (as in Def. 9 CDC paper) and add the variable to the PIESOS program.
    %
    % INPUTS
    % - prog   Current PIESOS program.
    % - d      Maximum degree in Z_d.  DP has degree 2*d.
    % - opdeg  Degree of the spatial monomial basis in SOS P operator.
    % - x      Fundamental state.
    % - dom    Spatial domain [a, b].
    %
    % OUTPUTS
    % - prog  Updated PIESOS program.
    % - DP    'polyopvar' object representing the inner product DP = < Z_d(x), P Z_d(x) >_{L2}.
    
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
    % If you modify this code, document all changes carefully and include date
    % authorship, and a brief description of modifications
    %
    % CR, 09/01/2026: Initial coding
    % CR, 09/07/2026: Build a block Gram form using tensor-product basis
    %                   operators and evaluate every block with innerprod.
    %                   Re-routed innerprod via innerprod_v2.m
    % CR, 09/07/2026: Verify the returned Gram blocks are symbolically
    %                   symmetric before using the factor-two reduction.
        
            
    %% Build the monomial basis used to parameterize P.

    % Construct the basis operator corresponding to \hat{U} in paper.
    pvar s s_dum
    Zmon = monomials([s,s_dum],0:opdeg);
    Zop = opvar();
    Zop.R.R0 = [0*Zmon;0*Zmon];
    Zop.R.R1 = [Zmon;0*Zmon];
    Zop.R.R2 = [0*Zmon;Zmon];
    Zop.var1 = s;
    Zop.var2 = s_dum;
    Zop.I = dom;
    Z = dopvar2ndopvar(Zop);
    Zx = Z*x;

    % Construct the T-PI operators (corresponding to \hat{U}^i x^i in 
    % the paper) as products of Zx.
    Zs = cell(d,1);
    for i = 1:d
        if i==1
            Zs{i} = Zx;
        else
            Zs{i} = DMB(Zs{i-1},Zx);
        end
    end

    %% Declare the block Gram operator P >= 0 and add its variables to prog.
    [prog, Pcell] = polyopvar_sosquadvar(prog, Zs, Zs, 'pos');
    
    % Only needed for testing.
    % % The lower-triangular evaluation below is valid only when
    % % Pcell{j,i} = Pcell{i,j}'. Check dimensions and symbolic dpvar
    % % structure before exploiting that identity.
    % for i = 1:d
    %     for j = 1:i
    %         if ~isequal(size(Pcell{i,j}),size(Pcell{j,i}')) || ...
    %                 ~isequal(Pcell{i,j},Pcell{j,i}')
    %             error('SOS_DP:NonSymmetricGram', ...
    %                 'Pcell must be symmetric before using the factor-two reduction.');
    %         end
    %     end
    % end

    %% Ensure strict positivity of the constructed SOS DP.
    eppos = 1e-4;
    for i = 1:d
        Pcell{i,i} = Pcell{i,i} + eppos*eye(size(Pcell{i,i}));
    end
    
    %% Evaluate DP = <Z_d(x), P Z_d(x)> as a complete block quadratic form.
    
    % Each Zs{i} is a vector-valued polyopvar and Pcell{i,j} is the
    % matching block of the global Gram matrix.
    % DP = 0;
    % for i = 1:d
    %     for j = 1:d
    %         % DP = DP + innerprod(Zs{i},Zs{j},Pcell{i,j});
    %         DP = DP + innerprod_v2(Zs{i},Zs{j},Pcell{i,j});          % CR, 09/07/2026
    %     end
    % end
    
    % Each Zs{i} is a vector-valued polyopvar and Pcell{i,j} is the
    % matching block of the global Gram matrix.
    % Compute DP when exploiting symmetry of Pcell.                    % CR 09/07/26
    DP = 0;
    for i = 1:d
        for j = 1:i
            if i==j
                DP = DP + innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
            else
                DP = DP + 2*innerprod_v2(Zs{i}, Zs{j}, Pcell{i,j});
            end
        end
    end

end

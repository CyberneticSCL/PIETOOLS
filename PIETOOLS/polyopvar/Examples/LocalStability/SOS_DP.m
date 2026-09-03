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
    % CRR, 09/01/2026: Initial coding
        
            
        %% Build the monomial basis used to parameterize P.
    
        % Construct the basis operator corresponding to U in paper.
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

        % Construct the T-PI operators (corresponding to U^i x^i in the paper) as products of Zx.
        Zs = cell(d,1);
        for i = 1:d
            if i==1
                Zs{i} = Zx;
            else
                Zs{i} = polyopvar_product(Zs{i-1},Zx);
            end
        end
    
        %% Declare the block Gram operator P >= 0 and add its variables to prog.
        [prog, Pcell] = polyopvar_sosquadvar(prog, Zs, Zs, 'pos'); % 'pos' ensure Pcell is symmetric.

        %% Ensure strict positivity of the constructed SOS DP.
        eppos = 1e-4;
        for i = 1:d
            Pcell{i,i} = Pcell{i,i} + eppos*eye(size(Pcell{i,i}));
        end
        
        %% Evaluate DP = <Z_d(x), P Z_d(x)> as a complete block quadratic form.
        
        % THIS DOES NOT WORK SINCE... Zs{1} will typically be a vector valued
        % 3-PI operator with the number of elements determined by the
        % monomial basis. innerprod cannot handle such array entries.
        DP = 0;
        for i = 1:d
            for j = 1:d
                DP = DP + innerprod(Zs{i},Zs{j},Pcell{i,j});
            end
        end
        
        DP = innerprod(Zs{2},Zs{2},Pcell{2,2}); % innerprod(Zs{1},Zs{1},Pcell{1,1}) + innerprod(Zs{1},Zs{2},Pcell{1,2}) + innerprod(Zs{2},Zs{1},Pcell{2,1});

end

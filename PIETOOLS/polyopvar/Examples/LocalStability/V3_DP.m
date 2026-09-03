function [prog, V3] = V3_DP(prog, d, opdeg, Top, x, dom)
    % [prog, V3] = V3_DP(prog, d, opdeg, x) Construct a degree-2*d distributed polynomial
    % V3 = < Z_d(x), P Z_d(T*x) >_{L2} (as in Sec. 7.5 Automatica paper) and add the variable to the PIESOS program.
    % In this function, P is constructed without any PSD and symmetry constraints. Positivity will be enforced when 
    % equating V3 to an SOS DP and the symmetry constraint will be enforced by...?
    %
    % INPUTS
    % - prog   Current PIESOS program.
    % - d      Maximum degree in Z_d.  DP has degree 2*d.
    % - opdeg  Degree of the spatial monomial basis in SOS P operator.
    % - Top    Inverse operator.
    % - x      Fundamental state.
    % - dom    Spatial domain [a, b].
    %
    % OUTPUTS
    % - prog  Updated PIESOS program.
    % - V3    'polyopvar' object representing the inner product DP = < Z_d(x), P Z_d(T*x) >_{L2}.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - V3_DP
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
    
        % Construct the basis Uhat.
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
        
        % Express as polyopvars which we can work with for constructing V3.
        Zx = Z*x;
        Tx = Top*x;
        ZTx = Z*Tx;
    
        % Construct the T-PI operators as products of Zx and ZTx.
        Zs1 = cell(d,1);
        Zs2 = cell(d,1);
        for i = 1:d
            if i==1
                Zs1{i} = Zx;
                Zs2{i} = ZTx;
            else
                Zs1{i} = Zs1{i-1}*Zx;
                Zs2{i} = Zs2{i-1}*ZTx;
            end
        end
    
    
        %% Declare the block Gram operator P and add its variables to prog.
        [prog, Pcell] = Polyopvar_sosquadvar(prog, Zs1, Zs2);
    
    
        %% Evaluate V3 = <Z_d(x), P Z_d(T*x)> as a complete block quadratic form.
        V3 = 0;
        for i = 1:d
            for j = 1:d
                V3 = V3 + innerprod(Zs1{i},Zs2{j},Pcell{i,j});
            end
        end

end

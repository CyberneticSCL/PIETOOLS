function [prog, V3] = V3_DP(prog, d, opdeg, Top, x, dom)
    % [prog, V3] = V3_DP(...) Construct a degree-2*d distributed polynomial
    % V3 = < Z_d(x), P Z_d(T*x) >_{L2} (as in Sec. 7.5 Automatica paper) and add the variable 
    % to the PIESOS program. In this function, P is constructed without any PSD constraints.
    % Positivity will be enforced when equating to an SOS DP and the symmetry 
    % constraint needed for the Lie derivative is enforced below.
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
    % CR, 09/01/2026: Initial coding
    % CR, 09/08/2026: Compose Z with Top before applying the result to x.
    %                   PIETOOLS does not support nopvar*polyopvar directly,
    %                   whereas (Z*Top)*x is the equivalent operator
    %                   composition Z(Top*x).
    % CR, 09/08/2026: Use an unconstrained Gram matrix and the vector-valued
    %                   inner-product implementation for the V3 blocks.
    % CR, 09/08/2026: Enforce the operator-adjoint condition through equality
    %                   of the corresponding scalar distributed-polynomial
    %                   forms.
    % DJ, 09/22/2026: Add monomials for multiplier operator in the basis
    %                   operator Zop.
    % CR, 09/28/2026: Removed symmetry constraints as these will be enforced 
    % when equating to an SOS_DP.
        
            
    %% Build the monomial basis used to parameterize P.

    % Construct the basis operator corresponding to \hat{U} in paper.
    % Currently just uses 2PI operators; however code is set to accept 3PI op. if d==1.
    % For d>1, code can only accept 2PI ops.
    pvar s s_dum
    Zmon     = monomials([s,s_dum],0:opdeg);
    Zop      = opvar();
    Zop.var1 = s;
    Zop.var2 = s_dum;
    Zop.I    = dom;    
    Zop.R.R0 = [0*Zmon;0*Zmon];
    Zop.R.R1 = [Zmon;0*Zmon];
    Zop.R.R2 = [0*Zmon;Zmon];

    
    % Express as tensopvar and polyopvars which we can work with for constructing V3.
    % Note that Z is a nopvar and Tx is a polyopvar. Applying Z directly to Tx invokes
    % an unsupported nopvar*polyopvar multiplication. Compose the two operators first, 
    % then apply the composed operator to the fundamental state. This gives 
    % ZTx = Z*(Top*x) = (Z*Top)*x.
    Z         = dopvar2ndopvar(Zop);
    ZTop      = Z*Top;                                                      % CR, 09/08/2026
    Zx        = Z*x;
    ZTx       = ZTop*x;                                                     % CR, 09/08/2026

    % Construct the T-PI operators (corresponding to \hat{U}^i x^i and 
    % (\hat{U} o T)^j x^j in the paper) as products of Zx and ZTx.
    Zs1 = cell(d,1);
    Zs2 = cell(d,1);
    for i = 1:d
        if i==1
            Zs1{i} = Zx;
            Zs2{i} = ZTx;
        else
            Zs1{i} = DMB(Zs1{i-1},Zx);
            Zs2{i} = DMB(Zs2{i-1},ZTx);
        end
    end


    %% Declare the block Gram operator P and add its variables to prog.
    % P in V3 is not constrained to be positive at this stage. The local helper 
    % represents this choice by the 'free' option.
    [prog, Pcell] = polyopvar_sosquadvar(prog, Zs1, Zs2, 'free');


    %% Evaluate V3 = <Z_d(x), P Z_d(Tx)> as a complete block quadratic form.

    % Each Zs1{i} and Zs2{j} are vector-valued polyopvars and Pcell{i,j} is the
    % matching block of the global Gram matrix. The symmetry constraint is naturally 
    % enforced when equating to an SOS DP, so no need to do it explicitly here.
    V3 = 0;
    for i = 1:d
        for j = 1:d
            V3 = V3 + innerprod_v2(Zs1{i},Zs2{j},Pcell{i,j});
            % if i >= j
            %     rhs = innerprod_v2(Zs2{i}, Zs1{j}, Pcell{j,i}');
            %     prog = piesos_eq(prog,Vij-rhs);
            % end
        end
    end

end

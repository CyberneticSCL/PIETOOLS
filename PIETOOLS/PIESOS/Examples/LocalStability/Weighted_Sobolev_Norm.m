function g = Weighted_Sobolev_Norm(alpha, Top, x)
    % G = WEIGHTED_SOBOLEV_NORM(ALPHA,TOP,X) returns the squared weighted
    % Sobolev norm (Sec. 7.1 in Automatica paper)
    %   g(x) = sum_j alpha(j+1)*||d^j/ds^j T*x||_{L2}^2
    % of the PDE state T*x, for the fundamental state x.
    %
    % INPUTS
    % - alpha:  m x (n+1) array of weights, with alpha(i,j+1) the weight of
    %           the jth derivative of component i of T*x. For m=1 (scalar
    %           PDE state), alpha can also be an array of n+1 weights.
    %           Derivatives are only taken of the components with nonzero
    %           weights, so that e.g. L2 states, without derivatives, can
    %           be included by setting the weights of their derivatives
    %           to zero;
    % - Top:    m x m' 'tensopmat' object with 'sopvar' operators,
    %           representing the map T from the fundamental state to the
    %           PDE state (PIE.T), in a single spatial variable;
    % - x:      m' x 1 'polyopvar' object representing the fundamental
    %           state (PIE.x);
    %
    % OUTPUTS
    % - g:      scalar 'polyopvar' object representing the squared
    %           weighted Sobolev norm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - Weighted_Sobolev_Norm
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
    % CR, 09/01/2026: Initial coding
    % CR, 09/07/2026: Route inner products through the local
    %                   MATLAB-compatible v2 implementation.
    % DJ, 10/06/2026: Directly avoid terms with zero factor.
    % CR, 10/12/2026: Changed from computing the weighted Sobolev ball (of
    %                   radius r) to computing just the weighted Sobolev
    %                   norm.
    % DJ, 10/15/2026: Updated to use new polyopvar functions.
    % DJ, 10/20/2026: Use '@sopvar/diff' function, instead of converting to
    %                   opvar.
    % DJ, 10/21/2026: Allow for vector-valued PDE states, with weights
    %                   alpha(i,j) for each component i.


    %n = size(alpha,2) - 1; % degree of PDE.
    % The map T is an m x m' block matrix of PI operators in a single
    % spatial variable, that we differentiate in this variable

    if ~isa(Top,'tensopmat')                                                % DJ, 10/21/2026
        error("Weighted_Sobolev_Norm: the map T must be a 'tensopmat'.")
    end 
    svar = Top.vars.out;                                                    % DJ, 10/20/2026

    % Component i of T*x gets the weights alpha(i,:): a vector alpha is
    % allowed only for m=1.
    m = size(Top,1);                                                        % DJ, 10/20/2026
    if m==1 && isvector(alpha)                                             
        alpha = reshape(alpha,1,[]);                                       
    end                                                                     
    if size(alpha,1)~=m               
        % ideally, each row of alpha would correspond to a vector-valued component...
        error("Weighted_Sobolev_Norm: alpha must have one row of weights per component of T*x.")  
    end                                                                  

    % g = r^2 - alpha(1)*innerprod(Top*x,Top*x);
    % g = r^2 - alpha(1)*innerprod_v2(Top*x,Top*x);                    % CR, 09/07/2026
    % g = r^2;
    g = 0;

    % Loop over (derivative order, component) pairs with nonzero weight, in
    % the order of the derivative as before; for m=1 this is the loop over
    % find(alpha) of the scalar case.
    for i=find(any(alpha,1))
        for c=find(alpha(:,i)')                                             % DJ, 10/20/2026
            % g = g - alpha(i)*innerprod_v2(Rop*x,Rop*x);                   % CR, 09/07/2026
            % R_{i-1} = d^{i-1}/ds^{i-1} T, as 'tensopmat' (T itself for i=1)
            % Component c only: the other components may be L2 states, that
            % cannot be differentiated.
            if m==1                                                         % DJ, 10/20/2026
                Rop = Top;
            else
                Rop = Top(c,:);     % ideally, we would allow for vector-valued component 
            end
            if i>1
                try                                                    
                    Rop = diff(Rop,svar,i-1);                                
                catch ME                                                   
                    error("Weighted_Sobolev_Norm: derivative of order %d of component %d of T*x is not defined (%s)",i-1,c,ME.message) 
                end                                                             
            end
            Rx = Rop*x;                   
            g = g + alpha(c,i)*innerprod(Rx,Rx);                            % DJ, 10/20/2026
        end
    end

end

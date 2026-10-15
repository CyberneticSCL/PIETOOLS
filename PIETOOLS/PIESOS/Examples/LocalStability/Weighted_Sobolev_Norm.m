function g = Weighted_Sobolev_Norm(alpha, Top, x)
    % G = WEIGHTED_SOBOLEV_NORM(ALPHA,TOP,X) returns the squared weighted
    % Sobolev norm (Sec. 7.1 in Automatica paper)
    %   g(x) = sum_j alpha(j+1)*||d^j/ds^j T*x||_{L2}^2
    % of the PDE state T*x, for the fundamental state x.
    %
    % INPUTS
    % - alpha:  array of n+1 weights, with alpha(j+1) the weight of the
    %           jth derivative;
    % - Top:    1 x 1 'tensopmat' object with a single 'sopvar' operator,
    %           representing the map T from the fundamental state to the
    %           PDE state (PIE.T);
    % - x:      scalar 'polyopvar' object representing the fundamental
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
    % DJ, 10/15/2026: Updated to use new polyopvar functions;


    %n = size(alpha,2) - 1; % degree of PDE.
    % The map is a single PI operator (scalar fundamental state, one
    % spatial variable): extract it as an 'opvar', to differentiate it
    if ~isa(Top,'tensopmat') || ~isequal(size(Top.ops),[1,1]) || ...
            ~isequal(size(Top.ops{1}.ops),[1,1]) || ~isa(Top.ops{1}.ops{1,1},'sopvar')
        error("Weighted_Sobolev_Norm: the map T must be a 1 x 1 'tensopmat' of a single 'sopvar' operator.")
    end
    Top_opvar = sopvar2opvar(Top.ops{1}.ops{1,1});                          % DJ, 10/14/2026

    % g = r^2 - alpha(1)*innerprod(Top*x,Top*x);
    % g = r^2 - alpha(1)*innerprod_v2(Top*x,Top*x);                    % CR, 09/07/2026
    % g = r^2;
    g = 0;

    for i=find(alpha(:)')                                                   % DJ, 10/06/2026
        % g = g - alpha(i)*innerprod(Rop*x,Rop*x);
        % g = g - alpha(i)*innerprod_v2(Rop*x,Rop*x);                  % CR, 09/07/2026
        % R_{i-1} = d^{i-1}/ds^{i-1} T, as 'tensopmat' (T itself for i=1)
        if i==1
            Rop = Top;
        else
            Rop = diff(Top_opvar,Top_opvar.var1,i-1,'pure');                % DJ, 10/14/2026
            Rop = tensopmat(sopvar2tensopvar(opvar2sopvar(Rop)));
        end
        Rx = Rop*x;
        g = g + alpha(i)*innerprod(Rx,Rx);                                  % DJ, 10/14/2026
    end

end

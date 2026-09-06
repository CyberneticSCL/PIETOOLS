function g = Weighted_Sobolev_Ball(r, alpha, Top, x)
    % g = weighted_sobolev_ball(...) Returns the weighted sobolev ball (Sec. 7.1 in Automatica paper).
    %
    % INPUTS
    % - r:         Radius of the local ball.
    % - alpha:     (n+1)-dim array containing parameters of weighted Sobolev ball.
    % - Top:       Map from fundamental to PDE state represented as an ndopvar object.
    % - x:         Fundamental state.
    
    %
    % OUTPUTS
    % - g:         Weighted sobolev ball of radius r.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - weighted_sobolev_ball
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
    
        %%
        n = size(alpha,2) - 1; % degree of PDE.
        Top_opvar = ndopvar2dopvar(Top);
        
        
        g = r^2 - alpha(1)*innerprod(Top*x,Top*x);
        
        for i=2:n+1
            Rop = dopvar2ndopvar(diff(Top_opvar,Top_opvar.var1,i-1,'pure'));
            g = g - alpha(i)*innerprod(Rop*x,Rop*x);
        end

end
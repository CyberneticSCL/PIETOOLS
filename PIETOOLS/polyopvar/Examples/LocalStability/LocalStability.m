function res = LocalStability(PDE, r, alpha, eppos, lambda, dist_degs, mon_degs, C)

    % res = Local_Stability(...) runs the local stability test (Thm. 1 in Automatica paper).
    %
    % INPUTS
    % - PDE:            The PDE to convert to a PIE.
    % - r:              Radius of the local ball.
    % - alpha:          (n+1)-dim array containing parameters of weighted Sobolev ball.
    % - eppos:          Treated as eppos^2, the lower bound on LF.
    % - lambda:         Exponential decay rate.
    % - dist_degs:      3d array containing degrees of distributed monomial basis for V, p1, p2 resptively.
    % - mon_degs:       3d array containing degrees of monomial basis for V, p1, p2 resptively.
    % - C:              Upper bound on LF (optional argument).
    % OUTPUTS
    % - res:            2d array containing C_sol (treated as C^2, upper bound on LF) and 
    %                   M_sol (multiplier of exponential upper bound on solution norm).
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % PIETOOLS - LocalStability.m
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


    %% Convert PDE to PIE.
    PIE = convert(PDE);
    Top = PIE.T; % Inverse map.
    f   = PIE.f; % Polynomial PIE.
    dom = PIE.dom; % Spatial domain [a,b].
    
    % Construct fundamental state.
    x = f.vartab;

    
    %% Initialize PIESOS program structure.
    
    if nargin < 8
        dpvar C % treated as C^2, the upper bound on the LF.
        prog = piesos_program(x,C);
        prog = piesos_setobj(prog,C); % Minimize C.
        fprintf(" --- Upper bound C to be minimised ---\n");
    else
        prog = piesos_program(x);
        fprintf(" --- Upper bound C is fixed ---\n");
    end
    
    
    %% Declare degrees of SOS distributed polynomials and define weighted Sobolev ball.
    
    % Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers (respectively).
    % Degree will be doubled when converted from quadratic to linear form.
    V_deg = dist_degs(1); p1_deg = dist_degs(2); p2_deg = dist_degs(3);
    
    % Declare monomial degrees in independent variables used to parametrize SOS LF and p1, p2 multipliers (respectively).
    V_mon = mon_degs(1); p1_mon = mon_degs(2); p2_mon = mon_degs(3);
    
    % Define weighted sobolev ball of radius r.
    g = Weighted_Sobolev_Ball(r, alpha, Top, x);
    fprintf(" --- Weighted Sobolev ball declared ---\n");

    % Define term related to lower and upper bounds.
    bound = Weighted_Sobolev_Ball(0.0, alpha, Top, x);
    

    %% Declare p1, p2 as SOS DPs.
    
    % p1_deg=2, p1_mon=4 takes approx. 2 mins to declare!
    [prog, p1] = SOS_DP(prog, p1_deg, p1_mon, x, dom);
    fprintf(" --- p1 declared ---\n");
    
    [prog, p2] = SOS_DP(prog, p2_deg, p2_mon, x, dom);
    fprintf(" --- p2 declared ---\n");
    

    %% Declare the LF (Not as SOS DP).
    [prog, V] = V3_DP(prog, V_deg, V_mon, Top, x, dom);
    fprintf(" --- LF declared and symmetry constraint imposed ---\n");
       

    %% Compute Lie derivative of the Lyapunov functional along the PIE.
    dV = Liediff(V,PIE);


    %% Define the lower bound on the LF and enforce constraint.
    
    % Define lower bound.
    V_low = V + eppos*bound; % bound term is already negated.
    
    % Enforce lower bound by defining and equating with new SOS DP.
    [prog, p3] = SOS_DP(prog, V_deg, V_mon, x, dom);
    prog = piesos_eq(prog, V_low-p3);
    fprintf(" --- Enforced lower bound equality ---\n");
    
    
    %% Define the upper bound on the LF and enforce constraint.

    % Define upper bound.
    V_up = -C*bound - V - p1*g; % bound term is already negated.
    
    % Enforce upper bound by defining and equating with new SOS DP.
    deg4 = max(V_deg,p1_deg+1);
    [prog, p4] = SOS_DP(prog, deg4, V_mon, x, dom);
    prog = piesos_eq(prog, V_up-p4);
    fprintf(" --- Enforced upper bound equality ---\n");
    

    %% Define constraint on the Lie derivative and enforce.
    
    % Define constraint.
    dV_con = -dV - 2*lambda*V - p2*g; 
    
    % Enforce constraint by defining and equating with new SOS DP.
    deg5 = max(V_deg,p2_deg+1);
    [prog, p5] = SOS_DP(prog, deg5, V_mon, x, dom);
    prog = piesos_eq(prog, dV_con-p5);
    fprintf(" --- Enforced Lie derivative equality ---\n");

    %% Solve the optimization program.
    
    sol_opts.simplify = true;
    prog_sol = piesos_solve(prog,sol_opts);


    % Extract the solution
    sol_info = prog_sol.solinfo.info;
    if sol_info.pinf || sol_info.dinf || sol_info.numerr || abs(sol_info.feasratio-1)>0.1
        res = [];
        fprintf(" --- PIESOS program was not solved ---\n")
    else
        if nargin < 8
            C_sol = piesos_getsol(prog_sol,C);
            C_sol = double(C_sol);
        else
            C_sol = C;
        end

        M_sol = sqrt(C_sol/eppos);
        res = [C_sol, M_sol];
        fprintf(" --- PIESOS program was solved ---\n")
    end

end

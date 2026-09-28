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
    % CR, 09/28/2026: Fixed bug with p1*g by creating polyopvar_times_v2.m
    % CR, 09/28/2026: Added positivity constraint to C.

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
        prog = piesos_ineq(prog,C-1e-6); % C>=0.
        prog = piesos_setobj(prog,C); % Minimize C.
        fprintf(" --- Upper bound C to be minimised ---\n");
    else
        prog = piesos_program(x);
        fprintf(" --- Upper bound C is fixed ---\n");
    end
    
    
    %% Unpack degrees of distributed polynomials and define weighted Sobolev ball.
    
    % Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers (respectively).
    % Degree will be doubled when converted from quadratic to linear form.
    V_deg = dist_degs(1); p1_deg = dist_degs(2); p2_deg = dist_degs(3);
    
    % Declare monomial degrees in independent variables used to parametrize SOS LF and p1, p2 multipliers (respectively).
    V_mon = mon_degs(1); p1_mon = mon_degs(2); p2_mon = mon_degs(3);
    
    % Define weighted sobolev ball of radius r.
    g     = Weighted_Sobolev_Ball(r, alpha, Top, x);
    g_deg = 2;
     
    % Define (negated) term related to lower and upper bounds.
    bound = Weighted_Sobolev_Ball(0.0, alpha, Top, x);
    
    fprintf(" --- Weighted Sobolev ball declared ---\n");


    %% Declare LF and enforce lower bound.
    [prog, V] = V3_DP(prog, V_deg, V_mon, Top, x, dom);

    % Define global lower bound on V.
    V_low = V + eppos*bound; % bound term is already negated.

    % Set lower bound by defining and equating with new SOS DP.
    [prog, sos1] = SOS_DP(prog, V_deg, V_mon, x, dom);
    prog = piesos_eq(prog, V_low-sos1);

    fprintf(" --- LF declared (without symmetry constraint) and lower bound equality set ---\n");

    %% Define the upper bound on the LF and enforce constraint.

    % Declare p1 as SOS DP.
    [prog, p1] = SOS_DP(prog, p1_deg, p1_mon, x, dom);

    % Define local upper bound on V.
    % V_up = -C*bound - V - p1*g; % bound term already negated.
    V_up = -C*bound - V - polyopvar_times_v2(p1,g); % bound term already negated.

    % Set upper bound by defining and equating with new SOS DP.
    deg = max(V_deg,p1_deg+1);
    g_mon = kernel_degree(g);
    g_mon = max(0,ceil((g_mon-1)/(2*g_deg)));
    mon = max(V_mon, p1_mon+g_mon);
    [prog, sos2] = SOS_DP(prog, deg, mon, x, dom);
    prog = piesos_eq(prog, V_up-sos2);
    fprintf(" --- p1 declared and LF upper bound equality set ---\n");


    %% Compute Lie derivative of the LF along the PIE and enforce constraint.
    
    % Declare p2 as SOS DP.
    % [prog, p2] = SOS_DP(prog, p2_deg, p2_mon, x, dom);

    % Compute and define upper bound on Lie derivative.
    % dV = Liediff(V,PIE);
    % dV_con = -dV - 2*lambda*V - p2*g; 
    % 
    % % Enforce constraint by defining and equating with new SOS DP.
    % deg = max(V_deg,p2_deg+1);
    % mon = max(V_mon, p2_mon+g_mon);
    % [prog, sos3] = SOS_DP(prog, deg, mon, x, dom);
    % prog = piesos_eq(prog, dV_con-sos3);
    % fprintf(" --- p2 declared and Lie derivative upper bound equality set ---\n");


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

function res = LocalStability(PDE, r, alpha, eppos, lambda, dist_degs, mon_degs, ...
                                psatz, active_conds)

    % res = Local_Stability(...) runs the local stability test (Thm. 1 in Automatica paper).
    %
    % INPUTS
    % - PDE:            The PDE to convert to a PIE.
    % - r:              Radius of the local ball.
    % - alpha:         (n+1)-dim array containing parameters of weighted Sobolev norm representing local ball.
    % - eppos:          Treated as eppos^2, the lower bound on LF.
    % - lambda:         Exponential decay rate.
    % - dist_degs:      array containing degrees of distributed monomial basis for V, p1, p2 resptively.
    % - mon_degs:       array containing degrees of monomial basis for V, p1, p2,
    %                   (and optionally sos3) resptively.
    % - psatz:          array containing indicators (0,1,2) which denote when to include psatz 
    %                   multipliers for sos1, sos2, sos3 resptively.
    % - active_conds:   Array containing logical values indicating which of the three SOS conditions to
    %                   enfore. All must be set to true for the complete stability test. Useful for 
    %                   determining appropriate degree values.
    % OUTPUTS
    % - res:            Cell containing C_sol, M_sol, V_sol, dV_sol, Pcell_sol.
    
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
    % CR, 10/05/2026: Replaced V3_DP.m and Liediff.m by joint function
    %                 V3_Liediff.m.
    % CR, 10/07/2026: Replaced kernel_degree by kernel_degree2 and added option
    %                 to pass in monomial degree of sos3 as mon_degs(4).
    % CR, 10/12/2026: Added active_conds as input for switching SOS constraints on/off.
    %                 Also replaced Weighted_Sobolev_Ball by Weighted_Sobolev_Norm.
    % CR, 10/14/2026: Added opts structure for constructing Pcell using
    %                 psatz multiplier.

    %% Convert PDE to PIE.
    PIE = convert(PDE);
    Top = PIE.T;    % Inverse map.
    f   = PIE.f;    % Polynomial PIE.
    dom = PIE.dom;  % Spatial domain [a,b].
    x   = f.vartab; % Fundamental state.

    
    %% Initialize PIESOS program structure.
    dpvar C % treated as C^2, the upper bound on the LF.
    prog = piesos_program(x,C);
    prog = piesos_ineq(prog,C-1e-6);
    prog = piesos_setobj(prog,C); % Minimize C.
    fprintf(" --- Upper bound C to be minimised ---\n");
    
    
    %% Unpack degrees of DPs and define weighted Sobolev ball.
    
    % Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers 
    % (respectively). Degree will be doubled when converted from quadratic 
    % to linear form.
    V_deg  = dist_degs(1);
    p1_deg = dist_degs(2);
    p2_deg = dist_degs(3);
    
    % Declare monomial degrees in independent variables used to parametrize
    % LF p1, p2 (respectively).
    V_mon    = mon_degs(1);
    p1_mon   = mon_degs(2);
    p2_mon   = mon_degs(3);

    % Declare psatz indicators for sos1, sos2, sos3 respectively.
    psatz1 = psatz(1);
    psatz2 = psatz(2);
    psatz3 = psatz(3);
     
    % Define term related to lower and upper bound on V.
    bound = Weighted_Sobolev_Norm(alpha, Top, x);
   
    % Define weighted sobolev ball of radius r.
    g = r^2 - bound;

    fprintf(" --- Weighted Sobolev ball and bound declared ---\n");


    %% Construct LF and its Lie derivative.
    % [prog, V, dV, Pcell] = V3_Liediff(prog, PIE, V_deg, V_mon);
    V  = -0.5 * innerprod(x, Top*x);
    dV = -innerprod(x, f);

    %% Define lower bound on the LF and enforce constraint.
    
    if active_conds(1)

        % Define global lower bound on V.
        V_low = V - eppos*bound;
    
        % Obtain kernel degree of V_low (not same as V_deg due to T) to use for 
        % defining sos1 of equal degree.
        V_low_mon = kernel_degree2(V_low);
    
        % Above are appropriate choices but computation requires limits to be placed 
        % on deg and mon. This could lead to lower order degrees in sos1 than V_low.
        sos1_mon = min(V_low_mon,6);
    
        fprintf(" --- sos1.deg = %d and sos1.mon = %d whilst V_low.deg = %d and V_low.mon = %d ---\n",V_deg, sos1_mon, V_deg, V_low_mon);
        
        % Set lower bound by defining and equating with new SOS DP.
        if psatz1 == 2
            [prog, sos1, ~, ~] = SOS_DP(prog, V_deg, sos1_mon, x, dom);
        elseif ~psatz1 || psatz1
            opts.psatz = psatz1;
            [prog, sos1, ~, ~] = SOS_DP(prog, V_deg, sos1_mon, x, dom, opts);
        else
            error("LocalStability: Invalid psatz entry!");
        end

        prog = piesos_eq(prog, V_low-sos1);
    
        fprintf(" --- LF lower bound equality set ---\n");

    end


    %% Define the upper bound on the LF and enforce constraint.
    
    if active_conds(2)

        % Declare p1 as SOS DP and define upper bound on V.
        if p1_mon>0
            [prog, p1, ~, ~] = SOS_DP(prog, p1_deg, p1_mon, x, dom);
            V_up = C*bound - V - polyopvar_times_v2(p1,g);
        elseif p1_mon==0
            V_up = C*bound - V;
        else
            error('LocalStability:p1_mon>=0.');
        end
    
        % Obtain kernel and distributed monomial degrees of V_up to use for defining
        % sos2 of equal degrees.
        V_up_deg = max(V_up.degmat); % degree of FDP in linear form - will always be even.
        V_up_deg = ceil(V_up_deg/2); % degree of FDP in quadratic form.
        V_up_mon = kernel_degree2(V_up);
    
        % Above are appropriate choices, but computation requires limits to be placed 
        % on deg and mon. This could lead to lower order degrees in sos2 than V_up.
        sos2_deg = min(V_up_deg,2);
        sos2_mon = min(V_up_mon,6);
    
        fprintf(" --- sos2.deg = %d and sos2.mon = %d whilst V_up.deg = %d and V_up.mon = %d ---\n",sos2_deg, sos2_mon, V_up_deg, V_up_mon);
    
        % Set upper bound by defining and equating with new SOS DP.
        if psatz2 == 2
            [prog, sos2, ~, ~] = SOS_DP(prog, V_up_deg, sos2_mon, x, dom);
        elseif ~psatz2 || psatz2
            opts.psatz = psatz2;
            [prog, sos2, ~, ~] = SOS_DP(prog, V_up_deg, sos2_mon, x, dom, opts);
        else
            error("LocalStability: Invalid psatz entry!");
        end
        
        prog = piesos_eq(prog, V_up-sos2);
    
        fprintf(" --- LF upper bound equality set ---\n");

    end


    %% Define upper bound on Lie derivative and enforce constraint.
    
    if active_conds(3)

        % Declare p2 as SOS DP and define local upper bound on dV.
        if p2_mon>0
            [prog, p2, ~, ~] = SOS_DP(prog, p2_deg, p2_mon, x, dom);
            dV_up = -dV - 2*lambda*V - p2*g; 
        elseif p2_mon==0
            dV_up = -dV - 2*lambda*V; 
        else
            error('LocalStability:p2_mon>=0.');
        end
    
        % Obtain kernel and distributed monomial degrees of V_up to use for defining
        % sos2 of equal degrees.
        dV_up_deg = max(dV_up.degmat); % degree of FDP in linear form - will always be even.
        dV_up_deg = ceil(dV_up_deg/2); % degree of FDP in quadratic form.
        dV_up_mon = kernel_degree2(dV_up);
    
        % Above are appropriate choices, but computation requires limits to be placed 
        % on deg and mon. This could lead to lower order degrees in sos3 than dV_up.
        sos3_deg = min(dV_up_deg,2);
        if size(mon_degs,2) == 4
            sos3_mon = mon_degs(4);
        else
            sos3_mon = min(dV_up_mon,10);
        end
    
        fprintf(" --- sos3.deg = %d and sos3.mon = %d whilst dV_up.deg = %d and dV_up.mon = %d ---\n",sos3_deg, sos3_mon, dV_up_deg, dV_up_mon);
    
        % Set upper bound by defining and equating with new SOS DP.
        if psatz3 == 2
            [prog, sos3, ~, ~] = SOS_DP(prog, sos3_deg, sos3_mon, x, dom);
        elseif ~psatz3 || psatz3
            opts.psatz = psatz3;
            [prog, sos3, ~, ~] = SOS_DP(prog, sos3_deg, sos3_mon, x, dom, opts);
        else
            error("LocalStability: Invalid psatz entry!");
        end     

        prog = piesos_eq(prog, dV_up-sos3);
    
        fprintf(" --- Lie derivative upper bound equality set ---\n");

    end


    %% Solve the optimization program.

    sol_opts.simplify = true;
    sol_opts.solver   = 'mosek'; % 'mosek', 'sdpt3'
    prog_sol = piesos_solve(prog,sol_opts);

    % Extract the solution
    sol_info = prog_sol.solinfo.info;
    if sol_info.pinf || sol_info.dinf || sol_info.numerr || abs(sol_info.feasratio-1)>0.1
        res = {};
        fprintf(" --- PIESOS program was not solved ---\n")
    else
        C_sol = piesos_getsol(prog_sol,C);
        C_sol = double(C_sol);
        M_sol     = sqrt(C_sol/eppos);
        V_sol     = piesos_getsol(prog_sol,V);
        dV_sol    = piesos_getsol(prog_sol,dV);
        res       = {C_sol, M_sol, V_sol, dV_sol};
        fprintf(" --- PIESOS program was solved ---\n")
    end

end

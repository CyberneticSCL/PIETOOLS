function res = LocalStability_v2(PDE, r, alpha, eppos, lambda, dist_degs, mon_degs, ...
                                psatz, active_conds, sos_degs, M)

    % RES = LOCALSTABILITY_V2(PDE,R,ALPHA,EPPOS,LAMBDA,DIST_DEGS,MON_DEGS,
    % PSATZ,ACTIVE_CONDS,SOS_DEGS,M) runs the local stability test
    % (Thm. 1 in Automatica paper) for the PDE, on the ball of radius r in the
    % weighted Sobolev norm with weights alpha. A Lyapunov functional V
    % and its Lie derivative dV along the PIE associated to the PDE are
    % declared using `piesos_LF`, and the constant C in the upper bound
    % on V is minimized subject to the conditions, declared using
    % `piesos_ineq`,
    %   sos1:  V - N_eppos >= 0,
    %   sos2:  C*N_M - V - p1*g >= 0,
    %   sos3:  -dV - 2*lambda*V - p2*g >= 0,
    % where N_w = sum_i sum_j w(i,j+1)*||d^j/ds^j (T*x)_i||^2 is the
    % squared weighted Sobolev norm with weights w, g = r^2 - N_alpha, and
    % p1, p2 are SOS multipliers.
    %
    % INPUTS
    % - PDE:            'pde_struct' object defining the PDE, to be
    %                   converted to a PIE. All PDE states must be scalar;
    % - r:              positive scalar, radius of the local ball;
    % - alpha:          m x (n+1) array of weights, with alpha(i,j+1) the
    %                   weight of the jth derivative of the ith PDE state
    %                   in the weighted Sobolev norm defining the local
    %                   ball, for m the number of states and n the maximal
    %                   order of the PDE. The states are ordered as they
    %                   were declared in the PDE. Derivatives of the
    %                   states with a zero weight are never taken, so that
    %                   states of order 0 can be included. For m=1, a
    %                   vector of n+1 weights can be used;
    % - eppos:          scalar or array of the same size as alpha,
    %                   specifying the weights of the lower bound on the
    %                   LF, V(x) >= N_eppos(x). If scalar, eppos*alpha is
    %                   used as the weights;
    % - lambda:         scalar, exponential decay rate;
    % - dist_degs:      1 x 3 array with the maximal degrees of the
    %                   distributed monomial bases of V, p1 and p2,
    %                   respectively. The basis of sos1, sos2 and sos3 is
    %                   selected by `piesos_ineq`;
    % - mon_degs:       1 x 3 or 1 x 4 array with the degrees of the
    %                   spatial monomials parameterizing the operators in
    %                   V, p1 and p2, and (optionally) sos3,
    %                   respectively. mon_degs(4) is the degree in each
    %                   factor of a distributed monomial;
    % - psatz:          (optional) 1 x 3 array of nonnegative integers
    %                   with the orders of the psatz multipliers to use
    %                   in sos1, sos2 and sos3, respectively, as in
    %                   `piesos_sosvar`. If psatz(i)=0, no multipliers
    %                   are used; if psatz(i)=q>0, multipliers
    %                   (s-a)^k*(b-s)^(q-k), k=0,...,q are used.
    %                   Defaults to [0,0,0];
    % - active_conds:   (optional) 1 x 3 logical array indicating which of
    %                   the conditions sos1, sos2 and sos3 to enforce.
    %                   All must be true for the complete stability
    %                   test; other choices are useful for determining
    %                   appropriate degrees. Defaults to true(1,3);
    % - sos_degs:       (optional) 1 x 3 array with the degrees of the
    %                   monomials in the operators of sos1, sos2 and
    %                   sos3, respectively, as opts.degs in
    %                   `piesos_ineq`. If sos_degs(i) is NaN, the degree
    %                   is selected from the kernels of the functional.
    %                   Takes precedence over mon_degs(4). Defaults to
    %                   NaN(1,3);
    % - M:              (optional) scalar or array of the same size as
    %                   alpha, specifying the weights of the upper bound
    %                   on the LF, V(x) <= C*N_M(x). If scalar, M*alpha
    %                   is used as the weights. Defaults to 1;
    %
    % OUTPUTS
    % - res:            1 x 5 cell {C_sol,M_sol,V_sol,dV_sol,rho} with the
    %                   optimal value C_sol of C, the gain M_sol, the
    %                   solved 'polyopvar' objects V_sol and dV_sol, and the
    %                   radius rho of the certified region of attraction.
    %                   An empty cell if the program is not solved. With
    %                   the norms ||.||_eppos and ||.||_M of the weights
    %                   eppos and M, the solution satisfies
    %                   ||x(t)||_eppos <= sqrt(C)*exp(-lambda*t)*
    %                   ||x(0)||_M. If eppos and M are both scalar, M_sol =
    %                   sqrt(C_sol*M/eppos) is the gain in the norm of
    %                   alpha, otherwise M_sol = sqrt(C_sol).
    %                   The radius rho = r*sqrt(gam/C_sol) is such that
    %                   every solution with ||x(0)||_M < rho remains in the
    %                   ball N_alpha <= r^2 and converges to zero, with
    %                   gam = min_{alpha(i,j)>0} eppos(i,j)/alpha(i,j) (so
    %                   that V >= N_eppos >= gam*N_alpha). Here, V(x(0)) <=
    %                   C*N_M(x(0)) < gam*r^2 <= min V on the boundary
    %                   N_alpha=r^2. If gam=0, then rho=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PIETOOLS - LocalStability_v2
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
    % CR, 09/28/2026: Fixed bug with p1*g by creating polyopvar_times_v2.m
    % CR, 09/28/2026: Added positivity constraint to C.
    % CR, 10/05/2026: Replaced V3_DP.m and Liediff.m by joint function
    %                 V3_Liediff.m.
    % CR, 10/07/2026: Replaced kernel_degree by kernel_degree2 and added
    %                   option to pass in monomial degree of sos3 as
    %                   mon_degs(4).
    % CR, 10/12/2026: Added active_conds as input for switching SOS
    %                   constraints on/off. Also replaced
    %                   Weighted_Sobolev_Ball by Weighted_Sobolev_Norm.
    % CR, 10/14/2026: Added opts structure for constructing Pcell using
    %                 psatz multiplier.
    % DJ, 10/14/2026: Alternative formulation using new `piesos_LF`
    %                   function.
    % DJ, 10/15/2026: Updated to use new polyopvar functions;
    % DJ, 10/21/2026: Add support for vector-valued case;
    % DJ, 10/21/2026: Certified radius rho as 5th output; progress times.

    t_start = tic;

    %% Set the default options.
    if nargin<8 || isempty(psatz)
        psatz = [0,0,0];
    end
    if nargin<9 || isempty(active_conds)
        active_conds = true(1,3);
    end
    if nargin<10 || isempty(sos_degs)
        sos_degs = NaN(1,3);
    end
    if nargin<11 || isempty(M)
        M = 1;
    end
    % Weights of the lower and upper bound on V: arrays like alpha, or scalar multiples of alpha
    is_scalar_bounds = isscalar(eppos) && isscalar(M);
    if isvector(alpha)
        alpha = reshape(alpha,1,[]);    % scalar PDE state: row of weights
    end
    if isscalar(eppos)  
        eppos_w = eppos*alpha;      
    else   
        eppos_w = eppos;    
    end
    if isscalar(M)
        M_w = M*alpha;              
    else
        M_w = M;     
    end
    if isvector(eppos_w) && size(alpha,1)==1
        eppos_w = reshape(eppos_w,1,[]);
    end
    if isvector(M_w) && size(alpha,1)==1
        M_w = reshape(M_w,1,[]);
    end
    if ~isequal(size(eppos_w),size(alpha)) || ~isequal(size(M_w),size(alpha))
        error('LocalStability: eppos and M must be scalars or of the same size as alpha!');
    end
    if any(psatz<0) || any(mod(psatz,1)~=0)
        error("LocalStability: Invalid psatz entry!");
    end

    %% Convert PDE to PIE.
    PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');                     % DJ, 10/14/2026
    Top = PIE.T;    % Inverse map.
    f   = PIE.f;    % Polynomial PIE.
    dom = PIE.dom;  % Spatial domain [a,b].
    x   = PIE.x;    % Fundamental state.

    % The rows of T follow the order of the PIE states, which differs from
    % the order in which the PDE states were declared: row i of T*x is PDE
    % state x_tab(i,1). Reorder the weights accordingly.
    if any(PIE.x_tab(:,2)~=1) || size(alpha,1)~=size(PIE.x_tab,1)
        error('LocalStability: alpha must have one row for each of the (scalar) PDE states!');
    end
    eppos_user = eppos_w;   % declaration order, used for the radius rho
    alpha_w = alpha(PIE.x_tab(:,1),:);
    eppos_w = eppos_w(PIE.x_tab(:,1),:);
    M_w     = M_w(PIE.x_tab(:,1),:);


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
    p1_deg = dist_degs(2);
    p2_deg = dist_degs(3);

    % Declare monomial degrees in independent variables used to parametrize
    % LF p1, p2 (respectively).
    p1_mon   = mon_degs(2);
    p2_mon   = mon_degs(3);

    % Declare psatz orders for sos1, sos2, sos3 respectively.
    psatz1 = psatz(1);
    psatz2 = psatz(2);
    psatz3 = psatz(3);

    % Define term related to lower and upper bound on V.
    bound = Weighted_Sobolev_Norm(alpha_w, Top, x);

    % Terms related to the lower and upper bound on V.
    bound_low = Weighted_Sobolev_Norm(eppos_w, Top, x);                     % DJ, 10/14/2026
    bound_up  = Weighted_Sobolev_Norm(M_w, Top, x);                         % DJ, 10/14/2026

    % Define weighted sobolev ball of radius r.
    g = r^2 - bound;

    fprintf(" --- Weighted Sobolev ball and bound declared ---\n"); fprintf("     (%.0f s elapsed)\n",toc(t_start));


    %% Construct LF and its Lie derivative.
    V_deg = dist_degs(1);
    V_mon = mon_degs(1);
    LF_opts.degs = V_mon;
    [prog, V, dV, Pcell] = piesos_LF(prog, PIE, 'pie', (1:V_deg)', LF_opts);

    %% Define lower bound on the LF and enforce constraint.

    if active_conds(1)

        % Define global lower bound on V.
        V_low = V - bound_low;

        % Set lower bound by enforcing V_low >= 0 as SOS DP: the basis of the SOS DP
        % and the degrees of its operators follow from the kernels of V_low.
        opts1 = struct();
        if psatz1>0
            opts1.psatz = psatz1;
        end
        if ~isnan(sos_degs(1))
            opts1.degs = sos_degs(1);
        end
        [prog, sos1, ~, ~, Z1] = piesos_ineq(prog, V_low, opts1);           % DJ, 10/14/2026

        fprintf(" --- sos1 has %d distributed monomials in its basis ---\n",size(Z1.Z,1));
        fprintf(" --- LF lower bound set ---\n"); fprintf("     (%.0f s elapsed)\n",toc(t_start));

    end


    %% Define the upper bound on the LF and enforce constraint.

    if active_conds(2)

        % Declare p1 as SOS DP and define upper bound on V.
        if p1_mon>0
            [prog, p1, ~, ~] = SOS_DP_v2(prog, p1_deg, p1_mon, x, dom);
            V_up = C*bound_up - V - times(p1,g);
        elseif p1_mon==0
            V_up = C*bound_up - V;
        else
            error('LocalStability:p1_mon>=0.');
        end

        % Set upper bound by enforcing V_up >= 0 as SOS DP, with the basis and
        % operator degrees following from the kernels of V_up.
        opts2 = struct();
        if psatz2>0
            opts2.psatz = psatz2;
        end
        if ~isnan(sos_degs(2))
            opts2.degs = sos_degs(2);
        end
        [prog, sos2, ~, ~, Z2] = piesos_ineq(prog, V_up, opts2);            % DJ, 10/14/2026

        fprintf(" --- sos2 has %d distributed monomials in its basis ---\n",size(Z2.Z,1));
        fprintf(" --- LF upper bound set ---\n"); fprintf("     (%.0f s elapsed)\n",toc(t_start));

    end


    %% Define upper bound on Lie derivative and enforce constraint.

    if active_conds(3)

        % Declare p2 as SOS DP and define local upper bound on dV.
        if p2_mon>0
            [prog, p2, ~, ~] = SOS_DP_v2(prog, p2_deg, p2_mon, x, dom);
            dV_up = -dV - 2*lambda*V - times(p2,g);
        elseif p2_mon==0
            dV_up = -dV - 2*lambda*V;
        else
            error('LocalStability:p2_mon>=0.');
        end

        % Set upper bound by enforcing dV_up >= 0 as SOS DP, with the basis and
        % operator degrees following from the kernels of dV_up, unless the
        % degree of the operators is specified as mon_degs(4).
        opts3 = struct();
        if psatz3>0
            opts3.psatz = psatz3;
        end
        if size(mon_degs,2) == 4
            opts3.degs = mon_degs(4);
        end
        if ~isnan(sos_degs(3))
            opts3.degs = sos_degs(3);
        end
        [prog, sos3, ~, ~, Z3] = piesos_ineq(prog, dV_up, opts3);           % DJ, 10/14/2026

        fprintf(" --- sos3 has %d distributed monomials in its basis ---\n",size(Z3.Z,1));
        fprintf(" --- Lie derivative upper bound set ---\n"); fprintf("     (%.0f s elapsed)\n",toc(t_start));

    end


    %% Solve the optimization program.

    sol_opts.simplify = true;
    sol_opts.solver   = 'mosek'; % 'sedumi', 'sdpt3'
    prog_sol = piesos_solve(prog,sol_opts);

    % Extract the solution
    sol_info = prog_sol.solinfo.info;
    if sol_info.pinf || sol_info.dinf || sol_info.numerr || abs(sol_info.feasratio-1)>0.1
        res = {};
        fprintf(" --- PIESOS program was not solved ---\n")
    else
        C_sol = piesos_getsol(prog_sol,C);
        C_sol = double(C_sol);
        if is_scalar_bounds
            M_sol = sqrt(C_sol*M/eppos);
        else
            M_sol = sqrt(C_sol);
        end
        V_sol     = piesos_getsol(prog_sol,V);
        dV_sol    = piesos_getsol(prog_sol,dV);
        % Radius of the certified region of attraction, in the norm N_M
        pos_alpha = alpha>0;
        gam       = min(eppos_user(pos_alpha)./alpha(pos_alpha));
        rho       = r*sqrt(max(gam,0)/C_sol);
        res       = {C_sol, M_sol, V_sol, dV_sol, rho};
        fprintf(" --- Certified radius rho = %g (in N_M norm) ---\n",rho);
        fprintf(" --- PIESOS program was solved ---\n")
    end

end

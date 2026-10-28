function res = LocalStability_wave(mu, a, k, rad, eppos, lambda, C_fix, sol_opts)

    % RES = LOCALSTABILITY_WAVE(MU,A,K,RAD,EPPOS,LAMBDA,C_FIX,SOL_OPTS)
    % runs a local stability test for the semilinear damped wave equation
    % with Robin boundary feedback,
    %   z_tt = z_ss - 2*mu*z_t + a*z + z^3,   s in [0,1],
    %   z(0) = 0,   z_s(1) = -k*z(1),
    % in the states z and y=z_t, on the energy ball
    %   N(x) = ||z_s||^2 + ||y||^2 <= rad^2,
    % using `piesos_LF` and `piesos_ineq`, as `LocalStability_v2`. In
    % contrast with `LocalStability_v2`, the multipliers p*g (g=rad^2-N)
    % are replaced by pointwise multipliers, and V gets a quartic term,
    %   V = V2 - cc/4*int z^4,   dV = dV2 - cc*int z^3*y,
    %   S(a_y,a_z) = int (a_y*y^2 + a_z*z^2)*(rad^2 - z^2) ds,
    % for V2 quadratic, cc>=0, and a_y, a_z>=0. S>=0 wherever |z|<=rad, in
    % particular on the ball, as |z|<=||z_s|| for z(0)=0. The conditions
    % are
    %   sos1:  V - eppos*N - S(p1y,p1z) >= 0,
    %   sos2:  C*N - V >= 0,
    %   sos3:  -dV - 2*lambda*V - S(q1,q2) >= 0,
    % with C minimized (or fixed).
    %
    % INPUTS
    % - mu,a,k:     scalars, damping, reaction coefficient and feedback
    %               gain of the wave equation (example: mu=2, a=5, k=3);
    % - rad:        positive scalar, radius of the energy ball;
    % - eppos:      (optional) positive scalar, V >= eppos*N on the ball.
    %               Defaults to 0.05;
    % - lambda:     (optional) nonnegative scalar, exponential decay
    %               rate. Defaults to 0.05;
    % - C_fix:      (optional) scalar. If nonempty, C is fixed to this
    %               value, and only the feasibility is tested (MOSEK
    %               stalls when minimizing C: numerr=2). Defaults to [];
    % - sol_opts:   (optional) struct with options passed to
    %               `piesos_solve`. Defaults to
    %               struct('solver','mosek','simplify',true);
    %
    % OUTPUTS
    % - res:        1 x 4 cell {C_sol,rho,V_sol,dV_sol} with the value
    %               C_sol of C, the radius rho=rad*sqrt(eppos/C_sol) such
    %               that every solution with ||(z_s,y)(0)||<rho remains in
    %               the ball and converges to zero, and the solved
    %               'polyopvar' objects V and dV. An empty cell if the
    %               program is not solved;
    %
    % NOTES
    % - The operators in V have degree 2, and those in sos1, sos2, sos3
    %   have degree 3, 2 and 4: sos1 is infeasible for degree 2, and sos3
    %   does not solve cleanly for degree 3. The program has about 380k x
    %   28k variables, and takes about 9 minutes and 8 GB (MOSEK).
    % - With the multiplier p*g, and V quadratic, sos3 is infeasible, as
    %   the cross term int (Px)_y*z^3 of dV (y*z^3) cannot be dominated:
    %   products of independent integrals, as in p*g, are not Gram terms
    %   of `piesos_ineq` (one primary variable), unlike pointwise terms.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PIETOOLS - LocalStability_wave
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
    % DJ, 10/28/2026: Initial coding;

    if nargin<5 || isempty(eppos),      eppos = 0.05;           end
    if nargin<6 || isempty(lambda),     lambda = 0.05;          end
    if nargin<7,                        C_fix = [];             end
    if nargin<8 || isempty(sol_opts),   sol_opts = struct('solver','mosek','simplify',true);    end

    sos_degs = [3,2,4];     % degrees of the operators in sos1, sos2, sos3
    V_mon = 2;              % degree of the operators in V2

    %% Convert the PDE to a PIE.
    pvar s t
    dom = [0,1];
    u1 = pde_var(s,dom);    % z
    u2 = pde_var(s,dom);    % y = z_t
    PDE = [diff(u1,t) == u2;
           diff(u2,t) == diff(u1,s,2) - 2*mu*u2 + a*u1 + u1^3;
           subs(u1,s,dom(1)) == 0;
           subs(diff(u1,s),s,dom(2)) + k*subs(u1,s,dom(2)) == 0];
    PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
    Top = PIE.T;    x = PIE.x;
    % The rows of T follow the order of the PIE states
    iz = find(PIE.x_tab(:,1)==1);       iy = find(PIE.x_tab(:,1)==2);
    z  = Top(iz,:)*x;                   y  = Top(iy,:)*x;

    %% Pointwise functionals, products at the same point s.
    z2 = z.*z;                          z3 = z2.*z;
    zz = innerprod(z,z);                yy = innerprod(y,y);
    z4 = innerprod(z2,z2);              y2z2 = innerprod(y.*z,y.*z);
    % Energy N = ||z_s||^2 + ||y||^2, with weights per PIE state
    alpha_w = [0,1,0; 1,0,0];           alpha_w = alpha_w(PIE.x_tab(:,1),:);
    N = Weighted_Sobolev_Norm(alpha_w,Top,x);
    % S(ay,az) = int (ay*y^2+az*z^2)*(rad^2-z^2) ds
    Sfun = @(ay,az) rad^2*(ay*yy + az*zz) - (ay*y2z2 + az*z4);

    %% Initialize the PIESOS program.
    dpvar C cc p1y p1z q1 q2
    prog = piesos_program(x,[C;cc;p1y;p1z;q1;q2]);
    prog = piesos_ineq(prog,C-1e-6);
    prog = piesos_setobj(prog,C);       % Minimize C.
    if ~isempty(C_fix)
        prog = soseq(prog,C-C_fix);
    end
    for v = {cc,p1y,p1z,q1,q2}
        prog = piesos_ineq(prog,v{1});
    end

    %% Lyapunov functional and its Lie derivative.
    [prog,V2,dV2] = piesos_LF(prog,PIE,'pie',1,struct('degs',V_mon));
    V  = V2 - (cc/4)*z4;
    dV = dV2 - cc*innerprod(z3,y);      % d/dt int z^4/4 = int z^3*y, as z_t = y
    fprintf(" --- LF and Lie derivative declared ---\n");

    %% Conditions sos1, sos2, sos3.
    prog = piesos_ineq(prog, V - eppos*N - Sfun(p1y,p1z), struct('degs',sos_degs(1)));
    prog = piesos_ineq(prog, C*N - V, struct('degs',sos_degs(2)));
    prog = piesos_ineq(prog, -dV - 2*lambda*V - Sfun(q1,q2), struct('degs',sos_degs(3)));

    %% Solve the program.
    prog_sol = piesos_solve(prog,sol_opts);
    sol_info = prog_sol.solinfo.info;
    if sol_info.pinf || sol_info.dinf || sol_info.numerr || abs(sol_info.feasratio-1)>0.1
        res = {};
        fprintf(" --- PIESOS program was not solved ---\n")
    else
        C_sol = double(piesos_getsol(prog_sol,C));
        rho   = rad*sqrt(eppos/C_sol);
        res   = {C_sol, rho, piesos_getsol(prog_sol,V), piesos_getsol(prog_sol,dV)};
        fprintf(" --- Certified radius rho = %g (in the energy norm) ---\n",rho);
        fprintf(" --- PIESOS program was solved ---\n")
    end

end

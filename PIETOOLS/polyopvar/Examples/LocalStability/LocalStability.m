function [g,p1] = LocalStability(PDE, r, alpha, eppos, lambda, dist_degs, mon_degs)

    % [C, M] = Local_Stability(...) runs the local stability test (Thm. 1 in Automatica paper).
    %
    % INPUTS
    % - PDE:            The PDE to convert to a PIE.
    % - r:              Radius of the local ball.
    % - alpha:          (n+1)-dim array containing parameters of weighted Sobolev ball.
    % - eppos:          Lower bound on SOS LF.
    % - lambda:         Exponential decay rate.
    % - dist_degs:      3d array containing degrees of distributed monomial basis for V, p1, p2 resptively.
    % - mon_degs:       3d array containing degrees of monomial basis for V, p1, p2 resptively.
    %
    % OUTPUTS
    % - C:              Upper bound on SOS LF.
    % - M:              Scalar multiplier of exponential upper bound on solution norm.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Convert PDE to PIE.
    PIE = convert(PDE);
    Top = PIE.T; % Inverse map.
    f   = PIE.f; % Polynomial PIE.
    dom = PIE.dom; % Spatial domain [a,b].
    
    % Construct fundamental state.
    x = f.vartab;
    % pvar s
    % x = polyopvar();
    % x.varname = f.varname;
    % x.varsize = 1;
    % x.degmat = 1;
    % x.C.ops = 1;
    % x.pvarname = f.pvarname;
    % x.dom = dom;
    % x.varmat = 1;

    
    %% Initialize PIESOS program structure.
    
    dpvar C % upper bound gamma on the LF.
    prog = piesos_program(x,C);
    prog = piesos_setobj(prog,C); % Minimize C.
    
    
    %% Declare degrees of SOS distributed polynomials and define weighted Sobolev ball.
    
    % Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers (respectively).
    % Degree will be doubled when converted from quadratic to linear form.
    V_deg = dist_degs(1); p1_deg = dist_degs(2); p2_deg = dist_degs(3);
    
    % Declare monomial degrees in independent variables used to parametrize SOS LF and p1, p2 multipliers (respectively).
    V_mon = mon_degs(1); p1_mon = mon_degs(2); p2_mon = mon_degs(3);
    
    % Define weighted sobolev ball of radius r.
    g = Weighted_Sobolev_Ball(r, alpha, Top, x);


    %% Declare p1, p2 as SOS DPs.
    [prog, p1] = SOS_DP(prog, p1_deg, p1_mon, x, dom);
    
    % [prog, p2] = SOS_DP(prog, p2_deg, p2_mon, x, dom);
    % 
    % 
    % %% Declare the LF as SOS DP.
    % [prog, V] = LF_V3(prog, V_deg, V_mon, Top, x, dom);
    % 
    % 
    % %% Compute Lie derivative of the Lyapunov functional along the PIE.
    % 
    % 
    % %% Equality and positivity conditions


end

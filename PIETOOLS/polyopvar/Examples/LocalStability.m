function [PIE, Top, f] = LocalStability(PDE, r, alpha, eppos, lambda, dist_degs, mon_degs)

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
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Convert PDE to PIE.
PIE = convert(PDE);
Top = PIE.T;
f = PIE.f;
v = f.vartab;
% Top_opvar = ndopvar2dopvar(Top);
% Rop = dopvar2ndopvar(diff(Top_opvar,Top_opvar.var1,1,'pure'));


%% Declare degrees of SOS distributed polynomials and define weighted Sobolev ball.


% Declare degrees of dist mon basis for SOS LF and p1, p2 multipliers (respectively).
% Degree will be doubled when converted from quadratic to linear form.
V_deg = dist_degs(1); p1_deg = dist_degs(2); p2_deg = dist_degs(3);

% Declare monomial degrees in independent variables used to parametrize SOS LF and p1, p2 multipliers (respectively).
V_mon = mon_degs(1); p1_mon = mon_degs(2); p2_mon = mon_degs(3);

% Set distributed monomial basis of SOS LF and p1, p2 multipliers (respectively).
V_Z  = dmonomials(v,(1:V_deg));
p1_Z = dmonomials(v,(1:p1_deg));
p2_Z = dmonomials(v,(1:p2_deg));

% % Set monomial basis of SOS LF and p1, p2 multipliers (respectively).
% V_U  = V_mon;
% p1_U = p1_mon;
% p2_U = p2_mon;


% Define weighted sobolev ball of radius r.
g = weighted_sobolev_ball(r, alpha, Top, v);


%% Initialize PIESOS program structure.

dpvar C % upper bound gamma on the LF.
prog = piesos_program(v,C);
prog = piesos_setobj(prog,C); % Minimize C.


%% Declare the Lyapunov functional.


%% Compute Lie derivative of the Lyapunov functional along the PIE.

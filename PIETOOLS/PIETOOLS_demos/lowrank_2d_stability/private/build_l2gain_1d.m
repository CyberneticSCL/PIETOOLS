function [prog,H] = build_l2gain_1d(PIE,st,dual)                            % CC, 09/23/2026
% BUILD_L2GAIN_1D  1-D Hinf-gain (L2-gain) LPI, returning OPERATOR HANDLES.
%
% Same LPI as executives/PIETOOLS_Hinf_gain.m (dual = false) and
% PIETOOLS_Hinf_gain_dual.m (dual = true), minus the solve.
%
% WHY THIS EXECUTIVE IS WORTH ADAPTING.  Three reasons, all structural:
%  (a) the example library publishes reference gamma values, so a certificate
%      can be scored on ACCURACY (gamma against the reference) instead of on a
%      binary gate -- which is what a benchmark for solver changes needs;
%  (b) gamma is a free scalar decision variable, so the program has Kf > 0.
%      restrict_solve and bm_resid both plumb a free block and neither has
%      ever been exercised with one (restrict_solve.m's own comment: "Kf = 0
%      in every program this package builds");
%  (c) on a face the restricted problem becomes min gamma s.t. M s = b,
%      S >= 0 -- a genuine tiny SDP with an objective.  The feasibility
%      restriction has no objective to pin its solution, which is why
%      restrict_solve has to take a sign-free determined solve and then
%      project onto the PSD cone; with an objective that hack is unnecessary.
%
% SOUNDNESS IS UNCHANGED, AND IT IS AN UPPER BOUND.  Restricting to a face can
% only lose feasible points, so any S >= 0 satisfying the restricted system
% gives a feasible point of the original LPI, hence a valid UPPER bound on the
% true optimal gamma.  A face that fails says nothing about the true gamma.
%
% THE LAYOUT TRAP THIS EXECUTIVE EXPOSES.  lpivar(prog,Top.dim,Qdeg) stores Qop
% as 'poly' entries in sos.var, INTERLEAVED between the 'sos' Gram blocks.
% MEASURED on Ex_Transport_Eq_with_Disturbance at 'light':
%   x(1) gam free | x(2..101) sos N=10 | x(102..144) poly free | x(145..433)
%   sos N=17 | x(434..497) sos N=8.
% raw_data's contiguous model returns Kf=1, Ns=[10 17 8] and so accounts for
% 454 of 497 coordinates.  Use pielr_rawdata, which reads the map from the
% program.  lpi_ineq(prog,gam) also produces an expr of type 'ineq' that
% raw_data would concatenate as an equality, pinning gamma instead of
% bounding it; pielr_rawdata separates those too.
%
% INPUT   PIE  pie_struct with T, A, Bw, Cz, Dzw (PIE.dim == 1)
%         st   1-D settings struct, as lpisettings(...) returns
%         dual optional, default false
% OUTPUT  prog unsolved LPI program
%         H    .Top .Aop .Bwop .Czop .Dzwop .Iw .Iz .Rop .Qop .Dop .Deop
%              .gam .dual .st .PIE

if nargin<3 || isempty(dual), dual = false; end
PIE = initialize(PIE);
Top = PIE.T;    Twop = PIE.Tw;
Aop = PIE.A;    Bwop = PIE.Bw;
Czop = PIE.Cz;  Dzwop = PIE.Dzw;

% The executive silently reroutes to the coercive variant when there is a
% disturbance at the boundary, which is a DIFFERENT LPI with a different
% residual expression.  Refuse instead of adapting the wrong one.
if ~(Twop==0)
    error('build_l2gain_1d:Tw', ...
        ['PIE.Tw is nonzero (disturbance at the boundary).  ' ...
         'PIETOOLS_Hinf_gain reroutes this case to PIETOOLS_Hinf_gain_coercive, ' ...
         'a different LPI; adapt that one rather than using this builder.']);
end
if isfield(st,'sosineq_on') && st.sosineq_on
    error('build_l2gain_1d:sosineq', ...
          'sosineq_on=1 has no Deop slack; the operator gate needs the equality route.');
end

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
% GAMFIX: build with a NUMERIC gamma, turning the optimisation into a pure
% feasibility test.  PIETOOLS_Hinf_gain documents this route itself ("a
% specific gain test specified by defining a specific desired value of gamma
% ... results in a feasibility test instead of an optimization problem").
%
% WHY IT MATTERS, measured by the parallel cuADMM session on this same
% executive and plant family: with the objective active the interior-point
% solve reaches rel_b 2.058e-07 (light) / 1.523e-07 (heavy), while bisection
% at fixed gamma on the same programs reaches 3.3e-08 to 7.8e-08 -- about an
% ORDER OF MAGNITUDE better equality residual for a gamma agreeing to 1.9e-05.
% With c = 0 the dual condition A'y + z = c is trivially satisfiable and the
% gap is ~0 for free, so the solver spends its effort on primal feasibility;
% with c nonzero the iterates trade feasibility against optimality.
%
% Substituting rather than pinning by an appended row: pinning would leave
% the gamma >= 0 cone block and the objective in the data while making that
% block redundant, and redundant equality rows are a known conditioning
% hazard in this stack.  Substitution removes the block instead.  Note b then
% DEPENDS on gamma (Dop carries -gam*Iw and -gam*Iz as constant blocks), so
% each fixed-gamma program is differently scaled -- harmless for a relative
% residual, misleading if absolute residuals are compared across a bracket.
if isfield(st,'gamfix') && ~isempty(st.gamfix)
    gam = st.gamfix;                 % a plain number: no decvar, no objective
else
    dpvar gam;
    prog = lpidecvar(prog, gam);
    prog = lpi_ineq(prog, gam);      % expr of type 'ineq'; see header
    prog = lpisetobj(prog, gam);
end

% ---- positive operator Rop, and the free operator Qop it couples to ------
[prog, R1op] = poslpivar(prog, Top.dim, st.dd1, st.options1);
if st.override1~=1
    [prog, P2op] = poslpivar(prog, Top.dim, st.dd12, st.options12);
    Rop = R1op + P2op;
else
    Rop = R1op;
end
Qdeg = get_lpivar_degs(Rop,Top);
[prog, Qop] = lpivar(prog,Top.dim,Qdeg);
if ~dual
    prog = lpi_eq(prog, Top'*Qop-Rop);        % PIETOOLS_Hinf_gain
else
    prog = lpi_eq(prog, Top*Qop-Rop);         % PIETOOLS_Hinf_gain_dual
end

% ---- the negativity block -----------------------------------------------
Iw = mat2opvar(eye(size(Bwop,2)), Bwop.dim(:,2), PIE.vars, PIE.dom);
Iz = mat2opvar(eye(size(Czop,1)), Czop.dim(:,1), PIE.vars, PIE.dom);
if ~dual
    Dop = [-gam*Iw,       Dzwop',     Bwop'*Qop;
            Dzwop,        -gam*Iz,    Czop;
            Qop'*Bwop,    Czop'       Aop'*Qop+Qop'*Aop];
else
    Dop = [-gam*Iz,       Dzwop,     Czop*Qop;
            Dzwop',       -gam*Iw,   Bwop';
            Qop'*Czop',   Bwop,      Qop'*Aop'+Aop*Qop];
end

[prog, De1op] = poslpivar(prog, Dop.dim, st.dd2, st.options2);
if st.override2~=1
    [prog, De2op] = poslpivar(prog, Dop.dim, st.dd3, st.options3);
    Deop = De1op + De2op;
else
    Deop = De1op;
end
prog = lpi_eq(prog,Deop+Dop,'symmetric');

H.Top=Top; H.Aop=Aop; H.Bwop=Bwop; H.Czop=Czop; H.Dzwop=Dzwop;
H.Iw=Iw;   H.Iz=Iz;   H.Rop=Rop;   H.Qop=Qop;   H.Dop=Dop;  H.Deop=Deop;
H.gam=gam; H.dual=dual; H.st=st;   H.PIE=PIE;
end

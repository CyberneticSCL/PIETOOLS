function A = pielr_lpi(lpi)                                                 % CC, 09/23/2026
% PIELR_LPI  Adapter registry: everything that is executive-specific about a
% low-rank certification, behind one struct of function handles.
%
% WHY AN ADAPTER LAYER.  The low-rank machinery (bm_setup, bm_lm2, bm_resid,
% restrict_solve, the rank ladder) sees only the SeDuMi triple and a rank
% profile -- nothing in it is aware of a spatial dimension or of which LPI is
% being certified.  What IS executive-specific is exactly three things:
%   (1) how the unsolved program is assembled;
%   (2) which operator expression must vanish, i.e. the executive's own
%       constraint lines repeated on the SUBSTITUTED operators;
%   (3) what to divide the residual by so the gate is scale-free and the
%       denominator cannot collapse.
% Putting those three behind one interface is what lets the same core certify
% stability and l2gain, in 1-D and 2-D, without a fork per case.  The package
% previously forked its whole core into tests_1d/ to reach 1-D, and two of the
% seven files have since drifted out of sync.
%
% (3) DESERVES ITS OWN NOTE, because it is the one place where a wrong choice
% silently changes what "certified" means.  The gate is
%       rel = max|coeff(Res)| / max|coeff(Nrm)|
% and Nrm must be an operator whose max coefficient is BOUNDED BELOW at every
% candidate, or rel reads 0/0 and improves for an artifactual reason.  This
% was measured in 1-D: on the Example 21 beam, dividing by max|Dop| gave
% max|Dop| = 2.4e-10 and manufactured a clean degree threshold that does not
% exist (tests_1d/opcheck.m, "Accept on relP, never on rel alone").  The 2-D
% checker opcheck_2d still divides by max|Qop| and has no such guard.  Each
% adapter therefore names its normaliser AND states why it cannot collapse.
%
% INTERFACE (fields of the returned struct)
%   .name      the lpi string, as lpiscript spells it
%   .dual      logical, whether this is the dual form
%   .has_obj   logical, whether the LPI carries an objective (gamma)
%   .build     [prog,H] = build(PIE,settings)   assembles, does NOT solve
%   .resid     S = resid(prog,H)   with prog.solinfo.RRx already set:
%                 S.Res  operator that must vanish
%                 S.Nrm  operator the residual is measured against
%                 S.why  one line: why max|Nrm| cannot collapse
%                 S.aux  struct of anything else worth reporting (e.g. gam)
%
% SUPPORTED.  'stability', 'stability-dual' in 1-D and 2-D.  'l2gain' and
% 'l2gain-dual' are declared here and dispatch to their own builders; the 2-D
% l2gain LPI is NOT the 1-D one re-dimensioned (PIETOOLS_Hinf_gain introduces
% a free lpivar operator Qop and couples it by Top'*Qop = Rop, while
% PIETOOLS_Hinf_gain_2D works with Pop directly), so they are separate
% transcriptions.

switch lower(char(lpi))
    case 'stability'
        A = stab_adapter(false);
    case 'stability-dual'
        A = stab_adapter(true);
    case 'l2gain'
        A = l2g_adapter(false);
    case 'l2gain-dual'
        A = l2g_adapter(true);
    case 'poincare'
        A = poincare_adapter();
    otherwise
        error('pielr_lpi:unknown', ...
            ['Unsupported LPI ''%s''.  This package currently certifies ' ...
             '''stability'', ''stability-dual'', ''l2gain'' and ''l2gain-dual''.\n' ...
             'The remaining lpiscript types (h2norm, well-posedness and the ' ...
             'four synthesis executives) need their own adapter: an assembly, ' ...
             'a residual expression and a non-collapsing normaliser.'],char(lpi));
end
end

% =========================================================================
function A = stab_adapter(dual)
if dual, A.name = 'stability-dual'; else, A.name = 'stability'; end
A.dual    = dual;
A.has_obj = false;
A.build   = @(PIE,st) stab_build(PIE,st,dual);
A.resid   = @(prog,H) stab_resid(prog,H);
end

function [prog,H] = stab_build(PIE,st,dual)
if pie_dim(PIE)==2
    [prog,H] = build_stab_2d_st2(PIE,stab_2d_settings(st),dual);
else
    [prog,H] = build_stab_1d(PIE,st,dual);
end
end

function s2d = stab_2d_settings(st)
% A 1-D settings struct carries the 2-D one under .settings_2d (lpisettings'
% own convention, and t1d_set's).  Accept either spelling so a caller can pass
% whichever it has.
if isfield(st,'settings_2d') && ~isempty(st.settings_2d)
    s2d = st.settings_2d;
else
    s2d = st;
end
end

function S = stab_resid(prog,H)
% Repeat the executive's own constraint lines on the SUBSTITUTED operators.
% Valid because the negativity operator is a LINEAR map of Pop and Pop is
% affine in the decision variables, so substituting then rebuilding equals
% rebuilding then substituting -- this tests the identity semantically
% instead of trusting the assembled equality rows.
Pop  = pielr_evalop(prog,H.Pop);
Deop = pielr_evalop(prog,H.Deop);
if isa(Pop,'opvar2d')
    % PIETOOLS_stability_2D lines 164-172 / _dual_2D lines 165-173
    if ~H.dual
        PTop  = Pop*H.Top;    APTop = H.Aop'*PTop;
        if H.epneg==0, Dop = APTop' + APTop;
        else,          Dop = APTop' + APTop + 2*H.epneg*(H.Top'*PTop);  end
    else
        TPop  = H.Top*Pop;    TPAop = TPop*H.Aop';
        if H.epneg==0, Dop = TPAop' + TPAop;
        else,          Dop = TPAop' + TPAop + 2*H.epneg*(TPop*H.Top');  end
    end
    Dop = clean_opvar(Dop,1e-12);      % same ztol the constraint was built at
else
    % PIETOOLS_PDEstability / _dual, one line each
    if ~H.dual
        Dop = H.Top'*Pop*H.Aop + H.Aop'*Pop*H.Top + H.epneg*H.Top'*Pop*H.Top;
    else
        Dop = H.Top*Pop*H.Aop' + H.Aop*Pop*H.Top' + H.epneg*H.Top*Pop*H.Top';
    end
end
S.Res = Dop + Deop;
S.Nrm = Pop;
S.why = ['Pop = poslpivar(..) + blkdiag(eppos*I,eppos2*I), so Pop >= eppos*I ' ...
         'by construction and max|Pop| is bounded below at every candidate.'];
% field assignment, NOT struct('Dop',Dop,...): with an opvar2d argument that
% call dispatches to opvar2d's own struct method and errors "Too many input
% arguments".  Same trap pielr_certify's norm_pie documents; it does not fire
% in 1-D, so it only shows up on the 2-D arm.
S.aux = struct();
S.aux.Dop = Dop;   S.aux.Pop = Pop;   S.aux.Deop = Deop;
end

% =========================================================================
function A = l2g_adapter(dual)
if dual, A.name = 'l2gain-dual'; else, A.name = 'l2gain'; end
A.dual    = dual;
A.has_obj = true;
A.build   = @(PIE,st) l2g_build(PIE,st,dual);
A.resid   = @(prog,H) l2g_resid(prog,H);
end

function [prog,H] = l2g_build(PIE,st,dual)
if pie_dim(PIE)==2
    error('pielr_lpi:l2gain2d', ...
        ['2-D l2gain is not yet adapted.  PIETOOLS_Hinf_gain_2D is a separate ' ...
         'formulation from the 1-D one (Pop directly, with Twop, rather than ' ...
         'a free lpivar Qop coupled by Top''*Qop = Rop), so it needs its own ' ...
         'transcription rather than a re-dimensioned copy of build_l2gain_1d.']);
end
[prog,H] = build_l2gain_1d(PIE,st,dual);
end

function S = l2g_resid(prog,H)
% PIETOOLS_Hinf_gain's own negativity block, rebuilt on substituted operators.
Qop  = pielr_evalop(prog,H.Qop);
Deop = pielr_evalop(prog,H.Deop);
gam  = pielr_evalop(prog,H.gam);
if ~isa(gam,'double'), gam = double(gam); end
if ~H.dual
    Dop = [-gam*H.Iw,        H.Dzwop',    H.Bwop'*Qop;
            H.Dzwop,        -gam*H.Iz,    H.Czop;
            Qop'*H.Bwop,     H.Czop'      H.Aop'*Qop+Qop'*H.Aop];
else
    Dop = [-gam*H.Iz,        H.Dzwop,    H.Czop*Qop;
            H.Dzwop',       -gam*H.Iw,   H.Bwop';
            Qop'*H.Czop',    H.Bwop,     Qop'*H.Aop'+H.Aop*Qop];
end
% BOTH equalities the builder imposes, not just the negativity one.  The
% coupling Top'*Qop = Rop is what ties the operator asserted PSD (Rop) to the
% one appearing in the Lyapunov inequality (Qop); without it in the gate a
% point can satisfy the negativity identity with the two unrelated.
Rop = pielr_evalop(prog,H.Rop);
if ~H.dual, Ceq = H.Top'*Qop - Rop;   % PIETOOLS_Hinf_gain
else,       Ceq = H.Top*Qop  - Rop;   % PIETOOLS_Hinf_gain_dual
end
S.Res = {Dop + Deop, Ceq};
% Normalisers.  Dop carries -gam*Iw and -gam*Iz on its diagonal so
% max|Dop| >= gam > 0; note there is NO eppos*I on Rop in this executive, so
% max|Rop| alone CAN collapse and must not be used.  For the coupling, the
% larger of the two sides is zero only when BOTH are, i.e. only at the trivial
% point -- which the negativity residual rejects independently here, because
% Dzwop and Czop enter Dop as standalone constant blocks.
nrm2 = max(pielr_maxop(H.Top'*Qop),pielr_maxop(Rop));
S.Nrm = {Dop, nrm2};              % a normaliser may be a plain scalar
S.why = ['negativity: max|Dop| >= gam > 0 from the -gam*I diagonal blocks. ' ...
         'coupling: max(|Top''Qop|,|Rop|), zero only at the trivial point, ' ...
         'which the negativity residual rejects on its own since Dzwop and ' ...
         'Czop are standalone constant blocks of Dop.'];
S.aux = struct();                          % see stab_resid: not struct(...)
S.aux.Dop = Dop;   S.aux.gam = gam;   S.aux.Qop = Qop;   S.aux.Deop = Deop;
end

% =========================================================================
function A = poincare_adapter()
A.name = 'poincare';
A.dual    = false;
A.has_obj = true;
A.build   = @(PIE,st) poincare_build(PIE,st);
A.resid   = @(prog,H) poincare_resid(prog,H);
end

function [prog,H] = poincare_build(PIE,st)
if pie_dim(PIE)==2
    error('pielr_lpi:poincare2d', ...
        ['2-D Poincare is not adapted.  lpi_ineq_2d is a separate routine ' ...
         'from lpi_ineq and would need its own transcription.']);
end
[prog,H] = build_poincare_1d(PIE,st);
end

function S = poincare_resid(prog,H)
% lpi_ineq's own constraint, rebuilt on the substituted operators:
%   Deop [+ De2op] - (gam*H1'H1 - H2'H2) = 0
Deop = pielr_evalop(prog,H.Deop);
gam  = pielr_evalop(prog,H.gam);
if ~isa(gam,'double'), gam = double(gam); end
Pop  = gam*(H.H1'*H.H1) - H.H2'*H.H2;
if isempty(H.De2op)
    Res = Deop - Pop;
else
    De2op = pielr_evalop(prog,H.De2op);
    Res   = Deop + De2op - Pop;
end
S.Res = Res;
% THE NORMALISER IS A CONSTANT OPERATOR, which is the cleanest case in the
% package: H2'*H2 carries no decision variables at all, so max|H2'H2| is a
% fixed positive number and cannot collapse for any candidate whatsoever --
% unlike max|Pop| (which depends on eppos and the solution) or max|Dop| (which
% is floored only by gamma).  Nothing here needs an argument about what keeps
% the denominator away from zero.
S.Nrm = H.H2'*H.H2;
S.why = ['H2''*H2 is a CONSTANT operator (no decision variables), so ' ...
         'max|H2''H2| is a fixed positive number independent of the candidate.'];
S.aux = struct();
S.aux.Dop = Pop;   S.aux.gam = gam;   S.aux.Deop = Deop;
% the reported quantity is the Poincare constant, sqrt(gam), against 1/pi
S.aux.poincare = sqrt(max(gam,0));
end

% =========================================================================
function d = pie_dim(PIE)
% PIE.dim is set by initialize() on a pie_struct; a raw struct of operators
% (which pielr_certify also accepts) is classified from the operator class.
if isstruct(PIE) && isfield(PIE,'dim') && ~isempty(PIE.dim)
    d = PIE.dim;   return
end
if isa(PIE,'pie_struct')
    d = PIE.dim;   return
end
if isa(PIE.T,'opvar2d') || isa(PIE.T,'dopvar2d'), d = 2; else, d = 1; end
end

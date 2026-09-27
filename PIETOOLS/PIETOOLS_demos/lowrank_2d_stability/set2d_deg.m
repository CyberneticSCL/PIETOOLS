function s2d = set2d_deg(Dup,dmult)
% set2d_deg(Dup,dmult) -- 2D settings with ONE basis knob, plus a MULTIPLIER knob.
%
% LF_deg is FIXED at settings_PIETOOLS_light_2D's values and the eq (negativity
% slack) degrees are LF_deg + Dup, which is exactly how the light settings file
% builds them (Dupx = Dupy = Dup2 = 3 there).  So Dup = 0,1,2,... sweeps the
% basis count of the negativity Gram at FIXED physics, FIXED state count and a
% FIXED Lyapunov block.  Dup = 3 reproduces settings_PIETOOLS_light_2D exactly
% (asserted against the stock struct when this was built).
%
% dmult (optional, default [] = leave light's value) sets LF_deg.d2{1,1} =
% dmult*ones(2,2).  At light's default, LF_deg.d2{1,1} = zeros(2,2), so the
% Lyapunov block's multiplier monomial Z2^oo = 1 is CONSTANT.  Raising dmult to
% >= 1 makes it non-constant, so the Lyapunov operator acquires a genuine
% non-constant multiplier while the physics and the negativity block are
% untouched.
%
% psatz is OFF everywhere (LF_use_psatz = 0, eq_use_psatz = [0;0]): poslpivar's
% psatz=1 REPLACES the plain term rather than adding to it, and enabling it
% made whole measured studies infeasible.  Keeping it off also keeps the block
% count at 2 (LF, eq) so per-block ranks map 1:1 onto operators.
%
% CC, 09/27/2026: eppos 1e-4/1e-6 -> 1e-2, matching pielr_settings.  The two
%   settings files in this package disagreed about one knob and the 2-D cells,
%   which use this one, kept the tiny-||b|| pathology the 1-D arm had already
%   shed.  See the comment at the assignment.  Invalidates banked 2-D numbers.
% CC, 09/27/2026: also see set2d_psatz4, which wraps this to add the four
%   linear face generators; it inherits this eppos.
if nargin<2, dmult = []; end
s2d = settings_PIETOOLS_light_2D();
if ~isempty(dmult)
    % the (1,1) cell is the multiplier direction in BOTH x and y; raising it is
    % what makes R22{1,1} non-constant
    s2d.LF_deg.d2{1,1} = dmult*ones(2,2);
end
dx = s2d.LF_deg.dx;   dy = s2d.LF_deg.dy;   d2 = s2d.LF_deg.d2;
s2d.eq_deg.dx = {Dup+dx{1};   Dup+dx{2};   Dup+dx{3}};
s2d.eq_deg.dy = {Dup+dy{1},   Dup+dy{2},   Dup+dy{3}};
s2d.eq_deg.d2 = cell(3,3);
for a = 1:3
    for b = 1:3
        s2d.eq_deg.d2{a,b} = Dup+d2{a,b};
    end
end
s2d.LF_use_psatz = 0;
s2d.eq_use_psatz = [0;0];
% eppos RAISED TO 1e-2 (CC, 09/27/2026).  b is proportional to eppos for the
% stability LPI, so the shipped 1e-4/1e-6 made ||b|| = 5.6e-06 on
% rd2d-deg3-f010 (m = 3456, 70 nonzeros in b).  At that scale the trivial
% point satisfies the equality rows to ~1e-9 and acceptance sits on the noise
% floor -- the same pathology diagnosed and fixed for 1-D, where
% pielr_settings pins eppos = eppos2 = 1e-2 and also writes 1e-2 into
% settings_2d.  set2d_deg is what the 2-D bench cells actually use, so the two
% settings files in this package disagreed about one knob and the 2-D arm kept
% the pathology the 1-D arm had shed.  Relative measures are invariant to this
% and absolute ones scale with it, so raising it lifts acceptance off the
% noise floor without flattering any ratio.
% INVALIDATES every banked 2-D number taken with this file: b, and hence every
% absolute residual, changes by ~1e4.
% s2d.eppos = [1e-4;1e-6;1e-6;1e-6];                                       % CC, 09/27/2026 (was)
s2d.eppos = 1e-2*ones(4,1);                                               % CC, 09/27/2026
s2d.epneg = 0;
s2d.use_sosineq = 0;
% MASTER-BRANCH PORTABILITY (measured, 09/19/2026): the shipped                % CC, 09/19/2026
% settings_PIETOOLS_light_2D reads use_sosineq FROM THE BASE WORKSPACE via     % CC, 09/19/2026
% evalin, defaulting to 1 -- and on that branch it never populates eq_opts /   % CC, 09/19/2026
% eq_deg_psatz / eq_opts_psatz.  Default them here so this function is         % CC, 09/19/2026
% deterministic on every PIETOOLS branch and whatever the caller's base        % CC, 09/19/2026
% workspace holds (a no-op where the settings file already set them).          % CC, 09/19/2026
if ~isfield(s2d,'eq_opts') || ~isfield(s2d.eq_opts,'exclude')                  % CC, 09/19/2026
    s2d.eq_opts = struct('psatz',0,'sep',zeros(1,6),'exclude',zeros(1,16));    % CC, 09/19/2026
end                                                                            % CC, 09/19/2026
if ~isfield(s2d,'eq_deg_psatz'),  s2d.eq_deg_psatz  = {s2d.eq_deg};  end       % CC, 09/19/2026
if ~isfield(s2d,'eq_opts_psatz'), s2d.eq_opts_psatz = {s2d.eq_opts}; end       % CC, 09/19/2026
end

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
s2d.eppos = [1e-4;1e-6;1e-6;1e-6];
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

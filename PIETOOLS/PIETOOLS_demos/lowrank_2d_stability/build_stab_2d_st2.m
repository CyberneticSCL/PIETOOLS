function [prog,H] = build_stab_2d_st2(PIE,s2d)
% build_stab_2d_st2 -- 2-D direct-form stability LPI, returning OPERATOR HANDLES.
%
% Same LPI as executives/2D/PIETOOLS_stability_2D.m, but
%   (a) takes the 2D settings STRUCT directly, so the 2D degrees can be varied
%       (the executive takes a stock settings NAME), and
%   (b) returns H, without which the solution cannot be checked at the operator
%       level.  The executive returns prog alone, so a 2D feasibility claim
%       could only rest on an equality residual -- worthless on these programs,
%       where norm(b) ~ 5e-6 means X = 0 already satisfies the equalities to
%       ~1e-9 relative while the operator residual is 1.
%
% The two poslpivar_2d calls request TWO outputs, never three.  Reason:
% poslpivar_2d assembles its Qmat return only when nargout>2 and does so by
% growing a dpvar inside a double loop over an object whose dvarname list is
% every decision variable of the block (q = N(N+1)/2 -- millions at large N).
% MEASURED: the third output costs 62% of build time and 75% of peak memory,
% scales n^2.96 against the problem's n^2.00, and the resulting program is
% bitwise identical without it.  H.Q is written by no code path here and read
% by none; every stock caller (all of executives/2D/*) also uses two outputs.
%
% The 2D executive builds the negativity operator as
%   Qop = (A'PT)' + A'PT [+ 2*epneg*T'PT]
% via PTop/APTop, with the poslpivar_2d slack Qeop.  opcheck_2d rebuilds Qop by
% exactly the lines below, so the two must not drift apart.
%
% INPUT
%   PIE  - PIE struct (2 spatial variables)
%   s2d  - 2D settings struct, e.g. settings_PIETOOLS_stripped_2D().  eppos and
%          epneg are defaulted exactly as the executive defaults them.
% OUTPUT
%   prog - unsolved lpiprogram
%   H    - .Top .Aop  PIE operators
%          .Pop       Lyapunov operator INCLUDING the eppos identity term
%          .Qop       negativity operator (executive's name; 1D analogue Dop)
%          .Deop      poslpivar_2d slack Qeop, so Qop + Deop == 0 is the equality
%          .Q         {Qmat} Gram handles from each poslpivar_2d call, in
%                     declaration order
%          .tags      a label per poslpivar_2d call, aligned with H.Q
%          .st .PIE   settings actually used, and the PIE

Top = PIE.T;    Aop = PIE.A;

% ---- eppos / epneg defaults, copied from PIETOOLS_stability_2D lines 81-95 ---
if ~isfield(s2d,'eppos')
    eppos = [1e-4; 1e-6; 1e-6; 1e-6];
else
    eppos = s2d.eppos;
    if numel(eppos)==1,  eppos = eppos*ones(4,1);  end %#ok<*ISCL> % executive's own spelling kept verbatim
end
if ~isfield(s2d,'epneg'),  epneg = 0;  else,  epneg = s2d.epneg;  end

% The inequality route declares no Qeop, so there is nothing for opcheck_2d to
% verify against; refuse rather than return an H that cannot be checked.
if isfield(s2d,'use_sosineq') && s2d.use_sosineq
    error('build_stab_2d_st2: use_sosineq=1 has no Qeop slack; opcheck_2d needs the equality route.');
end

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);

% ---- STEP 2: positive Lyapunov operator Pop -------------------------------
[prog,Pop] = poslpivar_2d(prog,Top.dim,s2d.LF_deg,s2d.LF_opts);              % 2 outputs: see header
H.Q = {};   H.tags = {'LF'};                                                  % Qmat not requested
up = s2d.LF_use_psatz;
dp = extract_psatz_deg_l(s2d.LF_deg_psatz,up);
op = extract_psatz_opts_l(s2d.LF_opts_psatz,up);
for j = 1:numel(up)
    if up(j)~=0
        [prog,P2op,QPj] = poslpivar_2d(prog,Top.dim,dp{j},op{j});
        Pop = Pop + P2op;
        H.Q{end+1} = QPj;   H.tags{end+1} = sprintf('LFpsatz%d',j);
    end
end
% strict positivity margin
if ~all(eppos==0)
    np_op = Pop.dim(:,1);
    Ip = blkdiag(eppos(1)*eye(np_op(1)),eppos(2)*eye(np_op(2)), ...
                 eppos(3)*eye(np_op(3)),eppos(4)*eye(np_op(4)));
    Pop = Pop + opvar2d(Ip,Pop.dim,PIE.dom,PIE.vars);
end

% ---- STEP 3: negativity operator Qop -------------------------------------
% opcheck_2d repeats these five lines on the SUBSTITUTED Pop.  Valid because
% Qop is a LINEAR map of Pop and Pop is affine in the decision variables, so
% substitution commutes with the operator algebra.
PTop  = Pop*Top;
APTop = Aop'*PTop;
if epneg==0
    Qop = APTop' + APTop;
else
    Qop = APTop' + APTop + 2*epneg*(Top'*PTop);
end
ztol = 1e-12;
Qop = clean_opvar(Qop,ztol);

% ---- STEP 4: negativity slack Qeop, and the equality ---------------------
% The executive's toggle/checkdeg_lpi_eq_2d loop is dead code (toggle = 0 at
% PIETOOLS_stability_2D line 194) and `prog = progQ` just commits this call.
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,ztol);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,eq_opts);                 % 2 outputs: see header
H.tags{end+1} = 'eq';                                                        % Qmat not requested
ue = s2d.eq_use_psatz;
de = extract_psatz_deg_l(s2d.eq_deg_psatz,ue);
oe = extract_psatz_opts_l(s2d.eq_opts_psatz,ue);
for j = 1:numel(ue)
    if ue(j)~=0
        oe{j}.exclude = oe{j}.exclude | eq_opts.exclude;
        oe{j}.sep     = oe{j}.sep     | eq_opts.sep;
        [prog,Qe2op,QDj] = poslpivar_2d(prog,Qop.dim,de{j},oe{j});
        Qeop = Qeop + Qe2op;
        H.Q{end+1} = QDj;   H.tags{end+1} = sprintf('eqpsatz%d',j);
    end
end
prog = lpi_eq_2d(prog,Qeop+Qop,'symmetric');

H.Top = Top;    H.Aop = Aop;    H.Pop = Pop;
H.Qop = Qop;    H.Deop = Qeop;
H.eppos = eppos;  H.epneg = epneg;
H.st = s2d;     H.PIE = PIE;
end

% ==== local copies of the executive's own local functions =================
% extract_psatz_deg / extract_psatz_opts are LOCAL to
% executives/2D/PIETOOLS_stability_2D.m (lines 258-309) and so are not on the
% path.  Kept local here too, verbatim, rather than relying on an external
% copy being addpath'd.  Lint pragmas aside, do not "improve" them: staying
% byte-comparable to the executive is the point.
%#ok<*AGROW>

function outcell = extract_psatz_deg_l(incell,use_psatz)
if all(use_psatz==0)
    outcell = {};
    return
end
if isa(incell,'struct')
    for j=1:length(use_psatz)
        outcell{j} = incell;
    end
elseif numel(incell)==1
    for j=1:length(use_psatz)
        outcell{j} = incell{1};
    end
elseif numel(incell)==length(use_psatz)
    outcell = incell;
elseif numel(incell)>=max(use_psatz)
    outcell = incell(use_psatz);
else
    error('For each element of ''use_psatz'', a ''deg'' field should be defined.')
end
end

function outcell = extract_psatz_opts_l(incell,use_psatz)
if all(use_psatz==0)
    outcell = {};
    return
end
if isa(incell,'struct')
    for j=1:length(use_psatz)
        outcell{j} = incell;
    end
elseif numel(incell)==1
    for j=1:length(use_psatz)
        outcell{j} = incell{1};
    end
elseif numel(incell)==length(use_psatz)
    outcell = incell;
elseif numel(incell)>=max(use_psatz)
    outcell = incell(use_psatz);
else
    error('For each element of ''use_psatz'', an ''opts'' field should be defined.')
end
end

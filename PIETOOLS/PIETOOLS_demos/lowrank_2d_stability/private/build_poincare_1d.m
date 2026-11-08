function [prog,H] = build_poincare_1d(PIE,st,~)                            % CC, 09/23/2026
% [PROG,H] = BUILD_POINCARE_1D(PIE,ST) assembles, without solving, the LPI
%
%       min gam   subject to   gam*(H1'*H1) - H2'*H2 >= 0,
%
% for the PIE of x_t = x_ss on [a,b] with x(a) = x(b) = 0 and output
% z = x_s, where H2 = T maps x_ss to x and H1 = C1 maps x_ss to x_s; sqrt(gam)
% bounds the Poincare constant. The inequality is posed by LPI_INEQ_SOP: the
% lift degree read from the support of the target, the weight degree raised
% by ST.poincare.dbump, and the psatz term ('auto') or none.
%
% INPUT
% - PIE:   'pie_struct' object of the 1-D PIE above, with 'copvar' (or
%          'opvar') operators T and C1;
% - st:    'struct' of settings; read are poincare.psatz (logical, default
%          true: the term with the weight (s-a)(b-s)), poincare.dbump
%          (nonnegative integer, default 0: added to the weight degree), and
%          gamfix (optional 'double', a fixed gam: a feasibility test).
% OUTPUT
% - prog:  'struct', the unsolved LPI program;
% - H:     'struct' with fields H1, H2 ('copvar'), gam ('dpvar', or 'double'
%          with ST.gamfix), Deop ('cdopvar', the positive operator Pop with
%          P - Pop = 0), De2op ([]), psatz, dw (the weight-degree raise
%          passed to LPI_INEQ_SOP), li ('struct', its info: degrees, terms),
%          st, PIE.
%
% CHANGE LOG
% CC, 09/23/2026: initial coding, the psatz branch of lpi_ineq on opvar, with
%   degrees from degbalance(P) plus the tier bump on every slot.
% CC, 09/24/2026: ST.gamfix.
% CC, 11/08/2026: rewritten on the container LPI layer.  opvar, poslpivar,
%   lpi_eq and degbalance left the path (f83486c8); the opvar body of
%   09/23/2026 and 09/24/2026 is deleted.  The degrees are now those of
%   LPI_INEQ_SOP (lift degree from the support, weight degree raised by
%   dbump), not degbalance(P) + dbump on every slot, and the psatz term sits
%   one weight degree lower: the program differs from the 09/23/2026 one.

% BEGIN CC, 11/08/2026: container build (replaces the opvar body, 09/23/2026
% and 09/24/2026)
% Tier options: psatz on/off and the weight-degree raise
pc = struct('psatz',1,'dbump',0);
if isfield(st,'poincare') && ~isempty(st.poincare)
    fn = fieldnames(st.poincare);
    for i = 1:numel(fn),  pc.(fn{i}) = st.poincare.(fn{i});  end
end
PIE  = exec_pie_sop(PIE);
H2op = PIE.T;       % (H2op*x_ss) = x
H1op = PIE.C1;      % (H1op*x_ss) = x_s
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
% gam: a decision variable with the objective min gam, or the fixed number
if isfield(st,'gamfix') && ~isempty(st.gamfix)
    gam = st.gamfix;
else
    [prog,gam] = lpidecvar(prog,'gam');
end
P = gam*(H1op'*H1op) - H2op'*H2op;
% P >= 0 as P - Pop = 0, Pop positive by construction.  'auto' is the base
% weight plus the (s-a)(b-s) term one weight degree lower, the interval
% Positivstellensatz form; 'none' the base weight alone, under which the
% LPI is infeasible for every degree and gam (proved 09/28/2026; BENCH.md).
dw = 1 + pc.dbump;                     % lpi_ineq_sop's default weight raise is 1
o  = struct('dw',dw);
if pc.psatz,  o.psatz = 'auto';  else,  o.psatz = 'none';  end
[prog,Pop,li] = lpi_ineq_sop(prog,P,o);
if ~(isfield(st,'gamfix') && ~isempty(st.gamfix))
    prog = lpisetobj(prog,gam);
end
% Handles for the operator gate (pielr_lpi poincare_resid)
H.H1 = H1op;    H.H2 = H2op;    H.gam = gam;
H.Deop = Pop;   H.De2op = [];
H.psatz = pc.psatz;     H.dw = dw;      H.li = li;
H.st = st;      H.PIE = PIE;
% END CC, 11/08/2026
end

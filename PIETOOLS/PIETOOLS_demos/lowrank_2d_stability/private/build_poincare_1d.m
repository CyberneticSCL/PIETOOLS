function [prog,H] = build_poincare_1d(PIE,st,~)                            % CC, 09/23/2026
% BUILD_POINCARE_1D  The Poincare-inequality LPI, returning OPERATOR HANDLES.
%
% Same LPI as PIETOOLS_demos/DEMO3_poincare_inequality.m: with
%   H2 = PIE.T   (H2*x_ss = x)      H1 = PIE.C1  (H1*x_ss = x_s)
% minimise gam subject to  gam*(H1'*H1) - H2'*H2 >= 0, and the Poincare
% constant is sqrt(gam).
%
% WHY THIS CASE IS WORTH ITS OWN ADAPTER.  It is the only benchmark here whose
% answer is ANALYTIC: on [0,1] with Dirichlet conditions at both ends the
% constant is exactly 1/pi = 0.3183099.  Every other case is scored against an
% interior-point solve of the same program, which is a reference and not a
% truth; this one is a truth.  It is therefore the right case to verify that a
% settings ladder actually buys accuracy rather than merely cost.
%
% THE DEGREES DO NOT COME FROM THE PRESET.  lpi_ineq sets them from
% degbalance(P), i.e. from the operator, so 'light' and 'veryheavy' give the
% identical program here.  The only knobs are the psatz flag and a bump on the
% balanced degree, which is why pielr_settings carries S.poincare separately
% from the stability degrees.
%
% Transcribed from lpi_ineq (the psatz branch), with the degree exposed:
%     [prog,Deop ] = poslpivar(prog,dim,d2,options2);
%     [prog,De2op] = poslpivar(prog,dim,d2,options3);   options3.psatz = 1
%     prog = lpi_eq(prog,Deop+De2op-P,'symmetric');
% The exclude logic is lpi_ineq's own, verbatim: it switches off the parameter
% blocks of P that are numerically zero, and dropping it inflates the program
% for nothing.
%
% INPUT   PIE  pie_struct with T and C1 (from convert of the Poincare PDE)
%         st   settings; only st.poincare.psatz and st.poincare.dbump are read
% OUTPUT  prog unsolved LPI program
%         H    .H1 .H2 .P .Deop .De2op .gam .psatz .st .PIE

pc = struct('psatz',1,'dbump',0);
if isfield(st,'poincare') && ~isempty(st.poincare)
    fn = fieldnames(st.poincare);
    for i=1:numel(fn), pc.(fn{i}) = st.poincare.(fn{i}); end
end

H2op = PIE.T;       % (H2op*x_ss) = x
H1op = PIE.C1;      % (H1op*x_ss) = x_s
dom  = H2op.I;      var1 = H2op.var1;

prog = lpiprogram(var1,dom);
% GAMFIX: numeric gamma, so the reference can be measured on the same kind of
% program the low-rank arm actually solves.  pielr_bisect_obj never optimises
% -- it pins gamma and tests feasibility -- so an OPTIMISING reference is not
% the same computation, and with score = rel/ipm_rel that difference sets the
% acceptance threshold.  See build_l2gain_1d for the residual measurements.
if isfield(st,'gamfix') && ~isempty(st.gamfix)
    gam = st.gamfix;                                                    % CC, 09/24/2026
else
    [prog,gam] = lpidecvar(prog,'gam');
end
P = gam*(H1op'*H1op) - H2op'*H2op;

dim = P.dim;
if any(dim(:,1)~=dim(:,2))
    error('build_poincare_1d:dim','P must be square to be sign definite.');
end
d2 = degbalance(P);
% the tier's bump, applied to every entry of degbalance's own answer
if pc.dbump > 0
    d2{1} = d2{1} + pc.dbump;
    d2{2} = d2{2} + pc.dbump;
    d2{3} = d2{3} + pc.dbump;
end

% lpi_ineq's exclude logic, verbatim: switch off the parameter blocks of P
% that are numerically zero
tol = 1e-14;
options2.exclude = [0,0,0,0];   options3.exclude = [0,0,0,0];
options2.sep = 0;               options3.sep = 0;
nm = {'P','R0','R1','R2'};
for i = 1:4
    if i==1, V = P.P; else, V = P.R.(nm{i}); end
    if isempty(V), continue, end
    if (isa(V,'double') && max(max(abs(V)))<tol) || ...
       (~isa(V,'double') && all(max(max(abs(V.C)))<tol))
        options2.exclude(i) = 1;   options3.exclude(i) = 1;
    end
end

if pc.psatz
    options3.psatz = 1;
    [prog, Deop ] = poslpivar(prog,dim,d2,options2);
    [prog, De2op] = poslpivar(prog,dim,d2,options3);
    prog = lpi_eq(prog,Deop+De2op-P,'symmetric');
else
    [prog, Deop ] = poslpivar(prog,dim,d2,options2);
    De2op = [];
    prog = lpi_eq(prog,Deop-P,'symmetric');
end
if ~(isfield(st,'gamfix') && ~isempty(st.gamfix))
    prog = lpisetobj(prog,gam);                                         % CC, 09/24/2026
end

H.H1 = H1op;   H.H2 = H2op;   H.gam = gam;
H.Deop = Deop; H.De2op = De2op;
H.psatz = pc.psatz;   H.d2 = d2;   H.st = st;   H.PIE = PIE;
end

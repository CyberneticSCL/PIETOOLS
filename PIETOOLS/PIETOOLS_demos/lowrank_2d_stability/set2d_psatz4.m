function s2d = set2d_psatz4(Dup,dmult)                                      % CC, 09/25/2026
% set2d_psatz4(Dup,dmult) -- set2d_deg plus the FOUR LINEAR PSATZ GENERATORS
% on the negativity constraint.  Same basis knobs, same physics, same
% Lyapunov block; the only change is eq_use_psatz.
%
% WHY THIS EXISTS.  set2d_deg turns the Positivstellensatz off everywhere,
% with the recorded reason that "poslpivar's psatz=1 REPLACES the plain term
% rather than adding to it, and enabling it made whole measured studies
% infeasible".  That reason is about psatz=1, the degree-4 PRODUCT of the four
% face inequalities, and it is sound: that multiplier vanishes on all four
% edges simultaneously, so it cannot correct a deficiency localised to one
% pair of edges.
%
% What has changed is that poslpivar_2d now carries psatz = 3,4,5,6, one
% NORMALISED LINEAR generator per face of the rectangle (commit 0cf1add1 on
% ndopvar).  Those are the useful form: separate multipliers, each degree 1.
%
% MEASURED by the parallel cuADMM session on this same plant (2-D
% reaction-diffusion, Dirichlet, searched Lyapunov operator, light settings,
% Mosek), reach as a fraction of lam* = 2*pi^2:
%     psatz off   certifies to ~0.25 lam* (SeDuMi); Mosek FAILS at every frac
%     linear4     certifies at rel_b ~1e-08 at 0.10, 0.50 AND 0.90 lam*
% so the generators EXTEND reach rather than being required to certify at all
% -- the earlier "reach 0.00 without them" figure was a Mosek artefact and has
% been corrected in poslpivar_2d's header.
%
% COST, same source: m 3456 -> 4492 (+30%), blocks [8 424] -> [8 424 424 424
% 424 424], nnz(At) 5.57e5 -> 3.66e6 (6.6x).  The nonzero count, not m, is
% what will be felt -- bm_setup's whitener forms Ssym*Ssym'.
%
% TWO TRAPS, both honoured below.
%  1. eq_opts_psatz must be CLONED, not built field by field.  The shipped
%     settings carry only as many entries as eq_use_psatz has, and assigning
%     e.g. eq_opts_psatz{3}.psatz = 5 creates a struct with no 'exclude'
%     field, whereupon poslpivar_2d dies with Unrecognized field name
%     "exclude".
%  2. eq_deg_psatz{j} must be the FULL eq_deg.  The usual convention of giving
%     psatz terms eq_deg-1 destroys the certificate outright here (measured
%     0.9999 -> 0.0000 reach), because the reduced blocks span no monomials
%     the base does not already have.  set2d_deg already defaults this way;
%     it is re-asserted here so the requirement is visible where it is used.
if nargin<1 || isempty(Dup),   Dup = 1;  end
if nargin<2,                   dmult = []; end
s2d = set2d_deg(Dup,dmult);

gens = [3;4;5;6];                       % the four faces of the rectangle
s2d.eq_use_psatz = gens;

% clone a known-good options struct per generator, then set only its psatz id
base_opts = s2d.eq_opts;
if isfield(s2d,'eq_opts_psatz') && ~isempty(s2d.eq_opts_psatz)
    base_opts = s2d.eq_opts_psatz{1};   % already a complete struct
end
s2d.eq_opts_psatz = cell(1,numel(gens));
s2d.eq_deg_psatz  = cell(1,numel(gens));
for j = 1:numel(gens)
    o = base_opts;                      % CLONE: keeps sep/exclude/etc intact
    o.psatz = gens(j);
    s2d.eq_opts_psatz{j} = o;
    s2d.eq_deg_psatz{j}  = s2d.eq_deg;  % FULL degree, not eq_deg-1
end
s2d.psatz_note = 'linear4: eq_use_psatz = [3;4;5;6], full eq_deg per term';
end

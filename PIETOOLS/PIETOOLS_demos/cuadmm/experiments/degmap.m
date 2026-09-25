% degmap.m -- does deg(Z'MZ) = 2*deg(Z) hold for the INTEGRAL cells?
%
% The whole degree algebra assumes the assembled operator's degree is twice the
% monomial-basis degree.  That is exact for a pointwise multiplier, but an
% integral cell involves int_a^s R(s,th) dth, and int_a^s th^k dth raises the
% s-degree by k+1 -- so the identity may simply not hold there, which would
% invalidate every balancing argument built on it (including my withdrawn one).
%
% IDENTIFICATION DESIGN: perturb ONE spec entry at a time from a baseline and
% record which assembled degree moves.  d = {d1, [dq dd dj], [dq dd dj]} where
% (verified at poslpivar.m:342-343) slot 1 is the QUADRATURE/integration
% variable and slot 2 the SURVIVING kernel variable.
% No solves -- pure structure.
cuadmm_path;
pvar s t
x = pde_var('state',1,s,[0,1]);
PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]));
DIM = PIE.T.dim;  VARS = PIE.vars;  DOM = PIE.dom;

base = {2,[1 2 3],[1 2 3]};
runs = { 'baseline'      , base
         'd1   +1'       , setc(base,1,3)
         'quad +1 (slot1)', setv(base,2,1,2)
         'surv +1 (slot2)', setv(base,2,2,3)
         'joint+1 (slot3)', setv(base,2,3,4) };

fprintf('DM variant|d1|quad|surv|joint|R0:degs|R1:degs,degth,tot|R2:degs,degth,tot\n');
for k = 1:size(runs,1)
    d = runs{k,2};
    prog = lpiprogram(VARS(:,1),VARS(:,2),DOM);
    o.psatz=0; o.exclude=[0 0 0 0]; o.sep=0;
    [~,P] = poslpivar(prog,DIM,d,o);
    [a0,~,t0] = pdeg(P.R.R0);
    [a1,b1,t1] = pdeg(P.R.R1);
    [a2,b2,t2] = pdeg(P.R.R2);
    fprintf('DM %s|%d|%d|%d|%d|%d|%d,%d,%d|%d,%d,%d\n', runs{k,1}, ...
        d{1}, d{2}(1), d{2}(2), d{2}(3), a0, a1,b1,t1, a2,b2,t2);
end

fprintf('DM --- pure scaling check: does deg scale as 2*d?\n');
for dv = 1:4
    d = {dv,[dv dv 2*dv],[dv dv 2*dv]};
    prog = lpiprogram(VARS(:,1),VARS(:,2),DOM);
    o.psatz=0; o.exclude=[0 0 0 0]; o.sep=0;
    [~,P] = poslpivar(prog,DIM,d,o);
    [a0,~,~]   = pdeg(P.R.R0);
    [a1,b1,t1] = pdeg(P.R.R1);
    fprintf('DM d=%d -> 2d=%d | R0 deg_s=%d | R1 deg_s=%d deg_th=%d total=%d\n', ...
            dv, 2*dv, a0, a1, b1, t1);
end

fprintf('DM --- with psatz (g = s(1-s), degree 2): how much does the multiplier add?\n');
for dv = 1:3
    d = {dv,[dv dv 2*dv],[dv dv 2*dv]};
    prog = lpiprogram(VARS(:,1),VARS(:,2),DOM);
    o.psatz=1; o.exclude=[0 0 0 0]; o.sep=0;
    [~,P] = poslpivar(prog,DIM,d,o);
    [a0,~,~]   = pdeg(P.R.R0);
    [a1,b1,t1] = pdeg(P.R.R1);
    fprintf('DM psatz d=%d | R0 deg_s=%d | R1 deg_s=%d deg_th=%d total=%d\n', ...
            dv, a0, a1, b1, t1);
end
fprintf('DMDONE\n');

function [ds,dt,tot] = pdeg(p)
% max degree in var1 (s), in var2 (th), and max TOTAL degree
ds=0; dt=0; tot=0;
if isempty(p), return; end
if isa(p,'double'), return; end
nm = p.varname; dm = full(p.degmat);
if isempty(dm), return; end
is = find(strcmp(nm,'s'),1);  it = find(~strcmp(nm,'s'));
if ~isempty(is), ds = max(dm(:,is)); end
if ~isempty(it), dt = max(max(dm(:,it))); end
tot = max(sum(dm,2));
end

function c = setc(c,i,v), c{i} = v; end
function c = setv(c,i,j,v), a = c{i}; a(j) = v; c{i} = a; c{3} = a; end

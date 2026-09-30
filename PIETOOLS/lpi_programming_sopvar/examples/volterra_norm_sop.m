function out = volterra_norm_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = VOLTERRA_NORM_SOP() is DEMO2 (PIETOOLS_demos/DEMO2_volterra_
% operator_norm.m) on the container path: an upper bound on the norm of the
% Volterra operator
%
%   (T x)(s) = int_0^s x(r) dr,   s in [0,1],   ||T|| = 2/pi = 0.63662,
%
% from the LPI  min gam  s.t.  gam - T'*T >= 0. The container form of
% DEMO2:68, lpi_ineq(prob, gam - Top'*Top), is 'gam - Tc'*Tc': a dpvar minus
% a copvar is gam*I - T'*T (Tier 1b). Container lpi_ineq is Tier 2, so the
% inequality is posed as its equality form with a positive slack: W from
% poscopvar at degree 1, the sum of a plain term and a product-psatz term
% (lpi_ineq with opts.psatz = 1 adds a psatz term to a plain one), and
% gam - T'*T - W = 0 by lpi_eq_cdopvar. The legacy DEMO2 LPI is solved
% alongside for comparison; the two positive cones need not coincide, so
% the bounds may differ, and both must lie above 2/pi.
%
% OUT: struct with the two bounds, sqrt(gam), and the MOSEK status and
% rel_b of the container solve.
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 0;  b = 1;
opvar Top;
Top.R.R1 = 1;   Top.I = [a,b];
% opvar2copvar needs the dummy variable named s_dum (opvar2sopvar).
Top.var1 = polynomial({'s'});   Top.var2 = polynomial({'s_dum'});
sopts = struct('solver','mosek','simplify',false);

% % % Legacy DEMO2 (lines 60-78), unchanged.
prob = lpiprogram(Top.vars,Top.I);
[prob,gam] = lpidecvar(prob,'gam');
opts.psatz = 1;
prob = lpi_ineq(prob,gam-Top'*Top,opts);
prob = lpisetobj(prob,gam);
prob = lpisolve(prob,sopts);
leg = sqrt(double(lpigetsol(prob,gam)));

% % % Container path.
Tc = opvar2copvar(Top);                     % L2[s] -> L2[s]
[sp,dm] = cx_space_list(Tc,'out');
prog = lpiprogram(Top.vars,Top.I);
[prog,gc] = lpidecvar(prog,'gam');
Km = gc - Tc'*Tc;                           % gam*I - T'*T, a cdopvar
[prog,W0] = poscopvar(prog,dm,sp,Top.I,1);                   % W0 >= 0
[prog,W1] = poscopvar(prog,dm,sp,Top.I,1,struct('psatz',1));  % >= 0 on [a,b]
W = W0 + W1;
prog = lpi_eq_cdopvar(prog,Km - W,'symmetric');
prog = lpisetobj(prog,gc);
prog = lpisolve(prog,sopts);
con = sqrt(double(lpigetsol(prog,gc)));
info = prog.solinfo.info;
rb = cx_resid(prog);

fprintf(['Volterra operator norm, exact 2/pi = %.6f: legacy DEMO2 bound %.6f, '...
         'container bound %.6f (MOSEK numerr %d, pinf %d, rel_b %.2e)\n'],...
        2/pi,leg,con,info.numerr,info.pinf,rb);
out = struct('exact',2/pi,'legacy',leg,'container',con,'numerr',info.numerr,...
    'pinf',info.pinf,'rel_b',rb);
if con < 2/pi - 1e-6 || leg < 2/pi - 1e-6 || info.numerr~=0 || rb > 1e-6
    error('volterra_norm_sop:check','A bound is below 2/pi or the container solve is not accepted.')
end
end

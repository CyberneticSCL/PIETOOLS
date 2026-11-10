% Test_polyopvar_times_v2.m -- verify polyopvar_times_v2 against `times`
% (the '.*' operator of the class 'polyopvar'), in the two regimes
% relevant to the "p1*g" product in LocalStability.m: p1 from SOS_DP.m
% (possibly involving a DMB/Def. 4-derived quartic term) times g from
% Weighted_Sobolev_Norm.m (a plain quadratic, non-Def. 4 term).
%
% PART 1 - plain quadratic p1 times plain quadratic g (no DMB involved):
%   agreement with `times`, plus a ground-truth check at a test state.
%
% PART 2 - p1 built with a DMB-derived (Def. 4) quartic term, as SOS_DP.m
%   produces for p1_deg>=2, times a plain quadratic g: agreement with
%   `times`, plus a ground-truth check at a test state.
%
% PART 3 - product of two functionals whose bases both carry a multiplier
%   component (nonzero Zop.R.R0, as in the d==1 basis of V3_Liediff.m):
%   `times` and polyopvar_times_v2 must not error, and the product must
%   equal the product of the values at a test state. This combination
%   does not arise in LocalStability.m: SOS_DP.m always zeros its R0
%   multiplier component, and Weighted_Sobolev_Norm.m uses no multiplier
%   basis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_polyopvar_times_v2
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ (CC), 10/15/2026: Initial coding;

clear; clear stateNameGenerator; close all; clc;
tol = 1e-8;

pvar s s_dum
dom = [0,1];
opdeg = 1;      % k=6 components, k^2=36 for the Kronecker product
Zmon = monomials([s,s_dum],0:opdeg);
Zop  = opvar();
Zop.var1 = s;   Zop.var2 = s_dum;   Zop.I = dom;
Zop.R.R0 = [0*Zmon;0*Zmon];
Zop.R.R1 = [Zmon;0*Zmon];
Zop.R.R2 = [0*Zmon;Zmon];

Z   = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop)));
x   = polyopvar('x1',1,s,dom);
Zx  = Z*x;
Zs2 = DMB(Zx,Zx);
n1 = size(Zx,1);
n2 = size(Zs2,1);
rng(21);

x_tst = rand_poly([1,1],s,3);

%% Part 1: plain quadratic p1, plain quadratic g.
p1a = innerprod_v2(Zx,Zx,rand(n1,n1));
ga  = innerprod_v2(Zx,Zx,rand(n1,n1));
pg_orig_a = times(p1a,ga);
pg_v2_a   = polyopvar_times_v2(p1a,ga);
pass1 = ops_match(pg_orig_a,pg_v2_a) && matches_ground_truth(pg_v2_a,p1a,ga,x_tst,tol);
fprintf(" --- Part 1 (plain quadratic x quadratic): %s ---\n",pass_str(pass1));

%% Part 2: DMB-derived quartic p1 (degrees 2,3,4), plain quadratic g.
Pcell11 = rand(n1,n1); Pcell11 = Pcell11+Pcell11';
Pcell12 = rand(n1,n2);
Pcell22 = rand(n2,n2); Pcell22 = Pcell22+Pcell22';
p1b = innerprod_v2(Zx,Zx,Pcell11) + 2*innerprod_v2(Zx,Zs2,Pcell12) + innerprod_v2(Zs2,Zs2,Pcell22);
gb  = innerprod_v2(Zx,Zx,rand(n1,n1));
pg_orig_b = times(p1b,gb);
pg_v2_b   = polyopvar_times_v2(p1b,gb);
pass2 = ops_match(pg_orig_b,pg_v2_b) && matches_ground_truth(pg_v2_b,p1b,gb,x_tst,tol);
fprintf(" --- Part 2 (DMB/Def.4 quartic x quadratic): %s ---\n",pass_str(pass2));

%% Part 3: multiplier terms in both factors.
opdeg0 = 1;
Zmon_m  = monomials([s,s_dum],0:opdeg0);
Zmon0_m = monomials(s,0:opdeg0);
Zop_m = opvar();
Zop_m.var1 = s;   Zop_m.var2 = s_dum;   Zop_m.I = dom;
Zop_m.R.R0 = [Zmon0_m;0*Zmon_m;0*Zmon_m];
Zop_m.R.R1 = [0*Zmon0_m;Zmon_m;0*Zmon_m];
Zop_m.R.R2 = [0*Zmon0_m;0*Zmon_m;Zmon_m];
Zx_m = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop_m)))*x;
nm = size(Zx_m,1);

Vm = innerprod_v2(Zx_m,Zx_m,rand(nm,nm));   % carries the special all-zero row
gm = innerprod_v2(Zx_m,[],rand(nm,nm));     % carries it too

try
    pg_orig_m = times(Vm,gm);
    pass3_crash = true;   % `times` of the class is expected NOT to crash here
catch
    pass3_crash = false;
end
try
    pg_v2_m = polyopvar_times_v2(Vm,gm);
    pass3_nocrash = true;
catch
    pass3_nocrash = false;
end
% The product must be the product of the values at x_tst 
pass3_val = false;
if pass3_crash && pass3_nocrash
    ref_m = double(eval_functional(Vm,x_tst))*double(eval_functional(gm,x_tst));
    val_m = double(eval_functional(pg_v2_m,x_tst));
    pass3_val = abs(val_m-ref_m)<tol*max(1,abs(ref_m));
end
fprintf(" --- Part 3 (multiplier terms in both factors): ---\n");
fprintf("       times of the class does not crash: %s\n", pass_str(pass3_crash));
fprintf("       polyopvar_times_v2 does not crash: %s\n", pass_str(pass3_nocrash));
fprintf("       the product equals the product of the values at x_tst: %s\n", pass_str(pass3_val));

all_pass = pass1 && pass2 && pass3_crash && pass3_nocrash && pass3_val;
fprintf("\n=== Test_polyopvar_times_v2: %s ===\n", pass_str(all_pass));


function tf = ops_match(C1,C2)
tf = polyopvar_equal(C1,C2,1e-12);
end

function tf = matches_ground_truth(C,p1,g,x_tst,tol)
p1_val = eval_functional(p1,x_tst);
g_val = eval_functional(g,x_tst);
ref_val = double(p1_val)*double(g_val);

C_val = eval_functional(C,x_tst);
tf = abs(double(C_val)-ref_val)<tol;
end

function Y = eval_functional(V,xval)                                        % DJ, 10/14/2026
% Y = EVAL_FUNCTIONAL(V,XVAL) evaluates the scalar functional V (polyopvar of
% 'intop' coefficients in a scalar state depending on s) at the polynomial
% state xval: the input of the operator of a monomial with d factors is
% kron(x(s_1),...,x(s_d)) (apply_op of 'intop'), summed over the monomials.
Zf = full(V.Z);
Y = 0;
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    d = nnz(Zf(k,:));
    v = 1;
    for j=1:d
        xj = xval;
        if d>1,     xj = subs(xval,pvar('s'),pvar(sprintf('s_%d',j)));     end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

function str = pass_str(tf)
if tf
    str = 'PASS';
else
    str = 'FAIL';
end
end

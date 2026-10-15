% Test_quad2lin_v2.m -- verify quad2lin_v2/quad2lin_term_v2 against
% `piesos_quadvar` and against direct evaluation, in the two regimes that
% arise in LocalStability.m:
%
% PART 1 - Regression check.
%   On the standard (non-Def. 4) basis Zx_std = dmonomials(x,(1:2)'),
%   quad2lin_v2(Pmat,Zop,Zx_std) must reproduce the functional V declared
%   by piesos_quadvar exactly (compared using polyopvar_equal).
%
% PART 2 - Def. 4 (vector-valued TDP) check.
%   Builds the same Zx and Zs2 = DMB(Zx,Zx) bases used by SOS_DP.m and
%   V3_Liediff.m to construct quartic functionals, and checks
%   innerprod_v2(.,.,.) [i.e. quad2lin_v2/quad2lin_term_v2] against a
%   direct single-integral evaluation
%       <Z(x),P*Z(x)>_{L2} = int_a^b Z(x)(s)'*P*Z(x)(s) ds,
%   with Z(x)(s) = kron(Zx(x)(s),Zx(x)(s)) for the degree-2 (quartic
%   functional) case, for a concrete test polynomial x(s). Both the
%   symmetric (deg2,deg2) pairing and the asymmetric (deg1,deg2) pairing
%   (in both orders) used by V3_Liediff/SOS_DP for d<=2 are checked.
%   Degrees above quartic (d>2) are not tested.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_quad2lin_v2
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
all_pass = true;

%% PART 1: regression check against piesos_quadvar (non-Def.4 basis).

pvar s1 s1_dum s1_dum_1 s1_dum_2 s1_dum_3 s1_dum_4
dom = [0,1];

vartab = polyopvar('x1',1,s1,dom);
Zx_std = dmonomials(vartab,(1:2)');

% V declared by piesos_quadvar, as the sum over the blocks of its Gram matrix
prog = piesos_program(vartab);
qopts.pos = true;       qopts.degs = 1;
[prog,Vx_orig,Pmat_std,Zop_std] = piesos_quadvar(prog,Zx_std,[],qopts);     % DJ, 10/14/2026

Vx_v2   = quad2lin_v2(Pmat_std,Zop_std,Zx_std);

pass1 = polyopvar_equal(Vx_orig,Vx_v2,1e-12);
fprintf(" --- Part 1 (regression vs. quad2lin, non-Def.4 basis): %s ---\n", pass_str(pass1));
all_pass = all_pass && pass1;


%% PART 2: Def.4 (vector-valued TDP) basis, as built by SOS_DP.m/V3_Liediff.m.

pvar s s_dum
% The Kronecker product of the vectors, as the pointwise product of the rows of the
% operators, has k^2 x k^2 terms for the weight in the quartic case: for opdeg 1,
% k=6 and the weight is 36 x 36
opdeg = 1;
Zmon  = monomials([s,s_dum],0:opdeg);
Zop   = opvar();
Zop.var1 = s;   Zop.var2 = s_dum;   Zop.I = dom;
Zop.R.R0 = [0*Zmon;0*Zmon];
Zop.R.R1 = [Zmon;0*Zmon];
Zop.R.R2 = [0*Zmon;Zmon];

Z  = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop)));
x  = polyopvar('x1',1,s,dom);
Zx = Z*x;
Zs2 = DMB(Zx,Zx);  % degree-2 Def.4 TDP basis: Zs2(x)(s) = kron(Zx(x)(s),Zx(x)(s))

n1 = size(Zx,1);
n2 = size(Zs2,1);

% Reference basis values, evaluated directly on a concrete test polynomial.
Zdop  = Zop;            % the PI operator (opvar), independent of the new classes
x_tst = rand_poly([1,1],s,3);
z1 = apply_opvar(Zdop,x_tst);   % n1 x 1 polynomial in s
z2 = kron(z1,z1);               % n2 x 1 polynomial in s

rng(0);

% --- Symmetric (deg2,deg2): the quartic-functional case used throughout
%     SOS_DP.m/V3_Liediff.m. ---
P22 = rand(n2,n2);
V22 = innerprod_v2(Zs2,Zs2,P22);
fval_ref22 = double(int(z2.'*P22*z2,s,dom(1),dom(2)));
fval_alt22 = double(eval_functional(V22,x_tst));
err22 = abs(fval_ref22-fval_alt22);
pass22 = err22<tol;
fprintf(" --- Part 2a (deg2,deg2) quartic functional: err=%.3g, %s ---\n",err22,pass_str(pass22));

% --- Asymmetric (deg1,deg2): cubic-functional pairing, left orientation. ---
P12 = rand(n1,n2);
V12 = innerprod_v2(Zx,Zs2,P12);
fval_ref12 = double(int(z1.'*P12*z2,s,dom(1),dom(2)));
fval_alt12 = double(eval_functional(V12,x_tst));
err12 = abs(fval_ref12-fval_alt12);
pass12 = err12<tol;
fprintf(" --- Part 2b (deg1,deg2) cubic functional: err=%.3g, %s ---\n",err12,pass_str(pass12));

% --- Asymmetric (deg2,deg1): cubic-functional pairing, right orientation. ---
P21 = rand(n2,n1);
V21 = innerprod_v2(Zs2,Zx,P21);
fval_ref21 = double(int(z2.'*P21*z1,s,dom(1),dom(2)));
fval_alt21 = double(eval_functional(V21,x_tst));
err21 = abs(fval_ref21-fval_alt21);
pass21 = err21<tol;
fprintf(" --- Part 2c (deg2,deg1) cubic functional: err=%.3g, %s ---\n",err21,pass_str(pass21));

all_pass = all_pass && pass22 && pass12 && pass21;

fprintf("\n=== Test_quad2lin_v2: %s ===\n", pass_str(all_pass));


function str = pass_str(tf)
if tf
    str = 'PASS';
else
    str = 'FAIL';
end
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
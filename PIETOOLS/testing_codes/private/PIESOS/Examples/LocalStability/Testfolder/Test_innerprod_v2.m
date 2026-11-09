% Test_innerprod_v2.m -- verify innerprod_v2 against `innerprod` on the
% "Z*x"-style bases used throughout the code (Zx, Tx, Rx, f, ...), all
% built by applying a tensor-PI operator to a 'polyopvar' state. A raw
% dmonomials-style multi-term identity basis is not a supported input to
% either innerprod or innerprod_v2.
%
% The numerical correctness for Def. 4 bases (vector-valued TDP, e.g.
% Zs2=DMB(Zx,Zx)) of the underlying quad2lin_v2/quad2lin_term_v2 pipeline
% is covered in Test_quad2lin_v2.m, via direct calls to innerprod_v2 on
% that basis; it is not repeated here.
%
% PART 1 - regression vs. innerprod, on the standard (non-Def. 4) basis:
%   1a: symmetric call, innerprod_v2(Zx,[],P);
%   1b: asymmetric call with explicit Z2, innerprod_v2(Zx,Zx,P) (path of
%       quad2lin_v2 with explicit ZopR, ZxR);
%   1c: mixed call with one raw ('double') and one 'polyopvar' argument,
%       innerprod_v2(const,Zx,P), as used e.g. for
%       innerprod(Zx,f,Pcell{1}).
%
% PART 2 - innerprod_v2-only convenience:
%   omitting P (or passing P=[]) must be equivalent to passing the
%   identity weight P=1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_innerprod_v2
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
all_pass = true;

pvar s s_dum
dom = [0,1];
opdeg = 2;
Zmon = monomials([s,s_dum],0:opdeg);
Zop  = opvar();
Zop.var1 = s;   Zop.var2 = s_dum;   Zop.I = dom;
Zop.R.R0 = [0*Zmon;0*Zmon];
Zop.R.R1 = [Zmon;0*Zmon];
Zop.R.R2 = [0*Zmon;Zmon];

Z  = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop)));
x  = polyopvar('x1',1,s,dom);
Zx = Z*x;   % genuine "Z*x" basis: type==[false,true], as used throughout.
n  = size(Zx,1);
rng(12);

%% Part 1a: symmetric call.
Pmat = rand(n,n);
V_orig = innerprod(Zx,Pmat,[]);
V_v2   = innerprod_v2(Zx,[],Pmat);
pass1a = ops_match(V_orig,V_v2);
fprintf(" --- Part 1a (symmetric) innerprod vs innerprod_v2: %s ---\n",pass_str(pass1a));

%% Part 1b: asymmetric call, explicit Z2 (5-arg quad2lin path).
Pmat2 = rand(n,n);
V_orig2 = innerprod(Zx,Pmat2,Zx);
V_v2_2  = innerprod_v2(Zx,Zx,Pmat2);
pass1b = ops_match(V_orig2,V_v2_2);
fprintf(" --- Part 1b (asymmetric, explicit Z2) innerprod vs innerprod_v2: %s ---\n",pass_str(pass1b));

%% Part 1c: mixed call, one raw 'double' argument.
Pmix = rand(1,n);
V_v2_3  = innerprod_v2(3,Zx,Pmix);
% innerprod of the class does not integrate for a raw double (V=(u1'*P)*u2 is the function
% 3*Pmix*Zx(s)), unlike innerprod_v2, the L2 inner product with the constant function 3:
% compare with the integral of this function at a test state, using the PI operator (opvar)
x_tst = 3*s^2-2*s+1;
z1    = apply_opvar(Zop,x_tst);
val_ref = double(int(3*Pmix*z1,s,dom(1),dom(2)));
val_v2  = double(eval_functional(V_v2_3,x_tst));
pass1c = abs(val_ref-val_v2)<1e-8*max(1,abs(val_ref));
fprintf(" --- Part 1c (mixed double/polyopvar) innerprod vs innerprod_v2: %s ---\n",pass_str(pass1c));

%% Part 2: empty/omitted P == identity weight, in innerprod_v2.
V_noP = innerprod_v2(Zx,[],[]);
V_P1  = innerprod_v2(Zx,[],1);
pass2 = ops_match(V_noP,V_P1);
fprintf(" --- Part 2 (empty P == identity weight, innerprod_v2 only): %s ---\n",pass_str(pass2));

all_pass = pass1a && pass1b && pass1c && pass2;
fprintf("\n=== Test_innerprod_v2: %s ===\n", pass_str(all_pass));


function tf = ops_match(V1,V2)
tf = polyopvar_equal(V1,V2,1e-12);
end

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
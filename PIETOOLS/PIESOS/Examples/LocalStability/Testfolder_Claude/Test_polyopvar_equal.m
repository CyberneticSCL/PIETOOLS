function Test_polyopvar_equal()
% TEST_POLYOPVAR_EQUAL verifies polyopvar_equal(F1,F2,tol), which tests
% whether two scalar 'polyopvar' objects (functionals) are equal, after
% canonicalization with clean_poly.
%
% Ground truth is by construction: F1 and F2 are built in different ways
% that are equal as functionals (different order of the sum, a product
% written in another order of the factors, a term added and subtracted),
% or differ by a known amount: F2=F1+delta*G for a functional G whose
% maximal kernel coefficient is known by hand, so that max_error must
% equal delta*max|G|.
%
% No SDPs are solved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_polyopvar_equal
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

clear stateNameGenerator
pvar s t
dom = [0,1];
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+x^2; subs(x,s,dom(1))==0; subs(x,s,dom(2))==0];
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;    xs = PIE.x;     f = PIE.f;

A = innerprod(xs,xs);               % delta term
B = innerprod(xs,Top*xs);           % kernel (s,theta) term
C = innerprod(Top*xs,Top*xs);       % 2 free variables

% (1) Equal functionals, built differently: order of the sum, a term added
% and subtracted, scalar factors
tol = 1e-12;
[ok,err] = polyopvar_equal(A+B+C,C+B+A,tol);
assert(ok && err==0,'(1a) sums in different order must be equal: error %g',err);
[ok,err] = polyopvar_equal(A+B,(A+B+C)-C,tol);
assert(ok && err<=tol,'(1b) term added and subtracted must be equal: error %g',err);
[ok,err] = polyopvar_equal(2*B,B+B,tol);
assert(ok && err<=tol,'(1c) 2B must equal B+B: error %g',err);
% symmetric product written with the factors in the other order: <x,Tx> and <Tx,x>
[ok,err] = polyopvar_equal(innerprod(xs,Top*xs),innerprod(Top*xs,xs),tol);
assert(ok && err<=tol,'(1d) <x,Tx> must equal <Tx,x> for the symmetric part: error %g',err);
% nonlinear term, the PIE f contains x o x: f minus itself
[ok,err] = polyopvar_equal(innerprod(xs,f),innerprod(xs,f),tol);
assert(ok && err==0,'(1e) a nonlinear functional must equal itself: error %g',err);

% (2) Differing functionals: F2 = F1 + delta*C, C = int (Tx)^2 has the
% kernel of Top*Top (known: maximal coefficient from the kernel of C itself)
delta = 1e-3;
Kc = kernels(C.C.ops{1,1});
cmax = 0;
for k=1:numel(Kc)
    cmax = max([cmax; abs(full(Kc{k}.coefficient(:)))]);
end
assert(cmax>0,'The test functional C must have a nonzero kernel.');
[ok,err] = polyopvar_equal(A+B,A+B+delta*C,tol);
assert(~ok,'(2a) differing functionals must not be equal.');
assert(abs(err-delta*cmax)<1e-12,'(2b) error %g, expected %g.',err,delta*cmax);
% within tolerance
[ok,~] = polyopvar_equal(A+B,A+B+delta*C,2*delta*cmax);
assert(ok,'(2c) functionals differing less than tol must be equal.');
% a difference only in the delta (multiplier) term
[ok,err] = polyopvar_equal(A+B,2*A+B,tol);
assert(~ok && abs(err-1)<1e-12,'(2d) delta term difference: error %g, expected 1.',err);
% a difference only in a degree-3 monomial
[ok,~] = polyopvar_equal(innerprod(xs,f),innerprod(xs,f)+delta*innerprod(xs,f),delta/2);
assert(~ok,'(2e) a difference in a nonlinear term must be detected.');

% (3) Inputs
threw = false;
try,    polyopvar_equal(A,1);   catch,  threw = true;   end
assert(threw,'(3a) an error is expected for a non-polyopvar input.');
threw = false;
try,    polyopvar_equal(A,B,-1);    catch,  threw = true;   end
assert(threw,'(3b) an error is expected for a negative tolerance.');
[ok,err] = polyopvar_equal(A,[A;B],tol);
assert(~ok && isinf(err),'(3c) different sizes must give logval=false, error=Inf.');

disp('Test_polyopvar_equal: all tests passed.')
end

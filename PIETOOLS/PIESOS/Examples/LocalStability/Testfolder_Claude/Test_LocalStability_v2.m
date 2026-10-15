function Test_LocalStability_v2()
% TEST_LOCALSTABILITY_V2 verifies LocalStability_v2(PDE,r,alpha,eppos,
% lambda,dist_degs,mon_degs,psatz,active_conds,sos_degs,M) on PDEs for
% which the outcome is known analytically, using convert_PIETOOLS_PDE_
% nonlinear, piesos_ineq and piesos_sosvar.
%
% LocalStability_v2 takes V = 1/2 <x,-T x> = 1/2 ||u_s||^2 for the
% Dirichlet PDEs below (x=u_ss, u=T*x), and minimizes C such that
%   (sos1)  V - eppos*B(x) >= 0,
%   (sos2)  C*B(x) - V - p1*g >= 0,           g = r^2 - B(x),
%   (sos3)  -dV - 2*lambda*V - p2*g >= 0,
% for B(x)=sum_j alpha_j*||d^j u/ds^j||^2. With alpha=[0,1,0],
% B=||u_s||^2=2V for these PDEs, so
%   (sos2) holds iff C>=1/2, and the optimal value is C=1/2,
%          M=sqrt(C/eppos);
%   (sos1) holds iff eppos<=1/2 (here 1/4);
%   (sos3) is -dV = ||u_ss||^2 - a*||u_s||^2 >= 0 for u_t=u_ss+a*u, true
%          iff a<=pi^2 (the multipliers p1, p2 can only relax the
%          conditions); for u_t=u_ss-u^3, -dV = ||x||^2 + 3||u*u_s||^2 >= 0
%          holds for all u.
% The checks:
% (1) heat, a=2, all conditions, no multipliers: C=1/2, M=sqrt(2), and
%     V_sol, dV_sol evaluated at the state of u=s(1-s) (x=-2) equal
%     1/2*int(1-2s)^2=1/6 and -||u_ss||^2+a||u_s||^2=-4+2/3, computed by
%     hand;
% (2) the same with the multipliers p1, p2, and with psatz order 1 in all
%     conditions: C=1/2 as well;
% (3) heat with a=12>pi^2, condition 3 only: infeasible (res empty);
% (4) cubic u_t=u_ss-u^3, all conditions with selected degrees, and the
%     call with the optional inputs omitted: C=1/2 (nonlinear PIE, `times`
%     and `innerprod` of nonlinear terms);
% (5) eppos and M as arrays (weights of the norms of the lower and upper
%     bound, like alpha), for the heat equation with alpha=[0,1,0] and
%     B=2V: with M=[0,2,0] (or the scalar M=2) V<=C*M*B holds iff C>=1/4;
%     with eppos=[0.1,0,0], an L2 lower bound, V-0.1*||T*x||^2>=0 holds,
%     and C=1/2 for M=[0,1,0]; M_sol=sqrt(C*M/eppos) for scalar eppos, M
%     and sqrt(C) otherwise;
% (6) errors for invalid psatz, and for eppos or M of another size than
%     alpha.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_LocalStability_v2
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

pvar s t
dom = [0,1];
ctol = 1e-4;
alpha = [0,1,0];
eppos = 0.25;

%% (1) Heat equation, a=2, all conditions
a = 2;
clear stateNameGenerator
x = pde_var(s,dom);
PDE_heat = [diff(x,t)==diff(x,s,2)+a*x; subs(x,s,dom(1))==0; subs(x,s,dom(2))==0];

[res,~] = run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[0,0,0],[1,1,1],[NaN,NaN,4]);
assert(~isempty(res),'(1a) the heat equation must be proven stable.');
check_C(res,1/2,eppos,ctol,'(1b)');
% V and dV at the state of u=s(1-s): x=u_ss=-2
Vval = double(eval_V(res{3},polynomial(-2)));
dVval = double(eval_V(res{4},polynomial(-2)));
assert(abs(Vval-1/6)<1e-8,'(1c) V(x)=%.10g, expected 1/6.',Vval);
assert(abs(dVval-(-4+a*(1/3)))<1e-8,'(1d) dV(x)=%.10g, expected %.10g.',dVval,-4+a/3);
fprintf('Passed (1): heat, a=2: C=%.6g, V=%.6g, dV=%.6g\n',res{1},Vval,dVval);

%% (2) With multipliers p1, p2, and with psatz multipliers
[res,~] = run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,1,1],[3,1,1],[0,0,0],[1,1,1],[NaN,NaN,4]);
assert(~isempty(res),'(2a) p1, p2 only relax the conditions: the heat equation must be proven stable.');
check_C(res,1/2,eppos,ctol,'(2b)');
[res,~] = run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[1,1,1],[1,1,1],[NaN,NaN,4]);
assert(~isempty(res),'(2c) psatz multipliers only relax the conditions: the heat equation must be proven stable.');
check_C(res,1/2,eppos,ctol,'(2d)');
fprintf('Passed (2): multipliers p1, p2 and psatz\n');

%% (3) Unstable heat equation, a=12>pi^2
clear stateNameGenerator
x = pde_var(s,dom);
PDE_unst = [diff(x,t)==diff(x,s,2)+12*x; subs(x,s,dom(1))==0; subs(x,s,dom(2))==0];
[res,~] = run_quiet(PDE_unst,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[0,0,0],[0,0,1],[NaN,NaN,4]);
assert(isempty(res),'(3a) the unstable heat equation (a=12>pi^2) cannot be proven stable: sos3 is infeasible.');
% the same PDE with a<pi^2 is feasible for these settings
[res,~] = run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[0,0,0],[0,0,1],[NaN,NaN,4]);
assert(~isempty(res),'(3b) the stable heat equation (a=2) must satisfy sos3 alone.');
fprintf('Passed (3): unstable heat equation is infeasible\n');

%% (4) Nonlinear PDE u_t = u_ss - u^3
clear stateNameGenerator
x = pde_var(s,dom);
PDE_cub = [diff(x,t)==diff(x,s,2)-x^3; subs(x,s,dom(1))==0; subs(x,s,dom(2))==0];
[res,~] = run_quiet(PDE_cub,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[0,0,0],[1,1,1],[NaN,NaN,4]);
assert(~isempty(res),'(4a) u_t=u_ss-u^3 must be proven stable.');
check_C(res,1/2,eppos,ctol,'(4b)');
% optional inputs omitted: no multipliers, all conditions, degrees selected
% from the kernels of the functionals
[res,~] = run_quiet(PDE_cub,1.0,alpha,eppos,0,[1,0,1],[3,0,0]);
assert(~isempty(res),'(4c) the call with the optional inputs omitted must prove u_t=u_ss-u^3 stable.');
check_C(res,1/2,eppos,ctol,'(4d)');
fprintf('Passed (4): cubic PDE: C=%.6g\n',res{1});

%% (5) eppos and M as arrays
[res,~] = run_quiet(PDE_heat,1.0,alpha,[0,0.25,0],0,[1,0,1],[3,0,0],[0,0,0],[1,1,0],[NaN,NaN,NaN],[0,2,0]);
assert(~isempty(res),'(5a) the heat equation must satisfy sos1 and sos2 for eppos=[0,1/4,0], M=[0,2,0].');
assert(abs(res{1}-1/4)<ctol,'(5a) C=%.8g, expected 1/4.',res{1});
assert(abs(res{2}-1/2)<10*ctol,'(5a) M=%.8g, expected sqrt(C)=1/2 for array bounds.',res{2});
% scalar eppos and M: M*alpha, and the gain in the norm of alpha
[res,~] = run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[0,0,0],[1,1,0],[NaN,NaN,NaN],2);
assert(~isempty(res),'(5b) the heat equation must satisfy sos1 and sos2 for the scalar M=2.');
check_C(res,1/4,eppos/2,ctol,'(5b)');   % M=sqrt(C*M/eppos)=sqrt(C/(eppos/2))
% an L2 lower bound with an H1 upper bound
[res,~] = run_quiet(PDE_heat,1.0,alpha,[0.1,0,0],0,[1,0,1],[3,0,0],[0,0,0],[1,1,0],[4,4,NaN],[0,1,0]);
assert(~isempty(res),'(5c) the heat equation must satisfy V>=0.1*||u||^2 and V<=C*||u_s||^2.');
assert(abs(res{1}-1/2)<ctol,'(5c) C=%.8g, expected 1/2.',res{1});
assert(abs(res{2}-sqrt(1/2))<10*ctol,'(5c) M=%.8g, expected sqrt(1/2).',res{2});
% array eppos, M omitted (M=1: M*alpha)
[res,~] = run_quiet(PDE_heat,1.0,alpha,[0,0.25,0],0,[1,0,1],[3,0,0],[0,0,0],[1,1,0]);
assert(~isempty(res) && abs(res{1}-1/2)<ctol && abs(res{2}-sqrt(1/2))<10*ctol,'(5d) C=1/2 and M=sqrt(C) expected for array eppos and M omitted.');
fprintf('Passed (5): eppos and M as arrays\n');

%% (6) Errors
threw = false;
try,    run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[-1,0,0],[1,0,0]);     catch,  threw = true;   end
assert(threw,'(6a) an error is expected for a negative psatz order.');
threw = false;
try,    run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[0.5,0,0],[1,0,0]);    catch,  threw = true;   end
assert(threw,'(6b) an error is expected for a non-integer psatz order.');
threw = false;
try,    run_quiet(PDE_heat,1.0,alpha,[1,2],0,[1,0,1],[3,0,0]);      catch,  threw = true;   end
assert(threw,'(6c) an error is expected for eppos of another size than alpha.');
threw = false;
try,    run_quiet(PDE_heat,1.0,alpha,eppos,0,[1,0,1],[3,0,0],[],[],[],[1,2]);   catch,  threw = true;   end
assert(threw,'(6d) an error is expected for M of another size than alpha.');

disp('Test_LocalStability_v2: all tests passed.')
end


function [res,txt] = run_quiet(PDE,varargin)
% Run LocalStability, suppressing output
txt = evalc('res = LocalStability_v2(PDE,varargin{:});');
end


function check_C(res,Cexp,eppos,tol,name)
% Check the optimal value C and M=sqrt(C/eppos)
assert(abs(res{1}-Cexp)<tol,'%s C=%.8g, expected %.8g.',name,res{1},Cexp);
assert(abs(res{2}-sqrt(Cexp/eppos))<10*tol,'%s M=%.8g, expected %.8g.',name,res{2},sqrt(Cexp/eppos));
end


function Y = eval_V(V,xval)
% Y = EVAL_V(V,XVAL) evaluates the scalar functional V (polyopvar of 'intop'
% coefficients, in a single scalar state depending on the variable s) at the
% state x=XVAL, polynomial in s. For a monomial of d factors, the dummy
% variables of the factors are s_1,...,s_d, so that the input to the
% operator is kron(x(s_1),...,x(s_d)) (as in apply_op of 'intop' and the
% evaluators of Test_piesos_LF).
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

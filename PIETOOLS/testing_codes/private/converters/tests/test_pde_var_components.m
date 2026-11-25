function results = test_pde_var_components()
% TEST_PDE_VAR_COMPONENTS Tests of numeric indexing x(idx) of a PDE
% variable x of size n, which returns the components I(idx,:)*x:
%   1. x(i) and x([i,j]) equal the products with rows of the identity,
%      for states with and without spatial variables, x(1:n) = x, and
%      x(1)^3 is a scalar power;
%   2. linear PDEs written with components give the same PIE as the same
%      PDEs written with matrices;
%   3. nonlinear terms x(1)^3 and x(1)*x(2): the PDE state and right-hand
%      side f of the PIE, evaluated at the fundamental state of a test
%      state, equal the values computed by hand;
%   4. indices outside 1..n, non-integer, empty or two indices give errors;
%   5. components of inputs and outputs, also in a PDE;
%   6. evaluation x(t-1,s), x(t,0), x(s1,0) and term indexing of sums of
%      terms and of systems are unchanged;
%   7. x(end), T(end), P(end) are the last component, term and equation,
%      and isequal compares 'pde_struct' objects;
%   8. rows of terms of different sizes, as in x+[x(2);x(1)], are added;
%   9. operations on components (diff, subs, int, delay, sums) equal those
%      on the products with rows of the identity, also in 2D;
%  10. end in other positions and for other objects; isequal is
%      structural, with polynomials compared by value;
%  11. stacked rows: commuted sums, terms moved over ==, stacked left-hand
%      sides of separate states and outputs, zero rows;
%  12. matrices times nonlinear terms (non-square, zero), x.*y, s.*x,
%      and errors for other factors;
%  13. PIEs from copvar operators (piess), and errors for other types.
%
% Run as: results = test_pde_var_components;
%
% OUTPUT
% - results: 'struct' array with fields name, passed and message, one
%            element per test.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_pde_var_components
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 11/24/2026

tests = {@case1_components, @case2_linear_pie, @case3_nonlinear_pie, ...
         @case4_errors, @case5_inputs_outputs, @case6_unchanged, ...
         @case7_end_isequal, @case8_rows_of_different_sizes, ...
         @case9_component_operations, @case10_end_isequal_more, ...
         @case11_stacked_rows, @case12_products, @case13_piess};
results = struct('name',cell(size(tests)),'passed',false,'message','');
for k = 1:numel(tests)
    results(k).name = func2str(tests{k});
    try
        tests{k}();
        results(k).passed = true;
        fprintf('%-30s passed\n',results(k).name);
    catch ME
        results(k).message = ME.message;
        fprintf('%-30s FAILED: %s\n',results(k).name,ME.message);
    end
end
if nargout==0
    clear results
end
end



function PIE = cv(P,opts)
% PIE = CV(P,OPTS) converts the PDE P to a PIE without printing.
% INPUT
% - P:      'pde_struct', the PDE;
% - opts:   (optional) cellstr, options of convert, {'silent'} by default;
% OUTPUT
% - PIE:    'pie_struct' or 'npie_struct', the PIE.
if nargin<2
    opts = {'silent'};
end
evalc('PIE = convert(P,opts);');
end

function e = poly_equal(A,B)
% E = POLY_EQUAL(A,B) checks that two polynomial arrays agree to 1e-10.
% INPUT
% - A, B:   'polynomial' or 'double' arrays of equal size;
% OUTPUT
% - e:      logical, true if all coefficients of A-B are below 1e-10.
D = polynomial(A-B);
e = isempty(D.coefficient) || full(max(abs(D.coefficient(:))))<1e-10;
end

function expect_error(f,msg_part,label)
% EXPECT_ERROR(F,MSG_PART,LABEL) checks that calling F throws an error
% whose message contains MSG_PART.
% INPUT
% - f:          function handle without arguments;
% - msg_part:   char, part of the expected error message;
% - label:      char, the call, for the failure message.
try
    tmp = f(); %#ok<NASGU>
catch ME
    assert(contains(ME.message,msg_part),'%s threw "%s", not the expected error.',label,ME.message);
    return
end
error('%s should throw an error.',label);
end


function case1_components()
% x(idx) equals I(idx,:)*x for a vector-valued state in 1 and 2 spatial
% variables and an ODE state; x(1:n) = x; x(1) = x for a scalar state.
pvar s s1 s2
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
assert(isequal(x(1),[1 0]*x),'x(1) differs from [1 0]*x.');
assert(isequal(x(2),[0 1]*x),'x(2) differs from [0 1]*x.');
assert(~isequal(x(1),x(2)) && ~isequal(x(1),x),'x(1) should differ from x(2) and from x.');
assert(isequal(x([2,1]),[0 1;1 0]*x),'x([2,1]) differs from [0 1;1 0]*x.');
assert(isequal(x([1,1]),[1 0;1 0]*x),'x([1,1]) differs from [1 0;1 0]*x.');
assert(isequal(x(1:2),x) && isequal(x([1;2]),x),'x(1:2) differs from x.');
xs = pde_var('state',1,s,[0,1]);
assert(isequal(xs(1),xs),'For a scalar state, x(1) differs from x.');
x3 = pde_var('state',3,[s1;s2],[0,1;0,1]);
E = eye(3);
assert(isequal(x3([3,1]),E([3,1],:)*x3),'For a 2D state, x([3,1]) differs from I([3,1],:)*x.');
X = pde_var('state',2);
assert(isequal(X(2),[0 1]*X),'For an ODE state, X(2) differs from [0 1]*X.');
% A component is a scalar, so its power is defined.
assert(isequal(x(1)^3,([1 0]*x)^3),'x(1)^3 differs from ([1 0]*x)^3.');
assert(isequal(x(1)*x(2),([1 0]*x)*([0 1]*x)),'x(1)*x(2) differs from ([1 0]*x)*([0 1]*x).');
end

function case2_linear_pie()
% A PDE written with components, x1_t = x1_ss + x2, x2_t = x2_ss with
% Dirichlet BCs, and an ODE-PDE system, give the same PIE as when written
% with matrices; the PIE differs from that without the coupling.
pvar s t
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
PIE1 = cv([diff(x,t)==[diff(x(1),s,2)+x(2); diff(x(2),s,2)];
           subs(x(1),s,0)==0; subs(x(1),s,1)==0; subs(x(2),s,0)==0; subs(x(2),s,1)==0]);
PIE2 = cv([diff(x,t)==diff(x,s,2)+[0 1;0 0]*x;
           subs([1 0]*x,s,0)==0; subs([1 0]*x,s,1)==0; subs([0 1]*x,s,0)==0; subs([0 1]*x,s,1)==0]);
assert(PIE1==PIE2,'With components, the PIE differs from that with matrices.');
PIE3 = cv([diff(x,t)==diff(x,s,2)+[0 1;0 0]*x; subs(x,s,0)==0; subs(x,s,1)==0]);
assert(PIE1==PIE3,'With components, the PIE differs from that with vector BCs.');
PIE0 = cv([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0]);
assert(~(PIE1==PIE0),'The coupling x2 in the equation of x1 is lost.');
% ODE-PDE: X1_t = -X1 + int x2, X2_t = -X2, x1_t = x1_ss + x2,
% x2_t = x2_ss + s*X1, with Robin BC x1_s(1) = -x1(1).
clear('stateNameGenerator');
X = pde_var('state',2);     x = pde_var('state',2,s,[0,1]);
PIE1 = cv([diff(X,t)==[-X(1)+int(x(2),s,[0,1]); -X(2)];
           diff(x,t)==[diff(x(1),s,2)+x(2); diff(x(2),s,2)+s*X(1)];
           subs(x(1),s,0)==0; subs(diff(x(1),s),s,1)==-subs(x(1),s,1);
           subs(x(2),s,0)==0; subs(x(2),s,1)==0]);
PIE2 = cv([diff(X,t)==-X+int([0 1;0 0]*x,s,[0,1]);
           diff(x,t)==diff(x,s,2)+[0 1;0 0]*x+[0 0;s 0]*X;
           [1 0]*subs(x,s,0)==0; [1 0]*subs(diff(x,s),s,1)==-[1 0]*subs(x,s,1);
           [0 1]*subs(x,s,0)==0; [0 1]*subs(x,s,1)==0]);
assert(PIE1==PIE2,'For the ODE-PDE system, the PIE with components differs from that with matrices.');
end

function case3_nonlinear_pie()
% For x_t = x_ss + N(x) with Dirichlet BCs, the PIE state is xf = x_ss,
% x_pde = T*xf, and f(xf) = xf + N(T*xf). At the test state
% x0 = [s(1-s); s^2(1-s)], check x_pde(xf) = x0 and f(xf) = xf + N(x0) for
% N(x) = -[1;0]*x1^3, 2*[0;1]*x1*x2, and the sum -[1;0]*x1^3+[0;1]*x1*x2,
% also written as one row per component and as a column of terms.
pvar s t
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
BCs = [subs(x,s,0)==0; subs(x,s,1)==0];
cases = {[diff(x,t)==diff(x,s,2)-[1;0]*x(1)^3; BCs], @(x1,x2) [-x1^3; 0];
         [diff(x,t)==diff(x,s,2)+[0;2]*(x(1)*x(2)); BCs], @(x1,x2) [0; 2*x1*x2];
         [diff(x,t)==diff(x,s,2)-[1;0]*x(1)^3+[0;1]*(x(1)*x(2)); BCs], @(x1,x2) [-x1^3; x1*x2];
         [diff(x,t)==[diff(x(1),s,2)-x(1)^3; diff(x(2),s,2)+x(1)*x(2)]; BCs], @(x1,x2) [-x1^3; x1*x2];
         [diff(x,t)==diff(x,s,2)+[-x(1)^3; x(1)*x(2)]; BCs], @(x1,x2) [-x1^3; x1*x2]};
for k=1:size(cases,1)
    N = cv(cases{k,1},{'silent','nonlinear'});
    assert(isa(N,'npie_struct'),'Case %d: the PIE should be nonlinear.',k);
    assert(isscalar(N.f.varname) && N.f.dims.in(1)==2,'Case %d: expected one PIE state of size 2.',k);
    v = polynomial(N.space.vars(1));
    x0 = [v*(1-v); v^2*(1-v)];      xf = diff(diff(x0,v),v);
    assert(poly_equal(apply_poly(N.x_pde,{xf}),x0),'Case %d: x_pde(xf) differs from the test state.',k);
    assert(poly_equal(apply_poly(N.f,{xf}),xf+cases{k,2}(x0(1),x0(2))),'Case %d: f(xf) differs from the hand value.',k);
end
end

function case4_errors()
% Indices outside 1..n, non-integer, empty, matrix or two indices.
pvar s
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
expect_error(@() x(3),'Component indices','x(3)');
expect_error(@() x(0),'Component indices','x(0)');
expect_error(@() x(-1),'Component indices','x(-1)');
expect_error(@() x(1.5),'Component indices','x(1.5)');
expect_error(@() x([]),'Component indices','x([])');
expect_error(@() x([1 2;2 1]),'Component indices','x([1 2;2 1])');
expect_error(@() x(1,1),'single index','x(1,1)');
expect_error(@() x(1,2),'single index','x(1,2)');
X = pde_var('state',1);
expect_error(@() X(2),'Component indices','X(2) for a scalar ODE state');
% Equations for single components are not supported.
pvar t
z = pde_var('output',2);
expect_error(@() diff(x(1),t)==diff(x(1),s,2),'single components','diff(x(1),t)==...');
expect_error(@() z(1)==int(x(1),s,[0,1]),'single components','z(1)==...');
expect_error(@() [diff(x(1),t);diff(x(2),t)]==[diff(x(1),s,2);diff(x(2),s,2)],'single components','[diff(x(1),t);diff(x(2),t)]==...');
end

function case5_inputs_outputs()
% Components of inputs and outputs, and a PDE with components of a
% distributed input, a finite-dimensional control and in an output
% equation, against the same PDE with matrices.
pvar s t
clear('stateNameGenerator');
w = pde_var('input',2,s,[0,1]);     u = pde_var('control',2);
z = pde_var('output',2);            y = pde_var('sense',1);
assert(isequal(w(2),[0 1]*w),'w(2) differs from [0 1]*w.');
assert(isequal(u([2,1]),[0 1;1 0]*u),'u([2,1]) differs from [0 1;1 0]*u.');
assert(isequal(z(1),[1 0]*z),'z(1) differs from [1 0]*z.');
assert(isequal(y(1),y),'For a scalar output, y(1) differs from y.');
expect_error(@() w(3),'Component indices','w(3)');
x = pde_var('state',2,s,[0,1]);
BCs = [subs(x,s,0)==0; subs(x,s,1)==0];
PIE1 = cv([diff(x,t)==[diff(x(1),s,2)+w(2); diff(x(2),s,2)+s*w(1)+u(2)];
           z==[int(x(1),s,[0,1]); subs(diff(x(2),s),s,1)]; BCs]);
PIE2 = cv([diff(x,t)==diff(x,s,2)+[0 1;s 0]*w+[0 0;0 1]*u;
           z==[1 0;0 0]*int(x,s,[0,1])+[0 0;0 1]*subs(diff(x,s),s,1); BCs]);
assert(PIE1==PIE2,'With components of inputs and outputs, the PIE differs from that with matrices.');
end

function case6_unchanged()
% Evaluation with a 'polynomial' argument, and (equation, term) indexing
% of objects other than a single PDE variable, are unchanged.
pvar s t s1 s2
clear('stateNameGenerator');
x = pde_var('state',1,s,[0,1]);
assert(isequal(x(t-1,s),subs(x,t,t-1)),'x(t-1,s) differs from subs(x,t,t-1).');
assert(isequal(x(t,0),subs(x,s,0)),'x(t,0) differs from subs(x,s,0).');
x2 = pde_var('state',2,[s1;s2],[0,1;0,1]);
assert(isequal(x2(s1,0),subs(x2,s2,0)),'x(s1,0) differs from subs(x,s2,0).');
% A sum of terms: T(2) is the second term, also as T(1,2).
T = diff(x,s,2)+3*x;
T2 = T(2);
assert(isscalar(T2.free{1}.term) && isequal(double(polynomial(T2.free{1}.term{1}.C)),3),'T(2) is not the second term of T.');
assert(isequal(T(1,2),T2),'T(1,2) differs from T(2).');
% A system: P(2) is the second equation.
P = [diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0];
P2 = P(2);
assert(isscalar(P2.free) && ~is_pde_var(P2),'P(2) is not a single equation.');
end

function case7_end_isequal()
% x(end) is the last component of a PDE variable, T(end) and T(1,end) the
% last term of a sum of terms, P(end) the last equation of a system; isequal
% compares 'pde_struct' objects, also with 'polynomial' values stored.
pvar s t
clear('stateNameGenerator');
x3 = pde_var('state',3,s,[0,1]);
assert(isequal(x3(end),x3(3)) && isequal(x3(end-1),x3(2)),'x(end) is not the last component.');
assert(isequal(x3([1,end]),x3([1,3])),'x([1,end]) differs from x([1,3]).');
x = pde_var('state',1,s,[0,1]);
T = diff(x,s,2)+3*x+s*x;
assert(isequal(T(end),T(3)) && isequal(T(1,end),T(3)) && ~isequal(T(end),T(1)),'T(end) is not the last term.');
P = [diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0];
assert(isequal(P(end),P(3)) && ~isequal(P(end),P(1)),'P(end) is not the last equation.');
assert(isequal(x3,x3) && isequal(x3,x3,x3),'isequal(x,x) should be true.');
assert(~isequal(x3,x) && ~isequal(x3,2) && ~isequal(2,x3) && ~isequal(x3(1),x3(2)),'isequal of different objects should be false.');
end

function case8_rows_of_different_sizes()
% Sums of columns of terms with rows of different sizes, x+[x(2);s*x(1)]
% and [X(1);X([2,3])]+[X([1,2]);X(3)], give the same PIE as with matrices.
pvar s t
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
BCs = [subs(x,s,0)==0; subs(x,s,1)==0];
PIE1 = cv([diff(x,t)==diff(x,s,2)+[x(2); s*x(1)]; BCs]);
PIE2 = cv([diff(x,t)==diff(x,s,2)+[0 1;s 0]*x; BCs]);
assert(PIE1==PIE2,'diff(x,s,2)+[x(2);s*x(1)] gives a different PIE than with a matrix.');
X = pde_var('state',3);
Q = [X(1);X([2,3])]+[X([1,2]);X(3)];
assert(isscalar(Q.free) && Q.free{1}.size==3,'The sum should be one row of size 3.');
assert(cv(diff(X,t)==Q)==cv(diff(X,t)==2*X),'[X(1);X([2,3])]+[X([1,2]);X(3)] gives a different PIE than 2*X.');
expect_error(@() X+[X(1);X(2)],'rows','X+[X(1);X(2)]');
end
function nl_check(P,Nfun,label)
% NL_CHECK(P,NFUN,LABEL) converts the nonlinear PDE P, x_t = x_ss + N(x)
% for one state x of size 2 with Dirichlet BCs, and checks that
% f(xf) = xf + N(x0) at the test state x0 = [s(1-s); s^2(1-s)], xf = x0_ss.
% INPUT
% - P:      'pde_struct', the PDE;
% - Nfun:   function handle, N as a function of (x1,x2,s);
% - label:  char, the case, for the failure message.
N = cv(P,{'silent','nonlinear'});
assert(isa(N,'npie_struct') && isscalar(N.f.varname) && N.f.dims.in(1)==2,'%s: expected a nonlinear PIE with one state of size 2.',label);
v = polynomial(N.space.vars(1));
x0 = [v*(1-v); v^2*(1-v)];      xf = diff(diff(x0,v),v);
assert(poly_equal(apply_poly(N.f,{xf}),xf+Nfun(x0(1),x0(2),v)),'%s: f(xf) differs from the hand value.',label);
end

function case9_component_operations()
% Derivatives, boundary values, integrals and delays of components equal
% those of the products with rows of the identity; sums of components,
% delayed components, BCs coupling components, and a 2D PDE give the PIE
% of the same PDE written with matrices.
pvar s t s1 s2
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
e1 = [1 0];     e2 = [0 1];
assert(isequal(x(:),x),'x(:) differs from x.');
assert(isequal(x(int8(2)),x(2)),'x(int8(2)) differs from x(2).');
assert(isequal(diff(x(2),s),diff(e2*x,s)),'diff(x(2),s) differs from diff([0 1]*x,s).');
assert(isequal(subs(x(1),s,0),subs(e1*x,s,0)),'subs(x(1),s,0) differs from subs([1 0]*x,s,0).');
assert(isequal(int(x(2),s,[0,1]),int(e2*x,s,[0,1])),'int(x(2),s,[0,1]) differs from int([0 1]*x,s,[0,1]).');
assert(isequal(subs(x(2),t,t-1),subs(e2*x,t,t-1)),'subs(x(2),t,t-1) differs from subs([0 1]*x,t,t-1).');
BCs = [subs(x,s,0)==0; subs(x,s,1)==0];
assert(cv([diff(x,t)==diff(x,s,2)+[x(1)+x(2); x(1)-2*x(2)]; BCs])==cv([diff(x,t)==diff(x,s,2)+[1 1;1 -2]*x; BCs]),...
        'Sums of components give a different PIE than [1 1;1 -2]*x.');
assert(cv([diff(x,t)==[diff(x(1),s,2)+subs(x(2),t,t-1); diff(x(2),s,2)]; BCs])==cv([diff(x,t)==diff(x,s,2)+[0 1;0 0]*subs(x,t,t-1); BCs]),...
        'A delayed component gives a different PIE than [0 1;0 0]*x(t-1).');
assert(cv([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(diff(x(1),s),s,1)==subs(x(2),s,1); subs(x(2),s,1)==0]) ...
        ==cv([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; e1*subs(diff(x,s),s,1)==e2*subs(x,s,1); e2*subs(x,s,1)==0]),...
        'BCs coupling components give a different PIE than with matrices.');
clear('stateNameGenerator');
y = pde_var('state',2,[s1;s2],[0,1;0,1]);
BC2 = [subs(y,s1,0)==0; subs(y,s1,1)==0; subs(y,s2,0)==0; subs(y,s2,1)==0];
assert(cv([diff(y,t)==[diff(y(1),s1,2)+diff(y(1),s2,2)+y(2); diff(y(2),s1,2)+diff(y(2),s2,2)]; BC2]) ...
        ==cv([diff(y,t)==diff(y,s1,2)+diff(y,s2,2)+[0 1;0 0]*y; BC2]),...
        'A 2D PDE with components gives a different PIE than with matrices.');
end

function case10_end_isequal_more()
% end for parts of components, ODE states, scalar states, inputs, sums with
% nonlinear terms and systems with outputs; isequal is structural, with
% polynomial coefficients compared by value.
pvar s t
clear('stateNameGenerator');
x3 = pde_var('state',3,s,[0,1]);
assert(isequal(x3(end-1:end),x3([2,3])),'x(end-1:end) differs from x([2,3]).');
X = pde_var('state',2);
assert(isequal(X(end),X(2)),'For an ODE state, X(end) differs from X(2).');
xs = pde_var('state',1,s,[0,1]);
assert(isequal(xs(end),xs),'For a scalar state, x(end) differs from x.');
w = pde_var('input',2,s,[0,1]);
assert(isequal(w(end),w(2)),'w(end) differs from w(2).');
T = diff(xs,s,2) - xs^3 + s*xs;
assert(isequal(T(end),T(3)) && ~isequal(T(end),T(2)),'T(end) is not the last term of a sum with a nonlinear term.');
z = pde_var('output',1);
P = [diff(xs,t)==diff(xs,s,2); z==int(xs,s,[0,1]); subs(xs,s,0)==0; subs(xs,s,1)==0];
assert(isequal(P(end),P(4)) && isequal(P(1,end),P(1,1)),'P(end) is not the last equation of a system with an output.');
y1 = pde_var('state',2,s,[0,1]);     y2 = pde_var('state',2,s,[0,1]);
assert(~isequal(y1,y2),'Different states of the same size should differ.');
assert(isequal((1+s)*xs,(s+1)*xs),'(1+s)*x should equal (s+1)*x.');
assert(~isequal((1+s)*xs,s*xs+xs),'One term and two terms are structurally different.');
assert(~isequal(xs,xs,y1),'isequal of three objects, one different, should be false.');
P2 = [diff(xs,t)==diff(xs,s,2); z==int(xs,s,[0,1]); subs(xs,s,0)==0; subs(xs,s,1)==0];
assert(isequal(P,P2) && isequal(initialize(P),initialize(P2)),'The same system declared twice should be equal.');
assert(~isequal(P,[diff(xs,t)==diff(xs,s,2); z==int(xs,s,[0,1]); subs(xs,s,1)==0; subs(xs,s,0)==0]),'BCs in another order are structurally different.');
end

function case11_stacked_rows()
% Stacked rows of components in any position of a sum or equation, stacked
% left-hand sides of separate states and outputs, and errors for rows of
% different total size.
pvar s t
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
BCs = [subs(x,s,0)==0; subs(x,s,1)==0];
PIEm = cv([diff(x,t)==diff(x,s,2)+[0 1;1 0]*x; BCs]);
assert(cv([diff(x,t)==[x(2);x(1)]+diff(x,s,2); BCs])==PIEm,'[x(2);x(1)]+diff(x,s,2) gives a different PIE.');
assert(cv([diff(x,t)-[x(2);x(1)]==diff(x,s,2); BCs])==PIEm,'diff(x,t)-[x(2);x(1)]==diff(x,s,2) gives a different PIE.');
assert(cv([diff(x,t)-diff(x,s,2)-[x(2);x(1)]==0; BCs])==PIEm,'All terms on the left-hand side give a different PIE.');
expect_error(@() x+x(1),'rows','x+x(1)');
expect_error(@() x+[x(1);x(2);x(1)],'rows','x+[x(1);x(2);x(1)]');
% A row of zeros, in a sum and in an equation.
PIEz = cv([diff(x,t)==diff(x,s,2)+[0 1;0 0]*x; BCs]);
assert(cv([diff(x,t)==diff(x,s,2)+[x(2);0]; BCs])==PIEz,'diff(x,s,2)+[x(2);0] gives a different PIE.');
assert(cv([diff(x,t)==[x(2);0]+diff(x,s,2); BCs])==PIEz,'[x(2);0]+diff(x,s,2) gives a different PIE.');
assert(cv([diff(x,t)==diff(x,s,2)+[0;x(1)]; BCs])==cv([diff(x,t)==diff(x,s,2)+[0 0;1 0]*x; BCs]),'diff(x,s,2)+[0;x(1)] gives a different PIE.');
expect_error(@() x+[x(1);x(2);0],'rows','x+[x(1);x(2);0]');
X = pde_var('state',2);
assert(cv(diff(X,t)==[X(2);0])==cv(diff(X,t)==[0 1;0 0]*X),'diff(X,t)==[X(2);0] gives a different PIE.');
clear('stateNameGenerator');
a = pde_var('state',1,s,[0,1]);     b = pde_var('state',1,s,[0,1]);
z1 = pde_var('output',1);           z2 = pde_var('output',1);
BCab = [subs(a,s,0)==0; subs(a,s,1)==0; subs(b,s,0)==0; subs(b,s,1)==0];
PIE1 = cv([[diff(a,t);diff(b,t)]==[diff(a,s,2)+b; diff(b,s,2)]; [z1;z2]==[int(a,s,[0,1]); int(b,s,[0,1])]; BCab]);
PIE2 = cv([diff(a,t)==diff(a,s,2)+b; diff(b,t)==diff(b,s,2); z1==int(a,s,[0,1]); z2==int(b,s,[0,1]); BCab]);
assert(PIE1==PIE2,'Stacked left-hand sides give a different PIE than separate equations.');
end

function case12_products()
% Non-square and zero matrices times nonlinear terms, elementwise products
% of components, polynomial coefficients; .* with scalars; errors for
% other factors.
pvar s t
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);
BCs = [subs(x,s,0)==0; subs(x,s,1)==0];
nl_check([diff(x,t)==[diff(x(1),s,2)+[1 0]*(x.*x); diff(x(2),s,2)]; BCs],@(x1,x2,v) [x1^2; 0],'[1 0]*(x.*x)');
nl_check([diff(x,t)==diff(x,s,2)+[2;0]*(x(2)*x(2)); BCs],@(x1,x2,v) [2*x2^2; 0],'[2;0]*(x(2)*x(2))');
nl_check([diff(x,t)==diff(x,s,2)+[0;1]*(x(1).*x(2)); BCs],@(x1,x2,v) [0; x1*x2],'[0;1]*(x(1).*x(2))');
nl_check([diff(x,t)==diff(x,s,2)+[1;0]*(s*x(1)^2); BCs],@(x1,x2,v) [v*x1^2; 0],'[1;0]*(s*x(1)^2)');
nl_check([diff(x,t)==diff(x,s,2)+[1 1;0 0]*(x.*x); BCs],@(x1,x2,v) [x1^2+x2^2; 0],'[1 1;0 0]*(x.*x)');
nl_check([diff(x,t)==diff(x,s,2)+[1;0]*(x(1)^2*x(2))-[0;1]*x(2)^3; BCs],@(x1,x2,v) [x1^2*x2; -x2^3],'x(1)^2*x(2) and x(2)^3');
assert(cv([diff(x,t)==diff(x,s,2)+zeros(2,1)*(x(1)*x(2)); BCs])==cv([diff(x,t)==diff(x,s,2); BCs]),...
        'A zero matrix times a nonlinear term should add nothing.');
assert(isequal(2.*x,2*x) && isequal(s.*x,s*x),'Scalar .* differs from scalar *.');
expect_error(@() [1 2].*x,'Elementwise','[1 2].*x');
expect_error(@() struct('P',1)*x,'first factor','struct*x');
expect_error(@() times(struct('P',1),x),'first factor','struct.*x');
end

function case13_piess()
% A PIE built by piess from the operators of a converted PIE has the same
% operators; matrices are accepted; other types give errors.
pvar s t
clear('stateNameGenerator');
x = pde_var('state',2,s,[0,1]);     w = pde_var('input',1);     z = pde_var('output',1);
PIE = cv([diff(x,t)==diff(x,s,2)+[1;0]*w; z==int(x(1),s,[0,1]); subs(x,s,0)==0; subs(x,s,1)==0]);
P1 = piess(PIE.T,PIE.A,PIE.Bw,PIE.Cz,PIE.Dzw);
assert(P1.T==PIE.T && P1.A==PIE.A && P1.Bw==PIE.Bw && P1.Cz==PIE.Cz && P1.Dzw==PIE.Dzw,'piess changed the operators.');
P2 = piess(eye(2),-eye(2));
assert(isa(P2,'pie_struct') && isa(P2.T,'copvar'),'piess with matrices should give a PIE with copvar operators.');
expect_error(@() piess(struct('a',1),1),'Operators Top','piess(struct,1)');
expect_error(@() initialize(setfield(pie_struct(),'T',struct('a',1))),'copvar','PIE.T a struct');
expect_error(@() terms2pde_struct([]),'legacy','terms2pde_struct');
end
function results = test_convert_PIETOOLS_PDE_static_constraint()
% TEST_CONVERT_PIETOOLS_PDE_STATIC_CONSTRAINT Static constraints 0=E(s).
% An equation 0 = E(s) that holds on the whole domain, such as the
% continuity equation of the Navier-Stokes equations,
%       0 = d/ds1 v2 + d/ds2 v3,
% is a static constraint on the state, not a boundary condition. It is
% declared together with the PDE, and 'initialize' stores it as an output
% flagged in PDE.constraint_ID. The converter then adds its operator (acting on the
% fundamental state) as extra rows of PIE.K:  K*xf = 0.
%
% Run as: results = test_convert_PIETOOLS_PDE_static_constraint;
%
% DB, 11/07/2026: Initial coding;
% DB, 11/08/2026: Add the tests of PDEs declared as a struct, without
%                   variables for the states, and with a boundary term on
%                   several states;

tests = {@constraint_1d, @navier_stokes_2d_dirichlet, ...
         @struct_declared_no_vars, @library_examples_no_vars};
results = struct('name',cell(size(tests)),'passed',false,'message','');
for k = 1:numel(tests)
    results(k).name = func2str(tests{k});
    try
        tests{k}();
        results(k).passed = true;
    catch ME
        results(k).message = ME.message;
    end
    fprintf('%-34s %s\n',results(k).name,ternary(results(k).passed,'passed',['FAILED: ',results(k).message]));
end
if nargout==0
    clear results
end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function constraint_1d()
% x1: Dirichlet heat equation, x2: ODE in s, with the constraint
% 0 = d/ds x1 + x2.
pvar t s
x1 = pde_var(1,s,[0 1]);    x2 = pde_var(1,s,[0 1]);
eqs = [diff(x1,t)==diff(x1,s,2); diff(x2,t)==-x2+x1;
       subs(x1,s,0)==0; subs(x1,s,1)==0;
       0==diff(x1,s)+x2];
PDE = initialize(eqs,true);
assert(numel(PDE.BC)==2,'A static constraint is not a boundary condition.');
%PDE.constraint_ID==1. The continuity equation is flagged in initialize.;
PIE = convert(PDE,'pie');
[~,dm] = copvar_space_list(PIE.K,'out');
assert(isequal(dm(:)',1),'K must have a single row of size 1.');
% The constraint is not an output of the PIE.
assert(isempty(PIE.Cz) || isempty(PIE.Cz.C),'The constraint must not remain an output.');
end

function navier_stokes_2d_dirichlet()
% Linearized rotating-flow equations on [-1,1]x[0,pi] with Dirichlet BCs on
% all sides, and the continuity equation of the velocity.
pvar t s1 s2
R0 = 0.5;   L = pi;     Re = 10;
dom = [-1 1; 0 L];
v1 = pde_var(1,[s1;s2],dom);   
v2 = pde_var(1,[s1;s2],dom); 
v3 = pde_var(1,[s1;s2],dom);
V = {v1,v2,v3};
lap = @(v) diff(v,s1,2)+diff(v,s2,2);
eqs = [diff(v1,t)==(1/Re)*lap(v1)+(R0-1)*v2;
       diff(v2,t)==(1/Re)*lap(v2)-R0*v1;
       diff(v3,t)==(1/Re)*lap(v3)];
for k = 1:3
    v = V{k};
    eqs = [eqs; 
        subs(v,s1,-1)== subs(v,s1,1); 
        subs(diff(v,s1),s1,-1)==subs(diff(v,s1),s1,1);
        subs(v,s2,0)==0;
        subs(v,s2,L)==0];
end
% Continuity, declared with the PDE
eqs = [eqs; 0==diff(v2,s1)+diff(v3,s2)];
%PDE.constraint_ID==1. The continuity equation is flagged in initialize.;
PDE = initialize(eqs);
PIE = convert(PDE,'pie');
[~,dm] = copvar_space_list(PIE.K,'out');
fprintf('K has %d block rows, of the sizes [%s].\n',numel(dm),num2str(dm(:)'));
assert(~isempty(PIE.K) && numel(dm)>=1,'K must contain the continuity row.');
display(PIE.K);
end

function struct_declared_no_vars()
% Two coupled heat equations, declared as a struct in the terms format: the
% states have no variables or domain (only the PDE has), and two boundary
% conditions are terms on several states, x=[1;2]. The static-constraint
% detection must skip these (they are boundary conditions), without error.
pvar s
PDE = struct();
PDE.vars = s;           PDE.dom = [0,1];
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.D = 2;
PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.D = 2;
PDE.BC{1}.term{1}.x = 1;    PDE.BC{1}.term{1}.loc = 0;
PDE.BC{2}.term{1}.x = 2;    PDE.BC{2}.term{1}.loc = 0;
% x1(1) = x2(1) and d/ds x1(1) + d/ds x2(1) = 0, each a term on 2 states
PDE.BC{3}.term{1}.x = [1;2];    PDE.BC{3}.term{1}.loc = 1;      PDE.BC{3}.term{1}.C = [1,-1];
PDE.BC{4}.term{1}.x = [1;2];    PDE.BC{4}.term{1}.loc = 1;      PDE.BC{4}.term{1}.D = 1;
PDE.BC{4}.term{1}.C = [1,1];
PIE = convert_PIETOOLS_PDE(PDE);
assert(~isempty(PIE.T) && ~isempty(PIE.A),'The PIE has not been constructed.');
end

function library_examples_no_vars()
% Library examples declared as structs without state variables, with boundary
% terms on several states (Timoshenko_Beam_2: BC{5}.term{1}.x=[2;4]).
exs = {@PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped, @PIETOOLS_PDE_Ex_Timoshenko_Beam_2};
for k = 1:numel(exs)
    [PDE_t,~] = exs{k}(0,{});
    PIE = convert_PIETOOLS_PDE(PDE_t);
    assert(~isempty(PIE.T),'The PIE of %s has not been constructed.',func2str(exs{k}));
end
end

function Test_piesos_program()
% TEST_PIESOS_PROGRAM tests, using assertions, the program utilities
%   prog = piesos_program(vartab,decvartab),
%   prog = piesos_setobj(prog,obj),
%   prog = piesos_solve(prog,opts),
%   sol  = piesos_getsol(prog,dvar),
% on programs with a handful of decision variables (function-level).
%
% Oracle, independent of the functions under test:
%  - piesos_program: the names passed in (spatial variables, state names,
%    decision variables);
%  - piesos_setobj: coefficients for the objective c2*2+c1, in the order of
%    prog.decvartable;
%  - piesos_solve: min c s.t. c>=1 has optimum 1; c>=1 and c<=0 is
%    infeasible;
%  - piesos_getsol: c1+c2=5, c1=3 gives c1=3, c2=2, so any linear
%    combination, in a dpvar, cell, intop, tensopmat or polyopvar
%    coefficient, is known in closed form.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_program
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

% NOTES
% - An equality-only program cannot be solved by SeDuMi ("Size b
%   mismatch"): solved_c1c2 and case_solve add a trivial PSD matrix with
%   sosposmatr;
% - piesos_setobj replaces the coefficient of a variable in a new
%   objective, and leaves those of variables not in it unchanged
%   (sossetobj behavior);
% - No distributed SOS constraints are tested beyond what is needed to
%   solve.

tests ={@case_program, @case_program_errors, @case_setobj, @case_setobj_errors, ...
         @case_solve, @case_solve_infeasible, @case_getsol_values, ...
         @case_getsol_passthrough, @case_getsol_operators, @case_getsol_errors};
nfail = 0;
for i=1:numel(tests)
    try
        tests{i}();
        fprintf('PASS  %s\n',func2str(tests{i}));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_piesos_program: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_program: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function x1 = state1()
x1 = polyopvar('x1',polynomial({'s'}),[0,1]);
end

function assert_err(f,needle)
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s".',err.message);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end

function prog = solved_c1c2()
% Program with c1+c2=5, c1=3, solved
dpvar c1 c2;
prog = piesos_program(state1(),[c1;c2]);
prog = soseq(prog,c1+c2-5);
prog = soseq(prog,c1-3);
[prog,~] = sosposmatr(prog,1);
evalc('prog = piesos_solve(prog);');
end

function v = val(prog,d)
v = double(piesos_getsol(prog,d));
end


%% ---- cases ------------------------------------------------------------
function case_program()
x1 = state1();
x2 = polyopvar('x2',polynomial({'s';'t'}),[0 1;0 1]);
prog = piesos_program([x1;x2]);
assert(isstruct(prog),'The program must be a struct.');
assert(isequal(sort(prog.pievartable(:)),{'x1';'x2'}),'State names must be in pievartable.');
assert(isequal(sort(prog.vartable(:)),{'s';'t'}),'Spatial variables of all states must be in vartable.');
assert(isempty(prog.decvartable),'No decision variables were declared.');
dpvar c1 c2;
prog = piesos_program(x1,[c1;c2]);
assert(isequal(sort(prog.decvartable(:)),{'c1';'c2'}),'Decision variables must be declared.');
assert(isequal(prog.vartable(:),{'s'}),'Only variable s.');
assert(isequal(prog.pievartable(:),{'x1'}),'State name x1.');
% a state without spatial variables
y = polyopvar('y');
prog = piesos_program(y);
assert(isempty(prog.vartable),'No spatial variables.');
assert(isequal(prog.pievartable(:),{'y'}),'State name y.');
end

function case_program_errors()
dpvar c1;
assert_err(@() piesos_program(3),'polyopvar');
assert_err(@() piesos_program(polynomial({'s'})),'polyopvar');
assert_err(@() piesos_program(state1(),5),'dpvar');
assert_err(@() piesos_program(state1(),polynomial({'s'})),'dpvar');
piesos_program(state1(),c1);
end

function case_setobj()
dpvar c1 c2;
prog = piesos_program(state1(),c1);
prog = piesos_setobj(prog,c1+2*c2);
assert(isequal(sort(prog.decvartable(:)),{'c1';'c2'}),'The new variable c2 must be declared.');
[~,idx] = ismember({'c1';'c2'},prog.decvartable(:));
cf = full(prog.objective(idx));
assert(isequal(cf(:),[1;2]),'Objective coefficients must be 1 for c1 and 2 for c2.');
prog2 = piesos_setobj(prog,3*c1);
assert(numel(prog2.decvartable)==2,'Declared variables must not be redeclared.');
[~,idx] = ismember('c1',prog2.decvartable(:));
assert(full(prog2.objective(idx))==3,'The coefficient of c1 in a new objective must replace the old one.');
% terms combined before use: c1+c1-c1
prog3 = piesos_setobj(prog,c1+c1-c1);
[~,idx] = ismember('c1',prog3.decvartable(:));
assert(full(prog3.objective(idx))==1,'Terms of the objective must be combined.');
end

function case_setobj_errors()
dpvar c1;
prog = piesos_program(state1(),c1);
s = polynomial({'s'});
assert_err(@() piesos_setobj(3,c1),'struct');
assert_err(@() piesos_setobj(prog,3),'dpvar');
assert_err(@() piesos_setobj(prog,s),'dpvar');
assert_err(@() piesos_setobj(prog,c1*s),'independent variables');
end

function case_solve()
% min c s.t. c>=1: the optimum is c=1
dpvar c;
prog = piesos_program(state1(),c);
prog = sosineq(prog,c-1);
prog = piesos_setobj(prog,c);
[prog,~] = sosposmatr(prog,1);
prog1 = prog;
evalc('prog = piesos_solve(prog);');
assert(prog.solinfo.info.pinf==0 && prog.solinfo.info.dinf==0,'Program must be feasible.');
assert(abs(val(prog,c)-1)<1e-5,'The optimum is c=1.');
% options are passed on: the same program with SeDuMi
evalc('progs = piesos_solve(prog1,struct(''solver'',''sedumi''));');
assert(progs.solinfo.info.pinf==0 && progs.solinfo.info.dinf==0,'Program must be feasible with SeDuMi.');
assert(abs(val(progs,c)-1)<1e-5,'The optimum is c=1 with SeDuMi.');
end

function case_solve_infeasible()
% c>=1 and c<=0
dpvar c;
prog = piesos_program(state1(),c);
prog = sosineq(prog,c-1);
prog = sosineq(prog,-c);
evalc('prog = piesos_solve(prog);');
info = prog.solinfo.info;
assert(info.pinf~=0 || info.dinf~=0,'The program must be infeasible.');
end

function case_getsol_values()
dpvar c1 c2;
prog = solved_c1c2();
assert(abs(val(prog,c1)-3)<1e-6 && abs(val(prog,c2)-2)<1e-6,'c1=3, c2=2 expected.');
assert(abs(val(prog,2*c1+c2)-8)<1e-6,'2*c1+c2=8 expected.');
sol = piesos_getsol(prog,[c1;c2-c1;c1-c1+4]);
assert(max(abs(double(sol)-[3;-1;4]))<1e-6,'Vector of expressions: [3;-1;4] expected.');
% cell: each element, with the size of the cell, and constants unchanged
C = piesos_getsol(prog,{c1,2*c2;[],3});
assert(isequal(size(C),[2,2]),'The size of the cell must be kept.');
assert(abs(double(C{1,1})-3)<1e-6 && abs(double(C{1,2})-4)<1e-6 && isempty(C{2,1}) && double(C{2,2})==3,...
    'Elements of the cell: 3, 4, [], 3 expected.');
C =piesos_getsol(prog,{{c1,c2}});
assert(abs(double(C{1}{1})-3)<1e-6 && abs(double(C{1}{2})-2)<1e-6,'Nested cell: 3, 2 expected.');
end

function case_getsol_passthrough()
prog = solved_c1c2();
x = polynomial({'s'});
assert(isequal(piesos_getsol(prog,4),4),'A double must be returned unchanged.');
assert(isempty(piesos_getsol(prog,[])),'An empty must be returned unchanged.');
p0 = x^2+1;
q = piesos_getsol(prog,p0);
assert(isa(q,'polynomial') && isequal(q.degmat,p0.degmat) && isequal(q.coeff,p0.coeff),'A polynomial must be returned unchanged.');
end

function case_getsol_operators()
dpvar c1 c2;
prog = solved_c1c2();
% intop with a decision variable in a constant kernel: the solved kernel is 2*c1+c2=8
K = intop(2*c1+c2);
Ks = piesos_getsol(prog,K);
assert(isa(Ks,'intop'),'The solution of an intop must be an intop.');
ker = kernels(Ks);
assert(numel(ker)==1 && abs(double(ker{1})-8)<1e-6,'Kernel of the solved functional: 8 expected.');
% the same for a tensopmat of intops
T = tensopmat({intop(c1-3),intop(c2+1)});
Ts = piesos_getsol(prog,T);
assert(isa(Ts,'tensopmat') && isequal(size(Ts),size(T)),'The solution of a tensopmat must be a tensopmat of the same size.');
k1 = kernels(Ts.ops{1});    k2 = kernels(Ts.ops{2});
assert(abs(double(k1{1}))<1e-6 && abs(double(k2{1})-3)<1e-6,'Solved blocks: 0 and 3 expected.');
% a polyopvar with a dpvar coefficient: only the coefficients are solved
x1 = state1();
V = x1*x1;
C = V.C;
C.ops = {c1+c2};
V.C = C;
V.chkval = true;
Vs = piesos_getsol(prog,V);
assert(isa(Vs,'polyopvar') && isequal(full(Vs.Z),full(V.Z)),'The monomials of a polyopvar must be kept.');
assert(abs(double(Vs.C.ops{1})-5)<1e-6,'Solved coefficient: 5 expected.');
end

function case_getsol_errors()
dpvar c1;
prog = solved_c1c2();
assert_err(@() piesos_getsol(3,c1),'struct');
assert_err(@() piesos_getsol(piesos_program(state1(),c1),c1),'has not been solved');
assert_err(@() piesos_getsol(prog,struct()),'should be specified as object');
assert_err(@() piesos_getsol(prog,ndopvar()),'ndopvar');
end

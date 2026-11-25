function Test_get_ineq_opts()
% TEST_GET_INEQ_OPTS tests, using assertions, `get_ineq_opts` and the
% option otype='auto' of `piesos_ineq`:
%   [Z,opts,kst] = get_ineq_opts(f,opts)
% returns the basis Z and the options of piesos_sosvar for f>=0, with
% opts.otype='auto' replaced by the smallest type of output of Zop (1
% pointwise, 2 tensor product, 3 both) whose terms have the structure of
% the terms of f.
%
% Ground truths are independent of the function under test:
%  - structure (case_structure, case_two_variables): the delta groups of
%    each functional written out by hand (e.g. int x^4: the four factors at
%    one point; (int x^2)^2: two pairs at two points; terms of u=T*x for the
%    PIE state x=u_ss: all factors integrated), and from them the type by
%    the rule of the NOTES of get_ineq_opts (one group: 1; two groups of
%    two: 2; also a group of three or more: 3);
%  - option 'auto' (case_auto): piesos_ineq with otype='auto' declares the
%    same program, Gram matrix, Zop and basis as with the numeric type of
%    the hand-derived structure;
%  - feasibility (case_feasibility): f certified with 'auto' where the
%    certificate is known (int x^4: Q=1 on x^2 pointwise; (int x^2)^2: P=1
%    on x(s_1)x(s_2); (int x^2)^2+0.1 int x^4: both), not certified for
%    -(int x^2)^2 < 0;
%  - the warning for a partial tie (case_partial_tie), (int x^4)(int x^2);
%  - the options that are not overwritten and the basis (case_options): a
%    specified otype, degs, exclude is kept; the basis of a homogeneous f
%    of degree 4 is x^2; piesos_ineq uses this basis; non-polyopvar input
%    is rejected.
% The programs of case_feasibility are solved by SeDuMi.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_get_ineq_opts
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
% DJ (CC), 11/20/2026: Initial coding;

tests = {@case_structure, @case_two_variables, @case_auto, ...
         @case_feasibility, @case_partial_tie, @case_options};
nfail = 0;
for i=1:numel(tests)
    tstart = tic;
    try
        tests{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(tests{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_get_ineq_opts: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_get_ineq_opts: all %d cases passed.\n',numel(tests));
end


%% ---- cases ------------------------------------------------------------

function case_structure()
% Scalar and vector states in one variable: groups by hand, type by rule
F = functionals();
% name, expected number of groups (largest per term), largest group, type
specs = {'q',1,2,1;  'x4',1,4,1;  'qq',2,2,2;  'qq_x4',2,4,3;  ...
         'yy_y4',2,4,3;  'p1q',1,2,1};
for c=1:size(specs,1)
    [ot,kst] = auto_otype(F.(specs{c,1}));
    assert(max(kst.ngrp)==specs{c,2},'%s: %d groups, expected %d.',specs{c,1},max(kst.ngrp),specs{c,2});
    assert(max(kst.gmax)==specs{c,3},'%s: largest group %d, expected %d.',specs{c,1},max(kst.gmax),specs{c,3});
    assert(ot==specs{c,4},'%s: otype %d, expected %d.',specs{c,1},ot,specs{c,4});
    assert(~kst.hastie,'%s: no partial tie expected.',specs{c,1});
end
% Energy of DEMO3 in the PIE state x=u_ss: u=T*x has no multiplier term,
% so all factors are integrated: no group, otype 1 (certified by otype 1)
[ot,kst] = auto_otype(demo3_energy());
assert(all(kst.ngrp==0) && all(kst.gmax==0),'DEMO3 energy: no delta group expected.');
assert(ot==1,'DEMO3 energy: otype %d, expected 1.',ot);
end

function case_two_variables()
% x(s1,s2) on [0,1]^2: int int x^4 ties the four factors in s1 and in s2
% (one group of four each); (int int x^2)^2 has two pairs in each variable
x = polyopvar('x',polynomial({'s1','s2'}),[0,1;0,1]);
q = innerprod(x,x);     x2 = x.*x;      x4 = innerprod(x2,x2);
[ot,kst] = auto_otype(x4);
assert(isequal(kst.ngrp,[1,1]) && isequal(kst.gmax,[4,4]),'2D int x^4: groups [1,1], sizes [4,4] expected.');
assert(ot==1,'2D int x^4: otype %d, expected 1.',ot);
[ot,kst] = auto_otype(q*q);
assert(isequal(kst.ngrp,[2,2]) && isequal(kst.gmax,[2,2]),'2D (int x^2)^2: groups [2,2], sizes [2,2] expected.');
assert(ot==2,'2D (int x^2)^2: otype %d, expected 2.',ot);
end

function case_auto()
% otype='auto' declares what the numeric type of the structure declares
F = functionals();
specs = {'x4',1; 'qq',2; 'qq_x4',3};
for c=1:size(specs,1)
    f = F.(specs{c,1});
    prog0 = piesos_program(F.x);
    [pA,gA,QA,ZopA,ZA] = piesos_ineq(prog0,f,struct('otype','auto'));
    [pN,gN,QN,ZopN,ZN] = piesos_ineq(prog0,f,struct('otype',specs{c,2}));
    assert(isequal(pA,pN) && isequal(gA,gN) && isequal(QA,QN) && isequal(ZopA,ZopN) && isequal(ZA,ZN), ...
        '%s: otype ''auto'' differs from otype %d.',specs{c,1},specs{c,2});
end
end

function case_feasibility()
% Known certificates, found with 'auto'; a negative functional is not
F = functionals();
ao = struct('otype','auto');
for nm = {'x4','qq','qq_x4'}
    assert(is_feasible(piesos_ineq(piesos_program(F.x),F.(nm{1}),ao)),'%s must be certified with otype ''auto''.',nm{1});
end
assert(~is_feasible(piesos_ineq(piesos_program(F.x),-1*F.qq,ao)),'-(int x^2)^2 must not be certified.');
end

function case_partial_tie()
% (int x^4)(int x^2): four factors at one point and two at another, which
% only a region tying two of three copies represents: warning, otype 3
F = functionals();
f = F.x4*F.q;
lastwarn('');
evalc('[ot,kst] = auto_otype(f);');
[~,wid] = lastwarn();
assert(strcmp(wid,'PIESOS:ineq:partialTie'),'(int x^4)(int x^2): partial-tie warning expected, got "%s".',wid);
assert(kst.hastie && ot==3,'(int x^4)(int x^2): hastie and otype 3 expected.');
end

function case_options()
% Specified options are kept; the basis is that of piesos_ineq
F = functionals();
f = F.qq_x4;
[Z,o] = get_ineq_opts(f,struct());
assert(~isfield(o,'otype'),'otype must not be set unless ''auto''.');
assert(isfield(o,'degs') && isfield(o,'exclude'),'degs and exclude must be set.');
% f is homogeneous of degree 4 in x: the basis is x^2 (half the degree)
assert(isequal(full(Z.Z),[1,1]),'Basis of (int x^2)^2 + 0.1 int x^4 must be x^2.');
spec = struct('otype',1,'degs',struct('in',[2;2],'out',1),'exclude',[false;true;true]);
[~,o2] = get_ineq_opts(f,spec);
assert(isequal(o2.otype,1) && isequal(o2.degs,spec.degs) && isequal(o2.exclude,spec.exclude), ...
    'Specified otype, degs, exclude must be kept.');
[~,~,~,~,Zi] = piesos_ineq(piesos_program(F.x),f,struct('otype',1));
assert(isequal(Zi,Z),'piesos_ineq must use the basis of get_ineq_opts.');
assert_err(@() get_ineq_opts('a'),'polyopvar');
end


%% ---- helpers ----------------------------------------------------------

function [ot,kst] = auto_otype(f)
% The type chosen for otype='auto', and the statistics of the kernels
[~,o,kst] = get_ineq_opts(f,struct('otype','auto'));
ot = o.otype;
end

function F = functionals()
% Functionals of x(s) on [0,1] and y(s) in R^2, with their delta groups:
%   q = int x^2 (one pair), x4 = int x^4 (four at one point),
%   qq = (int x^2)^2 (two pairs), qq_x4 = qq + 0.1 x4,
%   yy_y4 = (int |y|^2)^2 + 0.1 int |y|^4 (as qq_x4),
%   p1q = p1*<x,x>, p1 a quadratic functional of the lower and upper
%   integral rows of Zop (no group, all factors integrated) times one
%   pair: one group (certified by otype 3 only for some p1, see
%   Test_piesos_tensor_gram, case_certification)
pvar s
x = polyopvar('x',s,[0,1]);
F.x = x;
F.q = innerprod(x,x);
x2 = x.*x;
F.x4 = innerprod(x2,x2);
F.qq = F.q*F.q;
F.qq_x4 = F.qq + 0.1*F.x4;
y = polyopvar({'y'},2,{'s'},[0,1],true(1,1));
yy = otimes(y,y,2);     y2 = y.*y;
F.yy_y4 = innerprod(yy,yy) + 0.1*innerprod(y2,ones(2),y2);
prog = piesos_program(x);
[~,p1] = piesos_sosvar(prog,dmonomials(x,1),struct('degs',0,'exclude',[1;0;0],'otype',1));
F.p1q = p1*F.q;
end

function E = demo3_energy()
% Energy of PIESOS_DEMO3_energy_bound (Allen-Cahn, r=12), in the PIE state
clear stateNameGenerator
pvar t s
dom = [0,1];    r = 12;
u = pde_var(s,dom);
PDE = [diff(u,t)==diff(u,s,2)+r*u-u^3; subs(u,s,dom(1))==0; subs(u,s,dom(2))==0];
evalc('PIE = convert(PDE);');
x = PIE.x;  u = PIE.T*x;    u2 = u.^2;
E = -(u'*x) - r*(u'*u) + 0.5*(u2'*u2);
end

function [ok,prog_sol] = is_feasible(prog)
% Solve with SeDuMi; feasible if neither pinf nor dinf is reported. A
% program without PSD variable fails in SeDuMi: sosposmatr adds one
[prog,~] = sosposmatr(prog,1);
evalc('prog_sol = piesos_solve(prog,struct(''solver'',''sedumi'',''simplify'',true));');
info = prog_sol.solinfo.info;
ok = ~info.pinf && ~info.dinf;
end

function assert_err(fn,pattern)
% fn() must error with a message containing pattern
try
    fn();
catch err
    assert(contains(err.message,pattern),'Unexpected error message: "%s" (expected "%s").',err.message,pattern);
    return
end
error('No error thrown (expected "%s").',pattern);
end

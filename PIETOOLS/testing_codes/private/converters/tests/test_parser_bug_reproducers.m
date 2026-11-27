function results = test_parser_bug_reproducers()
% TEST_PARSER_BUG_REPRODUCERS Regression tests of the parser and the PIE
% structures, each checked against an independently declared reference
% rather than against the routine itself:
%   1. an integral over the variable that is the limit of an inner
%      integral;
%   2. reordering outputs that hold a static constraint;
%   3. expansion of higher-order temporal derivatives of a 'pie_struct'
%      with a map pde2pie.x.C with blocks outside the diagonal;
%   4. terms added to an equation;
%   5. a matrix times a nonlinear term;
%   6. a delayed boundary trace of a 2D state, fixed and distributed
%      delays;
%   7. concatenation of a PIE with an input on its own variable, and the
%      PIE built from its operators;
%   8. == of an 'npie_struct' and of the 'npie_struct' of a 'pie_struct'.
%
% Run as: results = test_parser_bug_reproducers;
%
% OUTPUT
% - results: 'struct' array with fields name, passed and message, one
%            element per test.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_parser_bug_reproducers
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
% Initial coding DJ - 11/23/2026

tests = {@case1_nested_integral, @case2_reorder_constraint, @case3_offdiagonal_map, ...
         @case4_terms_plus_equation, @case5_nonlinear_coefficient, @case6_delayed_trace, ...
         @case7_concatenation, @case8_npie_eq};
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



function PIE = cv(P)
% PIE = CV(P) converts the PDE P to a PIE without printing.
% INPUT
% - P:      'pde_struct', the PDE;
% OUTPUT
% - PIE:    'pie_struct' or 'npie_struct', the PIE.
evalc('PIE = convert_PIETOOLS_PDE(P,''silent'');');
end

function case1_nested_integral()
% int(int(x,s,[0,s]),s,[0,1]) = int_0^1 (1-th) x(th) dth, as a term, in
% the equation of a PDE state y(s), and in the equation of an ODE state.
pvar s t
clear('stateNameGenerator');
x = pde_var(s,[0,1]);   y = pde_var(s,[0,1]);   X = pde_var(1);
T = int(int(x,s,[0,s]),s,[0,1]);
tm = T.free{1}.term{1};
assert(isempty(T.free{1}.vars),'The result should not depend on s.');
assert(isequal(double(polynomial(tm.I{1})),[0,1]),'The result should be a full integral.');
assert(abs(double(subs(polynomial(tm.C),tm.loc(1),0.3))-0.7)<1e-12,'The kernel should be 1-th.');
base = [diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0; subs(y,s,0)==0; subs(y,s,1)==0];
assert(cv([diff(y,t)==diff(y,s,2)+T; base])==cv([diff(y,t)==diff(y,s,2)+int((1-s)*x,s,[0,1]); base]),...
        'In the equation of y(s), the PIE differs from that of int((1-s)*x,s,[0,1]).');
base = [diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0];
assert(cv([diff(X,t)==-X+T; base])==cv([diff(X,t)==-X+int((1-s)*x,s,[0,1]); base]),...
        'In an ODE equation, the PIE differs from that of int((1-s)*x,s,[0,1]).');
end

function case2_reorder_constraint()
% Reordering the outputs keeps the static constraint 0 = x1_s + x2
% flagged, also when it moves, and when it is the only output.
pvar s1 t
% The reproducer: the constraint keeps its position
clear('stateNameGenerator');
x1 = pde_var(s1,[0,1]); x2 = pde_var(s1,[0,1]); z = pde_var('output',1);
P = initialize([diff(x1,t)==diff(x2,s1); diff(x2,t)==diff(x1,s1); z==subs(x1,s1,1);
                subs(x1,s1,0)==0; subs(x2,s1,1)==0; diff(x1,s1)+x2==0],true);
PIE0 = cv(P);
P2 = reorder_comps(P,'z',true);
assert(isequal(P2.constraint_ID,2),'The constraint should have ID 2 after reordering.');
PIE2 = cv(P2);
assert(PIE2==PIE0 && size(PIE2.K,1)==1 && numel(PIE2.space.z.dims)==1,'Reordering the outputs changed the PIE.');
assert(cv(reorder_comps(P,'all',true))==PIE0,'Reordering all components changed the PIE.');
% A distributed output declared first: the constraint moves to the end
clear('stateNameGenerator');
x1 = pde_var(s1,[0,1]); x2 = pde_var(s1,[0,1]);
z1 = pde_var('output',1,s1,[0,1]); z2 = pde_var('output',1);
P = initialize([diff(x1,t)==diff(x2,s1); diff(x2,t)==diff(x1,s1); z1==x2; z2==subs(x1,s1,1);
                subs(x1,s1,0)==0; subs(x2,s1,1)==0; diff(x1,s1)+x2==0],true);
PIE0 = cv(P);
P2 = reorder_comps(P,'z',true);
is_con_terms = cellfun(@(zi) numel(zi.term)==2,P2.z);
assert(isequal(find(is_con_terms)',P2.constraint_ID),'constraint_ID does not point to the constraint.');
PIE2 = cv(P2);
assert(all(logical(PIE2.K==PIE0.K)) && all(logical(PIE2.A==PIE0.A)) && numel(PIE2.space.z.dims)==2,...
        'Reordering the outputs changed the constraint or the dynamics.');
% The constraint as the only output
clear('stateNameGenerator');
x1 = pde_var(s1,[0,1]); x2 = pde_var(s1,[0,1]);
Q = initialize([diff(x1,t)==diff(x2,s1); diff(x2,t)==diff(x1,s1); subs(x1,s1,0)==0; subs(x2,s1,1)==0; diff(x1,s1)+x2==0],true);
Q2 = reorder_comps(Q,'z',true);
assert(isequal(Q2.constraint_ID,1) && cv(Q2)==cv(Q),'A single constraint output is not kept.');
end

function case3_offdiagonal_map()
% Expansion of higher-order temporal derivatives of a 'pie_struct' rejects
% a map pde2pie.x.C with blocks
% outside the diagonal, and still expands a wave equation.
pvar s t
clear('stateNameGenerator');
x1 = pde_var(s,[0,1]); x2 = pde_var(s,[0,1]);
P = initialize(cv([diff(x1,t)==diff(x1,s,2)+x2; diff(x2,t)==diff(x2,s,2);
                   subs(x1,s,0)==0; subs(x1,s,1)==0; subs(x2,s,0)==0; subs(x2,s,1)==0]));
C = P.pde2pie.x.C;  Cc = C.C;  Cc{1,2} = Cc{1,1};
P.pde2pie.x.C = opvar(Cc,struct('vars',{C.vars},'dom',C.dom,'space_out',C.space_out,...
                    'space_in',C.space_in,'dim_out',C.dim_out,'dim_in',C.dim_in));
try
    expand_tderivatives(P,[2;1],true);
    error('test:noError','An off-diagonal map pde2pie.x.C was expanded without error.');
catch ME
    assert(strcmp(ME.identifier,'expand_tderivatives:mapOffDiagonal'),ME.message);
end
clear('stateNameGenerator');
x = pde_var(s,[0,1]);
PW = cv([diff(x,t,2)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0]);
assert(numel(PW.space.x.dims)==2,'The wave equation should have 2 PIE state components.');
end

function case4_terms_plus_equation()
% s*x + (x_t==x) adds s*x to the left-hand side, as s*x + x_t == x does;
% (x_t==x) + s*x adds it to the right-hand side.
pvar s t
clear('stateNameGenerator');
x = pde_var(s,[0,1]);
assert(cv(s*x+(diff(x,t)==x))==cv(s*x+diff(x,t)==x),'s*x+(x_t==x) differs from s*x+x_t==x.');
assert(~(cv(s*x+(diff(x,t)==x))==cv((diff(x,t)==x)+s*x)),'s*x+(x_t==x) should differ from (x_t==x)+s*x.');
end

function case5_nonlinear_coefficient()
% A matrix times a nonlinear term, [1;s]*(x*x) and [1 s]*(x2.*x2), equals
% the elementwise products of linear terms, compared by evaluation of the
% right-hand sides f of the PIEs.
pvar s t
clear('stateNameGenerator');
x = pde_var(s,[0,1]);   y = pde_var(2,s,[0,1]);
base = [diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0; subs(y,s,0)==0; subs(y,s,1)==0];
N1 = cv([diff(y,t)==diff(y,s,2)+[1;s]*(x*x); base]);
N2 = cv([diff(y,t)==diff(y,s,2)+([1;s]*x).*([1;1]*x); base]);
assert(f_equal(N1,N2),'[1;s]*(x*x) differs from ([1;s]*x).*([1;1]*x).');
clear('stateNameGenerator');
x2 = pde_var(2,s,[0,1]);   y = pde_var(s,[0,1]);
base = [diff(x2,t)==diff(x2,s,2); subs(x2,s,0)==0; subs(x2,s,1)==0; subs(y,s,0)==0; subs(y,s,1)==0];
N1 = cv([diff(y,t)==diff(y,s,2)+[1 s]*(x2.*x2); base]);
N2 = cv([diff(y,t)==diff(y,s,2)+([1 0]*x2).*([1 0]*x2)+([0 s]*x2).*([0 1]*x2); base]);
assert(f_equal(N1,N2),'[1 s]*(x2.*x2) differs from the sum of elementwise products.');
end

function e = f_equal(N1,N2)
% E = F_EQUAL(N1,N2) compares the right-hand sides f of two nonlinear PIEs
% with the same state, by evaluation at random quadratic values of the
% fundamental state.
% INPUT
% - N1, N2: 'npie_struct' objects;
% OUTPUT
% - e:      logical, true if the values agree to 1e-8.
v = polynomial(N1.space.vars(1));
e = true;
for trial=1:3
    xvals = cell(numel(N1.f.varname),1);
    for i=1:numel(xvals)
        c = randn(3,1);
        xvals{i} = (c(1)+c(2)*v+c(3)*v^2)*ones(N1.f.dims.in(i),1);
    end
    Y = polynomial(apply_poly(N1.f,xvals) - apply_poly(N2.f,xvals));
    e = e && (isempty(Y.coefficient) || full(max(abs(Y.coefficient(:))))<1e-8);
end
end

function case6_delayed_trace()
% Delayed terms, fixed and distributed, of 1D states and of boundary traces
% of 2D states, with and without integral: the PDE converts, and the term
% after expansion of the delays is the history state at the delay, integrated over
% the full domain of the other variable(s) where the PDE integrates.
for distributed = [false true]
    for k = 1:5
        if distributed && ismember(k,[2 5]), continue; end
        clear('stateNameGenerator');
        [PDE,iX,lbl] = delay_pde(k,distributed);
        cv(PDE);
        evalc('PDEd = expand_delays(initialize(PDE,true),true);');
        if numel(PDEd.x{iX}.term)>=2
            tm = PDEd.x{iX}.term{2};
        else
            iy = find(cellfun(@(xi) size(xi.vars,1)==1,PDEd.x),1);
            tm = PDEd.x{iy}.term{2};
        end
        dom = PDEd.x{tm.x}.dom;
        is_int = ~cellfun(@isempty,tm.I);
        if distributed
            assert(all(is_int),'%s (distributed): all variables should be integrated.',lbl);
        else
            is_const = arrayfun(@(j) isdouble(tm.loc(j)),1:size(dom,1));
            assert(sum(is_const)==1 && double(tm.loc(is_const))==dom(is_const,1) && ~is_int(is_const),...
                    '%s: the history state should be evaluated at the delay.',lbl);
            int_full = arrayfun(@(j) isequal(double(polynomial(tm.I{j})),dom(j,:)),find(is_int));
            assert(all(int_full) && (k==5 || all(is_int(~is_const))) && (k~=5 || ~any(is_int)),...
                    '%s: wrong integrals of the history state.',lbl);
        end
    end
end
end

function [PDE,iX,lbl] = delay_pde(k,distributed)
% [PDE,IX,LBL] = DELAY_PDE(K,DISTRIBUTED) declares delay test system K.
% INPUT
% - k:      'double', 1: int x(t-1,s) ds; 2: int s*x(t-1,s) ds;
%           3: int x(t-1,s1,0) ds1; 4: int x(t-1,0,s2) ds2; 5: x(t-1,0,s2)
%           in the equation of y(s2);
% - distributed:    logical, true for a distributed delay over [-1,0];
% OUTPUT
% - PDE:    'pde_struct', the system;
% - iX:     'double', the index of the ODE state X;
% - lbl:    char, a label of the system.
pvar s t s1 s2 tau
lbls = {'int x(t-1,s) ds','int s*x(t-1,s) ds','int x(t-1,s1,0) ds1','int x(t-1,0,s2) ds2','x(t-1,0,s2) in y(s2)'};
lbl = lbls{k};
dly = @(z) subs(z,t,t-1);
if distributed, dly = @(z) z; end
if k<=2
    x = pde_var(s,[0,1]); X = pde_var(1);
    if k==1, r = int(dly(x),s,[0,1]); else, r = int(s*dly(x),s,[0,1]); end
    PDE = [diff(X,t)==-X+r; diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0];
else
    % Dirichlet at s1=0 and s2=0: a delayed trace at a boundary where the
    % state is not fixed by a BC does not convert
    x = pde_var([s1;s2],[0,1;0,1]); X = pde_var(1);
    bcs = [subs(x,s1,0)==0; subs(diff(x,s1),s1,1)==0; subs(x,s2,0)==0; subs(diff(x,s2),s2,1)==0];
    if k==3
        PDE = [diff(X,t)==-X+int(dly(subs(x,s2,0)),s1,[0,1]); diff(x,t)==diff(x,s1,2)+diff(x,s2,2); bcs];
    elseif k==4
        PDE = [diff(X,t)==-X+int(dly(subs(x,s1,0)),s2,[0,1]); diff(x,t)==diff(x,s1,2)+diff(x,s2,2); bcs];
    else
        y = pde_var(s2,[0,1]);
        PDE = [diff(X,t)==-X; diff(y,t)==diff(y,s2,2)+dly(subs(x,s1,0)); diff(x,t)==diff(x,s1,2)+diff(x,s2,2);
               bcs; subs(y,s2,0)==0; subs(y,s2,1)==0];
    end
end
iX = find(cellfun(@(xi) isempty(xi.vars),PDE.x),1);
if distributed
    PDE.tau = [tau,1];
    PDE.x{iX}.term{2}.delay = tau;
end
end

function case7_concatenation()
% [PIE,PIE], [PIE;PIE] and the PIE built from the operators, for an ODE
% with an input w(s): the input
% depends on a variable on which no state depends.
pvar s t
clear('stateNameGenerator');
x = pde_var('state',1,[],[]); w = pde_var('in',1,s,[0,2]); z = pde_var('out',1);
PIE = cv([diff(x,t)==-x+int(w,s,[0,2]); z==x]);
Ph = initialize([PIE,PIE]);     Pv = initialize([PIE;PIE]);
assert(isequal(Ph.space.vars,PIE.space.vars) && isequal(size(Ph.Bw,'block'),[2 2]),'horzcat of the PIE is wrong.');
assert(isequal(Pv.space.vars,PIE.space.vars) && isequal(size(Pv.Bw,'block'),[2 1]),'vertcat of the PIE is wrong.');
assert(all(logical(Ph.Bw.C{1,1}==PIE.Bw.C{1,1})),'The (1,1) block of horzcat Bw differs from Bw.');
assert(piess(PIE.T,PIE.A,PIE.Bw,PIE.Cz,PIE.Dzw)==PIE,'piess of the operators differs from the PIE.');
end

function case8_npie_eq()
% The 'npie_struct' of the 'pie_struct' of a linear PDE equals the
% 'npie_struct' converted directly; a
% different PIE is not equal.
pvar s t
clear('stateNameGenerator');
x = pde_var(s,[0,1]); w = pde_var('in',1);
PDE = [diff(x,t)==diff(x,s,2)+w; subs(x,s,0)==0; subs(x,s,1)==0];
N = convert_PIETOOLS_PDE(PDE,{'silent','nonlinear'});
assert(npie_struct(cv(PDE))==N,'npie_struct(P) differs from N.');
PDE2 = [diff(x,t)==diff(x,s,2)+2*w; subs(x,s,0)==0; subs(x,s,1)==0];
assert(~(npie_struct(cv(PDE))==convert_PIETOOLS_PDE(PDE2,{'silent','nonlinear'})),'Different PIEs compare equal.');
end

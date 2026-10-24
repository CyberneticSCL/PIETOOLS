function Test_piesos_LF(only)
% TEST_PIESOS_LF(ONLY) tests, using assertions, piesos_LF,
%   [prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,xvars,degmat,opts),
% which declares a Lyapunov functional V(x) in the PIE state x, and
% computes its Lie derivative dV(x) along the PIE d/dt (T*x) = f(x), for
% u=T*x the PDE state,
%   (PDE states)  V(x) = <Zop*Z(u), Pmat*Zop*Z(u)>,
%   (PIE states)  V(x) = <Zop*Z(x), Pmat*Zop*Z(u)>, with the symmetry
%                 condition V(xL,T*xR) = V(xR,T*xL) imposed on Pmat.
%
% Oracle, independent of the function under test. The PIEs are converted
% from PDEs with Dirichlet conditions, so a polynomial PIE state x=u''
% (u(0)=u(1)=u''(0)=u''(1)=0) has the velocity xdot=(f(x))'' with
% T*xdot=f(x):
%  - V(x0) and dV(x0) are compared with direct polynomial arithmetic
%    (apply_op and integration, no composition of functionals):
%      PDE:  V = <Zop*Z(u0),P Zop*Z(u0)>,  dV = d/dh V(u0+h*f(x0)) at h=0,
%      PIE:  V = <Zop*Z(x0),P Zop*Z(T*x0)>, dV = d/dh V(x0+h*xdot0) at h=0,
%    with random decision variables for PDE states (V and dV are linear
%    in them), and for PIE states a feasible Pmat (symmetry condition)
%    solved with SeDuMi, with a random normalization. Covers the product
%    rule, doubling only the derivative of the right factors (the exact dV
%    varies both sides), the composition with T, and the symmetry
%    condition itself (a Pmat violating it gives another dV);
%  - states of different dependence, several states and subsets, sums of
%    functionals in PDE and PIE states, input formats, errors;
%  - PIEs with nonlinear terms (nested tensor-PI operators on substitution
%    into V, see piesos_LF): Burgers, cubic, and coupled nonlinear PDEs;
%  - opts.derivative=false and the nargout rule (case_skip_derivative):
%    program, V, Pmat, Zop equal those of the full call, and the calls of
%    subs of the derivative are gone, counted by the profiler;
%  - an ODE state coupled to a PDE state (case_ode_pde), compared with
%    polynomial arithmetic;
%  - vector-valued states (case_vector_states): V and dV equal those of
%    the equivalent scalar states with the same decision variables, also
%    for subsets of the states (rows of T and f selected by state).
%
% INPUTS
% - only:   (optional) cell of character arrays with the names of the
%           cases to run, e.g. Test_piesos_LF({'case_errors'}). Defaults
%           to running all cases;
%
% NOTES
% - Degree-2 cases use operators that only multiply or have an upper
%   integral;
% - check_mixed: the symmetry condition is a redundant equality system, so
%   SeDuMi needs a conic variable and more variables than equations
%   ("Size b mismatch" otherwise);
% - eval_V names dummy variables s_k_j as apply_op does; they are matched
%   by name;
% - fix_empty_params_* work around apply_sopvar failing on empty params
%   (as in Test_piesos_quadvar.m);
% - case_skip_derivative counts calls of polyopvar.subs with the profiler.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_LF
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
% DJ (CC), 10/21/2026: Add cases for an ODE state and vector-valued states;
% DJ (CC), 10/22/2026: apply_side, check_pde_blocks: Zop acts on the
%                        entries of the basis (rows of Z.C), not on its
%                        monomials;
% DJ (CC), 10/23/2026: check_pde_blocks: one block of Zop per state (vector
%                        variable of dmonomials), not per scalar entry;

tests ={@case_pde_options, @case_pde_degrees, @case_pde_coupled, ...
         @case_mixed, @case_mixed_coupled, @case_combined, ...
         @case_input_forms, @case_errors, @case_nonlinear, ...
         @case_skip_derivative, @case_ode_pde, @case_vector_states};
if nargin>0
    tests = tests(ismember(cellfun(@func2str,tests,'UniformOutput',false),only));
end
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
    error('Test_piesos_LF: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_LF: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function PIE = get_pie(kind)
% PIE of the nonlinear converter, for a PDE with Dirichlet conditions on [0,1]
persistent cache
if isempty(cache)
    cache = struct();
end
if isfield(cache,kind)
    PIE = cache.(kind);
    return
end
clear stateNameGenerator
pvar s t
dom = [0,1];
switch kind
    case 'heat'         % u_t = u_ss + 2*u
        x = pde_var(s,dom);
        PDE = [diff(x,t)==diff(x,s,2)+2*x;  subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];
    case 'coupled'      % u1_t = u1_ss + 2*u1 + 0.5*u2,  u2_t = u2_ss - u2 + 0.3*u1
        x1 = pde_var(s,dom);
        x2 = pde_var(s,dom);
        PDE = [diff(x1,t)==diff(x1,s,2)+2*x1+0.5*x2;
               diff(x2,t)==diff(x2,s,2)-x2+0.3*x1;
               subs(x1,s,dom(1))==0;  subs(x1,s,dom(2))==0;
               subs(x2,s,dom(1))==0;  subs(x2,s,dom(2))==0];
    case 'burgers'      % u_t = u_ss + 2*u - u*u_s
        x = pde_var(s,dom);
        PDE = [diff(x,t)==diff(x,s,2)+2*x-x*diff(x,s);  subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];
    case 'fisher'       % u_t = u_ss + 2*u - u^3
        x = pde_var(s,dom);
        PDE = [diff(x,t)==diff(x,s,2)+2*x-x*x*x;  subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];
    case 'coupled_nl'   % u1_t = u1_ss + 2*u1 + u1*u2,  u2_t = u2_ss - u2 - u1*u1_s
        x1 = pde_var(s,dom);
        x2 = pde_var(s,dom);
        PDE = [diff(x1,t)==diff(x1,s,2)+2*x1+x1*x2;
               diff(x2,t)==diff(x2,s,2)-x2-x1*diff(x1,s);
               subs(x1,s,dom(1))==0;  subs(x1,s,dom(2))==0;
               subs(x2,s,dom(1))==0;  subs(x2,s,dom(2))==0];
    case 'ode_pde'     % xo_t = xo + u_s(0),  u_t = u_ss,  u(0) = xo,  u(1) = -2*xo
        xo = pde_var();
        x = pde_var(s,dom);
        PDE = [diff(xo,t)==xo+subs(diff(x,s),s,dom(1));  diff(x,t)==diff(x,s,2);
               subs(x,s,dom(1))==xo;  subs(x,s,dom(2))==-2*xo];
    case 'sc3'          % three scalar states, u_t = u_ss + A3*u, Dirichlet
        x1 = pde_var(s,dom);
        x2 = pde_var(s,dom);
        x3 = pde_var(s,dom);
        A3 = [2 .5 .2; .3 -1 .4; .1 .6 -2];
        PDE = [diff(x1,t)==diff(x1,s,2)+A3(1,1)*x1+A3(1,2)*x2+A3(1,3)*x3;
               diff(x2,t)==diff(x2,s,2)+A3(2,1)*x1+A3(2,2)*x2+A3(2,3)*x3;
               diff(x3,t)==diff(x3,s,2)+A3(3,1)*x1+A3(3,2)*x2+A3(3,3)*x3;
               subs(x1,s,dom(1))==0;  subs(x1,s,dom(2))==0;
               subs(x2,s,dom(1))==0;  subs(x2,s,dom(2))==0;
               subs(x3,s,dom(1))==0;  subs(x3,s,dom(2))==0];
    case 'vec'          % vector state with 2 components: the same PDE as 'coupled'
        x = pde_var(2,s,dom);
        PDE = [diff(x,t)==diff(x,s,2)+[2 .5; .3 -1]*x;  subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];
    case 'vec_sc'       % vector state with 2 components, and a scalar state: the same PDE as 'sc3'
        x = pde_var(2,s,dom);
        y = pde_var(s,dom);
        A3 = [2 .5 .2; .3 -1 .4; .1 .6 -2];
        PDE = [diff(x,t)==diff(x,s,2)+A3(1:2,1:2)*x+A3(1:2,3)*y;
               diff(y,t)==diff(y,s,2)+A3(3,1:2)*x+A3(3,3)*y;
               subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0;
               subs(y,s,dom(1))==0;  subs(y,s,dom(2))==0];
    otherwise
        error('Unknown PIE.');
end
evalc('PIE = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
cache.(kind) = PIE;
end

function X0 = states_x0(n)
% Polynomial PIE states x_i with x_i(0)=x_i(1)=0 (so that f(x0) satisfies the
% boundary conditions, and xdot0 = (f(x0))'' is the state velocity)
sv = polynomial({'s'});     sv = sv(1);
X0 = {-12*sv+12*sv^2, 3*sv-3*sv^2, sv-sv^2};
X0 = X0(1:n);
end

function [vals,Tx,fx,u0,w0,xdot0] = pie_values(PIE,X0)
% PDE state u0=T*x0, its velocity w0=f(x0), and the PIE state velocity
% xdot0=(w0)'' (T*xdot0=w0 for w0 satisfying the Dirichlet conditions)
sv = polynomial({'s'});     sv = sv(1);
Tx = PIE.T*PIE.x;           fx = PIE.f;
u0 = eval_polyopvar(Tx,X0);
w0 = eval_polyopvar(fx,X0);
xdot0 = cell(1,numel(X0));
for i=1:numel(X0)
    xdot0{i} = diff(w0(i),sv,2);
end
vals = struct();
for i=1:numel(X0)
    vals.(PIE.varname_pie{i}) = X0{i};
end
end

function idcs = state_idcs(PIE,names,use_pde)
% Indices in the PIE of the states with the given names, sorted by name
if use_pde
    list = PIE.varname_pde(:);
else
    list = PIE.varname_pie(:);
end
[~,idcs] = ismember(names(:),list);
[~,srt] = sort(list(idcs));
idcs = idcs(srt);
end

function Zu = monomials_in(PIE,idcs,degmat)
% The monomials Z(u) in the PDE states idcs, as used by piesos_LF
xu = polyopvar(PIE.varname_pde(idcs),PIE.x_tab(idcs,2),{'s'},PIE.dom,logical(PIE.x_tab(idcs,3)));
Zu = dmonomials(xu,degmat);
end

function res = check_pde(PIE,X0,names,degmat,opts,seed)
% Declare V, dV for the PDE states 'names' and compare with polynomial arithmetic
pvar h
sv = polynomial({'s'});     sv = sv(1);
[vals,~,~,u0,w0] = pie_values(PIE,X0);
idcs = state_idcs(PIE,names,true);
prog = piesos_program(PIE.x);
[prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,names,degmat,opts); %#ok<ASGLU>
d0map = random_decvars(Pmat,seed);
Zu = monomials_in(PIE,idcs,degmat);
Xu = struct();
for i=1:numel(idcs)
    Xu.(PIE.varname_pde{idcs(i)}) = u0(idcs(i)) + h*w0(idcs(i));
end
Ph = block_matrix(cellfun(@(P) eval_dpvar_at(P,d0map),Pmat,'UniformOutput',false));
YL = apply_side(Zop,Zu,Xu);
Fh = int(YL.'*Ph*YL,sv,0,1);
res.Vtrue = double(subs(Fh,h,0));
res.dVtrue = double(subs(diff(Fh,h),h,0));
res.V = double(eval_dpvar_at(eval_V(V,vals),d0map));
res.dV = double(eval_dpvar_at(eval_V(dV,vals),d0map));
res.errV = abs(res.V-res.Vtrue)/max(1e-12,abs(res.Vtrue));
res.errdV = abs(res.dV-res.dVtrue)/max(1e-12,abs(res.dVtrue));
end

function res = check_mixed(PIE,X0,names_pde,degmat,opts,seed)
% Declare V, dV for the PIE states (associated to the PDE states 'names_pde'), solve for a
% feasible Pmat with a random normalization, and compare with polynomial arithmetic at x0
% with velocity xdot0
pvar h
sv = polynomial({'s'});     sv = sv(1);
[vals,~,~,u0,w0,xdot0] = pie_values(PIE,X0);
idcs = state_idcs(PIE,names_pde,true);
prog = piesos_program(PIE.x);
[prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,PIE.varname_pie(idcs),degmat,opts);
% Normalization: random linear combination of the decision variables equal to 1
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)]; %#ok<AGROW>
end
names = unique(names,'stable');
rng(seed);
dsum = dpvar(sparse([0;randn(numel(names),1)]),zeros(1,0),cell(0,1),names,[1 1]);
prog = soseq(prog,dsum-1);
% SeDuMi needs a conic variable and more variables than equations (see NOTES)
ncon = sum(cellfun(@(A) size(A,2),prog.expr.At));
[prog,~] = sosposmatr(prog,ceil(sqrt(ncon))+2);
evalc('[prog,info] = sossolve(prog);');
assert(info.pinf==0 && info.dinf==0,'The symmetry condition must be feasible with a nonzero Pmat.');
Psol = piesos_getsol(prog,Pmat);
Vsol = piesos_getsol(prog,V);
dVsol = piesos_getsol(prog,dV);
Zu = monomials_in(PIE,idcs,degmat);
Xl = struct();      Xr = struct();
for i=1:numel(idcs)
    Xl.(PIE.varname_pde{idcs(i)}) = X0{idcs(i)} + h*xdot0{idcs(i)};
    Xr.(PIE.varname_pde{idcs(i)}) = u0(idcs(i)) + h*w0(idcs(i));
end
Ph = block_matrix(Psol);
YL = apply_side(Zop,Zu,Xl);
YR = apply_side(Zop,Zu,Xr);
Fh = int(YL.'*Ph*YR,sv,0,1);
res.Vtrue = double(subs(Fh,h,0));
res.dVtrue = double(subs(diff(Fh,h),h,0));
res.V = double(eval_V(Vsol,vals));
res.dV = double(eval_V(dVsol,vals));
res.errV = abs(res.V-res.Vtrue)/max(1e-12,abs(res.Vtrue));
res.errdV = abs(res.dV-res.dVtrue)/max(1e-12,abs(res.dVtrue));
res.ndvars = numel(names);
end

function assert_res(res,label)
assert(abs(res.Vtrue)>1e-8 && abs(res.dVtrue)>1e-8,'%s: the values %g, %g are too small to test.',label,res.Vtrue,res.dVtrue);
assert(res.errV<=1e-8,'%s: V(x0)=%.10g differs from %.10g.',label,res.V,res.Vtrue);
assert(res.errdV<=1e-8,'%s: dV(x0)=%.10g differs from %.10g.',label,res.dV,res.dVtrue);
end

function assert_err(f,needle)
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s" (expected "%s").',err.message,needle);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end

function d0map = random_decvars(Pmat,seed)
% Random values for the decision variables in Pmat, by name
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)]; %#ok<AGROW>
end
names = unique(names,'stable');
rng(seed);
d0map = containers.Map(names,num2cell(0.5+rand(numel(names),1)));
end

function P = block_matrix(Pc)
% Matrix of the blocks in the cell Pc
rows = cell(size(Pc,1),1);
for i=1:size(Pc,1)
    rows{i} = horzcat(Pc{i,:});
end
P = vertcat(rows{:});
end

function Y = apply_side(Zop,Z,X)
% Zop applied to the entries of Z (the states in X): block i acts on the tensor
% product of the factors of the monomial in row i of Z.C, each with its
% coefficient (the selection of a component of a vector state) applied
Zf = full(Z.Z);
nr = size(Z.C.ops,1);
Xc = cell(1,nr);
for i=1:nr
    k = find(~cellfun(@isempty,Z.C.ops(i,:)));
    assert(isscalar(k),'Row %d of the basis does not have a single monomial.',i);
    idc = Zf(k,Zf(k,:)~=0);
    if isempty(idc)
        Xc{i} = {1};
    else
        blk = Z.C.ops{i,k};
        assert(size(blk.ops,1)==1 && size(blk.ops,2)==numel(idc),'Row %d of the basis has an unexpected coefficient.',i);
        Xc{i} = arrayfun(@(m) apply_sopvar(blk.ops{1,m},X.(Z.varname{idc(m)})),1:numel(idc),'UniformOutput',false);
    end
end
Y = apply_op(fix_empty_params_zop(Zop),Xc);
end

function Y = eval_V(V,vals)
% Evaluate the functional V (polyopvar of intop kernels) at the states vals.(name),
% polynomials in their own variables. Dummy variables of s_k are numbered over the factors
% depending on s_k, as s_k_j, if the monomial has more than 1 factor.
vnames = V.vars.in;
dep_state = double(V.depmat.in)>0;
Y = 0;
Zf = full(V.Z);
for k=1:size(V.Z,1)
    op = V.C.ops{1,k};
    if isempty(op), continue, end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    if d==0
        Y = Y + apply_op(op,polynomial(1));
    else
        slot = zeros(1,numel(vnames));
        v = 1;
        for j=1:d
            xj = vals.(V.varname{idcs(j)});
            for mm=1:numel(vnames)
                if dep_state(idcs(j),mm)
                    slot(mm) = slot(mm)+1;
                    if d>1
                        isv = strcmp(xj.varname,vnames{mm});
                        xj.varname(isv) = {[vnames{mm},'_',num2str(slot(mm))]};
                    end
                end
            end
            v = kron(v,xj);
        end
        Y = Y + apply_op(op,v);
    end
end
end

%% ---- workaround for apply_sopvar (as in Test_piesos_quadvar.m) ---------
function Zop2 = fix_empty_params_zop(Zop)
% monopmat leaves params{k}=[] for excluded multiplier/lower/upper combos, which
% apply_sopvar cannot handle: expand to zero sparse matrices
Zop2 = Zop;
[p,d] = size(Zop.ops);
for i=1:p
    for j=1:d
        if ~isempty(Zop.ops{i,j})
            Zop2.ops{i,j} = fix_empty_params_tensopvar(Zop.ops{i,j});
        end
    end
end
end

function T2 = fix_empty_params_tensopvar(T)
T2 = T;
[nt,nf] = size(T.ops);
for i=1:nt
    for j=1:nf
        entry = T.ops{i,j};
        if isa(entry,'sopvar')
            T2.ops{i,j} = fill_empty_params_sopvar(entry);
        elseif isa(entry,'tensopvar')
            T2.ops{i,j} = fix_empty_params_tensopvar(entry);
        end
    end
end
end

function Pop2 = fill_empty_params_sopvar(Pop)
params = Pop.params;
nL = prod(cellfun(@numel,Pop.ZL));
nR = prod(cellfun(@numel,Pop.ZR));
nrow = Pop.dims(1)*nL;      ncol = Pop.dims(2)*nR;
for k=1:numel(params)
    if isempty(params{k})
        params{k} = sparse(nrow,ncol);
    end
end
Pop2 = sopvar(params,Pop.vars,Pop.ZL,Pop.ZR,Pop.dom,Pop.dims);
end


%% ---- cases ------------------------------------------------------------
function case_pde_options()
% Quadratic LF in the PDE state of the heat equation: symmetric and SOS Pmat,
% psatz multipliers, excluded integrals, separability
PIE = get_pie('heat');      X0 = states_x0(1);
mult = struct('exclude',[false;true;true]);
cs = {};
cs{end+1} = {'sym',        struct('degs',struct('in',1,'out',1),'sym',true)};
cs{end+1} = {'pos',        struct('degs',struct('in',1,'out',1),'pos',true)};
cs{end+1} = {'pos, psatz 1', struct('degs',struct('in',1,'out',2),'pos',true,'psatz',1)};
cs{end+1} = {'exclude upper',struct('degs',struct('in',1,'out',1),'pos',true,'exclude',[false;false;true])};
cs{end+1} = {'exclude lower',struct('degs',struct('in',2,'out',1),'pos',true,'exclude',[false;true;false])};
cs{end+1} = {'sep',        struct('degs',struct('in',1,'out',1),'pos',true,'sep',true)};
cs{end+1} = {'multiplier only',struct('degs',struct('in',0,'out',2),'pos',true,'exclude',mult.exclude)};
cs{end+1} = {'default degs, sym',struct('sym',true)};
for c=1:numel(cs)
    res = check_pde(PIE,X0,{'x1'},1,cs{c}{2},c);
    assert_res(res,['PDE, ',cs{c}{1}]);
end
% opts.sym is the default for a functional that is not SOS
opts = struct('degs',struct('in',1,'out',1));
res = check_pde(PIE,X0,{'x1'},1,opts,20);
assert_res(res,'PDE, default (symmetric)');
end

function case_pde_degrees()
% Monomials of degree 2 (cumulative degree), and degrees [1;2], with operators that only multiply
PIE = get_pie('heat');      X0 = states_x0(1);
opts = struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[false;true;true]);
res = check_pde(PIE,X0,{'x1'},2,opts,1);
assert_res(res,'PDE, degree 2');
res = check_pde(PIE,X0,{'x1'},[1;2],opts,2);
assert_res(res,'PDE, degrees 1 and 2');
% Operators with an integral, degree 2: V is a polynomial of degree 4 (slow)
opts = struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[true;true;false]);
res = check_pde(PIE,X0,{'x1'},2,opts,3);
assert_res(res,'PDE, degree 2, upper integral');
end

function case_pde_coupled()
% Two coupled states: both, the second only, monomials mixing the states
PIE = get_pie('coupled');   X0 = states_x0(2);
opts = struct('degs',struct('in',1,'out',1),'pos',true);
res = check_pde(PIE,X0,{'x1','x2'},1,opts,1);
assert_res(res,'two states, pos');
opts = struct('degs',struct('in',0,'out',1),'sym',true);
res = check_pde(PIE,X0,{'x2','x1'},1,opts,2);
assert_res(res,'two states, sym, names unsorted');
opts = struct('degs',struct('in',1,'out',1),'pos',true);
res = check_pde(PIE,X0,{'x2'},1,opts,3);
assert_res(res,'second state only');
opts = struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[false;true;true]);
res = check_pde(PIE,X0,{'x1','x2'},[1 1],opts,4);
assert_res(res,'two states, monomial x1*x2');
end

function case_mixed()
% LF in the PIE state: Pmat with the symmetry condition
PIE = get_pie('heat');      X0 = states_x0(1);
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',1)),1);
assert_res(res,'mixed, in=1, out=1');
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',2,'out',2)),2);
assert_res(res,'mixed, in=2, out=2');
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',1),'exclude',[false;false;true]),3);
assert_res(res,'mixed, excluded upper');
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',2),'exclude',[false;true;false]),4);
assert_res(res,'mixed, excluded lower');
end

function case_mixed_coupled()
PIE = get_pie('coupled');   X0 = states_x0(2);
res = check_mixed(PIE,X0,{'x1','x2'},1,struct('degs',struct('in',1,'out',1)),1);
assert_res(res,'mixed, two states');
res = check_mixed(PIE,X0,{'x2'},1,struct('degs',struct('in',1,'out',1)),2);
assert_res(res,'mixed, second state only');
end

function case_combined()
% Functional in PDE and PIE states: the sum of the two functionals, declared in
% the same program
PIE = get_pie('coupled');
prog = piesos_program(PIE.x);
opts = struct('degs',struct('in',1,'out',1));
[prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,{'x1','x2_pie'},1,opts);
prog1 = piesos_program(PIE.x);
[prog1,V1,dV1,Pmat1,Zop1] = piesos_LF(prog1,PIE,{'x1'},1,opts);
[prog1,V2,dV2,Pmat2,Zop2] = piesos_LF(prog1,PIE,{'x2_pie'},1,opts);
assert(isequal(prog,prog1),'The program must be that of the separate declarations.');
assert(isequal(V,V1+V2) && isequal(dV,dV1+dV2),'V and dV must be the sums of the functionals.');
assert(iscell(Pmat) && isequal(size(Pmat),[2,1]) && isequal(size(Zop),[2,1]),'Pmat and Zop must be cells with the PDE and PIE parts.');
assert(isequal(Pmat{1},Pmat1) && isequal(Pmat{2},Pmat2) && isequal(Zop{1},Zop1) && isequal(Zop{2},Zop2),'Pmat and Zop must be those of the separate declarations.');
% degmat with a column for each variable
[~,V,dV] = piesos_LF(prog,PIE,{'x1','x2_pie'},[1 0;0 1],opts);
assert(isa(V,'polyopvar') && isa(dV,'polyopvar'));
end

function case_input_forms()
% The states as a row or column cellstr, or as polyopvar, and opts omitted
PIE = get_pie('coupled');
opts = struct('degs',struct('in',0,'out',0),'pos',true);
prog = piesos_program(PIE.x);
[~,Va,dVa] = piesos_LF(prog,PIE,{'x1','x2'},1,opts);
[~,Vb,dVb] = piesos_LF(prog,PIE,{'x1';'x2'},1,opts);
assert(isequal(Va,Vb) && isequal(dVa,dVb),'Row and column cellstr must be equivalent.');
xu = polyopvar({'x1';'x2'},[1;1],{'s'},[0 1],true(2,1));
[~,Vc,dVc] = piesos_LF(prog,PIE,xu,1,opts);
assert(isequal(Va,Vc) && isequal(dVa,dVc),'A polyopvar must be equivalent to the names of its states.');
% PIE states, specified by the polyopvar state of the PIE
opts = struct('degs',struct('in',0,'out',0));
[~,Vd,dVd] = piesos_LF(prog,PIE,PIE.x(2,:),1,opts);
[~,Ve,dVe] = piesos_LF(prog,PIE,{'x2_pie'},1,opts);
assert(isequal(Vd,Ve) && isequal(dVd,dVe),'A PIE state must be equivalent to its name.');
% opts omitted
[~,Vf] = piesos_LF(prog,PIE,{'x1'},1);
assert(isa(Vf,'polyopvar'),'V must be a polyopvar.');
% 'pde' and 'pie': all PDE or all PIE states, in any case or as string
opts = struct('degs',struct('in',0,'out',0),'pos',true);
[~,Vg,dVg] = piesos_LF(prog,PIE,'pde',1,opts);
[~,Vh,dVh] = piesos_LF(prog,PIE,PIE.varname_pde,1,opts);
assert(isequal(Vg,Vh) && isequal(dVg,dVh),'''pde'' must be equivalent to all PDE state names.');
[~,Vg,dVg] = piesos_LF(prog,PIE,"PDE",1,opts);
assert(isequal(Vg,Vh) && isequal(dVg,dVh),'Case and string must not matter.');
opts = struct('degs',struct('in',0,'out',0));
[~,Vi,dVi] = piesos_LF(prog,PIE,'pie',1,opts);
[~,Vj,dVj] = piesos_LF(prog,PIE,PIE.varname_pie,1,opts);
assert(isequal(Vi,Vj) && isequal(dVi,dVj),'''pie'' must be equivalent to all PIE state names.');
assert(numel(PIE.varname_pie)==2 && nnz(Vi.degmat(:))>0,'Vacuous: ''pie'' must involve both states.');
end

function case_errors()
PIE = get_pie('heat');
prog = piesos_program(PIE.x);
assert_err(@() piesos_LF(prog,PIE,{},1),'At least one state variable');
assert_err(@() piesos_LF(prog,PIE,{'y'},1),'must appear as state');
assert_err(@() piesos_LF(prog,PIE,5,1),'cellstr');
assert_err(@() piesos_LF(prog,PIE,'x1',1),'''pde'' or ''pie''');
assert_err(@() piesos_LF(prog,PIE,{'x1'},1,5),'''struct''');
assert_err(@() piesos_LF(prog,PIE,{'x1_pie'},1,struct('pos',true)),'positivity of mixed');
assert_err(@() piesos_LF(prog,PIE,{'x1'},0),'Constant');
% operators that differ on the left and right: V is not symmetric
exL = struct('exclude',[false;true;false]);
exR = struct('exclude',[false;false;true]);
assert_err(@() piesos_LF(prog,PIE,{'x1'},1,struct('L',exL,'R',exR)),'same operators');
% function space of provided state variables must match that of the PIE
xbad = polyopvar({'x1'},1,{'s'},[0 2],true);
assert_err(@() piesos_LF(prog,PIE,xbad,1),'Function spaces');
end

function case_skip_derivative()
% opts.derivative=false and fewer than 3 outputs: dV is not computed (returned as []), and
% prog, V, Pmat, Zop are those of the full call, for PDE states (with positivity), PIE states
% (with the symmetry condition), and PDE and PIE states together (recursive calls). The
% profiler counts the calls of subs, which the derivative makes for each factor
PIE = get_pie('coupled');
cs = {'PDE states',         {'x1'},         struct('degs',struct('in',1,'out',1),'pos',true)
      'PIE states',         {'x2_pie'},     struct('degs',struct('in',1,'out',1))
      'PDE and PIE states', {'x1','x2_pie'},struct('degs',struct('in',1,'out',1))};
for c=1:size(cs,1)
    lbl = cs{c,1};      names = cs{c,2};    opts = cs{c,3};
    prog = piesos_program(PIE.x);
    [prog1,V1,dV1,Pmat1,Zop1] = piesos_LF(prog,PIE,names,1,opts);
    assert(isa(dV1,'polyopvar'),'%s: dV must be computed by default.',lbl);
    optsT = opts;       optsT.derivative = true;
    [~,~,dV1b] = piesos_LF(prog,PIE,names,1,optsT);
    assert(isequal(dV1,dV1b),'%s: opts.derivative=true must be the default.',lbl);
    optsF = opts;       optsF.derivative = false;
    [prog2,V2,dV2,Pmat2,Zop2] = piesos_LF(prog,PIE,names,1,optsF);
    assert(isequal(dV2,[]),'%s: dV must be [] for opts.derivative=false.',lbl);
    assert(isequal(prog1,prog2) && isequal(V1,V2) && isequal(Pmat1,Pmat2) && isequal(Zop1,Zop2),...
            '%s: prog, V, Pmat, Zop must not depend on opts.derivative.',lbl);
    % nargout rule: fewer than 3 outputs skips the derivative
    [prog3,V3] = piesos_LF(prog,PIE,names,1,opts);
    prog4 = piesos_LF(prog,PIE,names,1,opts);
    assert(isequal(prog1,prog3) && isequal(V1,V3) && isequal(prog1,prog4),...
            '%s: prog, V must not depend on the number of outputs.',lbl);
    n_full = count_subs(prog,PIE,names,opts,5);
    n_opt = count_subs(prog,PIE,names,optsF,5);
    n_nargout = count_subs(prog,PIE,names,opts,2);
    assert(n_opt<n_full,'%s: opts.derivative=false must skip calls of subs (%d vs %d).',lbl,n_opt,n_full);
    assert(n_nargout==n_opt,'%s: 2 outputs must skip the derivative as well (%d vs %d).',lbl,n_nargout,n_opt);
end
% Outputs ignored with ~ count in nargout: the derivative is computed
[prog,~,dVt] = piesos_LF(piesos_program(PIE.x),PIE,{'x1'},1,struct('degs',struct('in',1,'out',1)));
assert(isa(dVt,'polyopvar'));
% The option must be a scalar boolean
prog = piesos_program(PIE.x);
assert_err(@() piesos_LF(prog,PIE,{'x1'},1,struct('derivative',[true,false])),'scalar boolean');
assert_err(@() piesos_LF(prog,PIE,{'x1'},1,struct('derivative','no')),'scalar boolean');
end

function n = count_subs(prog,PIE,names,opts,nout)
% Number of calls of the 'subs' method of 'polyopvar' made by piesos_LF with nout outputs
profile clear
profile on
cleanup = onCleanup(@() profile('off'));
if nout==2
    [~,~] = piesos_LF(prog,PIE,names,1,opts);
else
    [~,~,~,~,~] = piesos_LF(prog,PIE,names,1,opts);
end
profile off
info = profile('info');
fnames = {info.FunctionTable.FunctionName};
n = sum([info.FunctionTable(strcmp(fnames,'polyopvar.subs')).NumCalls]);
assert(n>0,'Vacuous: the profiler must see the calls of subs.');
end

function case_ode_pde()
% An ODE state coupled to a PDE state through the domain and the boundary (states of different
% spatial dependence), PDE-state functional: V and dV compared with polynomial arithmetic. The
% ODE state is a constant, the PDE-state x2 any polynomial (T*x satisfies the boundary conditions)
PIE = get_pie('ode_pde');
sv = polynomial({'s'});     sv = sv(1);
X0 = {polynomial(0.7), 2*sv-3*sv^2+sv^3};
cs = {};
cs{end+1} = {{'x1','x2'},  struct('degs',struct('in',1,'out',1),'pos',true)};
cs{end+1} = {{'x1','x2'},  struct('degs',struct('in',0,'out',1),'sym',true)};
cs{end+1} = {{'x2'},       struct('degs',struct('in',1,'out',1),'pos',true)};
cs{end+1} = {{'x1'},       struct('degs',struct('in',0,'out',0),'pos',true)};
cs{end+1} = {{'x1','x2'},  struct('degs',struct('in',1,'out',1),'pos',true,'exclude',[false;false;true])};
for c=1:numel(cs)
    res = check_pde(PIE,X0,cs{c}{1},1,cs{c}{2},c);
    assert_res(res,sprintf('ODE-PDE, case %d',c));
end
end

function vals = pie_vals(PIE,X0)
% Values of the PIE states by name
vals = struct();
for i=1:numel(X0)
    vals.(PIE.varname_pie{i}) = X0{i};
end
end

function res = check_equiv(PIEa,Xa,namesa,PIEb,Xb,namesb,degmat,opts,seed)
% V and dV declared in two equivalent PIEs (a: scalar states, b: with a vector-valued state) for
% the states namesa, namesb, with the same random decision variables and the same PIE state values
% (stacked as vectors in b). The scalar case is verified against polynomial arithmetic (check_pde)
va = pie_vals(PIEa,Xa);         vb = pie_vals(PIEb,Xb);
proga = piesos_program(PIEa.x);     progb = piesos_program(PIEb.x);
[~,Va,dVa,Pa] = piesos_LF(proga,PIEa,namesa,degmat,opts);
[~,Vb,dVb,Pb] = piesos_LF(progb,PIEb,namesb,degmat,opts);
% Pmat is partitioned differently (blocks per state), so match its values rather than the variables
Ma = block_matrix(Pa);              Mb = block_matrix(Pb);
assert(isequal(size(Ma),size(Mb)),'The declared Pmat must have the same size.');
da = random_decvars(Pa,seed);
db = match_decvars(eval_dpvar_at(Ma,da),Mb);
res.Va = double(eval_dpvar_at(eval_V(Va,va),da));       res.Vb = double(eval_dpvar_at(eval_V(Vb,vb),db));
res.dVa = double(eval_dpvar_at(eval_V(dVa,va),da));     res.dVb = double(eval_dpvar_at(eval_V(dVb,vb),db));
res.errV = abs(res.Va-res.Vb)/max(1e-12,abs(res.Va));
res.errdV = abs(res.dVa-res.dVb)/max(1e-12,abs(res.dVa));
end

function d0map = match_decvars(Pnum,Mb)
% Values of the decision variables of the matrix Mb ('dpvar') for which Mb equals the matrix of
% polynomials Pnum, by least squares on points in s (Mb is linear in its decision variables)
names = Mb.dvarname(:);         nd = numel(names);
pts = linspace(0.1,0.9,5);
f0 = eval_points(eval_dpvar_at(Mb,containers.Map(names,num2cell(zeros(nd,1)))),pts);
A = zeros(numel(f0),nd);
for k=1:nd
    e = zeros(nd,1);        e(k) = 1;
    A(:,k) = eval_points(eval_dpvar_at(Mb,containers.Map(names,num2cell(e))),pts) - f0;
end
target = eval_points(Pnum,pts);
d = pinv(A)*(target-f0);
assert(norm(A*d+f0-target)<=1e-8*max(1,norm(target)),'The matrices Pmat of the two representations cannot be matched.');
d0map = containers.Map(names,num2cell(d));
end

function v = eval_points(P,pts)
% Values of the matrix of polynomials P in s at the points pts, stacked
sv = polynomial({'s'});     sv = sv(1);
v = zeros(0,1);
for p=pts
    if isempty(P.varname)
        Pk = double(P);
    else
        Pk = double(subs(P,sv,p));
    end
    v = [v;Pk(:)]; %#ok<AGROW>
end
end

function assert_equiv(res,label)
assert(abs(res.Va)>1e-8 && abs(res.dVa)>1e-8,'%s: the values %g, %g are too small to test.',label,res.Va,res.dVa);
assert(res.errV<=1e-8,'%s: V(x0)=%.10g differs from %.10g of the scalar states.',label,res.Vb,res.Va);
assert(res.errdV<=1e-8,'%s: dV(x0)=%.10g differs from %.10g of the scalar states.',label,res.dVb,res.dVa);
end

function res = check_pde_blocks(PIE,X0,names,opts,seed)
% As check_pde, for the degree-1 monomials Z(u)=u of states that may be vector-valued: Zop and Pmat
% are those of the functional itself (one block per state), and V and dV are compared with
%   V = int Y'*P*Y,  Y = Zop*u,   dV = d/dh V(u0+h*w0) at h=0
% for the PDE states u0=T*x0 and their velocity w0=f(x0) (the rows of the states of size n are n rows).
% Zop and Pmat have one block per state (dmonomials has the states as variables): block i acts on
% state i, of size n, and has n rows
pvar h
sv = polynomial({'s'});     sv = sv(1);
[vals,~,~,u0,w0] = pie_values(PIE,X0);
idcs = state_idcs(PIE,names,true);
row0 = cumsum([0;PIE.x_tab(:,2)]);
prog = piesos_program(PIE.x);
[prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,names,1,opts); %#ok<ASGLU>
d0map = random_decvars(Pmat,seed);
uh = cell(numel(idcs),1);
for k=1:numel(idcs)
    r = row0(idcs(k))+1:row0(idcs(k)+1);
    uh{k} = u0(r)+h*w0(r);
end
% One input per block of Zop, the (vector-valued) state
Xc = cell(1,numel(uh));
for k=1:numel(uh)
    Xc{k} = {uh{k}};
end
Y = apply_op(fix_empty_params_zop(Zop),Xc);
Ph = block_matrix(cellfun(@(P) eval_dpvar_at(P,d0map),Pmat,'UniformOutput',false));
Fh = int(Y.'*Ph*Y,sv,0,1);
res.Vtrue = double(subs(Fh,h,0));
res.dVtrue = double(subs(diff(Fh,h),h,0));
res.V = double(eval_dpvar_at(eval_V(V,vals),d0map));
res.dV = double(eval_dpvar_at(eval_V(dV,vals),d0map));
res.errV = abs(res.V-res.Vtrue)/max(1e-12,abs(res.Vtrue));
res.errdV = abs(res.dV-res.dVtrue)/max(1e-12,abs(res.dVtrue));
end

function case_vector_states()
% Vector-valued states. PDE-state functionals with monomials Z(u)=u, compared with polynomial
% arithmetic (check_pde_blocks): a vector state, and a vector state with a scalar state, also for
% functionals in a subset of the states (rows of T and f selected by state). Functionals with only
% a multiplier (so that the order of the basis operators is the same) are compared with those of
% the equivalent scalar states, also for PIE states. All comparisons are run, failures reported together
fails = {};
% (i) one vector state of size 2: PDE states
PIE = get_pie('vec');           Xa = states_x0(2);          Xv = {[Xa{1};Xa{2}]};
cs = {'sym',            struct('degs',struct('in',1,'out',1),'sym',true)
      'pos',            struct('degs',struct('in',1,'out',1),'pos',true)
      'multiplier only',struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[false;true;true])
      'upper integral', struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[true;true;false])
      'psatz',          struct('degs',struct('in',1,'out',2),'pos',true,'psatz',1)};
for c=1:size(cs,1)
    fails = [fails,try_blocks(['vector state, ',cs{c,1}],PIE,Xv,{'x1'},cs{c,2},c)]; %#ok<AGROW>
end
% (ii) vector state of size 2 and a scalar state: all states, and each of the subsets
PIE = get_pie('vec_sc');        X3 = states_x0(3);          Xvs = {[X3{1};X3{2}],X3{3}};
opts = struct('degs',struct('in',1,'out',1),'pos',true);
cs = {'all states',{'x1','x2'}; 'vector state',{'x1'}; 'scalar state',{'x2'}};
for c=1:size(cs,1)
    fails = [fails,try_blocks(['vector and scalar state: ',cs{c,1}],PIE,Xvs,cs{c,2},opts,10+c)]; %#ok<AGROW>
end
% (iii) equivalence with scalar states, for functionals that only multiply, PDE and PIE states
opts = struct('degs',struct('in',0,'out',1),'exclude',[false;true;true]);
fails = [fails,try_equiv('equivalence, vector state, PDE states',get_pie('coupled'),Xa,{'x1','x2'},get_pie('vec'),Xv,{'x1'},1,opts,1)];
fails = [fails,try_equiv('equivalence, vector state, PIE states',get_pie('coupled'),Xa,{'x1_pie','x2_pie'},get_pie('vec'),Xv,{'x1_pie'},1,opts,2)];
PIEa = get_pie('sc3');          PIEb = get_pie('vec_sc');
fails = [fails,try_equiv('equivalence, vector and scalar state, PDE states',PIEa,X3,{'x1','x2','x3'},PIEb,Xvs,{'x1','x2'},1,opts,3)];
fails = [fails,try_equiv('equivalence, vector and scalar state, PIE state: scalar state',PIEa,X3,{'x3_pie'},PIEb,Xvs,{'x2_pie'},1,opts,4)];
fails = [fails,try_equiv('equivalence, vector and scalar state, PIE state: vector state',PIEa,X3,{'x1_pie','x2_pie'},PIEb,Xvs,{'x1_pie'},1,opts,5)];
assert(isempty(fails),'%d comparisons failed:\n%s',numel(fails),strjoin(fails,newline));
end

function msg = try_blocks(label,PIE,X0,names,opts,seed)
% Cell with the error message of check_pde_blocks/assert_res, empty if they pass
msg = {};
try
    res = check_pde_blocks(PIE,X0,names,opts,seed);
    assert_res(res,label);
catch err
    msg = {sprintf('%s: %s',label,err.message)};
end
end
function msg = try_equiv(label,PIEa,Xa,namesa,PIEb,Xb,namesb,degmat,opts,seed)
% Cell with the error message of check_equiv/assert_equiv, empty if they pass
msg = {};
try
    res = check_equiv(PIEa,Xa,namesa,PIEb,Xb,namesb,degmat,opts,seed);
    assert_equiv(res,label);
catch err
    msg = {sprintf('%s: %s',label,err.message)};
end
end

function case_nonlinear()
% PIEs with nonlinear terms: pointwise products of operators applied to the PIE state, which
% give nested tensor-PI operators when substituted into the functional. V and dV compared with
% polynomial arithmetic
PIE = get_pie('burgers');   X0 = states_x0(1);     % u_t = u_ss + 2u - u*u_s
cs = {};
cs{end+1} = {'multiplier',  struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[false;true;true])};
cs{end+1} = {'upper integral',struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[true;true;false])};
cs{end+1} = {'all integrals',struct('degs',struct('in',1,'out',1),'pos',true)};
cs{end+1} = {'sym',         struct('degs',struct('in',1,'out',2),'sym',true)};
for c=1:numel(cs)
    res = check_pde(PIE,X0,{'x1'},1,cs{c}{2},c);
    assert_res(res,['Burgers PDE state, ',cs{c}{1}]);
end
% PIE state, with the symmetry condition
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',1)),1);
assert_res(res,'Burgers PIE state');
% Cubic term
PIE = get_pie('fisher');
res = check_pde(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',1),'pos',true),1);
assert_res(res,'u^3, PDE state');
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',1)),2);
assert_res(res,'u^3, PIE state');
% Coupled states, with nonlinear terms in both (some rows of f lack some of the monomials)
PIE = get_pie('coupled_nl');    X0 = states_x0(2);
res = check_pde(PIE,X0,{'x1','x2'},1,struct('degs',struct('in',1,'out',1),'pos',true),1);
assert_res(res,'coupled nonlinear, PDE states');
res = check_mixed(PIE,X0,{'x1','x2'},1,struct('degs',struct('in',1,'out',1)),3);
assert_res(res,'coupled nonlinear, PIE states');
end
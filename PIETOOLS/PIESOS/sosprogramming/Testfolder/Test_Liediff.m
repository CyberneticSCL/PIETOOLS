function Test_Liediff(only)
% TEST_LIEDIFF(ONLY) tests, using assertions, that dV = Liediff(V,PIE)
% computes the Lie derivative of a functional V along the PIE
% d/dt (T*x) = f(x), for u=T*x the PDE state:
%   (PDE states)  V(u) = <Zop*Z(u), Pmat*Zop*Z(u)>, product rule over all
%                 factors;
%   (PIE states)  V(x,u) = <Zop*Z(x), Pmat*Zop*Z(u)>, with the symmetry
%                 condition V(xL,T*xR) = V(xR,T*xL) assumed: twice the
%                 product rule over the factors u.
%
% Oracle: the PIEs are converted from PDEs with Dirichlet conditions, so
% a polynomial PIE state x=u'' (u(0)=u(1)=u''(0)=u''(1)=0) has u0=T*x0,
% velocity w0=f(x0) and xdot0=(w0)'' with T*xdot0=w0. V and dV at x0 are
% compared with direct polynomial arithmetic (apply_op and integration,
% no composition of functionals):
%    PDE:  V = <Zop*Z(u0),P Zop*Z(u0)>,  dV = d/dh V(u0+h*w0) at h=0,
%    PIE:  V = <Zop*Z(x0),P Zop*Z(T*x0)>, dV = d/dh V(x0+h*xdot0) at h=0,
% with random decision variables for PDE states (general nonsymmetric
% Pmat), and for PIE states a feasible symmetric Pmat solved with SeDuMi,
% or a functional symmetric by construction.
% Covers the product rule, doubling for the factors u only, composition
% with T and f, nested operators of nonlinear f, monomials in different
% states, and a known function substituted for a factor. The errors are
% tested as well.
%
% INPUTS
% - only:   (optional) cell of character arrays with the names of the
%           cases to run, e.g. Test_Liediff({'case_mixed'}). Defaults to
%           running all cases;
%
% NOTES
% - PIE-state functionals V(x,u) assume V(xL,T*xR)=V(xR,T*xL); case_mixed
%   imposes it on Pmat with piesos_eq, and Liediff doubles the derivative
%   of the factors u only;
% - The degree-2 monomial cases use operators with only a multiplier part,
%   or only an upper integral part;
% - eval_V names dummy variables s_k_j as apply_op does; the dummy
%   variables are matched by name, so they must be renamed per factor;
% - fix_empty_params_* work around apply_sopvar failing on empty params
%   (as in Test_piesos_LF).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_Liediff
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
% DJ (CC), 10/22/2026: apply_side: Zop acts on the entries of the basis
%                        (rows of Z.C), not on its monomials;

tests ={@case_pde, @case_pde_degrees, @case_pde_coupled, @case_mixed, ...
         @case_mixed_known, @case_combined, @case_unused_states, @case_nonlinear, ...
         @case_multiplier, ...
         @case_errors};
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
    error('Test_Liediff: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_Liediff: all %d cases passed.\n',numel(tests));
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

function [valsx,valsu,u0,w0,xdot0] = pie_values(PIE,X0)
% PDE state u0=T*x0, its velocity w0=f(x0), and the PIE state velocity
% xdot0=(w0)'' (T*xdot0=w0 for w0 satisfying the Dirichlet conditions); the values of
% the PIE states and of the PDE states, by name
sv = polynomial({'s'});     sv = sv(1);
Tx = PIE.T*PIE.x;           fx = PIE.f;
u0 = eval_polyopvar(Tx,X0);
w0 = eval_polyopvar(fx,X0);
xdot0 = cell(1,numel(X0));
for i=1:numel(X0)
    xdot0{i} = diff(w0(i),sv,2);
end
valsx = struct();       valsu = struct();
for i=1:numel(X0)
    valsx.(PIE.varname_pie{i}) = X0{i};
    valsu.(PIE.varname_pde{i}) = u0(i);
end
end

function idcs = state_idcs(PIE,names)
% Indices in the PIE of the PDE states with the given names, sorted by name
list = PIE.varname_pde(:);
[~,idcs] = ismember(names(:),list);
[~,srt] = sort(list(idcs));
idcs = idcs(srt);
end

function Zu = monomials_in(PIE,idcs,degmat)
% The monomials Z(u) in the PDE states idcs
xu = polyopvar(PIE.varname_pde(idcs),PIE.x_tab(idcs,2),{'s'},PIE.dom,logical(PIE.x_tab(idcs,3)));
Zu = dmonomials(xu,degmat);
end

function res = check_pde(PIE,X0,names,degmat,opts,seed)
% Declare a functional V(u) in the PDE states 'names' with Pmat random, take the derivative
% with Liediff, and compare with polynomial arithmetic
pvar h
sv = polynomial({'s'});     sv = sv(1);
[valsx,valsu,u0,w0] = pie_values(PIE,X0);
idcs = state_idcs(PIE,names);
prog = piesos_program(PIE.x);
Zu = monomials_in(PIE,idcs,degmat);
[prog,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog,Zu,Zu,opts); %#ok<ASGLU>
dV = Liediff(V,PIE);
assert(all(ismember(dV.varname,PIE.varname_pie)),'dV must be expressed in the PIE states.');
d0map = random_decvars(Pmat,seed);
Xu = struct();
for i=1:numel(idcs)
    Xu.(PIE.varname_pde{idcs(i)}) = u0(idcs(i)) + h*w0(idcs(i));
end
Ph = block_matrix(cellfun(@(P) eval_dpvar_at(P,d0map),Pmat,'UniformOutput',false));
YL = apply_side(ZopL,Zu,Xu);
YR = apply_side(ZopR,Zu,Xu);
Fh = int(YL.'*Ph*YR,sv,0,1);
res.Vtrue = double(subs(Fh,h,0));
res.dVtrue = double(subs(diff(Fh,h),h,0));
res.V = double(eval_dpvar_at(eval_V(V,valsu),d0map));       % V in the PDE states
res.dV = double(eval_dpvar_at(eval_V(dV,valsx),d0map));     % dV in the PIE states
res.errV = abs(res.V-res.Vtrue)/max(1e-12,abs(res.Vtrue));
res.errdV = abs(res.dV-res.dVtrue)/max(1e-12,abs(res.dVtrue));
end

function res = check_mixed(PIE,X0,names,degmat,opts,seed)
% Declare a functional V(x,u) in the PIE states and the associated PDE states 'names', with
% the symmetry condition V(xL,T*xR)=V(xR,T*xL) imposed on Pmat, solve for a feasible Pmat with
% a random normalization, take the derivative with Liediff, and compare with polynomial
% arithmetic at the state x0 with velocity xdot0
pvar h
sv = polynomial({'s'});     sv = sv(1);
[valsx,valsu,u0,w0,xdot0] = pie_values(PIE,X0);
idcs = state_idcs(PIE,names);
names_pde = PIE.varname_pde(idcs);      names_x = PIE.varname_pie(idcs);
prog = piesos_program(PIE.x);
Zu = monomials_in(PIE,idcs,degmat);
Zx = subs(Zu,Zu.varname,names_x);
[prog,Vxu,Pmat,ZopL,ZopR] = piesos_quadvar(prog,Zx,Zu,opts);
assert(isequal(ZopL,ZopR),'The same operators must be used on the left and right.');
% Symmetry condition V(xL,T*xR) = V(xR,T*xL) for independent copies xL, xR of the PIE state
Top = PIE.T;
if ~isequal(idcs(:)',1:numel(PIE.varname_pde))
    Top = Top(idcs,:);
end
Tx = Top*PIE.x;
TxL = subs(Tx,Tx.varname,strcat(Tx.varname,'_L'));
TxR = subs(Tx,Tx.varname,strcat(Tx.varname,'_R'));
WLR = subs(subs(Vxu,names_x,strcat(names_x,'_L')),names_pde,TxR);
WRL = subs(subs(Vxu,names_x,strcat(names_x,'_R')),names_pde,TxL);
prog = piesos_eq(prog,WLR-WRL);
dV = Liediff(Vxu,PIE);
% Normalization: random linear combination of the decision variables equal to 1
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)]; %#ok<AGROW>
end
names = unique(names,'stable');
rng(seed);
dsum = dpvar(sparse([0;randn(numel(names),1)]),zeros(1,0),cell(0,1),names,[1 1]);
prog = soseq(prog,dsum-1);
% SeDuMi needs a conic variable, and more variables than equations (the symmetry condition
% is a redundant system: "Size b mismatch" otherwise)
ncon = sum(cellfun(@(A) size(A,2),prog.expr.At));
[prog,~] = sosposmatr(prog,ceil(sqrt(ncon))+2);
evalc('[prog,info] = sossolve(prog);');
assert(info.pinf==0 && info.dinf==0,'The symmetry condition must be feasible with a nonzero Pmat.');
Psol = piesos_getsol(prog,Pmat);
Vsol = piesos_getsol(prog,Vxu);
dVsol = piesos_getsol(prog,dV);
Xl = struct();      Xr = struct();
for i=1:numel(idcs)
    Xl.(PIE.varname_pde{idcs(i)}) = X0{idcs(i)} + h*xdot0{idcs(i)};
    Xr.(PIE.varname_pde{idcs(i)}) = u0(idcs(i)) + h*w0(idcs(i));
end
Ph = block_matrix(Psol);
YL = apply_side(ZopL,Zu,Xl);
YR = apply_side(ZopR,Zu,Xr);
Fh = int(YL.'*Ph*YR,sv,0,1);
res.Vtrue = double(subs(Fh,h,0));
res.dVtrue = double(subs(diff(Fh,h),h,0));
res.V = double(eval_V(Vsol,mergestruct(valsx,valsu)));      % V(x,u): both PIE and PDE states
res.dV = double(eval_V(dVsol,valsx));
res.errV = abs(res.V-res.Vtrue)/max(1e-12,abs(res.Vtrue));
res.errdV = abs(res.dV-res.dVtrue)/max(1e-12,abs(res.dVtrue));
end

function res = check_known(PIE,X0,V,names)
% Functional V(x,u) that is symmetric by construction, in the PIE states and the PDE states
% 'names', without decision variables: compare with polynomial arithmetic, for the exact
% derivative V(x0+h*xdot0,u0+h*w0) of the functional of independent x and u
pvar h
[valsx,valsu,u0,w0,xdot0] = pie_values(PIE,X0);
dV = Liediff(V,PIE);
res.V = double(eval_V(V,mergestruct(valsx,valsu)));
res.dV = double(eval_V(dV,valsx));
Vh = V_known(h,PIE,X0,u0,w0,xdot0,names);
res.Vtrue = double(subs(Vh,h,0));
res.dVtrue = double(subs(diff(Vh,h),h,0));
res.errV = abs(res.V-res.Vtrue)/max(1e-12,abs(res.Vtrue));
res.errdV = abs(res.dV-res.dVtrue)/max(1e-12,abs(res.dVtrue));
end

function Vh = V_known(h,PIE,X0,u0,w0,xdot0,names) %#ok<INUSL>
% sum_i int_0^1 x_i(s)*u_i(s) ds at (x0+h*xdot0, u0+h*w0), over the states in names
sv = polynomial({'s'});     sv = sv(1);
idcs = state_idcs(PIE,names);
Vh = 0;
for i=idcs(:)'
    Vh = Vh + int((X0{i}+h*xdot0{i})*(u0(i)+h*w0(i)),sv,0,1);
end
end

function out = mergestruct(a,b)
out = a;
fn = fieldnames(b);
for i=1:numel(fn)
    out.(fn{i}) = b.(fn{i});
end
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
% depending on s_k, as s_k_j, if the monomial has more than 1 factor. Zeros in Z are
% multipliers from R, with the dummy input 1.
vnames = V.vars.in;
dep_state = double(V.depmat.in)>0;
Y = 0;
Zf = full(V.Z);
for k=1:size(V.Z,1)
    op = V.C.ops{1,k};
    if isempty(op), continue, end
    pos = Zf(k,:);
    d = nnz(pos);
    if d==0
        Y = Y + apply_op(op,polynomial(1));
    else
        slot = zeros(1,numel(vnames));
        v = 1;
        for j=1:numel(pos)
            if pos(j)==0
                continue                % multiplier from R (or padding): input 1
            end
            xj = vals.(V.varname{pos(j)});
            for mm=1:numel(vnames)
                if dep_state(pos(j),mm)
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

%% ---- workaround for apply_sopvar (as in Test_piesos_LF.m) ---------------
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
function case_pde()
% Quadratic functional in the PDE state of the heat equation: symmetric, general
% (nonsymmetric Pmat: the product rule needs no symmetry) and SOS Pmat, psatz multipliers,
% excluded integrals
PIE = get_pie('heat');      X0 = states_x0(1);
mult = struct('exclude',[false;true;true]);
cs = {};
cs{end+1} = {'sym',        struct('degs',struct('in',1,'out',1),'sym',true)};
cs{end+1} = {'general',    struct('degs',struct('in',1,'out',1))};
cs{end+1} = {'pos',        struct('degs',struct('in',1,'out',1),'pos',true)};
cs{end+1} = {'pos, psatz 1', struct('degs',struct('in',1,'out',2),'pos',true,'psatz',1)};
cs{end+1} = {'exclude upper',struct('degs',struct('in',1,'out',1),'exclude',[false;false;true])};
cs{end+1} = {'exclude lower',struct('degs',struct('in',2,'out',1),'pos',true,'exclude',[false;true;false])};
cs{end+1} = {'sep',        struct('degs',struct('in',1,'out',1),'sep',true)};
cs{end+1} = {'multiplier only',struct('degs',struct('in',0,'out',2),'exclude',mult.exclude)};
cs{end+1} = {'default degs',struct()};
for c=1:numel(cs)
    res = check_pde(PIE,X0,{'x1'},1,cs{c}{2},c);
    assert_res(res,['PDE, ',cs{c}{1}]);
end
end

function case_pde_degrees()
% Monomials of degree 2 (cumulative degree), and degrees [1;2]: the product rule over 4 and
% over 2 and 4 factors, with operators that only multiply
PIE = get_pie('heat');      X0 = states_x0(1);
opts = struct('degs',struct('in',0,'out',1),'exclude',[false;true;true]);
res = check_pde(PIE,X0,{'x1'},2,opts,1);
assert_res(res,'PDE, degree 2');
res = check_pde(PIE,X0,{'x1'},[1;2],opts,2);
assert_res(res,'PDE, degrees 1 and 2');
% Operators with an upper integral, degree 2
opts = struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[true;true;false]);
res = check_pde(PIE,X0,{'x1'},2,opts,3);
assert_res(res,'PDE, degree 2, upper integral');
end

function case_pde_coupled()
% Two coupled states: both, the second only, monomials mixing the states
PIE = get_pie('coupled');   X0 = states_x0(2);
opts = struct('degs',struct('in',1,'out',1));
res = check_pde(PIE,X0,{'x1','x2'},1,opts,1);
assert_res(res,'two states, general');
opts = struct('degs',struct('in',0,'out',1),'sym',true);
res = check_pde(PIE,X0,{'x2','x1'},1,opts,2);
assert_res(res,'two states, sym, names unsorted');
opts = struct('degs',struct('in',1,'out',1),'pos',true);
res = check_pde(PIE,X0,{'x2'},1,opts,3);
assert_res(res,'second state only');
opts = struct('degs',struct('in',0,'out',1),'exclude',[false;true;true]);
res = check_pde(PIE,X0,{'x1','x2'},[1 1],opts,4);
assert_res(res,'two states, monomial x1*x2');
end

function case_mixed()
% Functional in PIE and PDE states, with the symmetry condition imposed on Pmat
PIE = get_pie('heat');      X0 = states_x0(1);
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',1)),1);
assert_res(res,'mixed, in=1, out=1');
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',2,'out',2)),2);
assert_res(res,'mixed, in=2, out=2');
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',1),'exclude',[false;false;true]),3);
assert_res(res,'mixed, excluded upper');
res = check_mixed(PIE,X0,{'x1'},1,struct('degs',struct('in',1,'out',2),'exclude',[false;true;false]),4);
assert_res(res,'mixed, excluded lower');
PIE = get_pie('coupled');   X0 = states_x0(2);
res = check_mixed(PIE,X0,{'x1','x2'},1,struct('degs',struct('in',1,'out',1)),1);
assert_res(res,'mixed, two states');
res = check_mixed(PIE,X0,{'x2'},1,struct('degs',struct('in',1,'out',1)),2);
assert_res(res,'mixed, second state only');
end

function case_mixed_known()
% V(x,u) = int x_i*u_i ds, symmetric since the Green's function T of the Dirichlet problem is
% self-adjoint, with no decision variables: dV = 2*int x*f(x) ds, as the exact derivative
PIE = get_pie('coupled');   X0 = states_x0(2);
xu = polyopvar(PIE.varname_pde,PIE.x_tab(:,2),{'s'},PIE.dom,logical(PIE.x_tab(:,3)));
V1 = innerprod(PIE.x(1,:),xu(1,:));
res = check_known(PIE,X0,V1,{'x1'});
assert_res(res,'int x1*u1');
% Two states, and a subset of the PDE states appearing in V
V12 = V1 + innerprod(PIE.x(2,:),xu(2,:));
res = check_known(PIE,X0,V12,{'x1','x2'});
assert_res(res,'int x1*u1 + x2*u2');
V2 = innerprod(PIE.x(2,:),xu(2,:));
res = check_known(PIE,X0,V2,{'x2'});
assert_res(res,'int x2*u2 (second state only)');
% Not in the order of the states
V21 = V2 + V1;
res = check_known(PIE,X0,V21,{'x1','x2'});
assert_res(res,'int x2*u2 + x1*u1');
end

function case_combined()
% Functional with terms in PDE states only, and terms in PDE and PIE states (the factor 2 is
% for the latter only), in the same polyopvar
PIE = get_pie('coupled');   X0 = states_x0(2);
pvar h
sv = polynomial({'s'});     sv = sv(1);
[valsx,valsu,u0,w0,xdot0] = pie_values(PIE,X0);
xu = polyopvar(PIE.varname_pde,PIE.x_tab(:,2),{'s'},PIE.dom,logical(PIE.x_tab(:,3)));
V1 = innerprod(xu(1,:),xu(1,:));                % int u1^2
V2 = innerprod(PIE.x(2,:),xu(2,:));             % int x2*u2
a = 0.7;    b = 1.3;
V = a*V1 + b*V2;
dV = Liediff(V,PIE);
Vh = a*int((u0(1)+h*w0(1))^2,sv,0,1) + b*int((X0{2}+h*xdot0{2})*(u0(2)+h*w0(2)),sv,0,1);
Vtrue = double(subs(Vh,h,0));       dVtrue = double(subs(diff(Vh,h),h,0));
res = struct('Vtrue',Vtrue,'dVtrue',dVtrue);
res.V = double(eval_V(V,mergestruct(valsx,valsu)));
res.dV = double(eval_V(dV,valsx));
res.errV = abs(res.V-Vtrue)/abs(Vtrue);     res.errdV = abs(res.dV-dVtrue)/abs(dVtrue);
assert_res(res,'u1^2 + x2*u2');
% Linearity
dV1 = Liediff(V1,PIE);      dV2 = Liediff(V2,PIE);
dVlin = a*dV1 + b*dV2;
assert(abs(double(eval_V(dVlin,valsx))-res.dV)<=1e-8*abs(res.dV),'Liediff must be linear in V.');
assert(isequal(sort(dV.varname),sort(unique([dV1.varname;dV2.varname]))),'dV must be in the PIE states of its terms.');
end

function case_unused_states()
% States in V.varname that appear in no monomial are ignored, as left by subs (the states
% x_L of V(x_L,u) after x_L -> T*x in piesos_LF), also if not a state of the PIE
PIE = get_pie('coupled');   X0 = states_x0(2);
pvar h
sv = polynomial({'s'});     sv = sv(1);
[valsx,~,u0,w0] = pie_values(PIE,X0);
xu = polyopvar(PIE.varname_pde,PIE.x_tab(:,2),{'s'},PIE.dom,logical(PIE.x_tab(:,3)));
V1 = innerprod(xu(1,:),xu(1,:));                % int u1^2
yv = polyopvar({'y'},1,{'s'},PIE.dom,true);
Vs = {common_states(V1,xu(2,:)), common_states(V1,PIE.x(2,:)), common_states(V1,yv)};
Vtrue = double(int(u0(1)^2,sv,0,1));
dVtrue = double(subs(diff(int((u0(1)+h*w0(1))^2,sv,0,1),h),h,0));
for i=1:numel(Vs)
    assert(numel(Vs{i}.varname)>numel(unique(nonzeros(full(Vs{i}.Z)))),'Vacuous: V must list an unused state.');
    dV = Liediff(Vs{i},PIE);
    res = struct('Vtrue',Vtrue,'dVtrue',dVtrue,'V',Vtrue,'dV',double(eval_V(dV,valsx)));
    res.errV = 0;       res.errdV = abs(res.dV-dVtrue)/abs(dVtrue);
    assert_res(res,sprintf('int u1^2 with an unused state (%d)',i));
end
end

function case_nonlinear()
% PIEs with nonlinear terms: pointwise products of operators applied to the PIE state, which
% give nested tensor-PI operators when substituted into the functional
PIE = get_pie('burgers');   X0 = states_x0(1);     % u_t = u_ss + 2u - u*u_s
cs = {};
cs{end+1} = {'multiplier',  struct('degs',struct('in',0,'out',1),'exclude',[false;true;true])};
cs{end+1} = {'upper integral',struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[true;true;false])};
cs{end+1} = {'all integrals',struct('degs',struct('in',1,'out',1))};
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
res = check_pde(PIE,X0,{'x1','x2'},1,struct('degs',struct('in',1,'out',1)),1);
assert_res(res,'coupled nonlinear, PDE states');
res = check_mixed(PIE,X0,{'x1','x2'},1,struct('degs',struct('in',1,'out',1)),3);
assert_res(res,'coupled nonlinear, PIE states');
end

function case_multiplier()
% Functional of a monomial of which one factor was replaced by a known function m(s), by subs:
% V = int u1*m*u2 ds, in which m is absorbed in the kernel (a functional has no zeros in Z), so
% that the derivative follows the product rule over the 2 remaining state factors
PIE = get_pie('coupled');   X0 = states_x0(2);
pvar h
sv = polynomial({'s'});     sv = sv(1);
[valsx,valsu,u0,w0] = pie_values(PIE,X0);
xu = polyopvar(PIE.varname_pde,PIE.x_tab(:,2),{'s'},PIE.dom,logical(PIE.x_tab(:,3)));
m = 1+2*sv;
% V = int u1*u1*u2 ds, with the second factor replaced by the function m: V = int u1*m*u2 ds
V3 = innerprod(xu(1,:).*xu(1,:),xu(2,:));
V3m = subs(V3,2,m);
assert(isequal(size(V3m.Z),[1,2]),'Vacuous: V must be a quadratic monomial u1*u2.');
dV = Liediff(V3m,PIE);
Vh = int((u0(1)+h*w0(1))*m*(u0(2)+h*w0(2)),sv,0,1);
Vtrue = double(subs(Vh,h,0));       dVtrue = double(subs(diff(Vh,h),h,0));
res = struct('Vtrue',Vtrue,'dVtrue',dVtrue);
res.V = double(eval_V(V3m,valsu));
res.dV = double(eval_V(dV,valsx));
res.errV = abs(res.V-Vtrue)/abs(Vtrue);     res.errdV = abs(res.dV-dVtrue)/abs(dVtrue);
assert_res(res,'int u1*m*u2');
end

function case_errors()
PIE = get_pie('coupled');
xu = polyopvar(PIE.varname_pde,PIE.x_tab(:,2),{'s'},PIE.dom,logical(PIE.x_tab(:,3)));
V1 = innerprod(xu(1,:),xu(1,:));
assert_err(@() Liediff(5,PIE),'polyopvar');
assert_err(@() Liediff([V1;V1],PIE),'scalar-valued');
assert_err(@() Liediff(V1,rmfield(PIE,'f')),'with fields');
assert_err(@() Liediff(V1,5),'with fields');
% states of V that are not states of the PIE
yv = polyopvar({'y'},1,{'s'},PIE.dom,true);
assert_err(@() Liediff(innerprod(yv,yv),PIE),'must appear as PDE or PIE state');
% a monomial in PIE states only
assert_err(@() Liediff(innerprod(PIE.x(1,:),PIE.x(1,:)),PIE),'PIE states only');
assert_err(@() Liediff(V1+innerprod(PIE.x(1,:),PIE.x(1,:)),PIE),'PIE states only');
% PDE and PIE states of the same name
PIE2 = PIE;     PIE2.varname_pie = PIE2.varname_pde;
assert_err(@() Liediff(V1,PIE2),'different names');
end

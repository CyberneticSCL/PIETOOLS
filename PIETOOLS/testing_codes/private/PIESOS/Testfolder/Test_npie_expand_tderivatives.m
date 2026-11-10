% Test_npie_expand_tderivatives.m -- tests expand_tderivatives for an
% 'npie_struct', checking the expanded PIE against the semantics of the
% original PIE. For (d/dt)^{p_i}(T*x)_i = f_i(x) (no inputs), the residual
%   R_i = (T*x^{(p_i)})_i - f_i(x),         x^{(k)} = (d/dt)^k x,
% at a fixed time t0 depends only on the values of the time derivatives
% x^{(k)}(t0,.), which can be chosen arbitrarily (random polynomials in the
% spatial variables). The expanded PIE, with the state x_{c,l}=x_c^{(l)}, has
% the residual R' = T'*(d/dt x') - f'(x'), which must equal R_i in the row of
% the equation of x_i, and be zero in the rows of the added equations
% (d/dt)x_{c,l} = x_{c,l+1}. The same is checked for the constraint K*x=0
% (for each derivative), the outputs gz, gy, and the PDE state
% X_pde*(T*x), computed from the operators of the original and the expanded
% PIE, without using the converters between the PIE structures.
%
% DJ (CC), 11/09/2026: Initial coding;

clear stateNameGenerator
rng(1);
ncases = 0;
tol = 1e-8;

% % Case 1: 1D, two coupled PDE components with a nonlinear term, the
% % boundary condition coupling the components, so that T is not diagonal
pvar s t
dom = [0,1];
u1 = pde_var(s,dom);        u2 = pde_var(s,dom);
PDE = [diff(u1,t)==diff(u1,s,2)-u1+u2-u1*u2;
       diff(u2,t)==-diff(u2,s)-u2+u1+u2*diff(u1,s);
       subs(u1,s,dom(1))==0;    subs(u1,s,dom(2))==subs(u2,s,dom(1));
       subs(u2,s,dom(2))==0];
evalc('PIE = convert(PDE,{''silent'',''nonlinear''});');
for tl = {[2;1],[1;3],[3;2],[2;2]}
    ncases = ncases + check_expansion(PIE,tl{1},sprintf('1D coupled, orders %s',mat2str(tl{1})),tol);
end

% % Case 1b: custom names of the components, kept, with the new components
% % named x<k>_pie after them; names for which this is not possible (the new
% % names would sort before y_pie, or x3_pie is already used) give a renaming
% % of all components with default names, with a warning
for nm = {{'u_pie';'v_pie'},{'x1';'x1_2'},{'a_pie';'a_pie2'},{'y_pie';'z_pie'},{'x1_pie';'x3_pie'}}
    PIEr = rename_x(PIE,nm{1});
    for tl = {[2;1],[3;2]}
        ncases = ncases + check_expansion(PIEr,tl{1},sprintf('1D, names %s, orders %s',strjoin(nm{1},'/'),mat2str(tl{1})),tol);
    end
end

% % Case 1c: 10 derivatives of one component: the added names x02_pie, ...,
% % x11_pie are zero-padded, so they sort after the old name x01_pie, but not
% % after x1_pie (renamed, with a warning)
clear stateNameGenerator
u = pde_var(s,dom);
PDE1 = [diff(u,t)==diff(u,s,2)-u*diff(u,s);   subs(u,s,dom(1))==0;   subs(u,s,dom(2))==0];
evalc('PIE1c = convert(PDE1,{''silent'',''nonlinear''});');
ncases = ncases + check_expansion(PIE1c,11,'1D, one component, order 11',tol);
ncases = ncases + check_expansion(rename_x(PIE1c,{'x01_pie'}),11,'1D, one component x01_pie, order 11',tol);

% % Case 2: 1D, with a finite-dimensional component
clear stateNameGenerator
u1 = pde_var(s,dom);        u2 = pde_var(s,dom);        x3 = pde_var('state',1,[],[]);
PDE = [diff(u1,t)==diff(u1,s,2)+u2+x3*u1;
       diff(u2,t)==-diff(u2,s)-u2+subs(u1,s,dom(2))*u2;
       diff(x3,t)==-x3+subs(u2,s,dom(2))+x3^2;
       subs(u1,s,dom(1))==0;    subs(u1,s,dom(2))==0;
       subs(u2,s,dom(1))==x3];
evalc('PIE = convert(PDE,{''silent'',''nonlinear''});');
for tl = {[2;1;1],[1;1;2],[3;1;2],[2;2;2]}
    ncases = ncases + check_expansion(PIE,tl{1},sprintf('1D with ODE state, orders %s',mat2str(tl{1})),tol);
end

% % Case 3: 2D, two components
clear stateNameGenerator
pvar s1 s2
dom2 = [0,1;0,2];
u1 = pde_var([s1,s2],dom2);        u2 = pde_var([s1,s2],dom2);
PDE = [diff(u1,t)==diff(u1,s1,2)+diff(u1,s2,2)+u2-u1*u2;
       diff(u2,t)==diff(u2,s1,2)+diff(u2,s2,2)+u1*diff(u1,s1)];
for u = [u1,u2]
    PDE = [PDE; subs(u,s1,0)==0; subs(u,s1,1)==0; subs(u,s2,0)==0; subs(u,s2,2)==0];
end
evalc('PIE = convert(PDE,{''silent'',''nonlinear''});');
for tl = {[2;1],[3;2]}
    ncases = ncases + check_expansion(PIE,tl{1},sprintf('2D, orders %s',mat2str(tl{1})),tol);
end

% % Case 4: 1D with a constraint K*x=0 (surplus boundary condition on u1)
clear stateNameGenerator
u1 = pde_var(s,dom);        u2 = pde_var(s,dom);
PDE = [diff(u1,t)==diff(u1,s)-u1*u2+u2;
       diff(u2,t)==diff(u2,s,2)-u2;
       subs(u1,s,dom(1))==0;    subs(u1,s,dom(2))==0;
       subs(u2,s,dom(1))==0;    subs(u2,s,dom(2))==0];
evalc('PIE_K = convert(PDE,{''silent'',''nonlinear''});');
if size(PIE_K.K.ops,1)==0
    error("Test_npie_expand_tderivatives: the PIE of case 4 has no constraint, the test is vacuous.")
end
for tl = {[2;1],[1;2],[3;2],[2;2]}
    ncases = ncases + check_expansion(PIE_K,tl{1},sprintf('1D with constraint, orders %s',mat2str(tl{1})),tol);
end

% % Case 5: 3D, one component
clear stateNameGenerator
pvar s3
dom3 = [0,1;0,1;0,1];
u = pde_var([s1,s2,s3],dom3);
PDE = [diff(u,t)==diff(u,s1,2)+diff(u,s2,2)+diff(u,s3,2)-u*diff(u,s1)];
for ss = [s1,s2,s3]
    PDE = [PDE; subs(u,ss,0)==0; subs(u,ss,1)==0];
end
evalc('PIE = convert(PDE,{''silent'',''nonlinear''});');
for tl = {2,3}
    ncases = ncases + check_expansion(PIE,tl{1},sprintf('3D, orders %s',mat2str(tl{1})),tol);
end

% % Case 6: inputs and outputs, and errors
clear stateNameGenerator
u1 = pde_var(s,dom);        u2 = pde_var(s,dom);
w = pde_var('in',1);        z = pde_var('out',1);
PDE = [diff(u1,t)==diff(u1,s,2)+u2-u1*diff(u1,s)+w;
       diff(u2,t)==diff(u2,s,2)-u2;
       z==int(u1,s,dom)+subs(u2,s,dom(2));
       subs(u1,s,dom(1))==0;    subs(u1,s,dom(2))==0;
       subs(u2,s,dom(1))==w;    subs(u2,s,dom(2))==0];
evalc('PIE = convert(PDE,{''silent'',''nonlinear''});');
% input in BC of u2 (boundary input), so only u1 may have a higher-order derivative
ncases = ncases + check_expansion(PIE,[2;1],'inputs, orders [2;1]',tol);
try
    expand_tderivatives(PIE,[1;2],true);
    error("Test_npie_expand_tderivatives: no error for a boundary input in an equation of higher order.")
catch ME
    if ~strcmp(ME.identifier,'expand_tderivatives:boundaryInputs'),     rethrow(ME),        end
end
ncases = ncases+1;
% constraint involving components with different numbers of derivatives
PIEc = PIE_K;
PIEc.K{1,2} = PIEc.K.ops{1,1};
try
    expand_tderivatives(PIEc,[2;1],true);
    error("Test_npie_expand_tderivatives: no error for a constraint with components of different orders.")
catch ME
    if ~strcmp(ME.identifier,'expand_tderivatives:constraintOrders'),   rethrow(ME),        end
end
ncases = ncases+1;
% map with off-diagonal block
PIEc = PIE_K;
PIEc.pde2pie.x.C{1,2} = PIEc.pde2pie.x.C.ops{1,1};
try
    expand_tderivatives(PIEc,[2;1],true);
    error("Test_npie_expand_tderivatives: no error for a map pde2pie.x.C with off-diagonal blocks.")
catch ME
    if ~strcmp(ME.identifier,'expand_tderivatives:mapOffDiagonal'),     rethrow(ME),        end
end
ncases = ncases+1;
% invalid inputs, and no expansion
for bad = {[2;1;1],[0;1],[1.5;1]}
    try
        expand_tderivatives(PIE_K,bad{1},true);
        error("Test_npie_expand_tderivatives: no error for the orders %s.",mat2str(bad{1}))
    catch ME
        if ~contains(ME.message,'positive integers'),   rethrow(ME),    end
    end
end
try
    expand_tderivatives(PIE_K);
    error("Test_npie_expand_tderivatives: no error for missing orders.")
catch ME
    if ~contains(ME.message,'Insufficient'),        rethrow(ME),    end
end
[PIE_same,tab] = expand_tderivatives(PIE_K,[1;1],true);
if ~isequal(tab,[1,0;2,0]) || ~(PIE_same==PIE_K)
    error("Test_npie_expand_tderivatives: PIE changed without higher-order derivatives.")
end
ncases = ncases+5;

% % Case 7: linear PIEs, compared with the expansion of the pie_struct, that
% % has its own implementation
for ex = 1:2
    clear stateNameGenerator
    if ex==1
        u1 = pde_var(s,dom);        u2 = pde_var(s,dom);
        PDE = [diff(u1,t)==diff(u1,s,2)-u1+u2;
               diff(u2,t)==-diff(u2,s)-u2+u1;
               subs(u1,s,dom(1))==0;    subs(u1,s,dom(2))==subs(u2,s,dom(1));
               subs(u2,s,dom(2))==0];
        tls = {[2;1],[1;3],[3;2],[2;2]};
    else
        u1 = pde_var([s1,s2],dom2);        u2 = pde_var([s1,s2],dom2);
        PDE = [diff(u1,t)==diff(u1,s1,2)+diff(u1,s2,2)+u2;
               diff(u2,t)==diff(u2,s1,2)+diff(u2,s2,2)+u1];
        for u = [u1,u2]
            PDE = [PDE; subs(u,s1,0)==0; subs(u,s1,1)==0; subs(u,s2,0)==0; subs(u,s2,2)==0];
        end
        tls = {[2;1],[3;2]};
    end
    evalc('PIEl = convert(PDE,''silent'');');
    PIEn = pie2npie_struct(PIEl);
    for tl = tls
        ncases = ncases + check_expansion(PIEn,tl{1},sprintf('linear example %d, orders %s',ex,mat2str(tl{1})),tol);
        % The derivatives follow the component they derive from in the
        % pie_struct, and come after the old components in the npie_struct:
        % compare after matching the components (tab2 row k2 is tab1 row p(k2))
        [E1,tab1] = expand_tderivatives(PIEn,tl{1},true);
        [Ep,tab2] = expand_tderivatives(PIEl,tl{1},true);
        E2 = pie2npie_struct(Ep);
        [is_match,p] = ismember(tab2,tab1,'rows');
        if ~all(is_match) || numel(unique(p))~=numel(p) || size(tab1,1)~=size(tab2,1)
            error("Test_npie_expand_tderivatives failed: the components of the npie_struct and the pie_struct expansion differ (example %d, orders %s).",ex,mat2str(tl{1}))
        end
        compare_with_pie(E1,E2,p,sprintf('linear example %d, orders %s',ex,mat2str(tl{1})),tol);
        % The result also converts to a pie_struct, with a component for each state
        if numel(npie2pie_struct(E1).space.x.dims)~=size(tab1,1)
            error("Test_npie_expand_tderivatives failed: the expanded npie_struct does not convert to a pie_struct (example %d, orders %s).",ex,mat2str(tl{1}))
        end
        ncases = ncases+1;
    end
end

fprintf('Test_npie_expand_tderivatives passed (%d cases).\n',ncases)



function ok = check_expansion(PIE0,tdiff_list,label,tol)
% Expand PIE0 and compare with the semantics of the original PIE. Returns 1.
tdiff_list = tdiff_list(:);
lastwarn('');
tic;    [PIE1,tab] = expand_tderivatives(PIE0,tdiff_list,true);     t_exp = toc;
nb = numel(PIE0.space.x.dims);      nnew = size(tab,1);
par = tab(:,1);     ord = tab(:,2);
host = accumarray(par,(1:nnew)',[],@max);
L = accumarray(par,ord,[],@max);

% % Structure of the new PIE
err = @(msg) error("Test_npie_expand_tderivatives failed (%s): %s",label,msg);
if ~PIE1.is_initialized,                                err("the PIE is not initialized."),     end
% The old components come first, followed by the derivatives ordered by
% component and order of the derivative
tab_exp = [(1:nb)',zeros(nb,1)];
for c=1:nb
    tab_exp = [tab_exp; repmat(c,L(c),1),(1:L(c))'];
end
if ~isequal(tab,tab_exp)
    err("tdiff_tab is not the old components followed by the derivatives."),                        end
if any(L(tdiff_list>1)<0) || any(L<0) || max(L)~=max(tdiff_list)-1
    err("the highest number of derivatives is not max(p_i)-1."),                                    end
if ~isequal(PIE1.space.x.dims,PIE0.space.x.dims(par)) || ~isequal(PIE1.space.x.depmat,PIE0.space.x.depmat(par,:))
    err("the function spaces of the new components are not those of their parents."),               end
% Names: the old components keep theirs, and the new component k is named
% x<k>_pie, unless these are not unique and sorted, in which case all
% components get the default names (and a warning)
names0 = PIE0.space.x.varname;
names_exp = names0;
for k=nb+1:nnew
    names_exp{k,1} = sprintf(['x%0',num2str(numel(num2str(nnew))),'d_pie'],k);
end
is_default = ~issorted(names_exp) || numel(unique(names_exp))<nnew;
if is_default
    names_exp = arrayfun(@(k) sprintf(['x%0',num2str(numel(num2str(nnew))),'d_pie'],k),(1:nnew)','UniformOutput',false);
end
if ~isequal(PIE1.space.x.varname,names_exp)
    err("the names of the new components are not those expected."),                                 end
[~,wid] = lastwarn;
if is_default~=strcmp(wid,'expand_tderivatives:namesNotSorted')
    err("the warning for names that cannot be sorted is not issued exactly when expected."),        end
if ~isequal(PIE1.f.varname,PIE1.space.x.varname) || ~isequal(PIE1.f.varname,PIE1.x.varname)
    err("the variables of f are not named as the new components."),                                 end
if size(PIE1.T.ops,1)~=nnew || size(PIE1.T.ops,2)~=nnew || size(PIE1.f.C.ops,1)~=nnew
    err("T and f do not have a row and column for each new component."),                            end
% maps between PDE and PIE
p0 = PIE0.pde2pie.x;    p1 = PIE1.pde2pie.x;
if ~isequal(p1.idx,p0.idx(par,:)) || ~isequal(p1.diff,p0.diff(par,:)) || ~isequal(p1.loc,p0.loc(par,:)) || ~isequal(p1.delay,p0.delay(par))
    err("pde2pie.x fields idx, diff, loc, delay are not those of the parents."),                    end
if ~isequal(p1.tdiff,p0.tdiff(par)+ord)
    err("pde2pie.x.tdiff is not that of the parent plus the order of the derivative."),             end
for k=1:nnew
    if ~isequal(p1.C.ops{k,k},p0.C.ops{par(k),par(k)})
        err("pde2pie.x.C block is not that of the parent."),                                        end
end

% % Semantics: random values of the time derivatives of the state at time t0
P = max(tdiff_list);
X = cell(P+1,1);
sp = PIE0.space;
for k=0:P
    X{k+1} = cell(1,nb);
    for c=1:nb
        X{k+1}{c} = randpoly(sp.vars(sp.x.depmat(c,:)),sp.x.dims(c));
    end
end
Xn = cell(1,nnew);      Xn1 = cell(1,nnew);
for k=1:nnew
    Xn{k} = X{ord(k)+1}{par(k)};            % x_{c,l} = x_c^{(l)}
    Xn1{k} = X{ord(k)+2}{par(k)};           % d/dt x_{c,l} = x_c^{(l+1)}
end
off0 = [0;cumsum(PIE0.T.dims.out)];         off1 = [0;cumsum(PIE1.T.dims.out)];
rows = @(Y,off,i) Y(off(i)+1:off(i+1));
% Residuals of the original and the expanded PIE, row by row
Yf0 = apply_poly(PIE0.f,X{1});
Yf1 = apply_poly(PIE1.f,Xn);
Yt1 = apply_op(PIE1.T,Xn1);
R0all = cell(1,nb);
for i=1:nb
    Yt0 = apply_op(PIE0.T,X{tdiff_list(i)+1});
    R0 = rows(Yt0,off0,i) - rows(Yf0,off0,i);
    R0all{i} = R0;
    R1 = rows(Yt1,off1,host(i)) - rows(Yf1,off1,host(i));
    if ~is_zero_poly(R0-R1,tol)
        err(sprintf("residual of equation %d differs between the original and the expanded PIE.",i)),    end
end
for k=reshape(setdiff(1:nnew,host),1,[])
    R1 = rows(Yt1,off1,k) - rows(Yf1,off1,k);
    if ~is_zero_poly(R1,tol)
        err(sprintf("the equation of component %d is not (d/dt)x_{c,l} = x_{c,l+1}.",k)),               end
end
% Outputs
for fname = {'gz','gy'}
    if size(PIE0.(fname{1}),1)>0
        if ~is_zero_poly(apply_poly(PIE0.(fname{1}),X{1}) - apply_poly(PIE1.(fname{1}),Xn),tol)
            err(['the output ',fname{1},' differs.']),                                                  end
    end
end
% The PDE state, X_pde*(T*x), with the values of x
if size(PIE0.pie2pde.x.ops,1)>0
    split = @(Y,off) arrayfun(@(i) rows(Y,off,i),1:numel(off)-1,'UniformOutput',false);
    xpde0 = apply_op(PIE0.pie2pde.x,split(apply_op(PIE0.T,X{1}),off0));
    xpde1 = apply_op(PIE1.pie2pde.x,split(apply_op(PIE1.T,Xn),off1));
    if ~is_zero_poly(xpde0-xpde1,tol)
        err("the PDE state X_pde*(T*x) differs."),                                                      end
    % The dependent properties x_pde = X_pde*T*x and f_pde = X_pde*f of the new
    % PIE. For the original PIE, d/dt x_pde = X_pde*(T*x^{(1)}), which equals
    % f_pde(x) up to X_pde times the residual of the equations of order 1
    Rcell = R0all;
    for i=find(tdiff_list'>1),   Rcell{i} = 0*R0all{i};     end
    dxpde0 = apply_op(PIE0.pie2pde.x,split(apply_op(PIE0.T,X{2}),off0)) - apply_op(PIE0.pie2pde.x,Rcell);
    % X_pde applied to the values f(x), as the operator product X_pde*f of the
    % property f_pde is not supported for X_pde other than the identity if f has
    % monomials of degree>1 (composition of a single-factor operator with a
    % two-factor operator), also for the original PIE if X_pde is not the identity
    fpde1 = apply_op(PIE1.pie2pde.x,split(Yf1,off1));
    if ~is_zero_poly(fpde1-dxpde0,tol)
        err("the dynamics X_pde*f(x) of the PDE state differ."),                                        end
    try
        fp = PIE1.f_pde;
        if ~is_zero_poly(apply_poly(fp,Xn)-fpde1,tol)
            err("the property f_pde differs from X_pde*f(x)."),                                         end
    catch ME
        if ~contains(ME.message,'Composition of tensor-PI operators'),  rethrow(ME),    end
    end
    if ~is_zero_poly(apply_poly(PIE1.x_pde,Xn)-xpde0,tol)
        err("the PDE state x_pde = X_pde*T*x differs."),                                                end
end
% The constraint K*x=0 holds for each derivative
nK = size(PIE0.K.ops,1);
if nK>0
    Kn = size(PIE1.K.ops,1);
    offK0 = [0;cumsum(PIE0.K.dims.out)];    offK1 = [0;cumsum(PIE1.K.dims.out)];
    % Number of derivatives of each constraint: maximum of L over the
    % components that it involves
    Lcon = zeros(nK,1);
    for r=1:nK
        for c=1:nb
            [~,z] = clean_op(PIE0.K.ops{r,c},1e-12);
            if ~z,  Lcon(r) = max(Lcon(r),L(c));    end
        end
    end
    if Kn~=sum(Lcon+1)
        err("the number of constraints is not the number of constraints times derivatives."),           end
    YK1 = apply_op(PIE1.K,Xn);
    q = 0;
    for l=0:max(Lcon)
        YK0 = apply_op(PIE0.K,X{l+1});
        for r=find(Lcon>=l)'
            q = q+1;
            if ~is_zero_poly(rows(YK0,offK0,r)-rows(YK1,offK1,q),tol)
                err(sprintf("constraint %d for derivative %d differs.",r,l)),                           end
        end
    end
end
% Operators on the inputs: rows of the equations move to the host, the rows
% of the added equations are zero, and the operators between inputs and
% outputs are not changed
for fname = {'Tw','Tu','Bw','Bu'}
    P0 = PIE0.(fname{1});       P1 = PIE1.(fname{1});
    for i=1:nb
        if ~isequal(P1.ops(host(i),:),P0.ops(i,:))
            err(['the row of ',fname{1},' did not move to the new equation.']),                         end
    end
    for k=reshape(setdiff(1:nnew,host),1,[])
        for j=1:size(P1.ops,2)
            [~,z] = clean_op(P1.ops{k,j},1e-12);
            if ~z,  err(['the row of ',fname{1},' for an added equation is not zero.']),    end
        end
    end
end
for fname = {'Dzw','Dzu','Dyw','Dyu'}
    if ~isequal(PIE1.(fname{1}),PIE0.(fname{1}))
        err(['the operator ',fname{1},' changed.']),                                                   end
end
if ~isequal(PIE1.pie2pde.w,PIE0.pie2pde.w) || ~isequal(PIE1.pie2pde.z,PIE0.pie2pde.z)
    err("the maps for the inputs and outputs changed."),                                               end
fprintf('  %s: ok (%d -> %d components, %.2f s)\n',label,nb,nnew,t_exp);
ok = 1;
end


function compare_with_pie(E1,E2,p,label,tol)
% Compare the expanded PIE E1 with the PIE E2 obtained from the expansion of a
% pie_struct, in which component k of E2 is component p(k) of E1: the operators
% T, f, K, X_pde applied to the same random values of the state, and the maps
% from the PDE. The components of E1 are renamed by the caller as needed.
err = @(msg) error("Test_npie_expand_tderivatives failed (%s): %s",label,msg);
n = numel(p);
if ~isequal(E1.space.x.dims(p),E2.space.x.dims) || ~isequal(E1.space.x.depmat(p,:),E2.space.x.depmat)
    err("the function spaces of the components differ from those of the pie_struct expansion."),    end
sp = E2.space;
V = cell(1,n);      V1 = cell(1,n);
for k=1:n,      V{k} = randpoly(sp.vars(sp.x.depmat(k,:)),sp.x.dims(k));     end
V1(p) = V;
off2 = [0;cumsum(E2.T.dims.out)];       off1 = [0;cumsum(E1.T.dims.out)];
rows = @(Y,off,i) Y(off(i)+1:off(i+1));
split = @(Y,off) arrayfun(@(i) rows(Y,off,i),1:numel(off)-1,'UniformOutput',false);
Y2 = apply_op(E2.T,V);          Y1 = apply_op(E1.T,V1);
F2 = apply_poly(E2.f,V);        F1 = apply_poly(E1.f,V1);
for k=1:n
    if ~is_zero_poly(rows(Y2,off2,k)-rows(Y1,off1,p(k)),tol)
        err("T*x differs from that of the pie_struct expansion."),                                  end
    if ~is_zero_poly(rows(F2,off2,k)-rows(F1,off1,p(k)),tol)
        err("f(x) differs from that of the pie_struct expansion."),                                 end
end
if size(E2.K.ops,1)>0
    if size(E1.K.ops,1)~=size(E2.K.ops,1) || ~is_zero_poly(apply_op(E2.K,V)-apply_op(E1.K,V1),tol)
        err("K*x differs from that of the pie_struct expansion."),                                  end
end
if ~is_zero_poly(apply_op(E2.pie2pde.x,split(Y2,off2))-apply_op(E1.pie2pde.x,split(Y1,off1)),tol)
    err("the PDE state differs from that of the pie_struct expansion."),                            end
for fname = {'idx','diff','tdiff','loc','delay'}
    f2 = E2.pde2pie.x.(fname{1});       f1 = E1.pde2pie.x.(fname{1});
    if ~isequal(f1(p,:),f2)
        err(['pde2pie.x.',fname{1},' differs from that of the pie_struct expansion.']),             end
end
end


function PIE = rename_x(PIE,names)
% Rename the state components of the PIE, and the variables of f, gz, gy
% (the new names must be sorted, as the old ones)
old = PIE.space.x.varname;
for fname = {'f','gz','gy'}
    P = PIE.(fname{1});
    [~,loc] = ismember(P.varname,old);
    P.varname = names(loc);
    PIE.(fname{1}) = P;
end
PIE.space.x.varname = names;
PIE = initialize(PIE);
end


function p = randpoly(names,d)
% Random d x 1 polynomial in the variables with the given names
c = cell(d,1);
for r=1:d
    pr = polynomial(randn);
    for k=1:numel(names)
        v = pvar(names{k});
        pr = pr + randn*v + randn*v^2;
    end
    if numel(names)>1
        pr = pr + randn*pvar(names{1})*pvar(names{2});
    end
    c{r} = pr;
end
p = vertcat(c{:});
end


function z = is_zero_poly(D,tol)
% True if the polynomial array D has all coefficients below tol
z = isempty(D.coefficient) || max(abs(full(D.coefficient(:))))<tol;
end

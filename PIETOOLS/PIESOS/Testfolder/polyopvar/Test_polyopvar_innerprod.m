% Test_polyopvar_innerprod.m -- verify innerprod.m, V = <u1,P*u2>,
% against an independent ground truth: evaluate u1, u2 and V at random
% polynomial state values (apply_poly / a local evaluator for the
% integral-operator coefficients of V) and compare with the integral of
% u1(x)'*P*u2(x) over the spatial domain, computed directly with
% polynomial integration. The ground truth does not use innerprod or any
% inverse of it.
%
% Covers: scalar and vector-valued (multi-block) arguments, 1-3 spatial
% variables, nonuniform domains, degree 1-3 arguments, stacked arguments
% of mixed monomial structure, weights of type double / polynomial / cell
% / dpvar (rectangular, nonsymmetric), and the one-argument form.
% innerprod(u1,P) and innerprod(u1,P,[]) are in Part 8, all-zero weights
% in Part 9, weights depending on space in Part 10.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - 'ntrial' (random trials per part) can be preset in the workspace;
%   setappdata(0,'innerprod_collect',true) logs all failures instead of
%   stopping at the first (counted in appdata 'innerprod_nfail').

tol = 1e-8;
if ~exist('ntrial','var')       % number of random trials per part, can be preset
    ntrial = 5;
end
setappdata(0,'innerprod_nfail',0);

%% Part 1: scalar states, scalar weight, 1-3 spatial variables
for nvars = 1:3
    for trial = 1:ntrial
        rng(20000+100*nvars+trial);
        vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
        dom = sort(randn(nvars,2),2);
        x = polyopvar('x',svars(nvars),dom);
        y = polyopvar('y',svars(nvars),dom);
        vals = struct('x',rand_poly([1,1],vars,2,4),'y',rand_poly([1,1],vars,2,4));
        a = randn;

        check(@() innerprod(x,x),     vals,x,x,1,dom,vars,tol,'scalar <x,x>');
        check(@() innerprod(x,y),     vals,x,y,1,dom,vars,tol,'scalar <x,y>');
        check(@() innerprod(x,a,y),   vals,x,y,a,dom,vars,tol,'scalar <x,a*y>');
        check(@() innerprod(x),       vals,x,x,1,dom,vars,tol,'one-argument <x>');
        check(@() innerprod(x,a,x),   vals,x,x,a,dom,vars,tol,'scalar <x,a*x>');
    end
end
disp('Test_polyopvar_innerprod Part 1 (scalar states, scalar weight) passed.')

%% Part 2: higher-degree arguments, scalar weight
for nvars = 1:2
    for trial = 1:ntrial
        rng(21000+100*nvars+trial);
        vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
        dom = sort(randn(nvars,2),2);
        x = polyopvar('x',svars(nvars),dom);
        y = polyopvar('y',svars(nvars),dom);
        vals = struct('x',rand_poly([1,1],vars,2,4),'y',rand_poly([1,1],vars,2,4));
        a = randn;

        check(@() innerprod(x,y.*y),      vals,x,y.*y,1,dom,vars,tol,'<x,y^2>');
        check(@() innerprod(x.*y,y.*y),   vals,x.*y,y.*y,1,dom,vars,tol,'<xy,y^2>');
        check(@() innerprod(x.*x,a,x.*y), vals,x.*x,x.*y,a,dom,vars,tol,'<x^2,a*xy>');
        check(@() innerprod(x+y,x-y),     vals,x+y,x-y,1,dom,vars,tol,'<x+y,x-y>');
        check(@() innerprod(x.*x.*y,x),   vals,x.*x.*y,x,1,dom,vars,tol,'<x^2 y,x>');
    end
end
disp('Test_polyopvar_innerprod Part 2 (higher degree) passed.')

%% Part 3: vector-valued, multi-block states with matrix weights
for trial = 1:2*ntrial
    rng(22000+trial);
    nvars = randi([1,3]);
    vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
    dom = sort(randn(nvars,2),2);
    p1 = randi([1,3]);      p2 = randi([1,3]);
    [x,~,dx] = rand_polyopvar(p1,3,vars,dom);
    [y,~,dy] = rand_polyopvar(p2,3,vars,dom);
    y.varname = cellfun(@(s)['y',s(2:end)],y.varname,'UniformOutput',false);
    n1 = sum(dx);   n2 = sum(dy);

    % State values by name
    vals = struct();
    for i=1:p1,  vals.(['x',num2str(i)]) = rand_poly([dx(i),1],vars,2,4);   end
    for i=1:p2,  vals.(['y',num2str(i)]) = rand_poly([dy(i),1],vars,2,4);   end
    vals_str = sprintf('(trial=%d,p1=%d,p2=%d,dx=%s,dy=%s,nvars=%d)',trial,p1,p2,mat2str(dx'),mat2str(dy'),nvars);

    % Rectangular weights as double, and as (constant) 'polynomial'
    P = randn(n1,n2);
    check(@() innerprod(x,P,y),vals,x,y,P,dom,vars,tol,['double rectangular weight ',vals_str]);
    Ppoly = polynomial(P);
    check(@() innerprod(x,Ppoly,y),vals,x,y,P,dom,vars,tol,['polynomial weight ',vals_str]);

    % Cell-partitioned weight
    Pcell = mat2cell(P,dx(:)',dy(:)');
    Pcell = reshape(Pcell,numel(dx),numel(dy));
    check(@() innerprod(x,Pcell,y),vals,x,y,P,dom,vars,tol,['cell weight ',vals_str]);

    % Same state vector on both sides, nonsymmetric weight
    if n1>=1
        Q = randn(n1,n1);
        check(@() innerprod(x,Q,x),vals,x,x,Q,dom,vars,tol,['<x,Qx> nonsymmetric ',vals_str]);
    end
    % One-argument form: <x,x> for multi-block x (symmetric branch)
    check(@() innerprod(x),vals,x,x,eye(n1),dom,vars,tol,['one-argument multi-block ',vals_str]);
    check(@() innerprod(x,y),vals,x,y,eye(n1),dom,vars,tol,['two-argument identity ',vals_str],isequal(dx,dy));
end
disp('Test_polyopvar_innerprod Part 3 (vector-valued, multi-block) passed.')

%% Part 4: vector-valued, higher degree
for trial = 1:ntrial
    rng(23000+trial);
    nvars = randi([1,2]);
    vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
    dom = sort(randn(nvars,2),2);
    p1 = randi([1,2]);
    [x,~,dx] = rand_polyopvar(p1,3,vars,dom);
    n1 = sum(dx);
    vals = struct();
    for i=1:p1,  vals.(['x',num2str(i)]) = rand_poly([dx(i),1],vars,2,4);   end
    vals_str = sprintf('(trial=%d,p1=%d,dx=%s,nvars=%d)',trial,p1,mat2str(dx'),nvars);
    Q = randn(n1,n1);
    check(@() innerprod(x,Q,x.*x),vals,x,x.*x,Q,dom,vars,tol,['<x,Qx^2> ',vals_str]);
    check(@() innerprod(x.*x,Q,x.*x),vals,x.*x,x.*x,Q,dom,vars,tol,['<x^2,Qx^2> ',vals_str]);
end
disp('Test_polyopvar_innerprod Part 4 (vector-valued, higher degree) passed.')

%% Part 5: stacked arguments of mixed monomial structure, scalar states
% Exercises products of Z1 (nZ1 monomials) and Z2 (nZ2 monomials) with
% nZ1~=nZ2, with a single nonzero in one of the Z's, and with a single
% multi-factor monomial (Z a row vector)
for nvars = 1:2
    for trial = 1:ntrial
        rng(25000+100*nvars+trial);
        vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
        dom = sort(randn(nvars,2),2);
        x = polyopvar('x',svars(nvars),dom);
        y = polyopvar('y',svars(nvars),dom);
        vals = struct('x',rand_poly([1,1],vars,2,4),'y',rand_poly([1,1],vars,2,4));
        xy = [x;y];     P2 = randn(2,2);   p1 = randn(1,2);   p2 = randn(2,1);

        check(@() innerprod(xy,P2,xy),            vals,xy,xy,P2,dom,vars,tol,'<[x;y],P[x;y]>');
        check(@() innerprod(xy,P2,[y;x]),         vals,xy,[y;x],P2,dom,vars,tol,'<[x;y],P[y;x]>');
        check(@() innerprod(xy,p2,x),             vals,xy,x,p2,dom,vars,tol,'<[x;y],p*x> (single nonzero in Z2)');
        check(@() innerprod(x,p1,xy),             vals,x,xy,p1,dom,vars,tol,'<x,p*[x;y]> (single nonzero in Z1)');
        check(@() innerprod(x.*y,p1,xy),          vals,x.*y,xy,p1,dom,vars,tol,'<xy,p*[x;y]> (Z1 a row vector)');
        check(@() innerprod(xy,p2,x.*y),          vals,xy,x.*y,p2,dom,vars,tol,'<[x;y],p*xy> (Z2 a row vector)');
        xyz = [x; x.*y];
        check(@() innerprod(xyz,P2,xyz),          vals,xyz,xyz,P2,dom,vars,tol,'<[x;xy],P[x;xy]> (mixed degree)');
        check(@() innerprod(xyz,P2,xy),           vals,xyz,xy,P2,dom,vars,tol,'<[x;xy],P[x;y]> (mixed degree)');
    end
end
disp('Test_polyopvar_innerprod Part 5 (stacked arguments) passed.')

%% Part 6: non-polyopvar arguments
rng(24000);
pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
vals = struct('x',rand_poly([1,1],{'s1','s2'},2,4));
a = randn;
V = innerprod(a,1,x);
chk = apply_poly(V,{vals.x}) - a*vals.x;
if max(max(abs(cleanpoly(chk,1e-10).C)))>tol
    error("Test_polyopvar_innerprod failed: double-left argument.")
end
V = innerprod(x,1,a);
chk = apply_poly(V,{vals.x}) - a*vals.x;
if max(max(abs(cleanpoly(chk,1e-10).C)))>tol
    error("Test_polyopvar_innerprod failed: double-right argument.")
end
disp('Test_polyopvar_innerprod Part 6 (constant arguments) passed.')

%% Part 7: dpvar weights, linear in decision variables
% V is linear in P, so for P = P0 + d1*P1 + d2*P2 the value of V must equal
% I0 + d1*I1 + d2*I2, with Ik the integral of u1'*Pk*u2 for numeric Pk.
dpvar d1 d2
for trial = 1:ntrial
    rng(26000+trial);
    nvars = randi([1,2]);
    vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
    dom = sort(randn(nvars,2),2);
    p1 = randi([1,3]);      p2 = randi([1,3]);
    [x,~,dx] = rand_polyopvar(p1,3,vars,dom);
    [y,~,dy] = rand_polyopvar(p2,3,vars,dom);
    y.varname = cellfun(@(s)['y',s(2:end)],y.varname,'UniformOutput',false);
    n1 = sum(dx);   n2 = sum(dy);
    vals = struct();
    for i=1:p1,  vals.(['x',num2str(i)]) = rand_poly([dx(i),1],vars,2,4);   end
    for i=1:p2,  vals.(['y',num2str(i)]) = rand_poly([dy(i),1],vars,2,4);   end
    u1v = eval_arg(x,vals);     u2v = eval_arg(y,vals);
    P0 = randn(n1,n2);   P1 = randn(n1,n2);   P2 = randn(n1,n2);
    Pd = P0 + d1*P1 + d2*P2;

    Vv = eval_V(innerprod(x,Pd,y),vals);
    ref = integrate(u1v'*P0*u2v,vars,dom) + d1*integrate(u1v'*P1*u2v,vars,dom) ...
            + d2*integrate(u1v'*P2*u2v,vars,dom);
    D = Vv-ref;
    err = full(D.C);
    if max(abs(err(:)))>tol*max(1,max(abs(full(ref.C(:)))))
        error("Test_polyopvar_innerprod failed: dpvar weight (trial=%d,dx=%s,dy=%s,nvars=%d).", ...
            trial,mat2str(dx'),mat2str(dy'),nvars)
    end

    % Scalar dpvar weight, multi-block arguments of matching block sizes
    y2 = x;     y2.varname = cellfun(@(s)['y',s(2:end)],x.varname,'UniformOutput',false);
    for i=1:p1,  vals.(['y',num2str(i)]) = rand_poly([dx(i),1],vars,2,4);   end
    for u2 = {x,y2}
        u2 = u2{1};     u2v = eval_arg(u2,vals);
        Vv = eval_V(innerprod(x,d1,u2),vals);
        D = Vv - d1*integrate(u1v'*u2v,vars,dom);
        err = full(D.C);
        if max(abs(err(:)))>tol*max(1,abs(integrate(u1v'*u2v,vars,dom)))
            error("Test_polyopvar_innerprod failed: scalar dpvar weight (trial=%d,dx=%s,nvars=%d).", ...
                trial,mat2str(dx'),nvars)
        end
    end
end
disp('Test_polyopvar_innerprod Part 7 (dpvar weights) passed.')

%% Part 8: u2 omitted or empty, innerprod(u1,P) = innerprod(u1,P,[]) = <u1,P*u1>
for trial = 1:ntrial
    rng(27000+trial);
    nvars = randi([1,2]);
    vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
    dom = sort(randn(nvars,2),2);
    p1 = randi([1,3]);
    [x,~,dx] = rand_polyopvar(p1,3,vars,dom);
    n1 = sum(dx);
    vals = struct();
    for i=1:p1,  vals.(['x',num2str(i)]) = rand_poly([dx(i),1],vars,2,4);   end
    vals_str = sprintf('(trial=%d,p1=%d,dx=%s,nvars=%d)',trial,p1,mat2str(dx'),nvars);
    Q = randn(n1,n1);       Qs = Q+Q';      a = randn;      Qcell = mat2cell(Q,dx(:)',dx(:)');

    check(@() innerprod(x,Q),         vals,x,x,Q,dom,vars,tol,['<x,Qx> via innerprod(x,Q) ',vals_str]);
    check(@() innerprod(x,Q,[]),      vals,x,x,Q,dom,vars,tol,['<x,Qx> via innerprod(x,Q,[]) ',vals_str]);
    check(@() innerprod(x,Qs),        vals,x,x,Qs,dom,vars,tol,['<x,Qsx> symmetric Qs ',vals_str]);
    check(@() innerprod(x,Qs,[]),     vals,x,x,Qs,dom,vars,tol,['<x,Qsx> symmetric Qs, u2=[] ',vals_str]);
    check(@() innerprod(x,Qcell),     vals,x,x,Q,dom,vars,tol,['<x,Qx> cell weight ',vals_str]);
    check(@() innerprod(x.*x,Qs),     vals,x.*x,x.*x,Qs,dom,vars,tol,['<x^2,Qsx^2> symmetric Qs ',vals_str]);
    check(@() innerprod(x,a,[]),      vals,x,x,a*eye(n1),dom,vars,tol,['<x,ax> scalar weight, u2=[] ',vals_str]);
end
disp('Test_polyopvar_innerprod Part 8 (u2 omitted or empty) passed.')

%% Part 9: all-zero weight returns the zero polynomial
for trial = 1:ntrial
    rng(28000+trial);
    nvars = randi([1,2]);
    vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
    dom = sort(randn(nvars,2),2);
    p1 = randi([1,3]);      p2 = randi([1,3]);
    [x,~,dx] = rand_polyopvar(p1,3,vars,dom);
    [y,~,dy] = rand_polyopvar(p2,3,vars,dom);
    y.varname = cellfun(@(s)['y',s(2:end)],y.varname,'UniformOutput',false);
    n1 = sum(dx);   n2 = sum(dy);
    vals = struct();
    for i=1:p1,  vals.(['x',num2str(i)]) = rand_poly([dx(i),1],vars,2,4);   end
    for i=1:p2,  vals.(['y',num2str(i)]) = rand_poly([dy(i),1],vars,2,4);   end
    vals_str = sprintf('(trial=%d,dx=%s,dy=%s,nvars=%d)',trial,mat2str(dx'),mat2str(dy'),nvars);

    check(@() innerprod(x,0,x),                 vals,x,x,0*eye(n1),dom,vars,tol,['scalar zero weight ',vals_str]);
    check(@() innerprod(x,zeros(n1,n2),y),      vals,x,y,zeros(n1,n2),dom,vars,tol,['zero matrix weight ',vals_str]);
    check(@() innerprod(x,zeros(n1,n1)),        vals,x,x,zeros(n1,n1),dom,vars,tol,['zero weight, u2 omitted ',vals_str]);
    check(@() innerprod(x,mat2cell(zeros(n1,n2),dx(:)',dy(:)'),y),vals,x,y,zeros(n1,n2),dom,vars,tol,['zero cell weight ',vals_str]);
    % The zero polynomial can be combined with other polynomials
    check(@() innerprod(x,0,x) + innerprod(x,x), vals,x,x,eye(n1),dom,vars,tol,['zero + <x,x> ',vals_str]);
end
disp('Test_polyopvar_innerprod Part 9 (zero weight) passed.')

%% Part 10: weights depending on the spatial variables
dpvar d1 d2
for trial = 1:ntrial
    rng(29000+trial);
    nvars = randi([1,3]);
    vars = arrayfun(@(i) ['s',num2str(i)],1:nvars,'UniformOutput',false);
    dom = sort(randn(nvars,2),2);
    p1 = randi([1,3]);      p2 = randi([1,3]);
    [x,~,dx] = rand_polyopvar(p1,3,vars,dom);
    [y,~,dy] = rand_polyopvar(p2,3,vars,dom);
    y.varname = cellfun(@(s)['y',s(2:end)],y.varname,'UniformOutput',false);
    n1 = sum(dx);   n2 = sum(dy);
    vals = struct();
    for i=1:p1,  vals.(['x',num2str(i)]) = rand_poly([dx(i),1],vars,2,4);   end
    for i=1:p2,  vals.(['y',num2str(i)]) = rand_poly([dy(i),1],vars,2,4);   end
    vals_str = sprintf('(trial=%d,dx=%s,dy=%s,nvars=%d)',trial,mat2str(dx'),mat2str(dy'),nvars);

    % Rectangular weights depending on all or some of the variables
    W = rand_poly([n1,n2],vars,2,4);
    check(@() innerprod(x,W,y),   vals,x,y,W,dom,vars,tol,['spatial weight, all variables ',vals_str]);
    Wp = rand_poly([n1,n2],vars(1),2,3);
    check(@() innerprod(x,Wp,y),  vals,x,y,Wp,dom,vars,tol,['spatial weight, first variable only ',vals_str]);
    Wc = polynomial(randn(n1,n2));
    check(@() innerprod(x,Wc,y),  vals,x,y,Wc,dom,vars,tol,['spatial weight, constant polynomial ',vals_str]);
    Q = rand_poly([n1,n1],vars,2,4);
    check(@() innerprod(x,Q),     vals,x,x,Q,dom,vars,tol,['spatial weight, u2 omitted ',vals_str]);
    check(@() innerprod(x,Q,x.*x),vals,x,x.*x,Q,dom,vars,tol,['spatial weight, degree 1 x 2 ',vals_str]);
    check(@() innerprod(x.*x,Q,x.*x),vals,x.*x,x.*x,Q,dom,vars,tol,['spatial weight, degree 2 x 2 ',vals_str]);

    % Scalar weight (all blocks)
    w = rand_poly([1,1],vars,2,4);
    check(@() innerprod(x,w,x),   vals,x,x,w*eye(n1),dom,vars,tol,['scalar spatial weight ',vals_str]);

    % Weight on a variable outside the output space is rejected
    pvar s9
    err_caught = false;
    try
        innerprod(x,s9*ones(n1,n2),y);
    catch
        err_caught = true;
    end
    if ~err_caught
        error("Test_polyopvar_innerprod failed: weight depending on a variable outside the output space was not rejected %s.",vals_str)
    end

    % dpvar weight depending on space: linear in the decision variables
    P0 = rand_poly([n1,n2],vars,2,4);   P1 = rand_poly([n1,n2],vars,2,4);   P2 = rand_poly([n1,n2],vars,2,4);
    Pd = P0 + d1*P1 + d2*P2;
    u1v = eval_arg(x,vals);     u2v = eval_arg(y,vals);
    Vv = eval_V(innerprod(x,Pd,y),vals);
    ref = integrate(u1v'*P0*u2v,vars,dom) + d1*integrate(u1v'*P1*u2v,vars,dom) ...
            + d2*integrate(u1v'*P2*u2v,vars,dom);
    D = Vv-ref;     err = full(D.C);
    if max(abs(err(:)))>tol*max(1,max(abs(full(ref.C(:)))))
        error("Test_polyopvar_innerprod failed: spatial dpvar weight %s.",vals_str)
    end
end
disp('Test_polyopvar_innerprod Part 10 (spatially varying weight) passed.')

nfail = getappdata(0,'innerprod_nfail');
if nfail>0
    error("Test_polyopvar_innerprod failed: %d case(s) failed (see FAIL/ERROR lines above).",nfail)
end
disp('Test_polyopvar_innerprod passed.')


%% Local functions
function pv = svars(nvars)
% Row of the first nvars spatial pvar objects s1,s2,s3
pvar s1 s2 s3
pv = [s1,s2,s3];
pv = pv(1:nvars);
end

function check(Vfun,vals,u1,u2,P,dom,vars,tol,msg,do_run)
% Compare V(x) with int u1(x)'*P*u2(x) ds, where V = Vfun() is the
% innerprod output. Set setappdata(0,'innerprod_collect',true) to log all
% failures instead of stopping at the first.
if nargin>=10 && ~do_run
    return
end
collect = isappdata(0,'innerprod_collect') && getappdata(0,'innerprod_collect');

% Ground truth: evaluate the arguments, multiply, integrate
u1v = eval_arg(u1,vals);
u2v = eval_arg(u2,vals);
if iscell(P)
    P = cell2mat(P);
end
integrand = u1v'*P*u2v;
ref = integrate(integrand,vars,dom);

% Evaluate V
try
    V = Vfun();
    Vv = double(eval_V(V,vals));
catch ME
    if collect
        fprintf(1,"ERROR: %s -- %s\n",msg,ME.message);
        setappdata(0,'innerprod_nfail',getappdata(0,'innerprod_nfail')+1);
        return
    end
    error("Test_polyopvar_innerprod failed: %s -- %s",msg,ME.message)
end
if abs(Vv-ref)>tol*max(1,abs(ref))
    if collect
        fprintf(1,"FAIL: %s (value %g, expected %g).\n",msg,Vv,ref);
        setappdata(0,'innerprod_nfail',getappdata(0,'innerprod_nfail')+1);
        return
    end
    error("Test_polyopvar_innerprod failed: %s (value %g, expected %g).",msg,Vv,ref)
end
end

function y = eval_arg(u,vals)
% Evaluate a polyopvar vector at the named state values
xv = cell(numel(u.varname),1);
for i=1:numel(xv)
    xv{i} = vals.(u.varname{i});
end
y = apply_poly(u,xv);
end

function ref = integrate(p,vars,dom)
% Integrate polynomial p over the box dom in the variables vars
sv = svars(numel(vars));
for i=1:numel(vars)
    v = sv(i);
    p = int_simple(p,v,dom(i,1),dom(i,2));
end
ref = double(p);
end

function Y = eval_V(V,vals)
% Evaluate a 1x1 polyopvar whose coefficients are integral operators
% acting on d-fold tensor products of the states
Y = 0;
nZ = size(V.Z,1);
Zf = full(V.Z);
for k=1:nZ
    op = V.C.ops{1,k};
    if isempty(op)
        continue
    end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    xv = cell(1,d);
    for j=1:d
        xv{j} = vals.(V.varname{idcs(j)});
    end
    if d==0
        Y = Y + apply_op(op,{1});
    elseif isa(op,'intop')
        % Variables of the jth factor are renamed to the jth dummy slot
        v = 1;
        for j=1:d
            xj = xv{j};
            if d>1
                for m=1:numel(op.vars.in)
                    isv = strcmp(xj.varname,op.vars.in{m});
                    xj.varname(isv) = {[op.vars.in{m},'_',num2str(j)]};
                end
            end
            v = kron(v,xj);
        end
        Y = Y + apply_op(op,v);
    else
        Y = Y + apply_op(op,xv);
    end
end
end







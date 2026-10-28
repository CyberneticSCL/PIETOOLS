% Test_polyopvar_class.m -- tests of the methods of the 'polyopvar' class
% that have no test of their own: num_factors, ispvar, common_states,
% common_vars, apply_poly, and plus/minus/uminus/mtimes/vertcat/horzcat
% with constants, numbers and mixed operands, including their error
% conditions. ('subs' is not tested here.)
%
% Expected values are computed independently of the class: a distributed
% polynomial P is evaluated at polynomial state values by eval_polyopvar
% (Testfolder/polyopvar, which does not use subs), and compared to the
% same expression in the 'polynomial' values themselves, e.g.
% (x.*y+2*z)(xv,yv,zv) = xv*yv+2*zv.
% Parts:
%   (1) num_factors: nfac, nops and is_empty for states, products,
%       constants, known factors, native blocks;
%   (2) ispvar: states, stacked states in any order, repeated states,
%       polynomials that are not states;
%   (3) common_states, common_vars: merged state and spatial variables,
%       depmat, Z, errors;
%   (4) plus, minus, uminus, mtimes with numbers, polynomials and zero,
%       and their errors;
%   (5) vertcat and horzcat with states, numbers and empty operands, and
%       their errors;
%   (6) apply_poly against eval_polyopvar, for scalar and vector-valued
%       states, and its errors.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - Adding a number to a function of spatial variables (x+a, a+x, x-a,
%   a-x) errors in plus.m, as a constant has no output variables.
%   chk_or_warn checks the value if the operation is supported, and
%   otherwise only warns.
% - Needs eval_polyopvar and rand_polyopvar (Testfolder/polyopvar) on the
%   path.

tol = 1e-8;
pvar s1 s2
polyopvar x(s1,s2) y(s1,s2) z(s1,s2)
polyopvar u(s1) w(s2)
vn = {'s1','s2'};
rp = @(vars,deg) rand_poly([1,1],vars,deg,5);       % random scalar polynomial in the variables vars (cellstr)
rv = @() struct('x',rp(vn,2),'y',rp(vn,2),'z',rp(vn,1));
ev = @(P,v) eval_polyopvar(P,{v.x;v.y;v.z});        % value of P in the states x,y,z (alphabetical order)

%% (1) num_factors
[nfac,nops,is_empty] = num_factors(x);
chk(isequal(full(nfac),1) && isequal(full(nops),1) && ~is_empty,'(1) state');
[nfac,nops,is_empty] = num_factors(x.*y);
chk(isequal(full(nfac),2) && isequal(full(nops),2) && ~is_empty,'(1) product of 2 states');
chk(isequal(full(num_factors(x.*y.*z)),3),'(1) product of 3 states');
chk(isequal(full(num_factors((x+y).*z)),[2;2]),'(1) (x+y)*z');
chk(isequal(full(num_factors(x.*x)),2),'(1) x*x');
c = constant2polyopvar(2);
[nfac,nops] = num_factors(c);
chk(isequal(full(nfac),1) && isequal(full(nops),1),'(1) constant: its only factor is a multiplier from R');
[nfac,nops] = num_factors(c.*x);
chk(isequal(full(nfac),2) && isequal(full(nops),2),'(1) c*x');
[nfac,nops] = num_factors(x.*c);
chk(isequal(full(nfac),2) && isequal(full(nops),2),'(1) x*c');
chk(isequal(full(num_factors([x;y])),[1;1]),'(1) stacked states');
% a function given by a block of type polynomial/double/dpvar acting on R: one factor, nops=0 (not a tensopvar)
vars_s = struct('in',{cell(1,0)},'out',{{'s1'}});       dom_s = struct('in',zeros(0,2),'out',[0,1]);
dpvar da
for k=1:3
    switch k
        case 1, blk = 1+s1;
        case 2, blk = 3;
        case 3, blk = da*s1+1;
    end
    f0 = mkdeg0(blk,vars_s,dom_s);
    [nfac,nops,is_empty] = num_factors(f0);
    chk(isequal(full(nfac),1) && isequal(full(nops),0) && ~is_empty,sprintf('(1) native block, type %d',k));
end

%% (2) ispvar
[tf,vname,vdims,vvars,vdom] = ispvar(x);
chk(tf && isequal(vname,{'x'}) && isequal(vdims,1) && isequal(vvars,{{'s1','s2'}}) && isequal(vdom,{[0,1;0,1]}),'(2) one state');
[tf,vname] = ispvar([z;x;y]);
chk(tf && isequal(vname,{'z';'x';'y'}),'(2) stacked states keep the order of the rows, not the sorted names');
[tf,vname] = ispvar([x;x;y]);
chk(tf && isequal(vname,{'x';'x';'y'}),'(2) repeated state');
[tf,vname,vdims,vvars,vdom] = ispvar([u;w]);
chk(tf && isequal(vname,{'u';'w'}) && isequal(vvars,{{'s1'};{'s2'}}) && isequal(vdom,{[0,1];[0,1]}),'(2) states with different spatial variables');
chk(iscolumn(vname),'(2) varname is a column');
[xv,xnm,xdm] = rand_polyopvar(3,[2;1;3],2);          % vector-valued states: one block per state
[tf,vname,vdims] = ispvar(xv);
chk(tf && isequal(vname,xnm) && isequal(vdims,xdm(:)),'(2) vector-valued states');
[tf,vname,vdims] = ispvar([xv;xv]);
chk(tf && isequal(vname,[xnm;xnm]) && isequal(vdims,[xdm(:);xdm(:)]),'(2) stacked vector-valued states');
bad = {x+y, 'x+y'; 2*x, '2*x'; x.*y, 'x*y'; constant2polyopvar(1), 'constant'; [x;2*y], '[x;2y]'; [x;x.*y], '[x;xy]'};
for k=1:size(bad,1)
    [tf,vname,vdims,vvars,vdom] = ispvar(bad{k,1});
    chk(~tf && isequal(size(vname),[0,1]) && isempty(vdims) && isempty(vvars) && isempty(vdom),['(2) not a state: ' bad{k,2}]);
end

%% (3) common_states and common_vars
% different states, same spatial variables
[A,B] = common_states(x,y);
chk(isequal(A.varname,{'x';'y'}) && isequal(B.varname,{'x';'y'}) && isequal(full(A.Z),1) && isequal(full(B.Z),2),'(3) merged states, Z');
chk(isequal(A.vars.in,vn) && isequal(B.vars.in,vn) && isequal(A.depmat.in,true(2,2)) && isequal(B.depmat.in,true(2,2)),'(3) merged spatial variables');
chk(isequal(A.dims.in,[1;1]) && isequal(A.vars.out,x.vars.out) && isequal(B.dom.out,y.dom.out),'(3) dims.in, output space untouched');
v = rv();
vals = struct('x',v.x,'y',v.y,'z',v.z);
chk_close(eval_polyopvar(A,{v.x;v.y}),v.x,tol,'(3) value of A after merging states');
chk_close(eval_polyopvar(B,{v.x;v.y}),v.y,tol,'(3) value of B after merging states');
% overlapping states: y is in both
P1 = x+y;       P2 = y.*z;
[A,B] = common_states(P1,P2);
chk(isequal(A.varname,{'x';'y';'z'}) && isequal(B.varname,{'x';'y';'z'}),'(3) merged state names');
chk_close(eval_polyopvar(A,{v.x;v.y;v.z}),v.x+v.y,tol,'(3) overlapping states, A');
chk_close(eval_polyopvar(B,{v.x;v.y;v.z}),v.y*v.z,tol,'(3) overlapping states, B');
[A,B] = common_states(x,x);
chk(isequal(A.varname,{'x'}) && isequal(full(A.Z),1) && isequal(full(B.Z),1),'(3) identical states');
% states in different spatial variables
[A,B] = common_states(u,w);
chk(isequal(A.vars.in,vn) && isequal(A.depmat.in,[true,false;false,true]) && isequal(B.depmat.in,[true,false;false,true]) && isequal(A.dom.in,[0,1;0,1]),'(3) spatial variables s1 and s2 merged');
% a varname stored as a row is returned as a column
Pr = [x;y];     Pr2 = Pr;       Pr2.varname = Pr.varname';
[A,B] = common_states(Pr2,Pr);
chk(iscolumn(A.varname) && isequal(A.varname,B.varname),'(3) row varname -> column');
% errors
xd = rand_polyopvar(2,[2;2],vn);                       % states of size 2, named as [x;y]
xd.varname = {'x';'y'};
expect_error(@() common_states([x;y],xd),'Dimensions of state vectors');
polyopvar q(s1)
q2 = q;     q2.dom.in = [0,2];      q2.dom.out = [0,2];
expect_error(@() common_states(q,q2),'same domain');
% common_vars: also the output space
[A,B] = common_vars(u,w);
chk(isequal(A.vars.out,vn) && isequal(B.vars.out,vn) && isequal(A.dom.out,B.dom.out) && isequal(size(A.depmat.out),size(B.depmat.out)),'(3) common_vars output variables');
chk(isequal(A.varname,B.varname) && isequal(A.vars.in,B.vars.in),'(3) common_vars states');
chk(isequal(A.C.vars.out,A.vars.out) && isequal(B.C.vars.out,B.vars.out),'(3) output space of C consistent');

%% (4) plus, minus, uminus, mtimes
for rep=1:5
    rng(10+rep);
    v = rv();
    a = randn;      b = rp(vn,1);
    chk_close(ev(x+y,v),v.x+v.y,tol,sprintf('(4.%d) x+y',rep));
    chk_close(ev(x-y,v),v.x-v.y,tol,sprintf('(4.%d) x-y',rep));
    chk_close(ev(-x,v),-v.x,tol,sprintf('(4.%d) -x',rep));
    chk_close(ev(a*x,v),a*v.x,tol,sprintf('(4.%d) a*x',rep));
    chk_close(ev(x*a,v),a*v.x,tol,sprintf('(4.%d) x*a',rep));
    chk_close(ev(a.*x+y.*z,v),a*v.x+v.y*v.z,tol,sprintf('(4.%d) a.*x+y.*z',rep));
    chk_close(ev(x.*y.*z-y,v),v.x*v.y*v.z-v.y,tol,sprintf('(4.%d) x.*y.*z-y',rep));
    chk_close(ev(x+0,v),v.x,tol,sprintf('(4.%d) x+0',rep));
    chk_close(ev(0+x,v),v.x,tol,sprintf('(4.%d) 0+x',rep));
    chk_close(ev((x+y)*z,v),(v.x+v.y)*v.z,tol,sprintf('(4.%d) (x+y)*z',rep));
    chk_close(ev(x*y,v),v.x*v.y,tol,sprintf('(4.%d) x*y (scalar mtimes)',rep));
    chk_close(ev(x.*x-x.*x+x,v),v.x,tol,sprintf('(4.%d) cancellation',rep));
    chk_close(ev(x.*(y+z),v),v.x*(v.y+v.z),tol,sprintf('(4.%d) x.*(y+z)',rep));
end
% number +- state: see NOTES
rng(50);    v = rv();       a = randn;
chk_or_warn(@() ev(x+a,v),v.x+a,tol,'x+a');
chk_or_warn(@() ev(a+x,v),v.x+a,tol,'a+x');
chk_or_warn(@() ev(x-a,v),v.x-a,tol,'x-a');
chk_or_warn(@() ev(a-x,v),a-v.x,tol,'a-x');
chk_close(ev(x+y+z-x,struct('x',2*s1,'y',s2,'z',1+s1)),s2+1+s1,tol,'(4) sum collects terms');
% the sum of a polynomial with itself is a single monomial
S = x+x;
chk(isequal(size(S.Z),[1,1]) && isequal(full(S.Z),1),'(4) x+x has one monomial');
% constants
cz = x-x;
v = rv();
chk_close(ev(cz+y,v),v.y,tol,'(4) zero polynomial plus y');
% errors
expect_error(@() x+'a','not supported');
expect_error(@() 'a'+x,'not supported');
expect_error(@() [x;y]+x,'Dimensions');
expect_error(@() [x;y]*[x;y],'Incorrect dimensions');
expect_error(@() [x;y;z]*[x;y],'Incorrect dimensions');
expect_error(@() x*[x;y],'Row dimensions of operators must match');      % scalar * vector of polynomials is not expanded

%% (5) vertcat and horzcat
for rep=1:3
    rng(100+rep);
    v = rv();       a = randn;
    V = [x;y];
    chk(isequal(size(V),[2,1]),sprintf('(5.%d) size of [x;y]',rep));
    chk_close(ev(V,v),[v.x;v.y],tol,sprintf('(5.%d) [x;y]',rep));
    chk_close(ev([x;a],v),[v.x;a],tol,sprintf('(5.%d) [x;a]',rep));
    chk_close(ev([a;x],v),[a;v.x],tol,sprintf('(5.%d) [a;x]',rep));
    chk_close(ev([x.*y;z;x+y],v),[v.x*v.y;v.z;v.x+v.y],tol,sprintf('(5.%d) [x*y;z;x+y]',rep));
    chk_close(ev([x;[];y],v),[v.x;v.y],tol,sprintf('(5.%d) [x;[];y]',rep));
    chk_close(ev([[];x],v),v.x,tol,sprintf('(5.%d) [[];x]',rep));
    chk_close(ev([V;z],v),[v.x;v.y;v.z],tol,sprintf('(5.%d) [[x;y];z]',rep));
    chk_close(ev([x;x],v),[v.x;v.x],tol,sprintf('(5.%d) [x;x]',rep));
end
chk(isequal(size([x;y;z]),[3,1]) && isequal(size(vertcat()),[0,0]) && isequal(size(horzcat()),[0,0]),'(5) sizes, empty calls');
chk(isequal(size(horzcat(x)),[1,1]),'(5) horzcat of one');
v = rv();
chk_close(ev(horzcat([],x),v),v.x,tol,'(5) horzcat([],x)');
% states with different spatial variables can be stacked
UW = [u;w];
chk(isequal(size(UW),[2,1]) && isequal(UW.vars.in,vn),'(5) [u;w]');
% errors
expect_error(@() [x,y],'not currently supported');
expect_error(@() [x,[y;z]],'Output function spaces|not currently supported');
expect_error(@() vertcat(x,'a'),'only supported with objects');
expect_error(@() horzcat(x,'a'),'only supported with objects');
Pbig = x;       Pbig.dom.in = [0,2;0,2];    Pbig.dom.out = [0,2;0,2];
expect_error(@() [x;Pbig],'domain|same');

%% (6) apply_poly against eval_polyopvar
for rep=1:3
    rng(200+rep);
    v = rv();
    Ps = {x, x+y, x.*y, x.*y+2*z, [x;y.*z], -x.*z+3*y, (x+y).*(y+z), [x;y;z]*1};
    for k=1:numel(Ps)
        P = Ps{k};
        chk_close(apply_poly(P,{v.x;v.y;v.z}),ev(P,v),tol,sprintf('(6.%d.%d) apply_poly',rep,k));
    end
end
% states not appearing in P need no value
P = x.*z;                                  % states x and z, values in the order of P.varname
chk_close(apply_poly(P,{2*s1;s2}),2*s1*s2,tol,'(6) x.*z');
P = clean_poly(x+y-y);                      % y cancels and its monomial is removed: no value needed for it
chk_close(apply_poly(P,{2*s1;[]}),2*s1,tol,'(6) empty value for a state not in P');
chk_close(apply_poly(P,{2*s1;s2^2}),2*s1,tol,'(6) value for a state not in P is ignored');
chk_close(apply_poly(x,{s1+1}),s1+1,tol,'(6) fewer values than states');
chk_close(apply_poly(x,{1.5;2;3}),polynomial(1.5),tol,'(6) double values');
% vector-valued states
[xx,xnm,xdm] = rand_polyopvar(2,[2;3],{'s1','s2'});
v1 = [rp(vn,1);rp(vn,1)];       v2 = [rp(vn,1);rp(vn,1);rp(vn,1)];
chk_close(apply_poly(xx,{v1;v2}),[v1;v2],tol,'(6) vector-valued states');
% a Kronecker output (otype(1)=2) renames the output variables s to s_1 in apply_op; apply_poly renames them back
Tx = x.C.ops{1,1};      Tx.otype = [2,2];       Tx = initialize(Tx);
xk = x;     xk.C = tensopmat({Tx});     xk.chkval = false;      xk = initialize(xk);
v = rv();
chk_close(apply_poly(xk,{v.x}),v.x,tol,'(6) Kronecker output, variables renamed back');
chk_close(apply_poly(2*xk+x,{v.x}),3*v.x,tol,'(6) Kronecker and pointwise output');
% degree-0 function with a 'tensopvar' coefficient with Kronecker output: apply_op names the output variable s1_1,
% apply_poly renames it back (subs is not involved without states)
vars_s1 = struct('in',{cell(1,0)},'out',{{'s1'}});      dom_s1 = struct('in',zeros(0,2),'out',[0,1]);
Tm = mat2tensopvar(1+s1,vars_s1,dom_s1);                % multiplier from R, pointwise output [1,2]
Tm.otype = [2,2];       Tm = initialize(Tm);
f0t = polyopvar();      f0t.varname = cell(0,1);        f0t.Z = sparse(1,0);
f0t.C = tensopmat({Tm});        f0t = initialize(f0t);
chk_close(apply_poly(f0t,{}),1+s1,tol,'(6) degree-0 function, variable renamed back');
chk_close(apply_poly(f0t+f0t,{}),2+2*s1,tol,'(6) degree-0 sum, variable renamed back');
% errors
expect_error(@() apply_poly(x.*y,{s1}),'Values must be specified for all states');

%% (7) apply_poly of functionals (intop) with monomials of degree d>1: the functional acts on the
% tensor product of the values (apply_op of 'intop' takes the cell of the factors, the copy s_j of
% a variable being that of the j-th factor). Compared with the integrals computed directly, not
% with eval_polyopvar, which cannot evaluate these monomials
I2 = @(p) int(int(p,s1,0,1),s2,0,1);
rng(300);
v = rv();       vx = v.x;   vy = v.y;   vz = v.z;
vals3 = {vx;vy;vz};
c3 = constant2polyopvar(3);
F1 = innerprod(x,y);                    I1 = I2(vx*vy);
F2 = innerprod(x,x)+3*innerprod(y,z);   I2v = I2(vx^2)+3*I2(vy*vz);
Fc = {F1,                    I1
      F2,                    I2v
      3*F1,                  3*I1
      F1.*c3,                3*I1
      c3.*F1,                3*I1
      (F1+2).*c3,            3*(I1+2)
      F1.*(F1+c3),           I1*(I1+3)
      F1.*F2,                I1*I2v
      F1+c3,                 I1+3};
for k=1:size(Fc,1)
    chk_close(apply_poly(Fc{k,1},vals3),Fc{k,2},tol,sprintf('(7.%d) functional',k));
end
% degree 3: <x.*y,z> = int vx*vy*vz
F3 = innerprod(x.*y,z);
chk(isequal(full(num_factors(F3)),3),'(7) degree 3 functional');
chk_close(apply_poly(F3,vals3),I2(vx*vy*vz),tol,'(7) degree 3');
% the operator applied to the cell of the factors, directly
op = F1.C.ops{1,1};
chk_close(apply_op(op,{vx,vy}),I1,tol,'(7) apply_op of intop with a cell of factors');
chk_close(apply_op(op,{vy,vx}),I1,tol,'(7) apply_op of intop, factors swapped (symmetric functional)');
% vector-valued states a, b of dimension 2: <a,M,b> = int a'*M*b with M not symmetric, so that the
% order of the Kronecker product of the factors (first factor slowest) is tested
M = [1,2;3,5];
va = [rp(vn,1);rp(vn,1)];       vb = [rp(vn,1);rp(vn,1)];
av = polyopvar({'a'},2,vn,[0,1;0,1],true(1,2));
bv = polyopvar({'b'},2,vn,[0,1;0,1],true(1,2));
chk_close(apply_poly(innerprod(av,M,bv),{va;vb}),I2(va.'*M*vb),tol,'(7) vector-valued states, asymmetric weight');
chk_close(apply_poly(innerprod(bv,M,av),{va;vb}),I2(vb.'*M*va),tol,'(7) vector-valued states, other order');
% spatially varying weight
wgt = 1+s1+s2^2;
chk_close(apply_poly(innerprod(x,wgt,y),vals3),I2(wgt*vx*vy),tol,'(7) weighted');
% a state in a single variable
pu = rp({'s1'},3);
Fu = innerprod(u,u);
chk_close(apply_poly(Fu,{pu}),int(pu^2,s1,0,1),tol,'(7) functional of a state in one variable');

fprintf('Test_polyopvar_class passed (%d cases).\n',counter(0))


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function n = counter(inc)
% number of checks performed
persistent c
if isempty(c),  c = 0;  end
c = c+inc;
n = c;
end

function chk(cond,label)
counter(1);
if ~cond
    error("Test_polyopvar_class failed %s.",label)
end
end

function expect_error(fun,msg)
counter(1);
try
    out = fun();            %#ok<NASGU>
catch ME
    if nargin>=2 && isempty(regexp(ME.message,msg,'once'))
        error("Test_polyopvar_class failed: error '%s' does not match '%s'.",ME.message,msg)
    end
    return
end
error("Test_polyopvar_class failed: no error for %s.",func2str(fun))
end

function chk_or_warn(fun,Yref,tol,label)
% check the value of fun() if it can be computed, else warn
try
    Y = fun();
catch ME
    warning('Test_polyopvar_class:unsupported',"%s errors: %s",label,ME.message);
    return
end
chk_close(Y,Yref,tol,label);
end

function chk_close(Y,Yref,tol,label)
counter(1);
Yref = polynomial(Yref);        Y = polynomial(Y);
if ~isequal(size(Y),size(Yref))
    error("Test_polyopvar_class failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
D = Y-Yref;
err = full(max([0;abs(nonzeros(D.coefficient))]));
scale = max(1,full(max([0;abs(nonzeros(Yref.coefficient))])));
if err>tol*scale
    error("Test_polyopvar_class failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient)) && ~contains(label,'zero')
    error("Test_polyopvar_class %s: reference value is zero; test would be vacuous.",label)
end
end

function f = mkdeg0(P,vars,dom)
% degree-0 polyopvar with the matrix P (or cell of blocks) of type polynomial/double/dpvar acting on R as coefficient
f = polyopvar();        f.varname = cell(0,1);      f.Z = sparse(1,0);
if ~iscell(P),     P = {P};    end
f.C = tensopmat(P,vars,dom);
f = initialize(f);
end

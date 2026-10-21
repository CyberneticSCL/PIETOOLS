% TEST_SOPVAR_DIFF
% Tests @sopvar/diff against the operator semantics: for a random 'sopvar'
% P and polynomial u, apply_sopvar(dP,u) must equal the derivative of
% apply_sopvar(P,u) computed with 'polynomial' calculus. apply_sopvar
% evaluates the kernels from the class definition and shares no code with
% diff for sopvar, so this is not a round-trip test. Also checks: errors
% for derivatives that are not PI operators, cancellation of boundary
% terms, stacked derivatives, non-output variables, and the 1D opvar 'diff'
% with option 'pure' as an independent implementation of the same Leibniz
% rule.
%
% Run pietools_path_update first.

% DJ (CC), 10/20/2026: Initial coding;

rng(7);
tol = 1e-9;
npass = 0;      nfail = 0;      n_nontriv = 0;

%% Random operators, all combinations of shared / input-only / output-only variables
cfgs = {struct('in',{{'s1'}},         'out',{{'s1'}}), ...
        struct('in',{{'s1'}},         'out',{{'s1','s2'}}), ...
        struct('in',{{'s1','s2'}},    'out',{{'s1'}}), ...
        struct('in',{{'s1','s2'}},    'out',{{'s1','s2'}}), ...
        struct('in',{{'s1','s2'}},    'out',{{'s1','s3'}}), ...
        struct('in',{{'s1','s2','s3'}},'out',{{'s1','s2','s3'}})};
dom = [0.5,2];
for c=1:numel(cfgs)
    vars = cfgs{c};
    vS3 = intersect(vars.in,vars.out);
    for trial=1:4
        m = randi(3);   n = randi(3);
        P0 = rand_sopvar([m,n],vars,dom,2,0.4);
        for v = vars.out
            vname = v{1};
            % Remove multiplier terms in direction vname: not a PI derivative otherwise
            P = zero_mult(P0,vname);
            x = rand_poly([n,2],polynomial(P.vars.in),2,5);
            [dP,is_zero] = diff(P,{vname},1);
            ok = check_deriv(P,dP,x,vname,1,tol);
            % count cases where the derivative is nonzero, so a trivial pass shows
            n_nontriv = n_nontriv + (norm_poly(diff(apply_sopvar(P,x),polynomial({vname})))>tol);
            if ~ok
                nfail = nfail+1;    fprintf('FAIL: cfg %d, trial %d, d/d%s\n',c,trial,vname);
            else
                npass = npass+1;
            end
            % Error when a multiplier term is present
            if ~isempty(intersect(vname,vS3))
                Pm = add_mult(P0,vname);
                if ~throws(@() diff(Pm,{vname},1),"cannot be expressed as a PI operator")
                    nfail = nfail+1;    fprintf('FAIL: no error, cfg %d, d/d%s\n',c,vname);
                else
                    npass = npass+1;
                end
            end
        end
    end
end

%% Kernel vanishing at t=s: second derivative is a multiplier
% P u = int_0^s (s-t) u(t) dt, d^2/ds^2 P u = u
vars = struct('in',{{'s1'}},'out',{{'s1'}});
dm = struct('in',[0,1],'out',[0,1]);
C2 = sparse([2,1],[1,2],[1,-1],2,2);           % ZL=ZR=[0;1]: s*1 - 1*t
P = sopvar({sparse(2,2);C2;sparse(2,2)},vars,{[0;1]},{[0;1]},dm,[1,1]);
[d2P,is_zero] = diff(P,{'s1'},2);
x = rand_poly([1,1],polynomial({'s1'}),3,4);
err = norm_poly(apply_sopvar(d2P,x)-x);
[npass,nfail] = tally(err<tol && ~is_zero,npass,nfail,'(s-t) kernel, d^2/ds^2');
if ~throws(@() diff(P,{'s1'},3),"cannot be expressed as a PI operator")
    nfail = nfail+1;    disp('FAIL: third derivative of (s-t) kernel should error');
else
    npass = npass+1;
end

%% Boundary terms that cancel up to roundoff: the second derivative is still a PI operator
% P u = int_0^s c(s-1)t u dt + int_s^1 c3*s(t-1) u dt, with c=0.1+0.2 and c3=0.3, so that the kernels
% are continuous at t=s up to roundoff (as for operators from a PDE conversion): no multiplier remains
% in d/ds P, so that d^2/ds^2 P is defined
cc = 0.1+0.2;       c3 = 0.3;
Cr2 = sparse([2,1],[2,2],[cc,-cc],2,2);             % ZL=ZR=[0;1]: c*(s*t - t)
Cr3 = sparse([2,2],[2,1],[c3,-c3],2,2);             % c3*(s*t - s)
Pr = sopvar({sparse(2,2);Cr2;Cr3},vars,{[0;1]},{[0;1]},dm,[1,1]);
d1Pr = diff(Pr,{'s1'},1);       d2Pr = diff(Pr,{'s1'},2);
xr = rand_poly([1,1],polynomial({'s1'}),3,4);
Prx = apply_sopvar(Pr,xr);
s1 = polynomial({'s1'});
[npass,nfail] = tally(nnz(d1Pr.params{1})==0 && norm_poly(apply_sopvar(d2Pr,xr)-diff(diff(Prx,s1),s1))<tol && norm_poly(apply_sopvar(d2Pr,xr))>0,npass,nfail,'boundary terms cancel up to roundoff');
% a boundary term that does not cancel is not removed
Pn = sopvar({sparse(2,2);Cr2;0.5*Cr3},vars,{[0;1]},{[0;1]},dm,[1,1]);
[npass,nfail] = tally(throws(@() diff(Pn,{'s1'},2),"cannot be expressed as a PI operator"),npass,nfail,'boundary terms that do not cancel');

%% Full integral (lower + upper with equal kernels): boundary terms cancel
% P u = int_0^1 (1 + s*t^2) u(t) dt, any number of derivatives is a PI operator
C = sparse([1,2],[1,3],[1,1],2,3);                % ZL=[0;1], ZR=[0;1;2]: 1 + s*t^2
P = sopvar({sparse(2,3);C;C},vars,{[0;1]},{[0;1;2]},dm,[1,1]);
for D=1:3
    dP = diff(P,{'s1'},D);
    x = rand_poly([1,1],polynomial({'s1'}),3,4);
    Px = apply_sopvar(P,x);
    for j=1:D, Px = diff(Px,polynomial({'s1'})); end
    err = norm_poly(apply_sopvar(dP,x)-Px);
    [npass,nfail] = tally(err<tol,npass,nfail,['full integral, order ',num2str(D)]);
end

%% Upper integral only: boundary term is -K_3(s,s) (lower integral absent)
% P u = int_s^1 (1+s*t) u dt, d/ds P u = -(1+s^2) u(s) + int_s^1 t u dt
Cu = sparse([1,2],[1,2],[1,1],2,2);                 % ZL=ZR=[0;1]: 1 + s*t
Pu = sopvar({sparse(2,2);sparse(2,2);Cu},vars,{[0;1]},{[0;1]},dm,[1,1]);
dPu = diff(Pu,{'s1'},1);
xu = rand_poly([1,1],polynomial({'s1'}),3,4);
[npass,nfail] = tally(norm_poly(apply_sopvar(dPu,xu)-diff(apply_sopvar(Pu,xu),polynomial({'s1'})))<tol && any(cellfun(@nnz,dPu.params)>0),npass,nfail,'upper integral only');

%% Boundary term when ZR has no degree 0: ZR must be enlarged
% 1x2 operator, ZL=[0;1], ZR=[1;2]; parameters are 2 x 4 (columns: matrix column outer, monomial inner),
% so a stale column stride for the enlarged ZR would misplace the second matrix column
C2z = sparse([2,1],[1,4],[1,2],2,4);        C3z = sparse([1,2],[2,3],[3,1],2,4);
Pz = sopvar({sparse(2,4);C2z;C3z},struct('in',{{'s1'}},'out',{{'s1'}}),{[0;1]},{[1;2]},struct('in',[0,1],'out',[0,1]),[1,2]);
dPz = diff(Pz,{'s1'},1);
xz = rand_poly([2,1],polynomial({'s1'}),3,4);
[npass,nfail] = tally(norm_poly(apply_sopvar(dPz,xz)-diff(apply_sopvar(Pz,xz),polynomial({'s1'})))<tol && any(dPz.ZR{1}==0),npass,nfail,'ZR enlarged by degree 0');

%% Stacked derivatives, mixed partials, repeated and scalar orders
vars = struct('in',{{'s1','s2'}},'out',{{'s1','s2'}});
dm = struct('in',[0,1;0,1],'out',[0,1;0,1]);
P0 = rand_sopvar([2,3],vars,dm,2,0.5);
P = zero_mult(zero_mult(P0,'s1'),'s2');
x = rand_poly([3,1],polynomial(P.vars.in),2,6);
% each row is one derivative; rows stack in the order of d
dmat = [1,0; 0,1; 1,1];
dP = diff(P,{'s1','s2'},dmat);
Px = apply_sopvar(P,x);     s1 = polynomial({'s1'});    s2 = polynomial({'s2'});
ref = [diff(Px,s1); diff(Px,s2); diff(diff(Px,s1),s2)];
[npass,nfail] = tally(norm_poly(apply_sopvar(dP,x)-ref)<tol && isequal(size(dP),[6,3]),npass,nfail,'stacked rows');
% repeated variable names add their orders (here 2 in s2, which is not a PI operator for this P); pvar input
[npass,nfail] = tally(throws(@() diff(P,[s2,s1,s2],[1,1,1]),"cannot be expressed as a PI operator"),npass,nfail,'repeated variable, second order error');
% orders in repeated names are added, [s1,s1] with [1,1] is d^2/ds1^2
Pfull = sopvar({sparse(2,3);C;C},struct('in',{{'s1'}},'out',{{'s1'}}),{[0;1]},{[0;1;2]},struct('in',[0,1],'out',[0,1]),[1,1]);
dP2 = diff(Pfull,[polynomial({'s1'}),polynomial({'s1'})],[1,1]);
dP2b = diff(Pfull,{'s1'},2);
xf = rand_poly([1,1],polynomial({'s1'}),3,4);
[npass,nfail] = tally(norm_poly(apply_sopvar(dP2,xf)-apply_sopvar(dP2b,xf))<tol && ...
    norm_poly(apply_sopvar(dP2,xf)-diff(apply_sopvar(Pfull,xf),polynomial({'s1'}),2))<tol,npass,nfail,'repeated variables, pvar input');
% scalar order
P1 = zero_mult(rand_sopvar([1,1],struct('in',{{'s1'}},'out',{{'s1'}}),[0,1],2,0.5),'s1');
x1 = rand_poly([1,1],polynomial({'s1'}),3,4);
dP = diff(P1,{'s1'},1);
[npass,nfail] = tally(norm_poly(apply_sopvar(dP,x1)-diff(apply_sopvar(P1,x1),polynomial({'s1'})))<tol,npass,nfail,'scalar order');
% default order 1, char variable name
dPa = diff(P1,'s1');
[npass,nfail] = tally(isequal(dPa.params,dP.params),npass,nfail,'default order, char name');

%% Non-output variables, zero order, invalid input
[dP,is_zero] = diff(P,{'s9'},1);
[npass,nfail] = tally(is_zero && all(cellfun(@nnz,dP.params(:))==0) && isequal(size(dP),[2,3]),npass,nfail,'non-output variable');
% row 1 differentiates a non-output variable, row 2 does not
dP = diff(P,{'s1','s9'},[1,1; 1,0]);
ref = [0*Px; diff(Px,s1)];
[npass,nfail] = tally(norm_poly(apply_sopvar(dP,x)-ref)<tol,npass,nfail,'mixed zero and nonzero rows');
[dP,is_zero] = diff(P,{'s1'},0);
[npass,nfail] = tally(~is_zero && isequal(dP.params,P.params),npass,nfail,'zero order');
[npass,nfail] = tally(throws(@() diff(P,{'s1'},-1),"nonnegative"),npass,nfail,'negative order');
[npass,nfail] = tally(throws(@() diff(1,P,1),"sopvar"),npass,nfail,'not a sopvar');
[npass,nfail] = tally(throws(@() diff(P,{'s1'},[1,1]),"each of the spatial"),npass,nfail,'order size');

%% Independent implementation: 1D opvar diff(P,s,'pure')
pvar s1 s1_dum
for trial=1:5
    nn = randi(3);      mm = randi(3);
    opvar opP;
    opP.I = [0.5,2];    opP.var1 = s1;      opP.var2 = s1_dum;
    opP.dim = [0,0;mm,nn];
    opP.R.R1 = rand_poly([mm,nn],[s1;s1_dum],2,5);
    opP.R.R2 = rand_poly([mm,nn],[s1;s1_dum],2,5);
    dopP = diff(opP,s1,1,'pure');
    sP = opvar2sopvar(opP);
    dsP = diff(sP,{'s1'},1);
    x = rand_poly([nn,1],s1,3,6);
    err = norm_poly(apply_opvar(dopP,x)-apply_sopvar(dsP,x));
    [npass,nfail] = tally(err<tol,npass,nfail,'vs opvar diff pure');
end

fprintf('\nTest_sopvar_diff: %d passed, %d failed (%d random cases with nonzero derivative).\n',npass,nfail,n_nontriv);
if nfail>0
    error('Test_sopvar_diff failed.');
end


%% Helpers
function ok = check_deriv(P,dP,x,vname,D,tol)
% d^D/dv apply(P,x) against apply(dP,x)
v = polynomial({vname});
Px = apply_sopvar(P,x);
for j=1:D, Px = diff(Px,v); end
ok = norm_poly(apply_sopvar(dP,x)-Px)<tol;
end

function e = norm_poly(p)
p = cleanpoly(polynomial(p),1e-12);
if isempty(p.C) || nnz(p.C)==0
    e = 0;
else
    e = max(abs(p.C(:)));
end
end

function [npass,nfail] = tally(ok,npass,nfail,name)
if ok
    npass = npass+1;
else
    nfail = nfail+1;    fprintf('FAIL: %s\n',name);
end
end

function tf = throws(f,msg)
tf = false;
try
    f();
catch ME
    tf = contains(ME.message,msg);
end
end

function P = zero_mult(P,vname)
% Zero all parameters with a multiplier (gamma==1) in direction vname
vS3 = intersect(P.vars.in,P.vars.out);
t = find(strcmp(vS3,vname));
if isempty(t), return, end
params = P.params;
for k=1:numel(params)
    g = gamma_of_cell(k,numel(vS3));
    if g(t)==1
        params{k} = sparse(size(params{k},1),size(params{k},2));
    end
end
P = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

function P = add_mult(P,vname)
% Make sure the multiplier parameters in direction vname are nonzero
vS3 = intersect(P.vars.in,P.vars.out);
t = find(strcmp(vS3,vname));
params = P.params;
for k=1:numel(params)
    g = gamma_of_cell(k,numel(vS3));
    if g(t)==1 && nnz(params{k})==0
        params{k}(1,1) = 1;                 % constant monomials: canonical
    end
end
P = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

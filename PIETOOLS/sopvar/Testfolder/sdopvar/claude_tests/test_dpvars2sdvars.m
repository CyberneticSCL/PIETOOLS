%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_DPVARS2SDVARS checks 'dpvars2sdvars' (sopvar/Testfolder/converters)
% against the definition in its header,
%
%   PARAMS{k} = (Im o ZL{1}(s1) o ... o ZL{M}(sM))^T unvec(A{k} + B{k}'*d)
%                   (In o ZR{1}(t1) o ... o ZR{N}(tN)),
%
% where d is the GLOBAL, unified decision-variable vector (dvars) shared
% across every parameter k, not each parameter's own local decision
% variables.
%
% Per case, against the semantics and never a round trip (CLAUDE.md S4):
% params are built directly ('dpvar' constructor, 'polynomial' constructor,
% or a plain double), never through 'sdvar2dpvar'/'dpvar2sdvar'. Ground
% truth per parameter is its own native evaluation (dpvar2poly+subs for
% 'dpvar', subs for 'polynomial', the value itself for 'double'); the value
% under test reconstructs the same point from (A{k},B{k},Z,dvars) via the
% header formula above, with ONE shared 'd' vector of length numel(dvars)
% indexed by name into each parameter's own (sub)set of decision variables.
%
% Defects this catches (fixed 10/05/2026, dpvars2sdvars.m header):
%  (1) 'isB' was never assigned in the 'dpvar' branch, so ANY dpvar
%      parameter with a decision variable crashed;
%  (2) 'dvars' was overwritten by whichever dpvar parameter ran last in the
%      loop, and B{k} was indexed against that parameter's own (local)
%      dvarname list -- so with >1 dpvar parameters on different decision
%      variables, different B{k} referred to incompatible 'd' vectors;
%  (3) the nargin>=3, chk_vars~=false branch (user-supplied Z) had four
%      copy/paste defects (an invalid dot-call, a negated isa test that
%      sent every dpvar parameter to the "must be polynomial/dpvar/
%      quadPoly" error, and two undefined-variable references);
%  (4) with chk_vars==false, a 'double' parameter crashed the same way,
%      since the (skipped) preprocessing was the only place 'double' got
%      wrapped as a constant 'polynomial';
%  (5) the docstring claimed 'dvars' is a ROW ("1 x K ... d = dvars'"); the
%      fixed union-building for (2) produces a COLUMN, matching 'sdopvar'
%      (class header: "d is column vector") and 'dpvar2sdvar' (P.dvarname,
%      fixed the same way).
% None of the five was caught before because nothing in the repository
% calls this function yet.
%
% Sweep: 1-3 parameters per trial, each 'double', 'polynomial' or 'dpvar'
% with independently random dims, native (possibly sparse/non-grid)
% monomial sets, and (for 'dpvar') a random subset of a shared 6-name pool
% of decision variables, so parameters overlap, partially overlap, or share
% no decision variables at all. Plus: the nargin>=3 branch with chk_vars
% true and false, mixing 'double'/'polynomial'/'dpvar' in one call; and that
% 'dvars' comes back a column with no duplicate names.
%
% Initial coding DJ, 10/05/2026
%

tol = 1e-9;         nsamp = 3;
nchk = 0;           fails = {};
dpool = arrayfun(@(i)sprintf('dd%d',i),1:6,'un',0);

%% NARGIN3: user-supplied Z, mixed double/polynomial/dpvar, both chk_vars
vars = struct('out',{{'s1'}},'in',{{'t1'}});
Zspec = struct('out',{{[0;1;2]}},'in',{{[0;1]}});
degmat_d = [1 0; 0 1];
Cd = sparse([1 2; 0.5 0]);
Dp = dpvar(Cd,degmat_d,{'s1';'t1'},{'d1'},[1,1]);
Pp = polynomial(sparse([1;2]),[1 0; 0 1],{'s1';'t1'},[1,1]);
params0 = {3.5, Pp, Dp};
for chk = [true,false]
    msg = '';
    try
        [A,B,Z,dvars] = dpvars2sdvars(params0,vars,Zspec,chk);
        if ~iscolumn(dvars)
            msg = 'dvars is not a column';
        else
            ZL = Z.out;     ZR = Z.in;
            s = 1.3;    t = -0.6;  d = randn(numel(dvars),1);
            zL = 1;     for i=1:numel(ZL), zL = kron(zL,s.^ZL{i}(:)); end
            zR = 1;     for i=1:numel(ZR), zR = kron(zR,t.^ZR{i}(:)); end
            for k = 1:3
                pk = params0{k};
                if isa(pk,'double')
                    mm=size(pk,1); nn=size(pk,2);  ref = pk;
                elseif isa(pk,'polynomial')
                    [mm,nn]=size(pk); ref = double(subs(pk,{'s1';'t1'},[s;t]));
                else
                    [mm,nn]=size(pk);
                    names=[pk.dvarname(:);{'s1';'t1'}];
                    [~,loc]=ismember(pk.dvarname(:),dvars(:));
                    ref = double(subs(dpvar2poly(pk),names,[d(loc);s;t]));
                end
                Cm = reshape(full(A{k}+B{k}'*d),mm*numel(zL),nn*numel(zR));
                got = kron(eye(mm),zL)'*Cm*kron(eye(nn),zR);
                err = norm(got-ref,'fro')/max(1,norm(ref,'fro'));
                if err>tol
                    msg = sprintf('chk_vars=%d k=%d: relative error %.2e',chk,k,err);
                    break
                end
            end
        end
    catch ME
        msg = ['error: ',ME.message];
    end
    nchk = nchk+1;
    if ~isempty(msg),   fails{end+1} = ['NARGIN3: ',msg];  end
end

%% VALUE: sweep over parameter count, type mix, variables and decision vars
ncase = 0;
for trial = 1:150
    ncase = ncase+1;
    rng(7000+trial);
    msg = check_trial(dpool,nsamp,tol);
    nchk = nchk+nsamp;
    if ~isempty(msg)
        fails{end+1} = sprintf('trial %d: %s',trial,msg);
    end
end

ntot = 2+ncase;
for k = 1:numel(fails),     fprintf('FAIL  %s\n',fails{k});     end
assert(isempty(fails),'test_dpvars2sdvars: %d of %d cases FAILED.',numel(fails),ntot);
fprintf('\ntest_dpvars2sdvars passed (%d cases, %d checks).\n',ntot,nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_trial(dpool,nsamp,tol)
% Empty msg on success. Errors are reported, not thrown, so one crashing
% case does not hide the others.
msg = '';
M = randi([0,2]);  N = randi([0,2]);
varnames1 = arrayfun(@(i)sprintf('s%d',i),1:M,'un',0);
varnames2 = arrayfun(@(i)sprintf('t%d',i),1:N,'un',0);
varnames_full = [varnames1,varnames2];
vars = struct('out',{varnames1},'in',{varnames2});
nparams = randi([1,3]);
params = cell(1,nparams);
truth = cell(1,nparams);
for k = 1:nparams
    m = randi(2);   n = randi(2);   typ = randi(3);
    nZ = randi([1,4]);
    degmat = randi([0,3],nZ,M+N);
    degmat = unique(degmat,'rows');
    nZ = size(degmat,1);
    if nZ==0,   degmat = zeros(0,M+N);     end
    if typ==1
        val = randn(m,n);
        params{k} = val;
        truth{k} = struct('type','double','m',m,'n',n,'val',val);
    elseif typ==2
        Cp = sprandn(nZ,m*n,0.6);
        p = polynomial(Cp,degmat,varnames_full(:),[m,n]);
        params{k} = p;
        truth{k} = struct('type','polynomial','m',m,'n',n,'obj',p);
    else
        ndk = randi([0,4]);
        dvn = dpool(randperm(6,ndk));
        Cd = sprandn((ndk+1)*m,nZ*n,0.6);
        for j = 1:ndk
            rr = (randi(m)-1)*(ndk+1) + 1 + j;
            cc = randi(max(nZ*n,1));
            Cd(rr,cc) = randn();
        end
        D = dpvar(Cd,degmat,varnames_full(:),dvn(:),[m,n]);
        params{k} = D;
        truth{k} = struct('type','dpvar','m',m,'n',n,'obj',D);
    end
end

try
    [A,B,Z,dvars] = dpvars2sdvars(params,vars);
catch ME
    msg = ['error: ',ME.message];   return
end
ndvars = numel(dvars);
if numel(unique(dvars))~=ndvars
    msg = 'dvars has duplicates';   return
end
if ~iscolumn(dvars)
    msg = 'dvars is not a column';   return
end
for k = 1:nparams
    if strcmp(truth{k}.type,'dpvar') && ~all(ismember(truth{k}.obj.dvarname,dvars))
        msg = sprintf('param %d: dvars union missing a local dvarname',k);   return
    end
end

ZL = Z.out;     ZR = Z.in;
for s_i = 1:nsamp
    d = randn(ndvars,1);
    s = rand_pt(M);     t = rand_pt(N);
    for k = 1:nparams
        m = truth{k}.m;     n = truth{k}.n;
        got = kernel_eval(A{k},B{k},ZL,ZR,[m,n],d,s,t);
        switch truth{k}.type
            case 'double'
                ref = truth{k}.val;
            case 'polynomial'
                ref = eval_poly(truth{k}.obj,varnames_full(:),[s;t]);
            case 'dpvar'
                Dobj = truth{k}.obj;
                names = [Dobj.dvarname(:);varnames_full(:)];
                [~,loc] = ismember(Dobj.dvarname(:),dvars(:));
                ref = eval_poly(dpvar2poly(Dobj),names,[d(loc);s;t]);
        end
        err = norm(got-ref,'fro')/max(1,norm(ref,'fro'));
        if err>tol
            msg = sprintf('param %d (%s) sample %d: relative error %.2e',...
                           k,truth{k}.type,s_i,err);
            return
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = kernel_eval(Ak,Bk,ZL,ZR,dims,d,s,t)
zL = 1;     for i=1:numel(ZL),  zL = kron(zL,s(i).^ZL{i}(:));  end
zR = 1;     for i=1:numel(ZR),  zR = kron(zR,t(i).^ZR{i}(:));  end
m = dims(1);    n = dims(2);
Cm = reshape(full(Ak + Bk'*d), m*numel(zL), n*numel(zR));
V = kron(eye(m),zL)'*Cm*kron(eye(n),zR);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = eval_poly(Pp,names,vals)
if isempty(Pp.varname)
    V = double(Pp);     return
end
[tf,loc] = ismember(Pp.varname(:),names(:));
assert(all(tf),'polynomial has a variable with no value');
V = double(subs(Pp,Pp.varname(:),vals(loc)));
end

function x = rand_pt(k)
x = (0.3+0.5*rand(k,1)+0.9*(rand(k,1)>0.5)).*sign(randn(k,1));
end

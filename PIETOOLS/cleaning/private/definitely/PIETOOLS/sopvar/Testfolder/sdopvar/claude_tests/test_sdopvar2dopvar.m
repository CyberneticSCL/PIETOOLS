%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SDOPVAR2DOPVAR checks 'sdopvar2dopvar' (sopvar/Testfolder/converters)
% against the semantics of the input 'sdopvar' object, never a round trip
% (CLAUDE.md S4): each block (P, Q1, Q2, or R={R0,R1,R2}) is built directly
% via the 'sdopvar' constructor (never through 'dopvar2sdopvar'), and the
% resulting 'dopvar' block's own 'dpvar2poly'+'subs' expansion is compared
% to 'pi_blk_kernels' (ANCHOR-style, as in test_sdvar2dpvar.m) evaluated on
% the ORIGINAL 'sdopvar' object. Random R0/R1/R2 content need not already be
% in canonical multiplier form; the 'sdopvar' constructor rewrites it (with
% a warning), and the reference is read back out from the (possibly
% rewritten) object itself, so the check holds regardless.
%
% Found and fixed during initial coding (see sdopvar2dopvar.m header): the
% L2->L2 branch passed the SAME physical
% variable name for both the output and (dummy) input roles when calling
% 'sdvar2dpvar' ('sdvar2dpvar(...,struct(''out'',{{vname}},''in'',
% {{vname}}),...)'). Unlike 'sdopvar.vars' (where the dummy is implicit and
% both sides legitimately share one name), 'sdvar2dpvar' builds an explicit
% two-variable 'dpvar' and needs DISTINCT names for the two roles; using the
% same name twice silently produced a dpvar with a duplicated degmat column,
% which later broke 'subs' ("entries of Old must be unique"). Fixed by using
% the dummy name ('s1_dum') for the 'in' role, matching how
% 'dopvar2sdopvar.m' already builds the dpvar-facing 'vars' for its own
% R-block conversion.
%
% Sweep: 40 random trials per block type (P, Q1, Q2, R), random dims (1-3),
% random (possibly sparse/non-grid) monomial degrees for Q1/Q2/R, and random
% decision-variable subsets of a shared 4-name pool per sub-block of R (so
% R0/R1/R2 overlap, partially overlap, or share no decision variables).
%
% Initial coding DJ, 10/05/2026
%

tol = 1e-8;
fails = {};     ntot = 0;
dpool = {'d1','d2','d3','d4'};

for trial = 1:40
    rng(9000+trial);
    ntot = ntot+1;
    msg = check_case1(dpool);
    if ~isempty(msg),  fails{end+1} = sprintf('Case1 trial %d: %s',trial,msg);  end

    ntot = ntot+1;
    msg = check_case2(dpool);
    if ~isempty(msg),  fails{end+1} = sprintf('Case2 trial %d: %s',trial,msg);  end

    ntot = ntot+1;
    msg = check_case3(dpool);
    if ~isempty(msg),  fails{end+1} = sprintf('Case3 trial %d: %s',trial,msg);  end

    ntot = ntot+1;
    msg = check_case4(dpool);
    if ~isempty(msg),  fails{end+1} = sprintf('Case4 trial %d: %s',trial,msg);  end
end

for k = 1:numel(fails),    fprintf('FAIL  %s\n',fails{k});    end
assert(isempty(fails),'test_sdopvar2dopvar: %d of %d cases FAILED.',numel(fails),ntot);
fprintf('\ntest_sdopvar2dopvar passed (%d cases).\n',ntot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case1(dpool)
% P block: R^p -> R^m, no spatial variable.
msg = '';
try
    m = randi(3);   p = randi(3);
    ndk = randi([0,3]);
    dvn = dpool(randperm(4,ndk));
    nc = m*p;
    A = sprandn(nc,1,0.8);     B = sprandn(ndk,nc,0.5);
    params = struct('A',{{A}},'B',{{B}});
    vars = struct('out',{cell(1,0)},'in',{cell(1,0)});
    dom = struct('out',zeros(0,2),'in',zeros(0,2));
    Ps = sdopvar(params,vars,dvn(:),cell(1,0),cell(1,0),dom,[m,p]);
    obj = sdopvar2dopvar(Ps);
    if ~isa(obj,'dopvar'),  msg = 'output is not a dopvar';  return;  end
    d = randn(ndk,1);
    [~,loc] = ismember(Ps.Zd,dvn(:));
    Kc = pi_blk_kernels(Ps,d(loc));
    if numel(Kc)~=1,    msg = 'expected 1 kernel';   return;     end
    ref = double(Kc{1});
    got = eval_poly(dpvar2poly(obj.P),dvn(:),d);
    msg = cmp(got,ref);
catch ME
    msg = ['error: ',ME.message];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case2(dpool)
% Q1 block: L2[s1] -> R^m.
msg = '';
try
    m = randi(3);   q = randi(3);
    ndk = randi([0,3]);
    dvn = dpool(randperm(4,ndk));
    degs = unique(randi([0,3],randi([1,4]),1));
    nc = m*q*numel(degs);
    A = sprandn(nc,1,0.8);     B = sprandn(ndk,nc,0.5);
    params = struct('A',{{A}},'B',{{B}});
    vars = struct('out',{cell(1,0)},'in',{{'s1'}});
    dom = struct('out',zeros(0,2),'in',[0,1]);
    Ps = sdopvar(params,vars,dvn(:),cell(1,0),{degs},dom,[m,q]);
    obj = sdopvar2dopvar(Ps);
    if ~isa(obj,'dopvar'),  msg = 'output is not a dopvar';  return;  end
    d = randn(ndk,1);
    [~,loc] = ismember(Ps.Zd,dvn(:));
    Kc = pi_blk_kernels(Ps,d(loc));
    if numel(Kc)~=1,    msg = 'expected 1 kernel';   return;     end
    s_val = rand_pt(1);
    ref = eval_poly(Kc{1},{'s1_dum'},s_val);   % pi_blk_kernels suffixes vars.in with '_dum'
    got = eval_poly(dpvar2poly(obj.Q1),[dvn(:);{'s1'}],[d;s_val]);
    msg = cmp(got,ref);
catch ME
    msg = ['error: ',ME.message];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case3(dpool)
% Q2 block: R^p -> L2[s1].
msg = '';
try
    n = randi(3);   p = randi(3);
    ndk = randi([0,3]);
    dvn = dpool(randperm(4,ndk));
    degs = unique(randi([0,3],randi([1,4]),1));
    nc = n*p*numel(degs);
    A = sprandn(nc,1,0.8);     B = sprandn(ndk,nc,0.5);
    params = struct('A',{{A}},'B',{{B}});
    vars = struct('out',{{'s1'}},'in',{cell(1,0)});
    dom = struct('out',[0,1],'in',zeros(0,2));
    Ps = sdopvar(params,vars,dvn(:),{degs},cell(1,0),dom,[n,p]);
    obj = sdopvar2dopvar(Ps);
    if ~isa(obj,'dopvar'),  msg = 'output is not a dopvar';  return;  end
    d = randn(ndk,1);
    [~,loc] = ismember(Ps.Zd,dvn(:));
    Kc = pi_blk_kernels(Ps,d(loc));
    if numel(Kc)~=1,    msg = 'expected 1 kernel';   return;     end
    s_val = rand_pt(1);
    ref = eval_poly(Kc{1},{'s1'},s_val);
    got = eval_poly(dpvar2poly(obj.Q2),[dvn(:);{'s1'}],[d;s_val]);
    msg = cmp(got,ref);
catch ME
    msg = ['error: ',ME.message];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case4(dpool)
% R block: L2[s1] -> L2[s1], via R0,R1,R2, each with an independently random
% subset of decision variables. Content need not start in canonical
% multiplier form (see file header).
msg = '';
try
    n = randi(3);   q = randi(3);
    degs = unique(randi([0,2],randi([2,4]),1));
    if ~ismember(0,degs),   degs = [0;degs];    end
    nZ = numel(degs);
    nc = n*nZ*q*nZ;

    ndk = [randi([0,2]),randi([0,2]),randi([0,2])];
    dv0 = dpool(randperm(4,ndk(1)));
    dv1 = dpool(randperm(4,ndk(2)));
    dv2 = dpool(randperm(4,ndk(3)));
    alldv = union(union(dv0,dv1,'stable'),dv2,'stable');

    % Every params.B{k} is sized against the SAME, shared decision-variable
    % list (sdopvar.Zd); a block using only a subset gets zero rows for the
    % others, placed here by 'ismember' into the unified list.
    A0 = sprandn(nc,1,0.8);    B0 = loc_sprandn(dv0,alldv,nc);
    A1 = sprandn(nc,1,0.8);    B1 = loc_sprandn(dv1,alldv,nc);
    A2 = sprandn(nc,1,0.8);    B2 = loc_sprandn(dv2,alldv,nc);
    params = struct('A',{{A0,A1,A2}},'B',{{B0,B1,B2}});
    vars = struct('out',{{'s1'}},'in',{{'s1'}});
    dom = struct('out',[0,1],'in',[0,1]);
    Ps = sdopvar(params,vars,alldv(:),{degs},{degs},dom,[n,q]);
    obj = sdopvar2dopvar(Ps);
    if ~isa(obj,'dopvar'),  msg = 'output is not a dopvar';  return;  end

    d = randn(numel(Ps.Zd),1);
    Kc = pi_blk_kernels(Ps,d);
    if numel(Kc)~=3,    msg = sprintf('expected 3 kernels, got %d',numel(Kc));  return;  end

    s_val = rand_pt(1);     t_val = rand_pt(1);
    objs = {obj.R.R0,obj.R.R1,obj.R.R2};
    for k = 1:3
        ref = eval_poly(Kc{k},{'s1';'s1_dum'},[s_val;t_val]);
        got = eval_poly(dpvar2poly(objs{k}),[Ps.Zd(:);{'s1';'s1_dum'}],[d;s_val;t_val]);
        msg = cmp(got,ref);
        if ~isempty(msg),   msg = sprintf('k=%d: %s',k,msg);   return;    end
    end
catch ME
    msg = ['error: ',ME.message];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = cmp(got,ref)
tol = 1e-8;
err = norm(got-ref,'fro')/max(1,norm(ref,'fro'));
if err>tol,     msg = sprintf('relative error %.2e',err);
else,           msg = '';
end
end

function V = eval_poly(Pp,names,vals)
% Substitute by name; values as a COLUMN (a row is read as several points).
if isempty(Pp.varname)
    V = double(Pp);     return
end
[tf,loc] = ismember(Pp.varname(:),names(:));
assert(all(tf),'polynomial has a variable with no value');
V = double(subs(Pp,Pp.varname(:),vals(loc)));
end

function x = rand_pt(k)
% k x 1 point with |x_i| in [0.3,0.8] or [1.2,1.7], away from 0 and +-1 so
% distinct degrees give distinct values.
x = (0.3+0.5*rand(k,1)+0.9*(rand(k,1)>0.5)).*sign(randn(k,1));
end

function B = loc_sprandn(dv_local,dv_global,nc)
% A random nnz(dv_local) x nc block, scattered into the GLOBAL (dv_global)
% row space via 'ismember' -- zero rows for decision variables this
% particular sub-block (R0/R1/R2) does not use.
nloc = numel(dv_local);
Bloc = sprandn(nloc,nc,0.5);
[~,rows] = ismember(dv_local(:),dv_global(:));
B = sparse(numel(dv_global),nc);
B(rows,:) = Bloc;
end

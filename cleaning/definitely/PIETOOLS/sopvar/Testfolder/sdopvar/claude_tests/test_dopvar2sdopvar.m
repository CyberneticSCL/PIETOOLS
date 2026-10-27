%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_DOPVAR2SDOPVAR checks 'dopvar2sdopvar' (sopvar/Testfolder/converters)
% against the semantics of the input 'dopvar' object, never a round trip
% (CLAUDE.md S4): each block of a 'dopvar' (P, Q1, Q2, or R = {R0,R1,R2}) is
% built directly from the 'dpvar' constructor, and the resulting 'sdopvar's
% kernel, read out by 'pi_blk_kernels' (ANCHOR-style, as in
% test_sdvar2dpvar.m), is compared to the ORIGINAL dpvar block's own
% 'dpvar2poly' expansion evaluated by 'subs' -- the identical semantic check
% used throughout this converter family. 'pi_blk_kernels' always names a
% kernel's input-side variable with a '_dum' suffix (even when that
% direction is not shared with the output space), so Q1's reference
% evaluates at 's1_dum', matching the kernel's own 'varname'.
%
% Defects this catches (fixed 10/05/2026, dopvar2sdopvar.m header):
%  (1)/(2) the 'L2 to R' (Q1) and 'R to L2' (Q2) branches called
%      'dpvar2sdvar(D,vars)' with 'D' undefined -- a hard crash on first use;
%  (3) all three non-R branches (P, Q1, Q2) built 'params.A'/'params.B' as
%      bare matrices, where 'sdopvar' and 'canonicalize_multiplier' require
%      cell arrays (one entry per alpha) -- confirmed against two other real
%      callers (mat2copvar_grid.m, and the already-fixed opvar2sopvar.m,
%      whose R->R branch hit and fixed the analogous "non-constant block"
%      defect on 09/21/2026, which is why the 'R to R' branch here now also
%      passes an explicit empty 'vars' instead of nargin==1, so a spatially-
%      dependent P errors cleanly instead of being silently mishandled.
% None of the three was caught before because nothing in the repository
% calls this function yet.
%
% Sweep: 60 random trials per block type (P, Q1, Q2, R), random dims (1-3),
% random (possibly sparse/non-grid) monomial degrees for Q1/Q2/R, and random
% decision-variable subsets of a shared 4-name pool per sub-block of R (so
% R0/R1/R2 overlap, partially overlap, or share no decision variables).
% Plus: the 'R to R' guard, on a P block that (invalidly) depends on a
% spatial variable; and that a plain 'opvar' input (not just 'dopvar') is
% accepted and converted correctly for all four block types.
%
% Initial coding DJ, 10/05/2026
% DJ, 10/05/2026: Add OPVAR_INPUT, covering dopvar2sdopvar.m's 'isa(Pop,
%                  ''opvar'')' branch across all four block types.

tol = 1e-8;
nfail = 0;      ntot = 0;       fails = {};
dpool = {'d1','d2','d3','d4'};

%% GUARD: a spatially-dependent P (R to R) block must error, not misbuild
msg = '';
degmat = [0;1];
C = sprandn(1,4,0.8);
D = dpvar(C,degmat,{'s1'},cell(0,1),[1,2]);
Pop = dopvar();     Pop.P = D;
try
    dopvar2sdopvar(Pop);
    msg = 'expected an error, got none';
catch ME
    if ~contains(ME.message,'s1')
        msg = ['wrong error: ',ME.message];
    end
end
ntot = ntot+1;
if ~isempty(msg),  fails{end+1} = ['GUARD: ',msg];  end

%% OPVAR_INPUT: dopvar2sdopvar also accepts a plain 'opvar' (converts via      % DJ, 10/05/2026
%  opvar2dopvar internally); exercise all four block types directly.         % DJ, 10/05/2026
msg = check_opvar_input();                                                   % DJ, 10/05/2026
ntot = ntot+1;                                                                % DJ, 10/05/2026
if ~isempty(msg),  fails{end+1} = ['OPVAR_INPUT: ',msg];  end                % DJ, 10/05/2026

%% Four block types, 60 random trials each
for trial = 1:60
    rng(8000+trial);
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
assert(isempty(fails),'test_dopvar2sdopvar: %d of %d cases FAILED.',numel(fails),ntot);
fprintf('\ntest_dopvar2sdopvar passed (%d cases).\n',ntot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case1(dpool)
% P block: R^p -> R^m, no spatial variable.
msg = '';
try
    m = randi(3);   p = randi(3);
    ndk = randi([0,3]);
    dvn = dpool(randperm(4,ndk));
    C = sprandn((ndk+1)*m,p,0.8);
    D = dpvar(C,zeros(1,0),cell(0,1),dvn(:),[m,p]);
    Pop = dopvar();     Pop.P = D;
    Psop = dopvar2sdopvar(Pop);
    d = randn(ndk,1);
    [~,loc] = ismember(Psop.Zd,dvn(:));
    Kc = pi_blk_kernels(Psop,d(loc));
    if numel(Kc)~=1,    msg = 'expected 1 kernel';   return;     end
    got = double(Kc{1});
    ref = eval_poly(dpvar2poly(D),dvn(:),d);
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
    C = sprandn((ndk+1)*m,q*numel(degs),0.8);
    D = dpvar(C,degs,{'s1'},dvn(:),[m,q]);
    Pop = dopvar();     Pop.Q1 = D;
    Psop = dopvar2sdopvar(Pop);
    d = randn(ndk,1);
    [~,loc] = ismember(Psop.Zd,dvn(:));
    Kc = pi_blk_kernels(Psop,d(loc));
    if numel(Kc)~=1,    msg = 'expected 1 kernel';   return;     end
    s_val = rand_pt(1);
    got = eval_poly(Kc{1},{'s1_dum'},s_val);     % pi_blk_kernels suffixes vars.in with '_dum'
    ref = eval_poly(dpvar2poly(D),[dvn(:);{'s1'}],[d;s_val]);
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
    C = sprandn((ndk+1)*n,p*numel(degs),0.8);
    D = dpvar(C,degs,{'s1'},dvn(:),[n,p]);
    Pop = dopvar();     Pop.Q2 = D;
    Psop = dopvar2sdopvar(Pop);
    d = randn(ndk,1);
    [~,loc] = ismember(Psop.Zd,dvn(:));
    Kc = pi_blk_kernels(Psop,d(loc));
    if numel(Kc)~=1,    msg = 'expected 1 kernel';   return;     end
    s_val = rand_pt(1);
    got = eval_poly(Kc{1},{'s1'},s_val);
    ref = eval_poly(dpvar2poly(D),[dvn(:);{'s1'}],[d;s_val]);
    msg = cmp(got,ref);
catch ME
    msg = ['error: ',ME.message];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case4(dpool)
% R block: L2[s1] -> L2[s1], via R0(s1), R1(s1,s1_dum), R2(s1,s1_dum), each
% with an independently random subset of decision variables.
msg = '';
try
    n = randi(3);   q = randi(3);
    degmat0 = unique(randi([0,3],randi([1,3]),1));
    degmat12 = unique(randi([0,2],randi([2,5]),2),'rows');

    nd0 = randi([0,2]);    dv0 = dpool(randperm(4,nd0));
    nd1 = randi([0,2]);    dv1 = dpool(randperm(4,nd1));
    nd2 = randi([0,2]);    dv2 = dpool(randperm(4,nd2));

    C0 = sprandn((nd0+1)*n,q*numel(degmat0),0.8);
    D0 = dpvar(C0,degmat0,{'s1'},dv0(:),[n,q]);
    C1 = sprandn((nd1+1)*n,q*size(degmat12,1),0.8);
    D1 = dpvar(C1,degmat12,{'s1','s1_dum'},dv1(:),[n,q]);
    C2 = sprandn((nd2+1)*n,q*size(degmat12,1),0.8);
    D2 = dpvar(C2,degmat12,{'s1','s1_dum'},dv2(:),[n,q]);

    Pop = dopvar();
    Pop.R.R0 = D0;  Pop.R.R1 = D1;  Pop.R.R2 = D2;
    Psop = dopvar2sdopvar(Pop);

    alldv = union(union(dv0,dv1,'stable'),dv2,'stable');
    d = randn(numel(alldv),1);
    [~,loc] = ismember(Psop.Zd,alldv(:));
    Kc = pi_blk_kernels(Psop,d(loc));
    if numel(Kc)~=3,    msg = 'expected 3 kernels';  return;     end

    s_val = rand_pt(1);     t_val = rand_pt(1);
    [~,l0] = ismember(dv0(:),alldv(:));    d0 = d(l0);
    [~,l1] = ismember(dv1(:),alldv(:));    d1v = d(l1);
    [~,l2] = ismember(dv2(:),alldv(:));    d2v = d(l2);
    refs = {eval_poly(dpvar2poly(D0),[dv0(:);{'s1'}],[d0;s_val]), ...
            eval_poly(dpvar2poly(D1),[dv1(:);{'s1';'s1_dum'}],[d1v;s_val;t_val]), ...
            eval_poly(dpvar2poly(D2),[dv2(:);{'s1';'s1_dum'}],[d2v;s_val;t_val])};
    for k = 1:3
        got = eval_poly(Kc{k},{'s1';'s1_dum'},[s_val;t_val]);
        msg = cmp(got,refs{k});
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_opvar_input()
% dopvar2sdopvar(Pop) also accepts a fixed 'opvar' Pop directly (converting
% it via opvar2dopvar internally); ground truth is the 'opvar' block itself.
msg = '';
try
    opvar P1;   P1.P = randn(2,3);
    Psop = dopvar2sdopvar(P1);
    Kc = pi_blk_kernels(Psop,[]);
    if ~isequal(double(Kc{1}),P1.P)
        msg = 'P block: value mismatch';  return
    end

    pvar s1 s1_dum;
    opvar P2;   P2.Q1 = [1+s1, s1^2; 2, 3-s1];
    Psop = dopvar2sdopvar(P2);
    Kc = pi_blk_kernels(Psop,[]);
    s_val = 1.3;
    got = double(subs(Kc{1},{'s1_dum'},s_val));
    ref = double(subs(polynomial(P2.Q1),{'s1'},s_val));
    if norm(got-ref,'fro')>1e-8,  msg = 'Q1 block: value mismatch';  return;  end

    opvar P3;   P3.Q2 = [1+s1, s1^2; 2, 3-s1];
    Psop = dopvar2sdopvar(P3);
    Kc = pi_blk_kernels(Psop,[]);
    s_val = -0.6;
    got = double(subs(Kc{1},{'s1'},s_val));
    ref = double(subs(polynomial(P3.Q2),{'s1'},s_val));
    if norm(got-ref,'fro')>1e-8,  msg = 'Q2 block: value mismatch';  return;  end

    opvar P4;
    P4.R.R0 = [1+s1, 0; 0, 2-s1];
    P4.R.R1 = [s1*s1_dum, 1; 0, s1+s1_dum];
    P4.R.R2 = [0, s1-s1_dum; 1, 0];
    Psop = dopvar2sdopvar(P4);
    Kc = pi_blk_kernels(Psop,[]);
    s_val = 0.4;    t_val = -0.8;
    refs = {double(subs(polynomial(P4.R.R0),{'s1'},s_val)), ...
            double(subs(polynomial(P4.R.R1),{'s1';'s1_dum'},[s_val;t_val])), ...
            double(subs(polynomial(P4.R.R2),{'s1';'s1_dum'},[s_val;t_val]))};
    for k = 1:3
        got = double(subs(Kc{k},{'s1';'s1_dum'},[s_val;t_val]));
        if norm(got-refs{k},'fro')>1e-8
            msg = sprintf('R block k=%d: value mismatch',k);   return
        end
    end
catch ME
    msg = ['error: ',ME.message];
end
end

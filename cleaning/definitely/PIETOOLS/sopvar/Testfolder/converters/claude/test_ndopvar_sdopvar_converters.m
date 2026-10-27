%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_NDOPVAR_SDOPVAR_CONVERTERS checks the pair of converters between the
% 'ndopvar' and 'sdopvar' classes.
%
% The essential point is that a round trip test is NOT sufficient here. The
% two converters used to be exact inverses of one another while both carried
% the same index-mapping error, so 'ndopvar -> sdopvar -> ndopvar' returned
% the original object bit for bit even though the intermediate sdopvar was a
% DIFFERENT operator. The existing Test_ndopvar_to_sdopvar_conversion only
% exercises that cancelling direction and is blind to it.
%
% Each direction is therefore checked against the operator itself. The
% kernels of each object are read straight off its own class definition, by
% PI_NDOPVAR_KERNELS and PI_SDOPVAR_KERNELS, evaluated at the same decision
% variable values, and compared as functions. Neither evaluator goes through
% a converter, so agreement means the conversion preserved the operator and
% not merely that two errors cancelled.
%
% Also checked: both round trips at the coefficient level, and the
% preconditions that sdopvar2ndopvar must reject, since an ndopvar carries a
% single degree per variable with the complete basis 0:deg and gives a
% multiplier direction no dummy monomials at all.
%
% MMP, 08/29/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(21);
tol = 1e-9;
npass = 0;

% Each case is {label, m, n, deg, nd}
cases = {
    {'1D m1 n1 deg0 nd1', 1, 1, 0,     1}
    {'1D m1 n1 deg1 nd1', 1, 1, 1,     1}
    {'1D m1 n1 deg1 nd3', 1, 1, 1,     3}
    {'1D m1 n1 deg2 nd5', 1, 1, 2,     5}
    {'1D m2 n2 deg1 nd3', 2, 2, 1,     3}
    {'1D m2 n3 deg2 nd4', 2, 3, 2,     4}
    {'1D m3 n1 deg1 nd2', 3, 1, 1,     2}
    {'2D m1 n1 deg1 nd3', 1, 1, [1;1], 3}
    {'2D m2 n2 deg1 nd2', 2, 2, [1;1], 2}
    {'2D m1 n2 deg[2;1]', 1, 2, [2;1], 4}
    };

fprintf('=== each direction against the operator itself ===\n');
for ic = 1:numel(cases)
    lbl = cases{ic}{1};
    m = cases{ic}{2};   n = cases{ic}{3};
    deg = cases{ic}{4}; nd = cases{ic}{5};
    N = numel(deg);

    vnames = arrayfun(@(k) sprintf('s%d',k),(1:N)','UniformOutput',false);
    var1 = polynomial(vnames);
    var2 = polynomial(strcat(vnames,'_dum'));
    dnames = arrayfun(@(k) sprintf('r%d',k),(1:nd)','UniformOutput',false);
    dom = repmat([0,1],N,1);

    P = rand_ndopvar([m,n],deg,dom,var1,var2,dnames);
    for ii = 1:numel(P.C)
        P.C{ii} = sparse(randn(size(P.C{ii})));     % dense, so nothing is vacuous
    end
    dv = randn(nd,1);
    dumnames = reshape(strcat(vnames,'_dum'),1,[]);

    Knd = pi_ndopvar_kernels(P,dv,dumnames);

    % ndopvar -> sdopvar must preserve the operator
    Ps = ndopvar2sdopvar(P);
    Ksd = pi_sdopvar_kernels(Ps,remap(dv,dnames,Ps.Zd));
    e_fwd = kernel_diff(Knd,Ksd,vnames,dumnames);
    if ~(e_fwd<tol)
        error('test_ndopvar_sdopvar_converters: ndopvar2sdopvar changed the operator by %.3g (%s).',...
              e_fwd,lbl);
    end

    % sdopvar -> ndopvar must preserve it too
    Pn = sdopvar2ndopvar(Ps);
    Kn2 = pi_ndopvar_kernels(Pn,remap(dv,dnames,Pn.dvarname),dumnames);
    e_bwd = kernel_diff(Knd,Kn2,vnames,dumnames);
    if ~(e_bwd<tol)
        error('test_ndopvar_sdopvar_converters: sdopvar2ndopvar changed the operator by %.3g (%s).',...
              e_bwd,lbl);
    end

    % and the coefficients must come back unchanged
    e_rt = 0;
    for ii = 1:numel(P.C)
        e_rt = max(e_rt,max(max(abs(full(P.C{ii}-Pn.C{ii})))));
    end
    if ~(e_rt<tol)
        error('test_ndopvar_sdopvar_converters: round trip changed C by %.3g (%s).',e_rt,lbl);
    end

    fprintf('  passed: %-19s nd->sd %.0e | sd->nd %.0e | C round trip %.0e\n',...
            lbl,e_fwd,e_bwd,e_rt);
    npass = npass+3;
end

%% Starting from an sdopvar instead, so that A and B are exercised directly
fprintf('\n=== sdopvar -> ndopvar -> sdopvar preserves A and B ===\n');
for cfg = {[1,2,1],[1,2,3],[2,2,4],[1,3,5],[2,3,2]}
    m = cfg{1}(1);  nZ = cfg{1}(2);  nd = cfg{1}(3);
    P = canonical_sdopvar(m,nZ,nd);
    Pb = ndopvar2sdopvar(sdopvar2ndopvar(P));
    [~,loc] = ismember(cellstr(string(Pb.Zd(:))),cellstr(string(P.Zd(:))));
    eA = 0;     eB = 0;
    for k = 1:numel(P.params.A)
        eA = max(eA,max(abs(full(P.params.A{k}-Pb.params.A{k}))));
        B1 = full(P.params.B{k});
        B2 = full(Pb.params.B{k});      B2 = B2(loc,:);
        eB = max(eB,max(abs(B1(:)-B2(:))));
    end
    if ~(eA<tol && eB<tol)
        error(['test_ndopvar_sdopvar_converters: sdopvar round trip changed A by %.3g ',...
               'and B by %.3g (m=%d nZ=%d nd=%d).'],eA,eB,m,nZ,nd);
    end
    fprintf('  passed: m=%d nZ=%d nd=%d   max|dA| %.0e  max|dB| %.0e\n',m,nZ,nd,eA,eB);
    npass = npass+1;
end

%% Preconditions that sdopvar2ndopvar must reject
fprintf('\n=== rejected inputs ===\n');
probes = { {'gapped left basis',       @() sdopvar2ndopvar(odd_sdopvar('gap'))}
           {'ZL differs from ZR',      @() sdopvar2ndopvar(odd_sdopvar('mismatch'))} };
for k = 1:numel(probes)
    raised = false;
    try
        probes{k}{2}();
    catch ME
        raised = true;
        msg = strtrim(ME.message);
        if numel(msg)>62, msg = [msg(1:62),'...']; end
        fprintf('  rejects : %-24s %s\n',probes{k}{1},msg);
    end
    if ~raised
        error('test_ndopvar_sdopvar_converters: expected %s to be rejected.',probes{k}{1});
    end
    npass = npass+1;
end

% % % The third precondition, that a multiplier parameter carry no dummy
% % % variable degree, is no longer something the converter has to reject: it
% % % is a class invariant, and the 'sdopvar' constructor folds any input that
% % % violates it into canonical form. What used to be a rejection probe is
% % % therefore an assertion that the fold happened. The guard inside
% % % 'sdopvar2ndopvar' is kept as defense in depth against a parameter set
% % % assigned directly to the properties, which bypasses the constructor.
fprintf('\n=== a non-canonical multiplier is folded, not rejected ===\n');
w = warning('off','sdopvar:noncanonicalMultiplier');
cleanup = onCleanup(@() warning(w));                                        %#ok<NASGU>
Pnc = odd_sdopvar('noncanon');
[tf,info] = is_canonical_multiplier(Pnc.params,Pnc.vars,Pnc.ZL,Pnc.ZR,Pnc.dims);
if ~tf
    error('test_ndopvar_sdopvar_converters: the constructor should have canonicalized. %s',info.message);
end
fprintf('  passed: constructor folded it, ZL grew to %s\n',mat2str(Pnc.ZL{1}(:)'));
npass = npass+1;

fprintf('\ntest_ndopvar_sdopvar_converters passed (%d checks).\n',npass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = kernel_diff(K1,K2,vnames,dumnames)
% Largest difference between two sets of kernels, sampled at a few points in
% the primary and dummy variables.

if numel(K1)~=numel(K2)
    error('Kernel sets have different numbers of parameters.');
end
pts = [0.31,0.77; 0.62,0.14; 0.05,0.93];
e = 0;
for j = 1:numel(K1)
    D = polynomial(K1{j}) - polynomial(K2{j});
    for r = 1:size(pts,1)
        Dr = D;
        for i = 1:numel(vnames)
            Dr = subs(Dr,polynomial(vnames(i)),pts(r,1)+0.07*i);
            Dr = subs(Dr,polynomial(dumnames(i)),pts(r,2)-0.05*i);
        end
        e = max(e,max(abs(double(Dr(:)))));
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dd = remap(dv,names_from,names_to)
% Re-express a decision variable assignment in another ordering.

nf = cellstr(string(names_from(:)));
nt = cellstr(string(names_to(:)));
[tf,loc] = ismember(nt,nf);
dd = zeros(numel(nt),1);
dd(tf) = dv(loc(tf));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = canonical_sdopvar(m,nZ,nd)
% A random sdopvar in the form an ndopvar can represent: complete bases
% 0:nZ-1 on both sides, and a multiplier parameter carrying no dummy degree.

ZL = {(0:nZ-1)'};   ZR = ZL;
NL = nZ;            NR = nZ;
len = (m*NL)*(m*NR);
Zd = arrayfun(@(k) sprintf('d%d',k),(1:nd)','UniformOutput',false);

A = cell(3,1);      B = cell(3,1);
for k = 1:3
    keep = true(NL,NR);
    if k==1
        keep(:) = false;    keep(:,1) = true;       % multiplier: no dummy degree
    end
    mask = repmat(keep,m,m);
    Ak = randn(m*NL,m*NR).*mask;
    A{k} = sparse(Ak(:));
    Bk = sparse(nd,len);
    for r = 1:nd
        Br = randn(m*NL,m*NR).*mask;
        Bk(r,:) = Br(:).';                                                  %#ok<SPRIX>
    end
    B{k} = Bk;
end

vv = struct();      vv.in = {'s1'};     vv.out = {'s1'};
dd = struct();      dd.in = [0,1];      dd.out = [0,1];
P = sdopvar(struct('A',{A},'B',{B}),vv,Zd,ZL,ZR,dd,[m,m]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = odd_sdopvar(kind)
% An sdopvar that an ndopvar cannot represent faithfully.

switch kind
    case 'gap'
        ZL = {[0;2]};       ZR = {[0;2]};
    case 'mismatch'
        ZL = {[0;1;2]};     ZR = {[0;1]};
    otherwise
        ZL = {[0;1]};       ZR = {[0;1]};
end
NL = numel(ZL{1});      NR = numel(ZR{1});
Zd = {'d1'};
A = cell(3,1);      B = cell(3,1);
for k = 1:3
    Ak = randn(NL,NR);
    if strcmp(kind,'noncanon') && k==1
        Ak(:,end) = 1;      % multiplier depending on the dummy variable
    elseif k==1
        Ak(:,2:end) = 0;
    end
    A{k} = sparse(Ak(:));
    B{k} = sparse(1,NL*NR);
end
vv = struct();      vv.in = {'s1'};     vv.out = {'s1'};
dd = struct();      dd.in = [0,1];      dd.out = [0,1];
P = sdopvar(struct('A',{A},'B',{B}),vv,Zd,ZL,ZR,dd,[1,1]);

end

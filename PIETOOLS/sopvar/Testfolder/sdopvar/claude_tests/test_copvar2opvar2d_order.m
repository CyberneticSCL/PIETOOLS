%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR2OPVAR2D_ORDER checks 'copvar2opvar2d' for every assignment of
% the registry variables to opvar2d's x and y roles (vars_xy), not only the
% sorted one. Before 09/28/2026 each block was converted by
% 'sopvar2opvar2d', which picks x and y BY ITS OWN RULE (sorted; a lone 's2'
% or 'y' is y), and the component of the same NAME was copied out. Where the
% two orientations differ the component was read from the wrong slot (left
% empty: the block was dropped) or, for R22, copied with its alpha axes
% swapped.
%
% Per case, against the semantics (CLAUDE.md S4), never a round trip:
%   ACTION  for every populated block C{i,j} and a random polynomial test
%           function x_j on s^j,
%             C{i,j}*x_j   from the SOPVAR DEFINITION (sopvar.m header;
%                          sopvar.pdf S4): kernel
%                          (I kron ZL(s2,s3))' C_alpha (I kron ZR(s3',s1)),
%                          S1 integrated out, alpha over S3 in vars.in order
%                          (canonical [S3,S1], S3 sorted);
%           equals the opvar2d component at slot R<out><in> applied to x_j
%           from the OPVAR2D DEFINITION (@opvar2d header, apply_opvar
%           formulas): cell axis d is direction d = var1(d), alpha 1/2/3 =
%           multiplier / int_a^s / int_s^b. The slot is found here from
%           vars_xy BY NAME, independently of the routine under test.
%           Integral kernels in every direction: a multiplier-only test sees
%           only R22{1,1}, where a transpose is invisible.
%   EMPTY   every slot no block maps to is empty or zero.
%   META    var1/var2 name vars_xy (and '_dum'), I(k,:) is the domain of
%           vars_xy{k}; dim rows/columns sit at the role-derived index.
%   OMIT    vars_xy omitted = vars_xy the sorted registry, byte for byte.
%   ROUND   (secondary only) opvar2d2copvar of the result acts as C by the
%           same sopvar evaluator. opvar2d2copvar is the inverse routine, so
%           this is never the sole evidence.
%
% Domains differ per variable (x in [0,1], y in [-1,2]) so a swapped I row
% changes every integral. Registries {a,b} and {y,z} put the lone-variable
% rule of sopvar2opvar2d against sorted order, so they fail the DEFAULT
% order too before 09/28/2026. A one-variable registry is converted with
% both roles named: a one-entry vars_xy gives a 1x1 var1, which opvar2d
% rejects (not addressed here).
%
% Initial coding MMP, 09/28/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','sopvar:noncanonicalMultiplier');
for f = {'copvar2opvar2d','sopvar2opvar2d','opvar2d2copvar','rand_copvar',...
         'pi_monom_vector'}
    nw = numel(which(f{1},'-all'));
    assert(nw==1,'test_copvar2opvar2d_order: %s resolves to %d files.',f{1},nw);
end
tol = 1e-9;

% {registry, domains (row k for registry{k}), spaces out, spaces in,
%  dim_out, dim_in, deg, list of vars_xy orders, label}
S4 = @(a,b) {{},{a},{b},{a,b}};     % R, L2[a], L2[b], L2[a,b]
CASES = {
 {'x','y'},   [0 1;-1 2], S4('x','y'),   S4('x','y'),   [2;1;2;1],[1;2;1;2],1, {{'x','y'},{'y','x'}},   'x,y deg 1';
 {'x','y'},   [0 1;-1 2], S4('x','y'),   S4('x','y'),   [1;2;1;2],[2;1;2;1],2, {{'x','y'},{'y','x'}},   'x,y deg 2';
 {'s1','s2'}, [0 1;-1 2], S4('s1','s2'), S4('s1','s2'), [2;1;1;2],[1;2;2;1],1, {{'s1','s2'},{'s2','s1'}},'s1,s2 deg 1';
 {'a','b'},   [0 1;-1 2], S4('a','b'),   S4('a','b'),   [1;2;1;2],[2;1;2;1],1, {{'a','b'},{'b','a'}},   'a,b (lone b)';
 {'y','z'},   [-1 2;0 1], S4('y','z'),   S4('y','z'),   [2;1;2;1],[1;2;1;2],1, {{'y','z'},{'z','y'}},   'y,z (lone y)';
 {'x','y'},   [0 1;-1 2], {{'y'},{'x','y'}}, {{},{'x'},{'y'}}, [2;1],[1;2;1],1, {{'x','y'},{'y','x'}}, 'subgrid';
 {'y'},       [-1 2],     {{},{'y'}},    {{},{'y'}},    [2;1],[1;2],1,         {{'x','y'},{'y','x'}},   'registry {y}';
 {'b'},       [-1 2],     {{},{'b'}},    {{},{'b'}},    [2;1],[1;2],1,         {{'a','b'},{'b','a'}},   'registry {b}'};

nchk = 0;   nfail = 0;
for c = 1:size(CASES,1)
    [reg,dm,so,si,dout,din,dg,orders,lbl] = CASES{c,:};
    rng(4200+c);
    sp = struct('out',{so},'in',{si},'vars',{reg});
    Pm = rand_copvar(sp,struct('out',dout,'in',din),dm,dg,0.7);
    % Test functions, one per input column of the container.
    X = cell(1,numel(si));
    for j = 1:numel(si)
        X{j} = testfun(Pm.vars(Pm.space_in(j,:)),Pm.dim_in(j));
    end
    % Container action per block, from the sopvar definition (computed once;
    % independent of vars_xy).
    Y = cell(size(Pm.C));
    for k = 1:numel(Pm.C)
        [i,j] = ind2sub(size(Pm.C),k);
        if ~isempty(Pm.C{k}),   Y{k} = act_sopvar(Pm.C{k},X{j});   end
    end
    for o = 1:numel(orders)
        vxy = orders{o};
        nchk = nchk+1;
        try
            Pop = copvar2opvar2d(Pm,vxy);
        catch ME
            % opvar2d's setter rejects a kernel in the wrong dummy, which is
            % how a swapped R22 with dummy dependence surfaced before.
            nfail = nfail+1;
            fprintf('  FAILED: %-14s vars_xy {%s}: conversion error: %s\n',...
                lbl,strjoin(vxy,','),ME.message);
            continue
        end
        [ok,msg] = check_case(Pm,Pop,vxy,X,Y,tol);
        if ok
            fprintf('  passed: %-14s vars_xy {%s}%s\n',lbl,strjoin(vxy,','),msg);
        else
            nfail = nfail+1;
            fprintf('  FAILED: %-14s vars_xy {%s}: %s\n',lbl,strjoin(vxy,','),msg);
        end
        % OMIT: default = sorted registry, byte for byte.
        if isequal(vxy,Pm.vars)
            nchk = nchk+1;
            if ~isequal(getByteStreamFromArray(copvar2opvar2d(Pm)),...
                        getByteStreamFromArray(Pop))
                nfail = nfail+1;
                fprintf('  FAILED: %-14s vars_xy omitted differs from sorted\n',lbl);
            end
        end
        % ROUND (secondary): only where opvar2d2copvar can take var1 (2 vars).
        if numel(vxy)==2
            try
                Pr = opvar2d2copvar(Pop);
                [ok,msg] = check_roundtrip(Pm,Pr,X,Y,tol);
            catch ME
                ok = false;     msg = ME.message;
            end
            nchk = nchk+1;
            if ~ok
                nfail = nfail+1;
                fprintf('  FAILED: %-14s vars_xy {%s} round trip: %s\n',lbl,strjoin(vxy,','),msg);
            end
        end
    end
end
if nfail>0
    error('test_copvar2opvar2d_order:failures','%d of %d checks failed.',nfail,nchk);
end
fprintf('test_copvar2opvar2d_order passed (%d checks).\n',nchk);


%% ========================== case check ================================
function [ok,msg] = check_case(Pm,Pop,vxy,X,Y,tol)
NM = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
ok = true;  msg = '';
% META: var1, var2, I and dim by role.
vn = pvar2varname(Pop.var1);    dn = pvar2varname(Pop.var2);
if ~isequal(vn(:)',vxy) || ~isequal(dn(:)',strcat(vxy,'_dum'))
    ok = false;     msg = [msg ' var1/var2 not vars_xy;'];
end
for k = 1:numel(vxy)
    if ~any(strcmp(Pm.vars,vxy{k})),    continue,   end   % unused role
    if ~isequal(Pop.I(k,:),Pm.dom(strcmp(Pm.vars,vxy{k}),:))
        ok = false;     msg = [msg sprintf(' I row %d;',k)];
    end
end
ro = role_index(Pm.space_out,Pm.vars,vxy);
ci = role_index(Pm.space_in ,Pm.vars,vxy);
d = zeros(4,2);     d(ro,1) = Pm.dim_out;   d(ci,2) = Pm.dim_in;
if ~isequal(Pop.dim,d)
    ok = false;     msg = [msg ' dim;'];
end
% ACTION per block, and which slots are targeted.
used = false(4,4);  emax = 0;   bad = {};
for k = 1:numel(Pm.C)
    if isempty(Pm.C{k}),    continue,   end
    [i,j] = ind2sub(size(Pm.C),k);
    used(ro(i),ci(j)) = true;
    slot = NM{ro(i),ci(j)};
    try
        Yo = act_opvar2d(Pop,slot,X{j});
        e = relerr(Y{k},Yo);
    catch ME
        e = Inf;    slot = [slot '(' ME.identifier ')'];
    end
    emax = max(emax,e);
    if ~(e<tol),    bad{end+1} = sprintf('%s %.2g',slot,e);   end %#ok<AGROW>
end
% EMPTY: untargeted slots carry nothing.
for k = find(~used(:))'
    if ~comp_zero(Pop.(NM{k}))
        bad{end+1} = [NM{k} ' nonempty'];   %#ok<AGROW>
    end
end
if ~isempty(bad)
    ok = false;     msg = [msg ' ' strjoin(bad,', ')];
end
if ok,  msg = sprintf(', %d blocks, max rel err %.1e',nnz(used),emax);  end
end

function [ok,msg] = check_roundtrip(Pm,Pr,X,Y,tol)
% Each block of Pr against the original's action, rows and columns matched
% by space (Pr is in opvar2d order, Pm in its own).
ok = true;  msg = '';   bad = {};
ro = match_rows(Pm.space_out,Pm.vars,Pr.space_out,Pr.vars);
ci = match_rows(Pm.space_in ,Pm.vars,Pr.space_in ,Pr.vars);
for k = 1:numel(Pm.C)
    [i,j] = ind2sub(size(Pm.C),k);
    B = Pr.C{ro(i),ci(j)};
    if isempty(Pm.C{k})
        if ~isempty(B) && relerr(act_sopvar(B,X{j}),0)>tol
            bad{end+1} = sprintf('(%d,%d) nonzero',i,j);  %#ok<AGROW>
        end
        continue
    end
    if isempty(B)
        e = relerr(Y{k},0);
    else
        e = relerr(Y{k},act_sopvar(B,X{j}));
    end
    if ~(e<tol),    bad{end+1} = sprintf('(%d,%d) %.2g',i,j,e);   end %#ok<AGROW>
end
if ~isempty(bad),   ok = false;     msg = strjoin(bad,', ');    end
end


%% ====================== the two definitions ===========================
function y = act_sopvar(B,x)
% (B x)(s2,s3) = sum_alpha int_S1 int_S3' I_alpha(s3-s3') K_alpha x(s3',s1),
% K_alpha = (I_m kron ZL(vars.out))' C_alpha (I_n kron ZR(vars.in')),
% coefficient rows (matrix row, monomial) monomial inner, columns likewise.
vin = B.vars.in;    vout = B.vars.out;  din = B.dom.in;
m = B.dims(1);      n = B.dims(2);
vdum = strcat(vin,'_dum');
ZL = pi_monom_vector(B.ZL,vout);
ZR = pi_monom_vector(B.ZR,vdum);
xd = polynomial(x);
for k = 1:numel(vin)
    xd = subs(xd,pv(vin{k}),pv(vdum{k}));
end
is3 = ismember(vin,vout);           % S3 in vars.in order, S1 the rest
k3 = find(is3);     k1 = find(~is3);
y = polynomial(zeros(m,size(x,2)));
for p = 1:numel(B.params)
    Cp = B.params{p};
    if isempty(Cp) || ~any(Cp(:)),  continue,   end
    Cp = full(Cp);
    K = kron(eye(m),ZL.')*Cp*kron(eye(n),ZR);
    t = K*xd;
    for k = k1
        t = int(t,pv(vdum{k}),din(k,1),din(k,2));
    end
    al = gamma_index(p,numel(k3));
    for l = 1:numel(k3)
        k = k3(l);      s = pv(vin{k});     sd = pv(vdum{k});
        switch al(l)
            case 1,     t = subs(t,sd,s);
            case 2,     t = int(t,sd,din(k,1),s);
            case 3,     t = int(t,sd,s,din(k,2));
        end
    end
    y = y + t;
end
end

function y = act_opvar2d(P,slot,x)
% Component R<o><i> of opvar2d P applied to x. Direction d is var1(d) on
% [I(d,1),I(d,2)], dummy var2(d). A direction in both spaces is a PI
% direction whose alpha is the cell index along cell axis d (Rxx/Rx2/R2x
% 3x1, Ryy/Ry2/R2y 1x3, R22 3x3); one in the input only is integrated out.
code = {'0',[];'x',1;'y',2;'2',[1 2]};
do = code{strcmp(code(:,1),slot(2)),2};
di = code{strcmp(code(:,1),slot(3)),2};
dp = reshape(intersect(di,do),1,[]);    % row: a 0x1 loop range runs once
df = reshape(setdiff(di,do),1,[]);
comp = P.(slot);
if ~iscell(comp),   comp = {comp};  end
y = 0;
for k = 1:numel(comp)
    Rk = comp{k};
    if isempty(Rk),     continue,   end
    sub = cell(1,2);    [sub{:}] = ind2sub(size(comp),k);
    xt = polynomial(x);
    for d = dp
        if sub{d}~=1,   xt = subs(xt,P.var1(d),P.var2(d));  end
    end
    t = Rk*xt;
    for d = dp
        switch sub{d}
            case 2,     t = int(t,P.var2(d),P.I(d,1),P.var1(d));
            case 3,     t = int(t,P.var2(d),P.var1(d),P.I(d,2));
        end
    end
    for d = df
        t = int(t,P.var1(d),P.I(d,1),P.I(d,2));
    end
    y = y + t;
end
end


%% ============================ helpers ================================
function idx = role_index(mask,vars,vxy)
% opvar2d index (1 R, 2 x, 3 y, 4 xy) of each space, x = vxy{1}, y = vxy{2}.
idx = ones(size(mask,1),1);
for i = 1:size(mask,1)
    sv = vars(mask(i,:));
    idx(i) = 1 + any(strcmp(sv,vxy{1})) + 2*(numel(vxy)>1 && any(strcmp(sv,vxy{2})));
end
end

function r = match_rows(maskA,varsA,maskB,varsB)
% Row of B over the same variable set as each row of A.
r = zeros(size(maskA,1),1);
for i = 1:size(maskA,1)
    sa = sort(varsA(maskA(i,:)));
    for k = 1:size(maskB,1)
        if isequal(sa,sort(varsB(maskB(k,:)))),  r(i) = k;  end
    end
end
end

function x = testfun(v,p)
% p x 1 polynomial of degree <= 2 in each variable of v, random coefficients.
if isempty(v)
    x = polynomial(randn(p,1));     return
end
Z = pi_monom_vector(repmat({(0:2)'},1,numel(v)),v);
x = randn(p,numel(Z))*Z;
end

function al = gamma_index(k,n3)
% Multi-index of linear index k of a 3 x ... x 3 cell, first axis fastest.
al = ones(1,n3);    k = k-1;
for i = 1:n3,   al(i) = mod(k,3)+1;     k = floor(k/3);     end
end

function p = pv(nm)
p = polynomial(1,1,{nm},[1,1]);
end

function tf = comp_zero(c)
if iscell(c)
    tf = all(cellfun(@comp_zero,c(:)));     return
end
if isempty(c),  tf = true;  return,     end
c = polynomial(c);
tf = isempty(c.coefficient) || ~any(c.coefficient(:));
end

function e = relerr(A,B)
% max |coef(A-B)| / max(1,max |coef(A)|).
A = polynomial(A);  B = polynomial(B);
D = A-B;
if isempty(D.coefficient) || ~any(D.coefficient(:)),   e = 0;  return,  end
sA = 1;
if ~isempty(A.coefficient),     sA = max(1,full(max(abs(A.coefficient(:)))));  end
e = full(max(abs(D.coefficient(:))))/sA;
end

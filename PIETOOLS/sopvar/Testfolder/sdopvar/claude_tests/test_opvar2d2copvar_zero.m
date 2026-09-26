%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_OPVAR2D2COPVAR_ZERO checks 'opvar2d2copvar' on opvar2d objects with
% zero components and zero-dimensional spaces. A zero component becomes a
% [] block; before 09/26/2026 a row or column holding only zeros was left
% with no block and the conversion errored (copvar:emptyRow/emptyColumn):
% D11 = 0 and Tw = 0 of a PIE, a B1/C1 with a zero part, eppos with zero
% entries. 'Test_copvar_converters' (section 5) covers dense operators only.
%
% Per case, against the semantics (CLAUDE.md S4), never a round trip:
%   SHAPE    one row per nonzero-dimensional output space, one column per
%            input space, opvar2d order; dims, variables and domains from
%            the opvar2d's own dim, var1 and I;
%   BLOCKS   a nonzero component is never []; a zero component is a block
%            only where its row or column has no nonzero component;
%   VALID    'verify' passes; copvar(C) rebuilds the same metadata from the
%            blocks alone;
%   KERNELS  every block's kernels (pi_blk_kernels) equal
%            (a) those of 'opvar2d2sopvar' on an opvar2d holding only that
%                component (places each component in its own cell), and
%            (b) the opvar2d DEFINITION: the stored component, input-only
%                directions renamed to their dummy, the alpha index of each
%                direction matched to the block's multi-index BY VARIABLE
%                NAME (pi_gamma_index over sorted S3).
%            A zero block must give zero kernels. In a multiplier direction
%            (gamma_k = 1) all sides are evaluated on the diagonal s_k' = s_k,
%            where the kernel is only determined up to moving degree between
%            s_k and s_k' (as in test_poscopvar_vs_poslpivar_2d's
%            sample_container).
%
% Zero components are identified from their coefficients, and for the
% hand-built cases also cross-checked against the list that was zeroed.
%
% Cases with var1 not in sorted name order (x named 's2') run check (b) as
% well: they test the sorted registry of the stated metadata, and R22's
% alpha axes, which opvar2d2sopvar swapped for such names until its
% 09/26/2026 fix (kernel error 0.83 against the definition by name).
%
% Initial coding MMP, 09/26/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','sopvar:noncanonicalMultiplier');
for f = {'opvar2d2copvar','opvar2d2sopvar','pi_blk_kernels','pi_gamma_index',...
         'rand_opvar2d','cx_plant'}
    nw = numel(which(f{1},'-all'));
    assert(nw==1,'test_opvar2d2copvar_zero: %s resolves to %d files.',f{1},nw);
end
NM = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
pvar s1 s2 s1_dum s2_dum
V1 = [s1;s2];   V2 = [s1_dum;s2_dum];
dom = [0 1;-1 2];               % unequal: a swapped registry row shows
dR  = [2 1;1 2;2 1;1 3];        % every space present, rectangular
tol = 1e-9;

% {opvar2d, label, zeroed components (cellstr), or [] if not declared}
CASES = cell(0,3);
% zero-dimensional spaces dropped
CASES(end+1,:) = {mk(dR,{},V1,V2,dom,1),                'dense rectangular', {}};
CASES(end+1,:) = {mk([0 1;2 0;0 0;1 2],{},V1,V2,dom,2), 'rows x,xy; cols R,xy', {}};
CASES(end+1,:) = {mk([2 3;0 0;0 0;0 0],{},V1,V2,dom,3), 'R only, empty registry', {}};
CASES(end+1,:) = {mk([0 0;0 0;1 2;0 0],{},V1,V2,dom,4), 'L2[y] only', {}};
CASES(end+1,:) = {mk([0 0;1 0;0 1;0 0],{},V1,V2,dom,5), 'Rxy only, L2[y] -> L2[x]', {}};
% each component zero in turn: no row or column emptied
for k = 1:16
    [i,j] = ind2sub([4,4],k);
    CASES(end+1,:) = {mk(dR,NM(i,j),V1,V2,dom,10+k),['zero ' NM{i,j}], NM(i,j)}; %#ok<SAGROW>
end
% rows emptied, columns emptied, both
for r = 1:4
    CASES(end+1,:) = {mk(dR,NM(r,:),V1,V2,dom,30+r),sprintf('zero row %d',r), NM(r,:)}; %#ok<SAGROW>
end
for c = 1:4
    CASES(end+1,:) = {mk(dR,NM(:,c)',V1,V2,dom,40+c),sprintf('zero column %d',c), NM(:,c)'}; %#ok<SAGROW>
end
for rc = [1 4; 2 3; 4 1; 3 3]'
    zc = union(NM(rc(1),:),NM(:,rc(2))');
    CASES(end+1,:) = {mk(dR,zc,V1,V2,dom,50+rc(1)),...
        sprintf('zero row %d and column %d',rc(1),rc(2)), zc};             %#ok<SAGROW>
end
CASES(end+1,:) = {mk(dR,NM(:)',V1,V2,dom,60),                 'fully zero', NM(:)'};
CASES(end+1,:) = {mk([2 3;0 0;0 0;0 0],{'R00'},V1,V2,dom,61), 'fully zero, R only', {'R00'}};
CASES(end+1,:) = {mk([0 0;0 0;0 0;2 1],{'R22'},V1,V2,dom,62), 'fully zero, L2[x,y] only', {'R22'}};
% eppos-like diagonal multiplier with zero entries: [1e-3, 1e-3, 0, 0]
Pe = opvar2d();     Pe.I = dom;     Pe.var1 = V1;   Pe.var2 = V2;
Pe.dim = [1 1;2 2;1 1;2 2];         % zero-filled by set.dim
Pe.R00 = 1e-3;      c = Pe.Rxx;     c{1} = polynomial(1e-3*eye(2));     Pe.Rxx = c;
CASES(end+1,:) = {Pe, 'eppos [1e-3 1e-3 0 0]', setdiff(NM(:)',{'R00','Rxx'})};
% var1 not in sorted name order: registry order and domains must follow
% the names, not var1's order (see NOT COVERED)
zc = union(NM(1,:),NM(:,4)');
CASES(end+1,:) = {mk(dR,zc,[s2;s1],[s2_dum;s1_dum],dom,70), 'reversed names, zero row 1, col 4', zc};
CASES(end+1,:) = {mk([0 0;1 2;0 0;0 0],{},[s2;s1],[s2_dum;s1_dum],dom,71), 'reversed names, L2[x] only', {}};
CASES(end+1,:) = {mk([1 1;0 0;2 1;0 0],{'R00','Ry0'},[s2;s1],[s2_dum;s1_dum],dom,72),...
    'reversed names, zero column R', {'R00','Ry0'}};
% the io2 plant: D11 = 0 and Tw = 0
PIE = cx_plant('io2');
for f = {'T','Tw','A','B1','C1','D11'}
    CASES(end+1,:) = {PIE.(f{1}), ['io2 PIE.' f{1}], []};                   %#ok<SAGROW>
end

nfail = 0;  nchk = 0;
for ic = 1:size(CASES,1)
    [Pop,lbl,zc] = deal(CASES{ic,:});
    try
        [n,nz] = check_case(Pop,zc,NM,tol);
        nchk = nchk+n;
        fprintf('  passed: %-28s dim %-26s %d zero block(s)\n',lbl,mat2str(Pop.dim),nz);
    catch ME
        nfail = nfail+1;
        fprintf('  FAILED: %-28s %s: %s\n',lbl,ME.identifier,ME.message);
    end
end
% Kept from 09/18/2026: an operator with no input space is rejected.
try
    opvar2d2copvar(PIE.D12);
    nfail = nfail+1;    fprintf('  FAILED: PIE.D12 (no input space) accepted\n');
catch ME
    if strcmp(ME.identifier,'opvar2d2copvar:empty')
        nchk = nchk+1;  fprintf('  passed: PIE.D12 (no input space) rejected\n');
    else
        nfail = nfail+1;    fprintf('  FAILED: PIE.D12: %s\n',ME.message);
    end
end
ntot = size(CASES,1)+1;
if nfail>0
    error('test_opvar2d2copvar_zero:failures','%d of %d cases failed.',nfail,ntot);
end
fprintf('test_opvar2d2copvar_zero passed (%d cases, %d checks).\n',ntot,nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,nzb] = check_case(Pop,zc,NM,tol)
% All checks on one opvar2d; errors on the first failure. n = checks made,
% nzb = blocks holding a zero component.
d = Pop.dim;
rows = d(:,1)>0;    cols = d(:,2)>0;
ro = find(rows)';   co = find(cols)';
vn = reshape(pvar2varname(Pop.var1),1,[]);
SP = {{}, vn(1), vn(2), vn};            % variables of R, L2[x], L2[y], L2[x,y]
defn = true;                            % check (b) for every case, see header

% Zero components, from the coefficients; cross-check the declared list.
Z = false(4,4);
for i = ro,     for j = co,     Z(i,j) = comp_zero(Pop.(NM{i,j}));   end,    end
if iscell(zc)
    Zd = false(4,4);
    for i = ro,     for j = co,     Zd(i,j) = ismember(NM{i,j},zc);  end,    end
    assert(isequal(Z,Zd),'test setup: zeroed components are not the declared ones');
end

Pm = opvar2d2copvar(Pop);

% SHAPE
assert(isequal(size(Pm),[numel(ro),numel(co)]),'grid is %s',mat2str(size(Pm)));
assert(isequal(Pm.dim_out(:),d(rows,1)) && isequal(Pm.dim_in(:),d(cols,2)),'dims');
for a = 1:numel(ro)
    assert(isempty(setxor(Pm.vars(Pm.space_out(a,:)),SP{ro(a)})),'row %d space',a);
end
for b = 1:numel(co)
    assert(isempty(setxor(Pm.vars(Pm.space_in(b,:)),SP{co(b)})),'column %d space',b);
end
for k = 1:numel(Pm.vars)
    assert(isequal(Pm.dom(k,:),Pop.I(strcmp(vn,Pm.vars{k}),:)),'domain of %s',Pm.vars{k});
end

% BLOCKS
nzb = 0;
for a = 1:numel(ro)
    for b = 1:numel(co)
        i = ro(a);  j = co(b);
        if ~Z(i,j)
            assert(~isempty(Pm.C{a,b}),'nonzero %s is []',NM{i,j});
        elseif ~isempty(Pm.C{a,b})
            nzb = nzb+1;
            assert(all(Z(i,co)) || all(Z(ro,j)),...
                'zero %s is a block though its row and column have nonzero ones',NM{i,j});
        end
    end
end

% VALID
v = verify(Pm);
assert(v.true==1,'verify: %s',strjoin(v.flags(:)',' | '));
assert(isequal(metadata(copvar(Pm.C)),metadata(Pm)),'copvar(C) does not rebuild the metadata');
n = 7;

% KERNELS
for a = 1:numel(ro)
    for b = 1:numel(co)
        B = Pm.C{a,b};
        if isempty(B),  continue,   end
        i = ro(a);  j = co(b);
        Kc = pi_blk_kernels(B,[]);
        % (a) the single-block converter on this component alone
        Pb = opvar2d();
        Pb.I = Pop.I;   Pb.var1 = Pop.var1;     Pb.var2 = Pop.var2;
        sel = zeros(4,2);   sel(i,1) = d(i,1);  sel(j,2) = d(j,2);
        Pb.dim = sel;       Pb.(NM{i,j}) = Pop.(NM{i,j});
        Kr = pi_blk_kernels(opvar2d2sopvar(Pb),[]);
        % (b) the definition
        comp = Pop.(NM{i,j});
        if ~iscell(comp),   comp = {comp};  end
        S3 = intersect(B.vars.in,B.vars.out);
        assert(isequal(size(Kr),size(Kc)) && numel(Kc)==3^numel(S3) ...
            && numel(comp)==numel(Kc),'%s: %d cells for %d shared variables',...
            NM{i,j},numel(Kc),numel(S3));
        ionly = setdiff(SP{j},SP{i});
        tx = find(strcmp(S3,vn{1}));    ty = find(strcmp(S3,vn{2}));
        for q = 1:numel(Kc)
            gam = pi_gamma_index(q,numel(S3));
            ax = 1;     ay = 1;
            if ~isempty(tx),    ax = gam(tx);   end
            if ~isempty(ty),    ay = gam(ty);   end
            Kd = polynomial(comp{ax,ay});
            if isempty(Kd),     Kd = polynomial(zeros(B.dims));   end
            for t = 1:numel(ionly)
                Kd = subs(Kd,polynomial(ionly(t)),polynomial({[ionly{t} '_dum']}));
            end
            Ka = Kc{q};     Kb = Kr{q};
            for t = find(gam==1)
                dv = polynomial({[S3{t} '_dum']});      pv = polynomial(S3(t));
                Ka = subs(Ka,dv,pv);    Kb = subs(Kb,dv,pv);    Kd = subs(Kd,dv,pv);
            end
            ea = perr(Ka,Kb);   eb = perr(Ka,Kd);
            assert(ea<tol,'%s cell %d differs from opvar2d2sopvar by %.2e',NM{i,j},q,ea);
            assert(eb<tol || ~defn,'%s cell %d differs from the definition by %.2e',NM{i,j},q,eb);
            n = n+1+defn;
        end
    end
end
end


function P = mk(d,zc,v1,v2,dom,seed)
% Random opvar2d, degree 1, with the listed components set to zero.
rng(seed);
P = rand_opvar2d(d,1,dom,v1,v2);
for k = 1:numel(zc)
    c = P.(zc{k});
    if iscell(c)
        for q = 1:numel(c),     c{q} = 0*c{q};  end
    else
        c = 0*c;
    end
    P.(zc{k}) = c;
end
end


function tf = comp_zero(c)
% True when every coefficient of a component (cells included) is zero.
if ~iscell(c),  c = {c};    end
tf = true;
for k = 1:numel(c)
    p = polynomial(c{k});
    if ~isempty(p) && ~isempty(p.coefficient) && any(p.coefficient(:))
        tf = false;     return
    end
end
end


function e = perr(A,B)
D = polynomial(A) - polynomial(B);
e = full(max(abs(D.coefficient(:))));
if isempty(e),  e = 0;  end
end
